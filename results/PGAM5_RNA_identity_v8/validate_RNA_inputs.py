"""Validate future RNA identity/known-mixture inputs; no protein fields.

This computes RNA agreement and cell-fraction recovery. Self-reported provenance
and calibration IDs need scientific review; the CLI never certifies TCGA use.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

SAMPLE_FIELDS = ['sample_id','patient_id','split','independent_patient_confirmed',
    'previously_used_for_model','RNA_platform','RNA_expression_unit','QC_status',
    'frozen_RNA_definition_id','independent_RNA_evidence_id']
IDENTITY_FIELDS = ['cell_id','sample_id','assigned_lineage','PGAM5_raw_counts',
    'independent_RNA_truth_available','independent_RNA_truth_macrophage',
    'independent_RNA_truth_PGAM5_detected','RNA_identity_method','QC_status']
MIXTURE_FIELDS = ['mixture_id','sample_id','preparation_replicate','mixture_kind',
    'target_cells','other_macrophage_cells','monocyte_cells','DC_cells','tumor_cells',
    'other_cells','total_cells','actual_target_fraction','estimated_target_fraction',
    'physical_mixture','bulk_platform','bulk_expression_unit','reference_expression_unit',
    'unit_calibration_evidence_id','QC_status']
COUNT_FIELDS = ['target_cells','other_macrophage_cells','monocyte_cells','DC_cells','tumor_cells','other_cells']


def boolean(series,field):
    values=series.astype(str).str.lower()
    assert values.isin(['true','false','1','0']).all(), f'Invalid or missing boolean: {field}'
    return values.isin(['true','1']).to_numpy()


def load(path,fields):
    d=pd.read_csv(path,sep='\t')
    assert set(fields)<=set(d.columns), f'Missing fields: {set(fields)-set(d.columns)}'
    return d


def validate(samples,identity,mix,output,computational_replay=False):
    output.mkdir(parents=True,exist_ok=True)
    status={'protein_validation_required':False,'RNA_scope':'RNA expression only',
        'usable_for_validated_TCGA_PGAM5_cell_abundance':False,
        'automatic_scientific_identity_or_unit_certification':False}
    if samples.empty or identity.empty or mix.empty:
        status.update(status='NOT_RUN_REQUIRED_RNA_INPUTS_MISSING',sample_rows=len(samples),
            identity_rows=len(identity),mixture_rows=len(mix),mixture_gates_passed=False)
        (output/'RNA_input_validation.json').write_text(json.dumps(status,indent=2)+'\n')
        return status
    assert samples.sample_id.notna().all() and samples.sample_id.is_unique
    assert samples.patient_id.notna().all() and samples.patient_id.astype(str).str.len().gt(0).all()
    assert samples.split.isin(['discovery','validation']).all()
    # A patient's samples cannot leak across discovery and validation.
    assert samples.groupby('patient_id').split.nunique().le(1).all(), 'Patient overlap across splits'
    independent=boolean(samples.independent_patient_confirmed,'independent_patient_confirmed')
    previously=boolean(samples.previously_used_for_model,'previously_used_for_model')
    samples=samples.copy(); samples['eligible_validation']=samples.split.eq('validation') & independent & ~previously & samples.QC_status.eq('PASS')
    assert identity[['sample_id','cell_id']].notna().all().all()
    assert not identity.duplicated(['sample_id','cell_id']).any(), 'Duplicate sample/cell key'
    assert identity.sample_id.isin(samples.sample_id).all(), 'Unknown identity sample ID'
    count=pd.to_numeric(identity.PGAM5_raw_counts,errors='raise').to_numpy(float)
    assert np.isfinite(count).all() and (count>=0).all() and (count==np.floor(count)).all()
    truth_available=boolean(identity.independent_RNA_truth_available,'independent_RNA_truth_available')
    truth_mac=boolean(identity.independent_RNA_truth_macrophage,'independent_RNA_truth_macrophage')
    truth_pg=boolean(identity.independent_RNA_truth_PGAM5_detected,'independent_RNA_truth_PGAM5_detected')
    # An independent method/evidence ID must be supplied for agreement assessment.
    join=identity.merge(samples,on='sample_id',validate='many_to_one',suffixes=('','_sample'))
    valid=truth_available & join.eligible_validation.to_numpy() & join.QC_status.eq('PASS').to_numpy()
    valid &= join.RNA_identity_method.fillna('').str.len().gt(0).to_numpy()
    valid &= join.independent_RNA_evidence_id.fillna('').str.len().gt(0).to_numpy()
    observed=identity.assigned_lineage.eq('Macrophage').to_numpy() & (count>0)
    truth=truth_mac & truth_pg
    rows=[]
    for patient in sorted(join.loc[valid,'patient_id'].unique()):
        ix=valid & join.patient_id.eq(patient).to_numpy()
        tp=int((observed[ix]&truth[ix]).sum()); fn=int((~observed[ix]&truth[ix]).sum())
        fp=int((observed[ix]&~truth[ix]).sum()); tn=int((~observed[ix]&~truth[ix]).sum())
        rows.append({'patient_id':patient,'cells':int(ix.sum()),'true_positive':tp,'false_negative':fn,
            'false_positive':fp,'true_negative':tn,'RNA_target_sensitivity':tp/(tp+fn) if tp+fn else np.nan,
            'RNA_target_specificity':tn/(tn+fp) if tn+fp else np.nan})
    pd.DataFrame(rows,columns=['patient_id','cells','true_positive','false_negative','false_positive','true_negative',
        'RNA_target_sensitivity','RNA_target_specificity']).to_csv(output/'RNA_identity_agreement_by_patient.csv',index=False)
    assert mix.mixture_id.notna().all() and mix.mixture_id.is_unique
    assert mix.sample_id.isin(samples.sample_id).all(), 'Unknown mixture sample ID'
    assert mix.mixture_kind.isin(['standard','dose_sensitivity','zero_target_competitor']).all()
    counts=mix[COUNT_FIELDS+['total_cells']].apply(pd.to_numeric,errors='raise').to_numpy(float)
    assert np.isfinite(counts).all() and (counts>=0).all() and (counts==np.floor(counts)).all()
    assert (counts[:,-1]>0).all() and np.array_equal(counts[:,:-1].sum(1),counts[:,-1]), 'Mixture cell totals disagree'
    actual=counts[:,0]/counts[:,-1]
    declared=pd.to_numeric(mix.actual_target_fraction,errors='raise').to_numpy(float)
    predicted=pd.to_numeric(mix.estimated_target_fraction,errors='raise').to_numpy(float)
    assert np.isfinite(declared).all() and np.allclose(actual,declared,atol=1e-10,rtol=0), 'Target cell fraction disagrees with counts'
    assert np.isfinite(predicted).all() and ((predicted>=0)&(predicted<=1)).all(), 'Predictions must be finite cell-fraction-scale coefficients'
    assert (actual[mix.mixture_kind.eq('zero_target_competitor')]==0).all()
    physical=boolean(mix.physical_mixture,'physical_mixture')
    m=mix.merge(samples,on='sample_id',validate='many_to_one',suffixes=('','_sample'))
    m['actual']=actual; m['predicted']=predicted
    evaluated=m.eligible_validation & m.QC_status.eq('PASS')
    if computational_replay:
        evaluated=m.QC_status.eq('PASS')
    assert evaluated.any(), 'No eligible validation mixture rows'
    rows=[]; gates=[]
    expected=np.array([0,.005,.01,.02,.05])
    for patient,d in m[evaluated].groupby('patient_id',sort=True):
        standard=d[d.mixture_kind.eq('standard')]
        # Average technical preparations within a patient/dose, not extra patients.
        average=standard.groupby('actual',sort=True).predicted.mean()
        mae=np.mean(np.abs(average.index.to_numpy()-average.to_numpy()))*100 if len(average) else np.nan
        correlation=pearsonr(average.index.to_numpy(),average.to_numpy()).statistic if len(average)>1 and average.nunique()>1 else np.nan
        zero=d[d.actual.eq(0)].predicted.to_numpy()*100
        p95=float(np.quantile(zero,.95)) if len(zero) else np.nan
        dose_cover=all(np.isclose(average.index.to_numpy(),dose,atol=1e-10,rtol=0).any() for dose in expected)
        challenge=bool(d.mixture_kind.eq('zero_target_competitor').any())
        rows.append({'patient_id':patient,'mixture_rows':len(d),'distinct_standard_doses':len(average),
            'standard_MAE_pp':mae,'standard_Pearson_r':correlation,'all_zero_P95_percent':p95,
            'dose_coverage_0_0.5_1_2_5_percent':dose_cover,'zero_competitor_challenge_present':challenge})
        for gate,value in [('MAE_le1pp',bool(mae<=1)),('r_ge0.7',bool(correlation>=.7)),
            ('zero_P95_le0.5percent',bool(p95<=.5)),('required_doses_present',dose_cover),('zero_competitor_challenge_present',challenge)]:
            gates.append({'patient_id':patient,'gate':gate,'passed':value})
    pd.DataFrame(rows).to_csv(output/'mixture_recovery_by_patient.csv',index=False)
    pd.DataFrame(gates).to_csv(output/'mixture_gates.csv',index=False)
    enough=len(rows)>=3
    numerical=bool(enough and all(g['passed'] for g in gates))
    metadata_physical=bool(np.all(physical[evaluated.to_numpy()]))
    calibration=bool(m.loc[evaluated,'unit_calibration_evidence_id'].fillna('').str.len().gt(0).all())
    frozen=bool(m.loc[evaluated,'frozen_RNA_definition_id'].fillna('').str.len().gt(0).all())
    status.update(status='COMPUTATIONAL_REPLAY_NOT_INDEPENDENT' if computational_replay else 'EVALUATED_REQUIRES_SCIENTIFIC_REVIEW',
        evaluated_patient_or_sample_labels=len(rows),
        eligible_validation_patients=m.loc[m.eligible_validation,'patient_id'].nunique(),
        eligible_independent_RNA_identity_patients=len(set(join.loc[valid,'patient_id'])),
        numerical_mixture_gates_passed=numerical,physical_mixtures_declared=metadata_physical,
        calibration_evidence_ids_present=calibration,frozen_RNA_definition_ids_present=frozen,
        declared_metadata_only_not_independently_verified=True,
        mixture_gates_passed=bool(numerical and metadata_physical and calibration and frozen and not computational_replay),
        note='Passing recovery does not establish a stable RNA state, verify external evidence, or calibrate TCGA by itself. Review independent RNA identity, reference specificity, raw evidence, and platform/RNA units.')
    (output/'RNA_input_validation.json').write_text(json.dumps(status,indent=2)+'\n')
    return status


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--samples',type=Path,required=True)
    parser.add_argument('--RNA-identity',dest='identity',type=Path,required=True)
    parser.add_argument('--mixtures',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--computational-replay',action='store_true',help='Evaluate historical/computational rows; never certifies independent or physical validation')
    args=parser.parse_args()
    result=validate(load(args.samples,SAMPLE_FIELDS),load(args.identity,IDENTITY_FIELDS),load(args.mixtures,MIXTURE_FIELDS),args.output_dir,args.computational_replay)
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
