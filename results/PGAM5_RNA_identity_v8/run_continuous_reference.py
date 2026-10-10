"""Add fixed continuous-ranking primary analysis, without rerunning finished DE."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import json
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from pooled_RNA_reference import reference,evaluate

R=Path(__file__).resolve().parent
P=json.loads((R/'protocol.json').read_text());S=Path(P['source_root'])
UNITS=['equalized_cell_fraction','library_RNA_contribution']


def main():
    split=pd.read_csv(R/'pooled_patient_sample_split.csv')
    ixparts=[];values=[]
    for cohort in P['primary_cohorts']:
        ind=pd.read_csv(S/f'{cohort}_group_stats_index.csv')
        membership=split[split.dataset.eq(cohort)].set_index('donor_or_sample').split
        ind['split']=ind.donor.map(membership)
        assert ind.split.notna().all()
        ixparts.append(ind);values.append(np.load(S/f'{cohort}_group_statistics.npz')['values'])
    index=pd.concat(ixparts,ignore_index=True);stats=np.concatenate(values)
    dimensions=[]
    for name in ['POOLED_INTERNAL_TRAINING','POOLED_ALL_TRAINING']:
        folder=R/name/'continuous_rank';folder.mkdir(exist_ok=True)
        de=pd.read_csv(R/name/'pooled_DE.csv.gz')
        chosen=np.flatnonzero(index.split.eq('training').to_numpy()) if name=='POOLED_INTERNAL_TRAINING' else np.arange(len(index))
        ind=index.iloc[chosen].reset_index(drop=True);v=stats[chosen]
        for unit in UNITS:
            info=reference(ind,v,de,folder,unit,None)
            dimensions.append({'model':name,'threshold_version':'continuous_rank','unit':unit,
                'macrophage_DE_cells':int(de.RNA_detected_cells.iloc[0]+de.RNA_undetected_cells.iloc[0]),
                'PGAM5_RNA_detected_DE_cells':int(de.RNA_detected_cells.iloc[0]),**info})
        if name=='POOLED_INTERNAL_TRAINING':
            meta=pd.read_csv(R/name/'log2FC_1/held_sample_mixture_metadata.csv')
            meta.to_csv(folder/'held_sample_mixture_metadata.csv',index=False)
            common=pd.Index(pd.read_csv(S/'common_genes.csv').gene)
            x=np.load(S/'all_mixture_expression.npz')['values'];queries={}
            for unit in UNITS:
                genes=np.load(folder/f'{unit}_model.npz')['genes']
                ids=meta.loc[meta.unit.eq(unit),'mixture_id'].to_numpy(dtype=int)
                queries[unit]=x[ids][:,common.get_indexer(genes)]
            np.savez_compressed(folder/'actual_held_sample_mixtures.npz',**queries)
    dim=pd.read_csv(R/'pooled_model_dimensions.csv');dim=dim[~dim.threshold_version.eq('continuous_rank')]
    pd.concat([dim,pd.DataFrame(dimensions)],ignore_index=True).to_csv(R/'pooled_model_dimensions.csv',index=False)
    with ProcessPoolExecutor(max_workers=2) as pool:
        new=pd.concat(list(pool.map(evaluate,[('continuous_rank',u) for u in UNITS])),ignore_index=True)
    predictions=pd.read_csv(R/'pooled_internal_mixture_predictions.csv.gz')
    predictions=pd.concat([predictions[~predictions.threshold_version.eq('continuous_rank')],new],ignore_index=True)
    predictions.to_csv(R/'pooled_internal_mixture_predictions.csv.gz',index=False)
    rows=[];gates=[];doses=[]
    for (level,unit),d in predictions.groupby(['threshold_version','unit'],sort=True):
        standard=d[d.scenario.eq('standard')];zero=d[d.truth.eq(0)]
        mae=float(np.abs(standard.truth-standard.predicted).mean()*100)
        corr=float(pearsonr(standard.truth,standard.predicted).statistic) if standard.predicted.nunique()>1 else np.nan
        p95=float(zero.predicted.quantile(.95)*100)
        patients=standard.loc[~standard.dataset.eq('GSE202642'),'donor_or_sample'].nunique()
        rows.append({'threshold_version':level,'unit':unit,'mixture_rows':len(d),'standard_rows':len(standard),
            'standard_patient_or_sample_labels':standard.donor_or_sample.nunique(),
            'standard_confirmed_patient_labels_excluding_GSE202642':patients,'standard_MAE_pp':mae,
            'standard_Pearson_r':corr,'all_zero_P95_percent':p95,'all_zero_max_percent':zero.predicted.max()*100,
            'physical_mixtures':False,'fresh_external_validation':False,'cross_cohort_replication_required':False,
            'usable_as_validated_TCGA_cell_abundance':False})
        for gate,passed in [('MAE_le1pp',mae<=1),('r_ge0.7',corr>=.7),('zero_P95_le0.5percent',p95<=.5),('internal_confirmed_patients_ge3',patients>=3)]:
            gates.append({'threshold_version':level,'unit':unit,'gate':gate,'passed':bool(passed)})
        for dose,part in standard.groupby('nominal_target_cell_fraction'):
            doses.append({'threshold_version':level,'unit':unit,'nominal_cell_fraction':dose,'mixtures':len(part),
                'patient_or_sample_labels':part.donor_or_sample.nunique(),'median_truth':part.truth.median(),
                'median_prediction':part.predicted.median(),'prediction_q25':part.predicted.quantile(.25),
                'prediction_q75':part.predicted.quantile(.75),'fraction_estimate_above0.5percent':part.predicted.gt(.005).mean()})
    pd.DataFrame(rows).to_csv(R/'pooled_internal_validation_summary.csv',index=False)
    pd.DataFrame(gates).to_csv(R/'pooled_internal_validation_gates.csv',index=False)
    pd.DataFrame(doses).to_csv(R/'pooled_internal_dose_recovery.csv',index=False)
    result=json.loads((R/'validation_result.json').read_text())
    result.update(primary_version='continuous_rank',FDR_hard_cutoff_required=False,log2FC_hard_cutoff_required=False,
        cross_cohort_replication_required=False,protein_validation_required=False,
        pooled_internal_mixture_gates_passed=all(x['passed'] for x in gates if x['threshold_version']=='continuous_rank'),
        pooled_state_gene_counts={x['threshold_version']+'_'+x['unit']:int(x['state_genes']) for x in pd.concat([dim,pd.DataFrame(dimensions)]).to_dict('records') if x['model']=='POOLED_ALL_TRAINING'})
    (R/'validation_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(pd.DataFrame(dimensions).to_string(index=False),flush=True)
    print(pd.DataFrame(rows).to_string(index=False),flush=True)


if __name__=='__main__':main()
