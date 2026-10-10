"""Apply the saved v8 solver to measured rows; freeze grouping before outcomes."""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
from guard_solver import fit

R = Path(__file__).resolve().parent
UNITS = ['library_RNA_contribution', 'equalized_cell_fraction']

def main():
    query = pd.read_csv(R/'TCGA_reference_union_CP10k.csv.gz', index_col=0)
    meta = pd.read_csv(R/'TCGA_primary_aliquot_selection.csv').set_index('patient').loc[query.columns]
    protocol = json.loads((R/'protocol.json').read_text())
    summaries = []
    for unit in UNITS:
        dest = R/unit
        dest.mkdir(exist_ok=True)
        modelpath = R/'reference'/unit/f'{unit}_model.npz'
        model = np.load(modelpath)
        genes, groups = model['genes'], model['groups']
        keep = pd.Index(genes).isin(query.index)
        measured = genes[keep]
        assert 'PGAM5' in measured and len(groups) == 20 and len(genes) == 395
        A, scales = model['A'][keep], model['scale'][keep]
        y = query.loc[measured].to_numpy(float)
        assert np.isfinite(y).all() and (y >= 0).all()
        pg = int(np.flatnonzero(measured == 'PGAM5')[0])
        target, mac = model['target_mask'], model['macrophage_mask']
        np.savez_compressed(dest/'used_reference_model.npz', A=A, scale=scales, genes=measured,
                            groups=groups, target_mask=target, macrophage_mask=mac)
        pd.DataFrame(A*scales[:,None], index=pd.Index(measured,name='gene'), columns=groups).to_csv(dest/'used_reference.tsv',sep='\t')
        pd.DataFrame({'gene':measured,'scale':scales}).to_csv(dest/'used_row_scales.csv',index=False)
        pd.DataFrame({'gene':genes[~keep]}).to_csv(dest/'unmeasured_reference_genes.csv',index=False)
        rows, cases = [], []
        for j, patient in enumerate(query.columns):
            b = y[:,j]/scales
            c, w, residual = fit(A,b,True,pg,target)
            row = {'patient':patient, 'RNA_sample':meta.loc[patient,'sample'],
                   'target_coefficient':float(c[target].sum()),
                   'total_macrophage_coefficient':float(c[mac].sum()),
                   'relative_unweighted_L2_residual':residual,
                   'target_PGAM5_scaled_contribution':float((A[pg]*target)@c),
                   'observed_PGAM5_scaled':float(b[pg]),
                   'usable_as_validated_cell_abundance':False}
            row.update({'component_'+str(g):float(v) for g,v in zip(groups,c)})
            rows.append(row)
            if j % 74 == 0:
                casefile = f'solver_case_{j:03d}.npz'
                np.savez_compressed(dest/casefile, A=A,b=b,c=c,weights=w,
                                    pg_index=pg,target_mask=target,groups=groups)
                cases.append({'patient':patient,'file':casefile})
            if (j+1)%50 == 0 or j+1 == len(query.columns): print(unit,j+1,'/',len(query.columns),flush=True)
        table = pd.DataFrame(rows)
        cut = float(table.target_coefficient.median())
        table['group'] = np.where(table.target_coefficient > cut,'High','Low')
        table.to_csv(dest/'patient_coefficients.csv',index=False)
        pd.DataFrame(cases).to_csv(dest/'solver_case_manifest.csv',index=False)
        summary = {'unit':unit, 'patients':len(table), 'original_reference_genes':len(genes),
            'used_genes':len(measured), 'gene_coverage':len(measured)/len(genes),
            'components':len(groups), 'global_median':cut,
            'high_patients':int(table.group.eq('High').sum()),'low_patients':int(table.group.eq('Low').sum()),
            'target_coefficient_gt_1e_8':int(table.target_coefficient.gt(1e-8).sum()),
            'target_groups':groups[target].tolist(), 'RNA_budget_max_violation':float((table.target_PGAM5_scaled_contribution-table.observed_PGAM5_scaled).max()),
            'candidate_state_genes_measured':sum(g in measured for g in json.loads((R/'reference'/unit/f'{unit}_definition.json').read_text())['candidate_state_genes']),
            'reference_validation_passed':False, 'platform_and_cell_RNA_calibrated':False,
            'output_unit':'Exploratory uncalibrated reference mixture coefficient, not a true cell fraction',
            'reference_model_sha256':hashlib.sha256(modelpath.read_bytes()).hexdigest()}
        summaries.append(summary)
        print(json.dumps(summary),flush=True)
    pd.DataFrame(summaries).to_csv(R/'deconvolution_summary.csv',index=False)
    (R/'frozen_grouping.json').write_text(json.dumps({'method':protocol['grouping'],
        'frozen_before_endpoint_analysis':True,'units':summaries},indent=2)+'\n')

if __name__ == '__main__': main()
