"""Exploratory TCGA coefficients; expressly NOT calibrated cell infiltration."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from reference_core import fit, T, N

R = Path(__file__).resolve().parent


def main():
    unit = 'library_RNA_contribution'
    ref = pd.read_csv(R/f'EXPLORATORY_reference_{unit}.tsv', sep='\t', index_col=0)
    scales = pd.read_csv(R/f'{unit}_row_scales.csv').set_index('gene').scale
    lock = json.loads((R/'locked_reference.json').read_text())
    source = np.load(R/'TCGA_primary_common_gene_counts.npz')
    counts = pd.DataFrame(source['counts'], index=source['genes'], columns=source['patients'])
    meta = pd.read_csv(R/'TCGA_primary_aliquot_selection.csv').set_index('patient').loc[counts.columns]
    measured = ref.index[ref.index.isin(counts.index)]
    assert 'PGAM5' in measured
    coverage = len(measured)/len(ref)
    missing = ref.index.difference(measured)
    pd.DataFrame({'gene': missing}).to_csv(R/'TCGA_missing_reference_genes.csv', index=False)
    query = counts.loc[measured].div(meta.full_library_total_counts, axis=1)*10000
    query.to_csv(R/'TCGA_used_reference_gene_CP10k.csv.gz')
    counts.loc[measured].to_csv(R/'TCGA_used_reference_gene_raw_counts.csv.gz')
    ref.loc[measured].to_csv(R/'TCGA_used_reference.tsv', sep='\t')
    scales.loc[measured].rename('scale').to_csv(R/'TCGA_used_row_scales.csv')
    A = ref.loc[measured].to_numpy()/scales.loc[measured].to_numpy()[:, None]
    model = {'matrix': A, 'features': np.arange(len(measured)), 'scale': scales.loc[measured].to_numpy()}
    target = ref.columns.get_loc(T)
    mac = [j for j, group in enumerate(ref.columns) if group in [T, N, 'Cycling_'+N]]
    without = np.flatnonzero(measured != 'PGAM5')
    without_model = {'matrix': A[without], 'features': without, 'scale': model['scale'][without]}
    coefficients, sensitivities, diagnostics, cases = [], [], [], []
    for j, patient in enumerate(query.columns):
        y = query[patient].to_numpy(float)
        c, residual, weights = fit(model, y, lock['solver'])
        no_pg, no_pg_residual, _ = fit(without_model, y, lock['solver'])
        coefficients.append(c)
        sensitivities.append(no_pg)
        diagnostics.append({'patient': patient, 'relative_weighted_residual': residual,
            'reference_feature_coverage': coverage, 'feature_coverage_ge95pct': coverage>=.95,
            'reference_validation_passed': False, 'bulk_platform_calibrated': False,
            'cell_RNA_content_calibrated': False, 'usable_as_cellular_infiltration': False})
        if j % 74 == 0:
            case = len(cases)
            np.savez_compressed(R/f'TCGA_solver_audit_case_{case}.npz', A=A, b=y/model['scale'],
                coefficients=c, weights=weights, solver=lock['solver'], target_index=target)
            cases.append({'case': case, 'patient': patient})
    table = pd.DataFrame(coefficients, index=query.columns, columns=ref.columns)
    table.index.name = 'patient'
    table.to_csv(R/'EXPLORATORY_TCGA_all_component_coefficients.csv')
    no_pg = pd.DataFrame(sensitivities, index=query.columns, columns=ref.columns)
    no_pg.index.name = 'patient'
    no_pg.to_csv(R/'TCGA_without_PGAM5_component_coefficients.csv')
    out = pd.DataFrame({'patient': table.index, 'RNA_sample': meta['sample'].to_numpy(),
        'exploratory_PGAM5_RNA_macrophage_coefficient': table.iloc[:, target].to_numpy(),
        'exploratory_total_macrophage_coefficient': table.iloc[:, mac].sum(axis=1).to_numpy(),
        'coefficient_without_PGAM5_anchor': no_pg.iloc[:, target].to_numpy(),
        'usable_as_cellular_infiltration': False})
    out['exploratory_target_within_macrophage_coefficient'] = out.exploratory_PGAM5_RNA_macrophage_coefficient/out.exploratory_total_macrophage_coefficient.replace(0, np.nan)
    out.to_csv(R/'EXPLORATORY_TCGA_coefficients_NOT_CELL_FRACTIONS.csv', index=False)
    pd.DataFrame(diagnostics).to_csv(R/'TCGA_fit_diagnostics.csv', index=False)
    pd.DataFrame(cases).to_csv(R/'TCGA_solver_audit_case_manifest.csv', index=False)
    rho = spearmanr(out.exploratory_PGAM5_RNA_macrophage_coefficient, out.coefficient_without_PGAM5_anchor).statistic
    result = {'patients': len(out), 'solver': lock['solver'], 'reference_genes': len(ref),
        'bulk_matched_genes': len(measured), 'coverage_fraction': coverage,
        'components': len(ref.columns), 'target_coefficients_gt1e_8': int(out.exploratory_PGAM5_RNA_macrophage_coefficient.gt(1e-8).sum()),
        'numerical_positive_tolerance': 1e-8,
        'median_exploratory_target_coefficient': float(out.exploratory_PGAM5_RNA_macrophage_coefficient.median()),
        'with_vs_without_PGAM5_anchor_Spearman': float(rho),
        'reference_validation_passed': False, 'valid_cellular_abundance_obtained': False,
        'output_unit': 'Exploratory relative RNA mixture coefficient on unmatched FeatureCounts/scRNA platforms; no valid conversion to cell fraction',
        'interpretation': 'Outputs are exploratory coefficients, not validated cellular infiltration. User-requested exploratory coefficient-survival associations are reported separately; significance cannot establish cell abundance or PGAM5 specificity.'}
    (R/'TCGA_trial_result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
