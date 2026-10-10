"""All-cell Spearman/Pearson correlations and a fixed rank permutation test."""
from pathlib import Path
import json
import platform
import importlib.metadata as md
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr, rankdata, beta
from statsmodels.stats.multitest import multipletests

R = Path(__file__).resolve().parent
B = 200000
SEED = 20261010


def main():
    protocol = {
        'primary': 'All 7448 existing primaryT macrophages; Spearman on log1p(CP10k), retaining all zeros; two-sided asymptotic p and absolute-rho Monte Carlo permutation p',
        'secondary_family': ['Pearson_all_cells_log1p_CP10k', 'Spearman_all_cells_raw_counts',
                             'Pearson_all_cells_raw_counts'],
        'multiple_testing': 'One designated primary gene-pair Spearman test; BH over the three secondary asymptotic p-values. No selection of results by significance.',
        'permutation': {'shuffles': B, 'seed': SEED, 'statistic': 'absolute Spearman rho',
                        'formula': '(extreme_count+1)/(B+1)',
                        'null': 'Shuffle MYO19 ranks across all cells; cells assumed exchangeable',
                        'uncertainty': 'Clopper-Pearson 95% binomial interval for the tail probability; Monte Carlo estimate, not an exact enumeration'},
        'zero_handling': 'No PGAM5>0, MYO19>0 or co-detection filter; zeros remain numeric zeros on all scales and in every test',
        'normalization': 'raw_count / stored full-library total_counts * 10000, then natural log1p; no matching, covariate regression or partial correlation',
        'scope': 'Previous primary HCC T sample scope and latest conservative macrophage labels; all cycling macrophages included; existing source QC retained',
        'limitations': 'Cell-level tests; same-patient cells are not independent patient replicates. Sparse counts, tied zeros, normalization and unequal sample cell counts limit biological interpretation. No causal/protein-interaction inference.'
    }
    (R / 'protocol.json').write_text(json.dumps(protocol, indent=2) + '\n')
    d = pd.read_csv(R / 'macrophage_expression.csv.gz')
    assert len(d) == 7448 and d.cell_id.is_unique
    x = d.PGAM5_log1p_CP10k.to_numpy()
    y = d.MYO19_log1p_CP10k.to_numpy()
    specifications = [
        ('Spearman_all_cells_log1p_CP10k', x, y, spearmanr, 'primary'),
        ('Pearson_all_cells_log1p_CP10k', x, y, pearsonr, 'secondary'),
        ('Spearman_all_cells_raw_counts', d.PGAM5_raw_count.to_numpy(),
         d.MYO19_raw_count.to_numpy(), spearmanr, 'secondary'),
        ('Pearson_all_cells_raw_counts', d.PGAM5_raw_count.to_numpy(),
         d.MYO19_raw_count.to_numpy(), pearsonr, 'secondary')]
    rows = []
    for name, xx, yy, method, family in specifications:
        z = method(xx, yy)
        rows.append({'analysis': name, 'unit': 'cell', 'n': len(xx),
                     'correlation': float(z.statistic), 'asymptotic_p': float(z.pvalue),
                     'family': family, 'secondary_BH_q': None,
                     'permutation_p': None, 'permutations': None,
                     'extreme_permutations': None})
    primary = rows[0]
    rx = rankdata(x) - (len(x) + 1) / 2
    ry = rankdata(y) - (len(y) + 1) / 2
    denominator = np.linalg.norm(rx) * np.linalg.norm(ry)
    assert np.isclose(np.dot(rx, ry) / denominator, primary['correlation'])
    rng = np.random.default_rng(SEED)
    extreme = 0
    for _ in range(B):
        statistic = np.dot(rx, rng.permutation(ry)) / denominator
        extreme += abs(statistic) >= abs(primary['correlation']) - 1e-12
    primary.update(permutation_p=(extreme + 1) / (B + 1),
                   permutations=B, extreme_permutations=int(extreme))
    interval = [float(beta.ppf(.025, extreme, B - extreme + 1)) if extreme else 0.,
                float(beta.ppf(.975, extreme + 1, B - extreme)) if extreme < B else 1.]
    qs = multipletests([z['asymptotic_p'] for z in rows[1:]], method='fdr_bh')[1]
    for row, q in zip(rows[1:], qs):
        row['secondary_BH_q'] = float(q)
    statistics = pd.DataFrame(rows)
    statistics.to_csv(R / 'correlation_statistics.csv', index=False)
    pg = d.PGAM5_raw_count.gt(0)
    my = d.MYO19_raw_count.gt(0)
    counts = {'cells': len(d), 'PGAM5_detected': int(pg.sum()), 'MYO19_detected': int(my.sum()),
              'both_detected': int((pg & my).sum()), 'neither_detected': int((~pg & ~my).sum()),
              'PGAM5_only': int((pg & ~my).sum()), 'MYO19_only': int((~pg & my).sum()),
              'cycling_macrophages': int(d.cycling_flag.sum())}
    (R / 'detection_summary.json').write_text(json.dumps(counts, indent=2) + '\n')
    summary = d.groupby('sample', sort=True).agg(
        macrophages=('cell_id', 'size'), patient=('patient', 'first'),
        PGAM5_detected=('PGAM5_detected', 'sum'), MYO19_detected=('MYO19_detected', 'sum'),
        PGAM5_mean_CP10k=('PGAM5_CP10k', 'mean'), MYO19_mean_CP10k=('MYO19_CP10k', 'mean'))
    summary['both_detected'] = d.assign(both=pg & my).groupby('sample')['both'].sum()
    summary['neither_detected'] = d.assign(neither=~pg & ~my).groupby('sample')['neither'].sum()
    summary.to_csv(R / 'sample_summary.csv')
    d.groupby(['PGAM5_raw_count', 'MYO19_raw_count'], sort=True).size().rename(
        'cells').reset_index().to_csv(R / 'raw_count_pairs.csv', index=False)
    versions = {'python': platform.python_version(), **{k: md.version(k) for k in
                ['numpy', 'pandas', 'scipy', 'statsmodels', 'matplotlib', 'anndata', 'h5py']}}
    (R / 'environment_versions.json').write_text(json.dumps(versions, indent=2) + '\n')
    (R / 'requirements.txt').write_text('\n'.join(k + '==' + v for k, v in versions.items()
                                                 if k != 'python') + '\n')
    result = {'dataset': 'GSE149614', 'scope': '10 primary HCC tumor T samples',
              'n_cells': len(d), 'n_samples': d['sample'].nunique(), 'n_patients': d.patient.nunique(),
              'all_zero_values_retained': True, 'detection_counts': counts,
              'Spearman_rho': primary['correlation'], 'Spearman_asymptotic_p': primary['asymptotic_p'],
              'Spearman_permutation_p': primary['permutation_p'], 'permutation_extremes': int(extreme),
              'permutation_tail_probability_95pct_CI': interval,
              'Pearson_r': rows[1]['correlation'], 'Pearson_p': rows[1]['asymptotic_p'],
              'interpretation': 'Very weak positive rank association in pooled cells; normalized Pearson is not significant. P-values are cell-level, not evidence of independent patient-level association or causal/protein interaction.'}
    (R / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    from plot_correlation import render
    render(R)
    print(statistics.to_string(index=False))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
