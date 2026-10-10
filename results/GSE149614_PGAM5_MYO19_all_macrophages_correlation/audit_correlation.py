"""Independently audit full-cell membership, HDF5 CSR counts, scales and statistics."""
from pathlib import Path
import argparse
import hashlib
import json
import h5py
import numpy as np
import pandas as pd
from scipy.stats import rankdata, t, beta

R = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--h5ad', type=Path, default=Path(
        '/workspace/scratch/PGAM5_myeloid_reference_v2/GSE149614_harmonized_counts_QC.h5ad'))
    parser.add_argument('--saved-results-only', action='store_true', help='Audit portable outputs without the source h5ad')
    args = parser.parse_args()
    d = pd.read_csv(R / 'macrophage_expression.csv.gz')
    s = pd.read_csv(R / 'correlation_statistics.csv').set_index('analysis')
    checks = []

    def check(name, condition):
        assert condition, name
        checks.append(name)

    def close(x, y):
        return np.allclose(x, y, rtol=1e-10, atol=1e-12, equal_nan=True)

    check('Exactly 7448 unique macrophages, all primary T and Tumor',
          len(d) == 7448 and d.cell_id.is_unique and d.harmonized_group.eq('Macrophage').all()
          and d['sample'].str.endswith('T').all() and d.site.eq('Tumor').all())
    check('All 10 samples/patient labels included; one patient label per sample',
          d['sample'].nunique() == 10 and d.patient.nunique() == 10
          and d.groupby('sample').patient.nunique().eq(1).all())
    check('Full-library normalization denominators positive',
          np.isfinite(d.total_counts).all() and d.total_counts.gt(0).all())
    for gene in ['PGAM5', 'MYO19']:
        raw = d[gene + '_raw_count'].to_numpy()
        cp = raw / d.total_counts.to_numpy() * 10000
        check(gene + ' raw integer counts nonnegative and finite',
              np.isfinite(raw).all() and (raw >= 0).all() and (raw == np.floor(raw)).all())
        check(gene + ' CP10k/log1p independently reconstructed',
              close(cp, d[gene + '_CP10k']) and close(np.log1p(cp), d[gene + '_log1p_CP10k']))
        check(gene + ' all raw zeros preserved exactly on both transformed scales',
              np.array_equal(raw == 0, d[gene + '_CP10k'].eq(0))
              and np.array_equal(raw == 0, d[gene + '_log1p_CP10k'].eq(0)))
        check(gene + ' detection flag is raw_count>0',
              np.array_equal(raw > 0, d[gene + '_detected']))
    check('Detection counts 197,311,24; neither 6964; only 173,287',
          d.PGAM5_detected.sum() == 197 and d.MYO19_detected.sum() == 311
          and (d.PGAM5_detected & d.MYO19_detected).sum() == 24
          and (~d.PGAM5_detected & ~d.MYO19_detected).sum() == 6964
          and (d.PGAM5_detected & ~d.MYO19_detected).sum() == 173
          and (~d.PGAM5_detected & d.MYO19_detected).sum() == 287)
    check('All 520 cycling_flag-positive macrophages retained', int(d.cycling_flag.sum()) == 520)
    if not args.saved_results_only:
        with h5py.File(args.h5ad, 'r') as f:
            def column(node):
                if isinstance(node, h5py.Group):
                    if 'values' in node:
                        assert not node['mask'][:].any()
                        return column(node['values'])
                    categories = column(node['categories'])
                    codes = node['codes'][:]
                    assert (codes >= 0).all()
                    return categories[codes]
                return node.asstr()[:] if h5py.check_string_dtype(node.dtype) is not None else node[:]

            check('All source QC cells are GSE149614 primaryT tumors',
                  (column(f['obs/dataset']) == 'GSE149614').all()
                  and all(str(z).endswith('T') for z in column(f['obs/sample']))
                  and (column(f['obs/site']) == 'Tumor').all())
            mask = column(f['obs/harmonized_group']) == 'Macrophage'
            ids = column(f['obs/cell_id'])[mask]
            check('Full 7448 source macrophage membership equals saved cell IDs',
                  mask.sum() == 7448 and len(set(ids)) == 7448 and set(ids) == set(d.cell_id))
            aligned = d.set_index('cell_id').loc[ids]
            for col in ['sample', 'patient', 'donor_id', 'total_counts', 'cycling_flag', 'fine_group']:
                check('Source metadata matches: ' + col,
                      np.array_equal(column(f['obs/' + col])[mask], aligned[col].to_numpy()))
            genes = column(f['var'][f['var'].attrs['_index']])
            indices = [np.flatnonzero(genes == gene) for gene in ['PGAM5', 'MYO19']]
            check('Exactly one symbol-merged column per gene', all(len(z) == 1 for z in indices))
            pointer = f['X/indptr'][:]
            columns = f['X/indices'][:]
            values = f['X/data'][:]
            recovered = []
            for index in indices:
                locations = np.flatnonzero(columns == index[0])
                rows = np.searchsorted(pointer, locations, side='right') - 1
                counts = np.bincount(rows, weights=values[locations], minlength=len(mask))
                recovered.append(counts[mask])
            check('All 14896 raw gene-cell counts independently recovered from source HDF5 CSR',
                  np.array_equal(np.column_stack(recovered),
                                 aligned[['PGAM5_raw_count', 'MYO19_raw_count']].to_numpy()))
            check('Source PGAM5 count annotation matches',
                  np.array_equal(column(f['obs/PGAM5_counts'])[mask], aligned.PGAM5_raw_count))

    for suffix in ['log1p_CP10k', 'raw_counts']:
        column_suffix = 'raw_count' if suffix == 'raw_counts' else suffix
        original_x = d['PGAM5_' + column_suffix].to_numpy()
        original_y = d['MYO19_' + column_suffix].to_numpy()
        for method in ['Spearman', 'Pearson']:
            x = rankdata(original_x) if method == 'Spearman' else original_x.copy()
            y = rankdata(original_y) if method == 'Spearman' else original_y.copy()
            x = x - x.mean()
            y = y - y.mean()
            n = len(x)
            r = np.sum(x * y) / np.sqrt(np.sum(x * x) * np.sum(y * y))
            if method == 'Spearman':
                p = 2 * t.sf(abs(r) * np.sqrt((n - 2) / ((1 + r) * (1 - r))), n - 2)
            else:
                p = 2 * beta(n / 2 - 1, n / 2 - 1, loc=-1, scale=2).cdf(-abs(r))
            row = s.loc[method + '_all_cells_' + suffix]
            check(method + ' ' + suffix + ': n=7448 and manual coefficient/analytic p',
                  int(row.n) == 7448 and close(r, row.correlation)
                  and close(np.log(p), np.log(row.asymptotic_p)))
    secondary = s[s.family.eq('secondary')]
    pvalues = secondary.asymptotic_p.to_numpy()
    order = np.argsort(pvalues)
    adjusted = np.empty(3)
    adjusted[order] = np.minimum(1, np.minimum.accumulate(
        (pvalues[order] * 3 / np.arange(1, 4))[::-1])[::-1])
    check('Exactly one primary and three secondary tests; BH independently checked',
          s.family.eq('primary').sum() == 1 and len(secondary) == 3
          and close(adjusted, secondary.secondary_BH_q))
    primary = s.loc['Spearman_all_cells_log1p_CP10k']
    check('Fixed200000 permutation p uses plus-one formula; valid extreme count',
          int(primary.permutations) == 200000 and 0 <= primary.extreme_permutations <= 200000
          and close(primary.permutation_p, (primary.extreme_permutations + 1) / 200001))
    summary = pd.read_csv(R / 'sample_summary.csv').set_index('sample')
    grouped = d.groupby('sample')
    check('All 10 per-sample counts and means independently reconciled',
          len(summary) == 10 and summary.macrophages.sum() == 7448
          and summary.macrophages.eq(grouped.size()).all()
          and close(summary.PGAM5_mean_CP10k, grouped.PGAM5_CP10k.sum() / grouped.size())
          and close(summary.MYO19_mean_CP10k, grouped.MYO19_CP10k.sum() / grouped.size())
          and summary.neither_detected.sum() == 6964)
    pairs = pd.read_csv(R / 'raw_count_pairs.csv')
    check('Raw count pair frequencies reconstruct all cells, including (0,0)',
          pairs.cells.sum() == 7448 and pairs.loc[pairs.PGAM5_raw_count.eq(0)
          & pairs.MYO19_raw_count.eq(0), 'cells'].iloc[0] == 6964)
    provenance = json.loads((R / 'source_provenance.json').read_text())
    check('Expression file SHA256 matches recorded provenance',
          hashlib.sha256((R / 'macrophage_expression.csv.gz').read_bytes()).hexdigest()
          == provenance['extracted_expression_sha256'])
    result = {'status': 'PASS', 'checks_passed': len(checks), 'checks': checks,
              'source_HDF5_counts_audited': not args.saved_results_only,
              'scope': 'Implementation and data retention checks, not independent biological validation or proof of cell/patient exchangeability'}
    filename = 'audit_saved_results.json' if args.saved_results_only else 'audit.json'
    (R / filename).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
