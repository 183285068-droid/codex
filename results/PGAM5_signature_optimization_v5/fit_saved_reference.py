"""Portable inference on a full gene-symbol x sample matrix, with honest unit flags."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from reference_core import fit, T

R = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--bulk', type=Path, required=True, help='TSV, exact gene symbols in first column, samples in remaining columns')
    p.add_argument('--unit', choices=['TPM','CPM','raw_counts','CP10k'], required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    bulk = pd.read_csv(args.bulk, sep='\t', index_col=0)
    assert bulk.columns.is_unique
    bulk = bulk.groupby(level=0, sort=True).sum()
    assert np.isfinite(bulk.to_numpy()).all() and (bulk.to_numpy() >= 0).all()
    if args.unit == 'raw_counts':
        assert bulk.sum().gt(0).all()
        bulk = bulk.div(bulk.sum(axis=0), axis=1)*10000
    elif args.unit in ['TPM','CPM']:
        bulk = bulk/100  # Already normalized to the complete library; no subset re-normalization.
    ref = pd.read_csv(R/'EXPLORATORY_reference_library_RNA_contribution.tsv', sep='\t', index_col=0)
    scales = pd.read_csv(R/'library_RNA_contribution_row_scales.csv').set_index('gene').scale
    lock = json.loads((R/'locked_reference.json').read_text())
    genes = ref.index[ref.index.isin(bulk.index)]
    assert 'PGAM5' in genes
    assert len(genes) > 0
    scale = scales.loc[genes].to_numpy()
    model = {'matrix': ref.loc[genes].to_numpy()/scale[:, None], 'features': np.arange(len(genes)), 'scale': scale}
    result = []
    for sample in bulk.columns:
        c, residual, _ = fit(model, bulk.loc[genes, sample].to_numpy(), lock['solver'])
        result.append({'sample': sample, **dict(zip(ref.columns, c)), 'relative_residual': residual,
            'reference_feature_coverage': len(genes)/len(ref),
            'reference_validation_passed': False, 'usable_as_cellular_infiltration': False})
    pd.DataFrame(result).to_csv(args.output, index=False)
    print('Exploratory coefficients saved. This reference failed validation; outputs are not cellular infiltration.')


if __name__ == '__main__':
    main()
