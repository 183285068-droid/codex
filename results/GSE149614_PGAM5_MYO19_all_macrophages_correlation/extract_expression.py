"""Extract ALL existing GSE149614 primary-tumor macrophages, retaining zeros."""
from pathlib import Path
import argparse
import hashlib
import json
import anndata as ad
import numpy as np
import pandas as pd

R = Path(__file__).resolve().parent
GENES = ['PGAM5', 'MYO19']


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--h5ad', type=Path, default=Path(
        '/workspace/scratch/PGAM5_myeloid_reference_v2/GSE149614_harmonized_counts_QC.h5ad'))
    parser.add_argument('--annotation', type=Path, default=Path(
        '/workspace/codex/results/PGAM5_RNA_annotation_deconvolution_v4/GSE149614_all_cell_annotation.csv.gz'))
    args = parser.parse_args()
    a = ad.read_h5ad(args.h5ad, backed='r')
    assert a.obs.dataset.astype(str).eq('GSE149614').all()
    assert a.obs['sample'].astype(str).str.endswith('T').all()
    assert a.obs.site.astype(str).eq('Tumor').all()
    mask = a.obs.harmonized_group.astype(str).eq('Macrophage').to_numpy()
    assert int(mask.sum()) == 7448
    assert a.var_names.is_unique
    indices = a.var_names.get_indexer(GENES)
    assert (indices >= 0).all()
    b = a[mask, indices].to_memory()
    raw = b.X.toarray().astype(float)
    assert np.isfinite(raw).all() and (raw >= 0).all() and (raw == np.floor(raw)).all()
    assert np.array_equal(raw[:, 0], b.obs.PGAM5_counts.to_numpy())
    cols = ['cell_id', 'sample', 'patient', 'donor_id', 'site', 'stage', 'virus',
            'total_counts', 'fine_group', 'cycling_flag', 'harmonized_group']
    d = b.obs[cols].reset_index(drop=True).copy()
    assert d.cell_id.is_unique and (d.total_counts > 0).all()
    assert d['sample'].nunique() == 10 and d.patient.nunique() == 10
    for j, gene in enumerate(GENES):
        d[gene + '_raw_count'] = raw[:, j].astype(int)
        d[gene + '_CP10k'] = raw[:, j] / d.total_counts.to_numpy() * 10000
        d[gene + '_log1p_CP10k'] = np.log1p(d[gene + '_CP10k'])
        d[gene + '_detected'] = raw[:, j] > 0
    annotation = pd.read_csv(args.annotation)
    annotation = annotation[annotation.harmonized_group.eq('Macrophage')]
    assert annotation.cell_id.is_unique and set(annotation.cell_id) == set(d.cell_id)
    v = annotation.set_index('cell_id').loc[d.cell_id]
    assert np.array_equal(v.PGAM5_counts.to_numpy(), d.PGAM5_raw_count.to_numpy())
    d['PGAM5_macrophage_annotation'] = v.PGAM5_macrophage_annotation.to_numpy()
    assert len(d) == int(mask.sum())  # No expression-dependent cell filtering.
    d.to_csv(R / 'macrophage_expression.csv.gz', index=False)
    provenance = {
        'dataset': 'GSE149614',
        'scope': 'All existing QC macrophages from the 10 primary HCC tumor samples ending T',
        'annotation_source_commit': 'a7bf04f34197b3b888a59fc74bbbbfb187c99769',
        'all_primaryT_QC_cells': a.n_obs,
        'macrophages': len(d),
        'samples': sorted(d['sample'].unique().tolist()),
        'patients': sorted(d.patient.unique().tolist()),
        'selection': 'harmonized_group == Macrophage only; no filter on either gene; cycling macrophages retained',
        'expression': 'Symbol-merged raw integer X counts; CP10k uses previously saved FULL-library total_counts; natural log1p; no depth matching/regression',
        'source_hashes': {str(p): {'bytes': p.stat().st_size, 'sha256': sha256(p)}
                          for p in [args.h5ad, args.annotation]},
        'extracted_expression_sha256': sha256(R / 'macrophage_expression.csv.gz')
    }
    a.file.close()
    (R / 'source_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps(provenance, indent=2))


if __name__ == '__main__':
    main()
