"""Check all published plotting summaries against original count matrices."""
from pathlib import Path
import hashlib
import json
import anndata as ad
import pandas as pd
import numpy as np
from scipy import sparse

R = Path(__file__).resolve().parent
reports = []
provenance = json.loads((R/'source_provenance.json').read_text())
for name, expected in provenance['source_sha256'].items():
    with (R.parent/name).open('rb') as f:
        assert hashlib.file_digest(f, 'sha256').hexdigest() == expected, name


def check(key, raw, genes, obs, cluster_key):
    values = pd.read_csv(R/'inputs'/(key+'_values.csv'), dtype={'cluster': str})
    selected = values.gene.drop_duplicates().tolist()
    idx = genes.get_indexer(selected)
    assert idx.min() >= 0
    raw = raw[:, idx].astype(np.float64)
    log = raw.multiply((10000/obs.total_counts.to_numpy())[:, None]).tocsr()
    log.data = np.log1p(log.data)
    max_error = 0
    for cluster, rows in values.groupby('cluster', sort=False):
        mask = obs[cluster_key].astype(str).eq(cluster).to_numpy()
        assert mask.sum() == rows.cells.iloc[0]
        expected = rows.set_index('gene').loc[selected]
        frac = np.asarray((raw[mask] > 0).mean(axis=0)).ravel()
        mean = np.asarray(log[mask].mean(axis=0)).ravel()
        np.testing.assert_allclose(frac, expected.detection_fraction, atol=1e-12, rtol=1e-12)
        # Original all-cell evidence used float32 normalized expression.
        np.testing.assert_allclose(mean, expected.mean_log1p_CP10k, atol=2e-5, rtol=2e-5)
        max_error = max(max_error, float(np.max(np.abs(mean-expected.mean_log1p_CP10k))))
    reports.append({'dataset_version': key, 'cells': len(obs), 'points': len(values),
        'raw_detection_fraction_agreement': True, 'full_library_log_normalization_agreement': True,
        'maximum_abs_mean_difference': max_error})


p = R.parent/'GSE202642_HCC_allcell_annotation'
a = ad.read_h5ad(p/'HCC_annotation_marker_subset.h5ad')
check('allcell_C0_C20', a.layers['counts'], a.var_names, a.obs, 'cluster')
p = R.parent/'GSE202642_macrophage_reclustering'
a = ad.read_h5ad(p/'macrophages_marker_subset.h5ad')
check('macrophage_M0_M8', a.layers['counts'], a.var_names, a.obs, 'mac_cluster')
p = R.parent/'GSE202642_PGAM5_subpopulation_assessment'
raw = sparse.load_npz(p/'portable_raw_feature_counts.npz')
genes = pd.Index(pd.read_csv(p/'portable_feature_genes.csv').gene)
obs = pd.read_csv(p/'cell_annotations.csv.gz', index_col=0, dtype={'primary_cluster': str})
check('macrophage_C0_C14_legacy', raw, genes, obs, 'primary_cluster')
(R/'independent_value_audit.json').write_text(json.dumps({'source_hashes_unchanged': True,
    'all_points_verified': reports}, indent=2)+'\n')
print(json.dumps(reports, indent=2))
