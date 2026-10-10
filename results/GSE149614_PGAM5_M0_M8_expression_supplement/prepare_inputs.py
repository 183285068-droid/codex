from pathlib import Path
import hashlib
import json

import anndata as ad
import numpy as np
import pandas as pd

R = Path(__file__).resolve().parent
P = R.parent / 'GSE149614_primaryT_hierarchical_clustering'
I = R / 'inputs'
I.mkdir(exist_ok=True)
a = ad.read_h5ad(P / 'macrophage_marker_subset.h5ad')
raw = a.layers['counts'][:, a.var_names.get_loc('PGAM5')].toarray().ravel()
assert np.array_equal(raw, a.obs.PGAM5_counts)
coords = pd.read_csv(P / 'macrophage_UMAP_coordinates.csv.gz', index_col=0)
assert a.obs_names.equals(coords.index)
np.testing.assert_array_equal(a.obsm['X_umap'], coords.to_numpy().astype(a.obsm['X_umap'].dtype))
cells = a.obs[['mac_cluster', 'Sample', 'total_counts']].copy()
cells.index.name = 'cell_id'
cells['PGAM5_raw_counts'] = raw
cells['UMAP1'] = a.obsm['X_umap'][:, 0]
cells['UMAP2'] = a.obsm['X_umap'][:, 1]
cells.to_csv(I / 'cell_expression_coordinates.csv.gz')
pd.read_csv(P / 'macrophage_cluster_annotations.csv').to_csv(I / 'cluster_annotations.csv', index=False)
pd.read_csv(P / 'PGAM5_by_macrophage_cluster.csv').to_csv(I / 'previous_cluster_statistics.csv', index=False)
names = ['macrophage_marker_subset.h5ad', 'macrophage_UMAP_coordinates.csv.gz', 'PGAM5_by_macrophage_cluster.csv', 'macrophage_cluster_annotations.csv']
hashes = {}
for name in names:
    with (P / name).open('rb') as f:
        hashes[name] = hashlib.file_digest(f, 'sha256').hexdigest()
(R / 'source_provenance.json').write_text(json.dumps({'dataset': 'GSE149614', 'scope': 'HCC01T-HCC10T', 'source_directory': P.name, 'source_commit': '8b0722be93ea99ef7878ee2c83616d6f4d2e8f15', 'source_sha256': hashes, 'cells': len(cells), 'clusters': 9}, indent=2))
print('Prepared verified cell expression and original UMAP coordinates:', len(cells))
