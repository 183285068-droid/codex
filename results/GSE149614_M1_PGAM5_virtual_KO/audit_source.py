"""Check the target against the previous all-macrophage object and saved hashes."""
from pathlib import Path
import hashlib
import json
import anndata as ad
import numpy as np
import pandas as pd

R=Path(__file__).resolve().parent
summary=json.loads((R/'input_summary.json').read_text())
source=R.parent/'GSE149614_primaryT_hierarchical_clustering/macrophages_annotated_full.h5ad'
with source.open('rb') as f:
    assert hashlib.file_digest(f,'sha256').hexdigest()==summary['source_sha256']
b=ad.read_h5ad(source)
original=b[b.obs.mac_cluster.astype(str).eq('M1')]
a=ad.read_h5ad(R/'inputs/M1_raw.h5ad')
assert original.obs_names.equals(a.obs_names) and original.var_names.equals(a.var_names)
assert (original.layers['counts']!=a.X).nnz==0
context=pd.read_csv(R/'inputs/macrophage_context.csv.gz',index_col=0)
assert context.index.equals(b.obs_names)
np.testing.assert_array_equal(context.mac_cluster.to_numpy(),b.obs.mac_cluster.astype(str))
np.testing.assert_array_equal(context[['UMAP1','UMAP2']].to_numpy().astype(b.obsm['X_umap'].dtype),b.obsm['X_umap'])
x=pd.read_csv(R/'inputs/network_raw_counts.csv.gz',index_col=0)
selection=pd.read_csv(R/'gene_selection.csv').set_index('gene')
assert set(selection.index[selection.network_selected])==set(x.index)
assert selection.loc[x.index.difference(['PGAM5']),'eligible'].all()
assert a.obs.Sample.nunique()==10 and a.obs.PGAM5_counts.gt(0).sum()==19
for name,expected in summary['input_sha256'].items():
    assert hashlib.sha256((R/name).read_bytes()).hexdigest()==expected
    assert hashlib.sha256((R/'unsparsified'/name).read_bytes()).hexdigest()==expected
vendor=json.loads((R/'vendor/provenance.json').read_text())
for name,expected in vendor['exact_runtime_python_source_sha256'].items():
    assert hashlib.sha256((R/'vendor'/name).read_bytes()).hexdigest()==expected
    assert hashlib.sha256((R/'unsparsified/vendor'/name).read_bytes()).hexdigest()==expected
report={'original_source_hash_matches':True,'all_target_cell_IDs_and_raw_counts_match_source':True,
        'original_UMAP_coordinates_unchanged':True,'gene_selection_table_matches_input':True,
        'both_models_use_identical_inputs':True,'runtime_source_snapshots_match_hashes':True,
        'target_cells':len(a),'PGAM5_detected':19,'PGAM5_zeros':189,'samples':10}
(R/'source_audit_summary.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
