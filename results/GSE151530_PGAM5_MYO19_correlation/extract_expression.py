"""Extract two-gene raw counts for the previously annotated HCC macrophages."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
import anndata as ad

R = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('--h5ad', type=Path, default=Path('/workspace/scratch/PGAM5_myeloid_reference_v2/GSE151530_harmonized_counts_QC.h5ad'))
parser.add_argument('--annotation', type=Path, default=Path('/workspace/codex/results/PGAM5_RNA_annotation_deconvolution_v4/GSE151530_all_cell_annotation.csv.gz'))
args = parser.parse_args()
a = ad.read_h5ad(args.h5ad, backed='r')
assert a.obs.dataset.astype(str).eq('GSE151530').all()
mask = a.obs.harmonized_group.astype(str).eq('Macrophage').to_numpy()
assert int(mask.sum()) == 3493
gene_indices = a.var_names.get_indexer(['PGAM5','MYO19'])
assert (gene_indices >= 0).all() and a.var_names.is_unique
b = a[mask,gene_indices].to_memory()
raw = b.X.toarray().astype(np.float64)
assert np.isfinite(raw).all() and (raw >= 0).all() and (raw == np.floor(raw)).all()
assert np.array_equal(raw[:,0], b.obs.PGAM5_counts.to_numpy())
assert b.obs.diagnosis.astype(str).eq('Hepatocellular carcinoma').all()
cols = ['cell_id','S_ID','Sample','donor_id','diagnosis','matrix_column','total_counts','fine_group','cycling_flag']
d = b.obs[cols].reset_index(drop=True).copy()
assert d.cell_id.is_unique and d.matrix_column.is_unique and (d.total_counts > 0).all()
for j, gene in enumerate(['PGAM5','MYO19']):
    d[gene+'_raw_count'] = raw[:,j].astype(int)
    d[gene+'_CP10k'] = raw[:,j]/d.total_counts.to_numpy()*10000
    d[gene+'_log1p_CP10k'] = np.log1p(d[gene+'_CP10k'])
    d[gene+'_detected'] = raw[:,j] > 0
assert np.allclose(d.PGAM5_log1p_CP10k, b.obs.PGAM5_log1p_10k, atol=1e-6)
annotation = pd.read_csv(args.annotation)
annotation = annotation[annotation.harmonized_group.eq('Macrophage')]
assert annotation.matrix_column.is_unique and set(annotation.matrix_column) == set(d.matrix_column)
v = annotation.set_index('matrix_column').loc[d.matrix_column]
assert np.array_equal(v.cell_id.to_numpy(),d.cell_id.to_numpy())
assert np.array_equal(v.PGAM5_counts.to_numpy(),d.PGAM5_raw_count.to_numpy())
d['PGAM5_macrophage_annotation'] = v.PGAM5_macrophage_annotation.to_numpy()
d.to_csv(R/'macrophage_expression.csv.gz',index=False)
a.file.close()
provenance = {'dataset':'GSE151530','scope':'Existing HCC QC macrophages, harmonized_group==Macrophage; no ICC or other lineage','annotation_source_commit':'a7bf04f34197b3b888a59fc74bbbbfb187c99769','all_HCC_QC_cells':44776,'macrophages':len(d),'sample_libraries':d.Sample.nunique(),'donor_proxy_labels':d.donor_id.nunique(),'definition':'All macrophages retained, including cycling and PGAM5/MYO19 undetected; no selection by gene positivity for primary analysis','expression':'X is harmonized symbol-merged raw integer counts; CP10k uses stored full-library total_counts; natural log1p. No depth matching/regression','source_hashes':{str(p):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in [args.h5ad,args.annotation]},'extracted_expression_sha256':hashlib.sha256((R/'macrophage_expression.csv.gz').read_bytes()).hexdigest()}
(R/'source_provenance.json').write_text(json.dumps(provenance,indent=2))
print(json.dumps(provenance,indent=2))
