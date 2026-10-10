from pathlib import Path
import scanpy as sc,pandas as pd,numpy as np,json
from sklearn.metrics import adjusted_rand_score
R=Path(__file__).resolve().parent;meta=pd.read_csv(R/'inputs/all_cell_metadata.csv.gz');sel=meta.diagnosis.eq('Hepatocellular carcinoma')&meta.n_genes.ge(500)&meta.pct_mt.lt(20);expected=pd.Index(meta.loc[sel,'S_ID']+'_'+meta.loc[sel,'Cell']);checks=[];parent=pd.read_csv(R/'allcell_cell_annotations.csv.gz',index_col=0,dtype={'cluster':str})
def check(name,condition):checks.append(dict(check=name,passed=bool(condition)));assert condition,name
for fn,key,coord,summary in [('allcell_marker_subset.h5ad','cluster','allcell_UMAP_coordinates.csv.gz','allcell_summary.json'),('macrophage_marker_subset.h5ad','mac_cluster','macrophage_UMAP_coordinates.csv.gz','macrophage_summary.json')]:
 a=sc.read_h5ad(R/fn);e=expected if key=='cluster' else pd.Index(parent.index[parent.author_type.eq('TAMs')&parent.cluster.isin(['9','11'])]);check(key+'_cell_selection_exact',a.obs_names.equals(e));u=pd.read_csv(R/coord,index_col=0);check(key+'_saved_coordinates_match',np.array_equal(a.obsm['X_umap'],u.to_numpy().astype(a.obsm['X_umap'].dtype)));j=a.var_names.get_loc('PGAM5');raw=a.layers['counts'][:,j].toarray().ravel();check(key+'_PGAM5_raw_metadata_match',np.array_equal(raw,a.obs.PGAM5_counts));values=np.log1p(raw/a.obs.total_counts.to_numpy()*1e4);check(key+'_PGAM5_normalized_full_library',np.allclose(a.X[:,j].toarray().ravel(),values,rtol=1e-5,atol=1e-7));check(key+'_finite_UMAP',np.isfinite(a.obsm['X_umap']).all());recon=a.copy();sc.tl.leiden(recon,resolution=.8 if key=='cluster' else .6,key_added='audit_cluster',random_state=151530,flavor='igraph',n_iterations=2,directed=False);check(key+'_graph_reconstruction_ARI1',adjusted_rand_score(a.obs[key],recon.obs.audit_cluster)==1)
 if key=='mac_cluster':
  s=pd.read_csv(R/'PGAM5_by_macrophage_cluster.csv');cp=raw/a.obs.total_counts.to_numpy()*1e4
  for row in s.itertuples():
   ix=a.obs.mac_cluster.astype(str).eq(row.cluster).to_numpy();check(row.cluster+'_PGAM5_mean_and_detection',np.isclose(cp[ix].mean(),row.PGAM5_mean_CP10k_all_cells) and (raw[ix]>0).sum()==row.PGAM5_detected)
  check('macro_stats_cells_sum',s.cells.sum()==len(a))
full_counts_available=(R/'HCC_raw_QC.h5ad').exists()
if full_counts_available:
 raw=sc.read_h5ad(R/'HCC_raw_QC.h5ad',backed='r');tot=np.asarray(raw.X[:].sum(1)).ravel();check('allcell_original_full_counts_sum_match',np.array_equal(tot,raw.obs.total_counts));raw.file.close()
pd.DataFrame(checks).to_csv(R/'independent_audit_checks.csv',index=False);(R/'audit_summary.json').write_text(json.dumps(dict(checks=len(checks),passed=sum(x['passed'] for x in checks),full_counts_available=full_counts_available,skipped=[] if full_counts_available else ['Full raw HCC counts audit requires a complete full-cell rerun.']),indent=2));print('AUDIT',len(checks),'passed',flush=True)
