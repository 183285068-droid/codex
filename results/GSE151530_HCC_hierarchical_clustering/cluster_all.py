from pathlib import Path
import json,gzip,hashlib,shutil
import numpy as np,pandas as pd,scanpy as sc,anndata as ad,harmonypy
from scipy import sparse
from scipy.io import mmread
from threadpoolctl import threadpool_limits
from sklearn.metrics import adjusted_rand_score
R=Path(__file__).resolve().parent;I=R/'inputs';I.mkdir(exist_ok=True)
matrix=R/'GSE151530_matrix.mtx.gz'
if not matrix.exists():raise FileNotFoundError('Run python download_matrix.py or copy the official matrix to this directory.')
meta=pd.read_csv(I/'all_cell_metadata.csv.gz');genes=pd.read_csv(I/'GSE151530_genes.tsv.gz',sep='\t',header=None,names=['ensembl','gene']);barcodes=pd.read_csv(I/'GSE151530_barcodes.tsv.gz',header=None)[0];assert meta.Cell.equals(barcodes);assert hashlib.sha256(matrix.read_bytes()).hexdigest()=='50dad60e779b3332350e2552ed4260369ccba8db4f73918bbe2a9f6549db6403'
print('Loading raw counts',flush=True)
with threadpool_limits(limits=4):
 with gzip.open(matrix,'rb') as f:x=mmread(f).T.tocsr().astype(np.int32)
assert x.shape==(56721,18667);assert np.array_equal(np.asarray(x.sum(1)).ravel(),meta.total_counts);unique=pd.Index(pd.unique(genes.gene));mer=sparse.csr_matrix((np.ones(len(genes),dtype=np.int32),(np.arange(len(genes)),unique.get_indexer(genes.gene))),shape=(len(genes),len(unique)));x=(x@mer).tocsr();sel=meta.diagnosis.eq('Hepatocellular carcinoma')&meta.n_genes.ge(500)&meta.pct_mt.lt(20);o=meta.loc[sel].copy();o.index=o.S_ID+'_'+o.Cell;o.index.name='cell_id';assert o.index.is_unique
raw=x[sel.to_numpy()];a=ad.AnnData(raw,obs=o,var=pd.DataFrame(index=unique));a.obs['Sample']=a.obs.Sample.astype('category');a.obs['author_type']=a.obs.Type;a.layers['counts']=raw.copy();assert len(a)==44776;np.testing.assert_array_equal(raw[:,a.var_names.get_loc('PGAM5')].toarray().ravel(),a.obs.PGAM5_counts);a.write_h5ad(R/'HCC_raw_QC.h5ad',compression='gzip');del x,raw
sc.settings.n_jobs=4;sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a);sc.pp.highly_variable_genes(a,n_top_genes=3000,flavor='seurat',batch_key='Sample');a.var.loc['PGAM5','highly_variable']=False;a.var.to_csv(R/'HVG_selection.csv.gz');v=a[:,a.var.highly_variable].copy();v.layers.clear();sc.pp.scale(v,max_value=10);sc.tl.pca(v,n_comps=50,random_state=151530);a.obsm['X_pca']=v.obsm['X_pca'].copy();del v;print('PCA done',flush=True)
sc.pp.neighbors(a,n_neighbors=20,n_pcs=40,key_added='unintegrated',random_state=151530);sc.tl.umap(a,neighbors_key='unintegrated',random_state=151530);a.obsm['X_umap_unintegrated']=a.obsm['X_umap'].copy()
with threadpool_limits(limits=4):h=harmonypy.run_harmony(a.obsm['X_pca'][:,:40],a.obs,'Sample',random_state=151530,max_iter_harmony=20)
z=np.asarray(h.Z_corr);a.obsm['X_pca_harmony']=z.T if z.shape[0]==40 else z;sc.pp.neighbors(a,n_neighbors=20,use_rep='X_pca_harmony',random_state=151530);sc.tl.umap(a,random_state=151530)
for key,seed in [('cluster',151530),('seed42_cluster',42)]:sc.tl.leiden(a,resolution=.8,key_added=key,random_state=seed,flavor='igraph',n_iterations=2,directed=False)
print('Clusters',a.obs.cluster.value_counts().to_dict(),flush=True);sc.tl.rank_genes_groups(a,groupby='cluster',method='wilcoxon',tie_correct=True,n_genes=40,use_raw=False);sc.get.rank_genes_groups_df(a,group=None).to_csv(R/'allcell_top40_markers.csv.gz',index=False)
a.obs.to_csv(R/'allcell_metadata.csv.gz');pd.DataFrame(a.obsm['X_umap'],index=a.obs_names,columns=['UMAP1','UMAP2']).to_csv(R/'allcell_UMAP_coordinates.csv.gz');a.obs.groupby(['cluster','author_type'],observed=True).size().unstack(fill_value=0).to_csv(R/'allcell_cluster_author_counts.csv');a.obs.groupby(['cluster','Sample'],observed=True).size().unstack(fill_value=0).to_csv(R/'allcell_cluster_sample_counts.csv')
panels=json.loads((R/'annotation_panels.json').read_text());extra=['PGAM5','FOLR2','CD163','TREM2','SPP1','CD5L','MARCO','APOE','IL7R','GZMK','FGFBP2','FOXP3','CTLA4','IL2RA','LYZ'];gs=[g for g in dict.fromkeys(sum(panels.values(),[])+extra) if g in a.var_names];rows=[]
for c in a.obs.cluster.cat.categories:
 ix=a.obs.cluster.eq(c).to_numpy()
 for g in gs:
  j=a.var_names.get_loc(g);vals=a.X[ix,j].toarray().ravel();counts=a.layers['counts'][ix,j].toarray().ravel();rows.append(dict(cluster=c,gene=g,mean_log1p_CP10k=float(vals.mean()),detection_percent=float((counts>0).mean()*100)))
pd.DataFrame(rows).to_csv(R/'allcell_marker_evidence.csv',index=False);(R/'annotation_panels.json').write_text(json.dumps(panels,indent=2));a.write_h5ad(R/'HCC_allcell_clustered.h5ad',compression='gzip');summary=dict(HCC_before_QC=50023,HCC_after_QC=len(a),samples=int(a.obs.Sample.nunique()),patient_proxies=int(a.obs.patient_proxy.nunique()),allcell_clusters=int(a.obs.cluster.nunique()),seed_ARI=adjusted_rand_score(a.obs.cluster,a.obs.seed42_cluster),genes=a.n_vars,QC='n_genes>=500,pct_mt<20',HVG='3000 batch-aware Seurat, PGAM5 not used in geometry',matrix_SHA256='50dad60e779b3332350e2552ed4260369ccba8db4f73918bbe2a9f6549db6403');(R/'allcell_summary.json').write_text(json.dumps(summary,indent=2));print('DONE',summary,flush=True)
