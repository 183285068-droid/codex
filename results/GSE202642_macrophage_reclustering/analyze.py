from pathlib import Path
import scanpy as sc,numpy as np,pandas as pd,json,harmonypy
from sklearn.metrics import adjusted_rand_score
from threadpoolctl import threadpool_limits
R=Path(__file__).resolve().parent
src=R/'inputs/macrophages_raw_full.h5ad'
a=sc.read_h5ad(src);a.obs['sample_name']=a.obs.sample_name.astype('category');assert len(a)==11536
if 'counts' not in a.layers:a.layers['counts']=a.X.copy()
assert np.array_equal(np.asarray(a.layers['counts'].sum(1)).ravel(),a.obs.total_counts)
sc.settings.n_jobs=4
sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a);sc.pp.highly_variable_genes(a,n_top_genes=2500,flavor='seurat',batch_key='sample_name');a.var.loc['PGAM5','highly_variable']=False
v=a[:,a.var.highly_variable].copy();v.layers.clear();sc.pp.scale(v,max_value=10);sc.tl.pca(v,n_comps=40,random_state=202642);a.obsm['X_pca']=v.obsm['X_pca'];del v
sc.pp.neighbors(a,n_neighbors=20,n_pcs=30,random_state=202642,key_added='unintegrated');sc.tl.umap(a,neighbors_key='unintegrated',random_state=202642);a.obsm['X_umap_unintegrated']=a.obsm['X_umap'].copy()
with threadpool_limits(limits=4):h=harmonypy.run_harmony(a.obsm['X_pca'][:,:30],a.obs,'sample_name',random_state=202642,max_iter_harmony=20)
hz=np.asarray(h.Z_corr);a.obsm['X_pca_harmony']=hz.T if hz.shape[0]==30 else hz
sc.pp.neighbors(a,n_neighbors=20,use_rep='X_pca_harmony',random_state=202642);sc.tl.umap(a,random_state=202642)
for res in [.4,.6,.8]:sc.tl.leiden(a,resolution=res,key_added='leiden_'+str(res),random_state=202642,flavor='igraph',n_iterations=2,directed=False)
a.obs['mac_cluster']=a.obs['leiden_0.6'].map(lambda x:'M'+str(x)).astype('category');sc.tl.leiden(a,resolution=.6,key_added='seed2',random_state=42,flavor='igraph',n_iterations=2,directed=False)
print('Clusters',a.obs.mac_cluster.value_counts().to_dict(),flush=True)
sc.tl.rank_genes_groups(a,groupby='mac_cluster',method='wilcoxon',tie_correct=True,n_genes=50,use_raw=False);sc.get.rank_genes_groups_df(a,group=None).to_csv(R/'top50_markers.csv.gz',index=False)
a.write_h5ad(R/'macrophages_full.h5ad',compression='gzip');a.obs.to_csv(R/'cell_annotations.csv.gz');pd.DataFrame(a.obsm['X_umap'],index=a.obs_names,columns=['UMAP1','UMAP2']).to_csv(R/'UMAP_coordinates.csv.gz')
rows=[]
genes=['C1QA','C1QB','CSF1R','CD68','FOLR2','CD163','TREM2','SPP1','APOE','MARCO','CD5L','HLA-DRA','CD74','IL1B','CXCL8','FCN1','VCAN','S100A8','S100A9','MKI67','TOP2A','ISG15','IFIT1','MT1G','MT2A','LYVE1','MRC1','NKG7','CD3D','CD1C','FCER1A','PGAM5']
for c in a.obs.mac_cluster.cat.categories:
 ix=a.obs.mac_cluster.eq(c).to_numpy()
 for g in genes:
  j=a.var_names.get_loc(g);x=a.X[ix,j].toarray().ravel();raw=a.layers['counts'][ix,j].toarray().ravel();rows.append(dict(cluster=c,gene=g,mean_log1p_CP10k=x.mean(),detection_percent=(raw>0).mean()*100))
pd.DataFrame(rows).to_csv(R/'marker_evidence.csv',index=False)
(R/'summary.json').write_text(json.dumps(dict(cells=len(a),genes=a.n_vars,clusters=int(a.obs.mac_cluster.nunique()),seed_ARI=adjusted_rand_score(a.obs.mac_cluster,a.obs.seed2),resolution=.6,resolution_sensitivity={str(r):int(a.obs['leiden_'+str(r)].nunique()) for r in [.4,.6,.8]},parents=['C0','C2','C16'],sample_count=int(a.obs.sample_name.nunique())),indent=2));print('DONE',flush=True)
