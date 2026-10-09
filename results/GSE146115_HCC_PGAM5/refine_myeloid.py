from pathlib import Path
import pandas as pd,numpy as np,scanpy as sc,anndata as ad
R=Path(__file__).resolve().parent
a=ad.read_h5ad(R/'all_tumor_counts_QC.h5ad');o=pd.read_csv(R/'all_tumor_QC_cluster_metadata.csv.gz',index_col=0,dtype={'leiden':str})
assert list(a.obs_names)==list(o.index)
a.obs['global_cluster']=o.leiden.to_numpy()
a=a[a.obs.global_cluster.eq('4').to_numpy()].copy()
sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
sc.pp.highly_variable_genes(a,n_top_genes=1500,flavor='seurat');a.var.loc[a.var.gene.eq('PGAM5')|a.var.ambiguous_symbol,'highly_variable']=False
sc.pp.pca(a,n_comps=15,random_state=0);sc.pp.neighbors(a,n_neighbors=10,n_pcs=15,random_state=0)
sc.tl.leiden(a,resolution=.6,random_state=0,flavor='igraph',n_iterations=2,directed=False,key_added='myeloid_cluster')
rows=[]
for cluster in a.obs.myeloid_cluster.cat.categories:
 sel=a.obs.myeloid_cluster.eq(cluster).to_numpy();row={'cluster':str(cluster),'cells':int(sel.sum())}
 for gene in ['C1QA','C1QB','C1QC','CD68','CSF1R','LST1','TYROBP','FCN1','S100A8','S100A9','CD1C','FCER1A','CD3D','MS4A1','ALB','SPP1']:
  ix=np.flatnonzero(a.var.gene.eq(gene));v=a.X[sel][:,ix].toarray().sum(1)
  row['fraction_'+gene]=float((v>0).mean());row['mean_log1p_'+gene]=float(v.mean())
 rows.append(row)
q=pd.DataFrame(rows);q.to_csv(R/'myeloid_cluster_marker_QA.csv',index=False)
a.obs.to_csv(R/'myeloid_cluster_metadata.csv.gz')
print(q[['cluster','cells']+['fraction_'+g for g in ['C1QA','C1QB','C1QC','CD68','CSF1R','TYROBP','FCN1','CD1C','CD3D','ALB']]].round(3).to_string(index=False))
