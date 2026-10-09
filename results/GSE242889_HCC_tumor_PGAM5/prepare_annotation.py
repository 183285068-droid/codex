from pathlib import Path
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,json
from scipy.io import mmread
from threadpoolctl import threadpool_limits
R=Path(__file__).resolve().parent;sc.settings.n_jobs=4
parts=[];genes=[]
for d in sorted((R/'counts').iterdir()):
 g=pd.read_csv(d/'genes.tsv',sep='\t',header=None,names=['ensembl','gene']);b=pd.read_csv(d/'barcodes.tsv',sep='\t',header=None)[0]
 assert g.ensembl.is_unique and b.is_unique
 with threadpool_limits(limits=4):x=mmread(d/'matrix.mtx').tocsr().T.tocsr().astype(np.int32)
 assert x.shape==(len(b),len(g)) and (x.data>0).all()
 sample=d.name.split('_')[0];o=pd.DataFrame({'barcode':b.values,'Sample':sample,'patient':sample[:-1]},index=sample+'_'+b.astype(str))
 o.index.name='cell_id'
 o['total_counts']=np.asarray(x.sum(1)).ravel();o['n_genes']=np.diff(x.indptr)
 o['pct_mt']=np.asarray(x[:,g.gene.str.startswith('MT-').to_numpy()].sum(1)).ravel()/o.total_counts.to_numpy()*100
 o['passes_QC']=o.n_genes.ge(500)&o.pct_mt.lt(20)
 g.index=g.ensembl;g.index.name='gene_ID';genes.append(g)
 parts.append(ad.AnnData(x,obs=o,var=g))
 print(sample,'cells',len(o),'QC',o.passes_QC.sum(),flush=True)
all_o=pd.concat([a.obs for a in parts]);all_o.to_csv(R/'all_tumor_cell_metadata_before_QC.csv.gz')
a=ad.concat(parts,join='outer',merge='first')
g=pd.concat(genes).drop_duplicates('ensembl').set_index('ensembl',drop=False);a.var=g.reindex(a.var_names)
assert a.var.gene.notna().all() and a.obs_names.is_unique
a=a[a.obs.passes_QC.to_numpy()].copy();pg=np.flatnonzero(a.var.gene.eq('PGAM5'));assert len(pg)==1
a.obs['PGAM5_counts']=a.X[:,pg[0]].toarray().ravel();a.obs['PGAM5_detected']=a.obs.PGAM5_counts.gt(0)
a.write_h5ad(R/'all_tumor_counts_QC.h5ad',compression='gzip')
sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
sc.pp.highly_variable_genes(a,n_top_genes=2000,flavor='seurat')
a.var.loc[a.var.gene.eq('PGAM5'),'highly_variable']=False
sc.pp.pca(a,n_comps=30,random_state=0);sc.pp.neighbors(a,n_neighbors=15,n_pcs=30,random_state=0)
sc.tl.leiden(a,resolution=1,random_state=0,flavor='igraph',n_iterations=2,directed=False)
markers=['C1QA','C1QB','C1QC','CD68','CSF1R','LST1','TYROBP','FCER1G','AIF1','CD163','MSR1','SPP1','FCN1','S100A8','S100A9','CD1C','FCER1A','CLEC9A','LAMP3','GZMB','JCHAIN','TPSAB1','CD3D','TRAC','NKG7','GNLY','MS4A1','CD79A','ALB','APOA1','TTR','EPCAM','KRT19','PECAM1','VWF','COL1A1','DCN']
rows=[]
for cluster in a.obs.leiden.cat.categories:
 sel=a.obs.leiden.eq(cluster).to_numpy();row={'cluster':str(cluster),'cells':int(sel.sum())}
 for gene in markers:
  ix=np.flatnonzero(a.var.gene.eq(gene))
  if not len(ix):continue
  v=a.X[sel][:,ix].toarray().sum(1);row['fraction_'+gene]=float((v>0).mean());row['mean_log1p_'+gene]=float(v.mean())
 rows.append(row)
qa=pd.DataFrame(rows);qa.to_csv(R/'global_cluster_marker_QA.csv',index=False)
a.obs.to_csv(R/'all_tumor_QC_cluster_metadata.csv.gz')
a.write_h5ad(R/'all_tumor_normalized_clustered.h5ad',compression='gzip')
cols=['cluster','cells']+['fraction_'+g for g in ['C1QA','C1QB','CD68','CSF1R','LST1','TYROBP','FCN1','CD1C','CD3D','GNLY','MS4A1','ALB','EPCAM','PECAM1','COL1A1']]
print(qa[cols].round(2).to_string(index=False),flush=True)
