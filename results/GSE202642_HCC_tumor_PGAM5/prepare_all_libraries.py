from pathlib import Path
import gzip,json
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
from scipy.io import mmread
from threadpoolctl import threadpool_limits
R=Path(__file__).resolve().parent;sc.settings.n_jobs=4
g=pd.read_csv(R/'GSE202642_features.tsv.gz',sep='\t',header=None,names=['ensembl','gene','feature_type'])
b=pd.read_csv(R/'GSE202642_barcodes.tsv.gz',header=None)[0]
assert g.ensembl.is_unique and b.is_unique and g.feature_type.eq('Gene Expression').all()
g.index=g.ensembl;g.index.name='gene_ID'
o=pd.DataFrame({'barcode':b.to_numpy(),'library_suffix':b.str.extract(r'-(\d+)$')[0].to_numpy()},index=b.to_numpy());o.index.name='cell_id'
assert set(o.library_suffix)==set(map(str,range(1,12)))
o['Sample']='library_'+o.library_suffix
# No suffix-to-tissue assumption is embedded here. All libraries are prepared pending verified mapping.
with threadpool_limits(limits=4):
 with gzip.open(R/'GSE202642_matrix.mtx.gz','rb') as f:x=mmread(f).tocsr().T.tocsr().astype(np.int32)
assert x.shape==(len(b),len(g)) and (x.data>0).all()
o['total_counts']=np.asarray(x.sum(1)).ravel();o['n_genes']=np.diff(x.indptr)
mt=g.gene.str.startswith('MT-').to_numpy();assert mt.sum()>=13
o['pct_mt']=np.asarray(x[:,mt].sum(1)).ravel()/o.total_counts.to_numpy()*100
o['passes_QC']=o.n_genes.ge(500)&o.pct_mt.lt(20)
pg=np.flatnonzero(g.gene.eq('PGAM5'));assert len(pg)==1
o['PGAM5_counts']=x[:,pg[0]].toarray().ravel();o['PGAM5_detected']=o.PGAM5_counts.gt(0)
o.to_csv(R/'all_library_metadata_before_QC.csv.gz')
o.groupby('library_suffix',observed=True).agg(cells=('barcode','size'),QC_cells=('passes_QC','sum')).to_csv(R/'library_cell_counts.csv')
a=ad.AnnData(x,obs=o,var=g);a=a[a.obs.passes_QC.to_numpy()].copy()
a.write_h5ad(R/'all_library_counts_QC.h5ad',compression='gzip')
print('All library cells',len(o),'QC',len(a),flush=True)
sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
sc.pp.highly_variable_genes(a,n_top_genes=2000,flavor='seurat');a.var.loc[a.var.gene.eq('PGAM5'),'highly_variable']=False
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
q=pd.DataFrame(rows);q.to_csv(R/'global_cluster_marker_QA.csv',index=False)
a.obs.to_csv(R/'all_library_QC_cluster_metadata.csv.gz')
a.write_h5ad(R/'all_library_normalized_clustered.h5ad',compression='gzip')
print(q[['cluster','cells']+['fraction_'+g for g in ['C1QA','C1QB','C1QC','CD68','CSF1R','LST1','TYROBP','FCN1','CD1C','CD3D','MS4A1','ALB','EPCAM']]].round(2).to_string(index=False),flush=True)
(R/'preparation_summary.json').write_text(json.dumps({'source_cells':len(o),'QC_cells':len(a),'libraries':11,'genes':len(g),'mt_gene_rows':int(mt.sum()),'scope':'All libraries only for preparation; tissue suffix mapping is not assumed','QC':'n_genes>=500,pct_mt<20','clustering':'2000 Seurat HVGs minus PGAM5,30PC,15NN,igraph Leiden1,seed0'},indent=2))
