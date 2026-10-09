from pathlib import Path
import json,re
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
from scipy.sparse import csr_matrix
R=Path(__file__).resolve().parent;sc.settings.n_jobs=4
c=pd.read_csv(R/'GSE146115_HCC1-2-5-9_count_with_ERCC.txt.gz',sep='\t',index_col=0)
assert c.columns.is_unique
spike=c.index.str.match(r'^ERCC-\d+$')
pd.DataFrame({'gene':c.index[spike],'reason':'external_ERCC_spikein'}).to_csv(R/'excluded_spikeins.csv',index=False)
c=c.loc[~spike]
assert np.isfinite(c.to_numpy()).all() and (c.to_numpy()>=0).all() and (c.to_numpy()%1==0).all()
# Preserve source rows with unique IDs; do not guess identities for date-corrupted symbols.
g=pd.DataFrame({'gene':c.index.to_numpy()},index=[f'row_{i+1:05d}' for i in range(len(c))]);g.index.name='gene_ID'
g['ambiguous_symbol']=g.gene.duplicated(False)|g.gene.str.match(r'^\d{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$')
g[g.ambiguous_symbol].to_csv(R/'ambiguous_source_gene_symbols.csv')
x=csr_matrix(c.to_numpy(dtype=np.int32).T);o=pd.DataFrame(index=c.columns);o.index.name='cell_id'
o['Sample']=o.index.str.extract(r'^(HCC\d+)',expand=False);assert set(o.Sample)=={'HCC1','HCC2','HCC5','HCC9'}
o['total_counts']=np.asarray(x.sum(1)).ravel();o['n_genes']=np.diff(x.indptr)
mt=g.gene.str.startswith('MT-').to_numpy();assert not mt.any(),'Review mitochondrial naming/QC availability'
mt_proxy=g.gene.isin(['RNR1','RNR2']).to_numpy();assert mt_proxy.sum()==2
o['pct_mt_proxy']=np.asarray(x[:,mt_proxy].sum(1)).ravel()/o.total_counts.to_numpy()*100
o['pct_mt']=o.pct_mt_proxy;o['mitochondrial_QC_available']=False;o['passes_QC']=o.n_genes.ge(500)&o.pct_mt_proxy.lt(20)
o['PGAM5_counts']=x[:,np.flatnonzero(g.gene.eq('PGAM5'))[0]].toarray().ravel();o['PGAM5_detected']=o.PGAM5_counts.gt(0)
# Check zero-padded processed-header mappings against actual count columns.
m=pd.read_csv(R/'GSE146115_raw_to_prcessed_data_column.txt.gz',sep='\t',skiprows=1)
valid=m['processed data column headers'].fillna('').str.match(r'^HCC\d+_COL\d+_ROW\d+$')
m=m[valid].copy()
m['cell_id']=m['processed data column headers'].map(lambda z:re.sub(r'_COL0*(\d+)_ROW0*(\d+)',lambda k:f'_COL{int(k[1])}_ROW{int(k[2])}',z))
assert set(m.cell_id)==set(o.index)
m.to_csv(R/'source_cell_header_mapping.csv',index=False)
o.to_csv(R/'all_cell_metadata_before_QC.csv.gz')
a=ad.AnnData(x,obs=o,var=g);a=a[a.obs.passes_QC.to_numpy()].copy()
a.write_h5ad(R/'all_tumor_counts_QC.h5ad',compression='gzip')
sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
sc.pp.highly_variable_genes(a,n_top_genes=2000,flavor='seurat');a.var.loc[a.var.gene.eq('PGAM5')|a.var.ambiguous_symbol,'highly_variable']=False
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
a.obs.to_csv(R/'all_tumor_QC_cluster_metadata.csv.gz');a.write_h5ad(R/'all_tumor_normalized_clustered.h5ad',compression='gzip')
s={'source_cells':len(o),'QC_cells':len(a),'human_gene_rows':len(g),'excluded_ERCC_rows':int(spike.sum()),'ambiguous_gene_rows':int(g.ambiguous_symbol.sum()),'patients':4,'mitochondrial_QC_available':False,'QC':'n_genes>=500, RNR1+RNR2 mitochondrial rRNA fraction <20%. Full mitochondrial fraction is unavailable; proxy is a lower bound, not assumed zero','count_type':'Fluidigm C1 HTSeq read counts (not UMI)','annotation':'No author cell type labels supplied; marker-inferred annotation','HVG':'2000 Seurat HVGs minus PGAM5 and ambiguous symbols'}
(R/'preparation_summary.json').write_text(json.dumps(s,indent=2))
print(json.dumps(s,indent=2),flush=True)
print(q[['cluster','cells']+['fraction_'+g for g in ['C1QA','C1QB','C1QC','CD68','CSF1R','LST1','TYROBP','FCN1','CD1C','CD3D','MS4A1','ALB','EPCAM']]].round(2).to_string(index=False),flush=True)
