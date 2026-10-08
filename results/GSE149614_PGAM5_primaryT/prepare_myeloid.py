from pathlib import Path
import gzip,json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
R=Path('/workspace/scratch/GSE149614');meta=pd.read_csv(R/'GSE149614_HCC.metadata.updated.txt.gz',sep='\t')
assert meta.Cell.is_unique
assert meta.loc[meta['sample'].str.endswith('T'),'site'].eq('Tumor').all()
obs=meta[meta.celltype.eq('Myeloid')].copy().set_index('Cell',drop=False)
rows=[];genes=[]
print('Streaming counts for',len(obs),'myeloid cells',flush=True)
with gzip.open(R/'GSE149614_HCC.scRNAseq.S71915.count.txt.gz','rb') as f:
 header=f.readline().strip().decode().split('\t')
 assert len(header)==len(meta) and len(set(header))==len(header) and set(header)==set(meta.Cell)
 indices=pd.Index(header).get_indexer(obs.index);assert (indices>=0).all()
 for k,line in enumerate(f):
  gene,sep,values=line.rstrip(b'\r\n').partition(b'\t');assert sep
  vec=np.fromstring(values,dtype=np.int32,sep='\t');assert len(vec)==len(header) and (vec>=0).all()
  genes.append(gene.decode());rows.append(sparse.csr_matrix(vec[indices].reshape(1,-1)))
  if (k+1)%3000==0:print('genes parsed',k+1,flush=True)
x=sparse.vstack(rows,format='csr').T.tocsr();del rows
var=pd.DataFrame({'gene_symbol':genes},index=genes)
a=ad.AnnData(x,obs=obs,var=var);a.var_names_make_unique()
a.obs['total_counts']=np.asarray(x.sum(1)).ravel();a.obs['n_genes']=np.diff(x.indptr)
mt=a.var.gene_symbol.str.startswith('MT-').to_numpy();a.obs['pct_mt']=np.asarray(x[:,mt].sum(1)).ravel()/a.obs.total_counts.to_numpy()*100
pg=np.flatnonzero(a.var.gene_symbol.eq('PGAM5'));assert len(pg)==1
a.obs['PGAM5_counts']=x[:,pg[0]].toarray().ravel();a.obs['PGAM5_detected']=a.obs.PGAM5_counts.gt(0)
a.write_h5ad(R/'all_sites_myeloid_counts.h5ad',compression='gzip')
a.obs.to_csv(R/'all_sites_myeloid_metadata.csv.gz')
# Population-level original-cluster marker summaries for auditable macrophage annotation.
markers=['LST1','TYROBP','FCER1G','CSF1R','CD68','CD14','CD163','MRC1','C1QA','C1QB','C1QC','APOE','SPP1','MMP9','MARCO','MSR1','VSIG4','FCN1','VCAN','S100A8','S100A9','CCR2','SELL','CD1C','FCER1A','CLEC10A','CD1E','CLEC9A','XCR1','CADM1','LAMP3','CCR7','FSCN1','CCL19','GZMB','JCHAIN','TCF4','IL3RA','CLEC4C','CD3D','CD3E','TRAC','MS4A1','EPCAM','KRT19','TPSAB1','TPSB2','CPA3','KIT','FCGR3B','CSF3R','CXCR2']
norm=x.astype(float).multiply(1e4/a.obs.total_counts.to_numpy()[:,None]).tocsr();norm.data=np.log1p(norm.data)
records=[]
for cluster in sorted(a.obs['res.3'].unique()):
 mask=a.obs['res.3'].eq(cluster).to_numpy()
 for g in markers:
  idx=np.flatnonzero(a.var.gene_symbol.eq(g))
  if len(idx):
   c=np.asarray(x[mask][:,idx].sum(1)).ravel();v=np.asarray(norm[mask][:,idx].sum(1)).ravel()
   records.append({'cluster':int(cluster),'gene':g,'cells':int(mask.sum()),'fraction':float((c>0).mean()),'mean_log1p_10k':float(v.mean())})
pd.DataFrame(records).to_csv(R/'original_myeloid_cluster_markers.csv',index=False)
summary={'source_cells':len(meta),'primary_T_samples':int(meta.loc[meta['sample'].str.endswith('T'),'sample'].nunique()),'primary_T_cells':int(meta['sample'].str.endswith('T').sum()),'primary_T_myeloid_cells':int((meta['sample'].str.endswith('T')&meta.celltype.eq('Myeloid')).sum()),'all_sites_myeloid_cells':len(obs),'genes':len(genes),'myeloid_nnz':int(x.nnz)}
(R/'source_summary.json').write_text(json.dumps(summary,indent=2));print(summary,flush=True)
