from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
R=Path(__file__).resolve().parent;a=pd.read_csv(R/'cell_annotations.csv.gz',index_col=0,dtype={'cluster':str});d=pd.read_csv(R/'cluster_annotations.csv');named=list(dict.fromkeys('/'.join(d.subtype_EN.str.split().str[0]).split('/')));genes=pd.read_csv(R/'marker_genes.csv').gene.tolist();raw=sparse.load_npz(R/'marker_raw_counts.npz');full=ad.read_h5ad(R/'HCC_raw_counts_QC.h5ad',backed='r');assert full.obs_names.equals(a.index);missing=[g for g in named if g in full.var_names and g not in genes]
if missing:
 extra=full[:,missing].X;raw=sparse.hstack([raw,extra],format='csr');genes+=missing;sparse.save_npz(R/'marker_raw_counts.npz',raw);pd.Series(genes,name='gene').to_csv(R/'marker_genes.csv',index=False)
full.file.close();vals=raw.astype(float).multiply(10000/a.total_counts.to_numpy()[:,None]).tocsr();vals.data=np.log1p(vals.data);rows=[]
for c in sorted(a.cluster.unique(),key=int):
 ix=np.flatnonzero(a.cluster.eq(c));fra=np.asarray((raw[ix]>0).mean(0)).ravel();means=np.asarray(vals[ix].mean(0)).ravel()
 for j,g in enumerate(genes):rows.append({'cluster':c,'gene':g,'cells':len(ix),'fraction':fra[j],'mean_log1p_CP10k':means[j]})
pd.DataFrame(rows).to_csv(R/'cluster_marker_evidence.csv',index=False)
compact=ad.read_h5ad(R/'HCC_annotation_marker_subset.h5ad');new=ad.AnnData(vals.astype(np.float32),obs=compact.obs.copy(),var=pd.DataFrame(index=pd.Index(genes,name='gene')));new.layers['counts']=raw
for key in compact.obsm:new.obsm[key]=compact.obsm[key].copy()
new.uns['scope']=compact.uns['scope'];new.write_h5ad(R/'HCC_annotation_marker_subset.h5ad',compression='gzip');print('Named markers added',missing,'total',len(genes),flush=True)
