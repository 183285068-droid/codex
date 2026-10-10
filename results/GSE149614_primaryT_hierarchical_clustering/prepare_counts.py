from pathlib import Path
import gzip,json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
R=Path(__file__).resolve().parent;I=R/'inputs';I.mkdir(exist_ok=True)
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--source-dir',type=Path,default=I);S=parser.parse_args().source_dir
names=['GSE149614_HCC.metadata.updated.txt.gz','GSE149614_HCC.scRNAseq.S71915.count.txt.gz'];expected=['4c547861aec44be917ed2340645d6b30009fff4c77853e878fd61c037bc6e8c8','6010d634ad9a22f30b3331d2dfbbdd7934b52aae54e25182d93c69c451ab46bf']
for name,h in zip(names,expected):
 with open(S/name,'rb') as f:assert hashlib.file_digest(f,'sha256').hexdigest()==h
meta=pd.read_csv(S/names[0],sep='\t');assert meta.Cell.is_unique
obs=meta[meta['sample'].str.endswith('T')].copy().set_index('Cell',drop=False);assert obs.site.eq('Tumor').all() and sorted(obs['sample'].unique())==[f'HCC{i:02}T' for i in range(1,11)]
obs.index.name='cell_id';obs['Sample']=obs['sample'];obs['patient_proxy']=obs.patient;obs['author_type']=obs.celltype;obs['diagnosis']='Hepatocellular carcinoma';rows=[];genes=[]
with gzip.open(S/names[1],'rb') as f:
 header=f.readline().strip().decode().split('\t');assert len(header)==len(meta) and set(header)==set(meta.Cell)
 indices=pd.Index(header).get_indexer(obs.index);assert (indices>=0).all()
 for k,line in enumerate(f):
  gene,sep,values=line.rstrip(b'\r\n').partition(b'\t');assert sep
  vec=np.fromstring(values,dtype=np.int32,sep='\t');assert len(vec)==len(header) and (vec>=0).all();genes.append(gene.decode());rows.append(sparse.csr_matrix(vec[indices].reshape(1,-1)))
  if (k+1)%3000==0:print('genes parsed',k+1,flush=True)
x=sparse.vstack(rows,format='csr').T.tocsr();del rows;unique=pd.Index(genes).unique();mer=sparse.csr_matrix((np.ones(len(genes),dtype=np.int32),(np.arange(len(genes)),unique.get_indexer(genes))),shape=(len(genes),len(unique)));x=(x@mer).tocsr();x.eliminate_zeros()
obs['total_counts']=np.asarray(x.sum(1)).ravel();obs['n_genes']=np.diff(x.indptr);mt=unique.str.startswith('MT-');obs['pct_mt']=np.asarray(x[:,mt].sum(1)).ravel()/obs.total_counts.to_numpy()*100;obs['PGAM5_counts']=x[:,unique.get_loc('PGAM5')].toarray().ravel();obs['PGAM5_detected']=obs.PGAM5_counts.gt(0)
obs.to_csv(I/'all_cell_metadata.csv.gz');sel=obs.n_genes.ge(500)&obs.pct_mt.lt(20);a=ad.AnnData(x[sel.to_numpy()],obs=obs.loc[sel].copy(),var=pd.DataFrame(index=unique));a.layers['counts']=a.X.copy();a.obs['Sample']=a.obs.Sample.astype('category');a.write_h5ad(R/'HCC_raw_QC.h5ad',compression='gzip');obs.groupby('Sample').agg(before=('Cell','size')).join(a.obs.groupby('Sample',observed=True).size().rename('after')).to_csv(R/'sample_QC_counts.csv');(R/'source_summary.json').write_text(json.dumps(dict(source_cells=len(meta),primary_T_samples=10,primary_T_cells=len(obs),after_QC=len(a),source_genes=len(genes),unique_genes=len(unique),input_hashes=dict(zip(names,expected))),indent=2));print('PREPARED',a,flush=True)
