from pathlib import Path
import sys,gzip,json,os
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.io import mmread
from threadpoolctl import threadpool_limits
R=Path(__file__).resolve().parent;S=Path(os.environ.get('PGAM5_SOURCE_ROOT','/workspace/scratch'));name=sys.argv[1];P=S/name
if name=='GSE151530':
 g=pd.read_csv(P/'GSE151530_genes.tsv.gz',sep='\t',header=None,names=['ensembl','gene']);o=pd.read_csv(P/'all_cell_metadata.csv.gz');sel=o.diagnosis.eq('Hepatocellular carcinoma').to_numpy()
 with threadpool_limits(limits=4):
  with gzip.open(P/'GSE151530_matrix.mtx.gz','rb') as f:x=mmread(f).tocsr()
 assert x.shape==(len(g),len(o));x=x[:,sel].T.tocsr().astype(np.int32);o=o.loc[sel].copy().reset_index(drop=True)
 o['donor']=o.patient_proxy.astype(str);o['cell_id']=o.Cell
 coarse=o.Type.map({'TAMs':'TAM','Malignant cells':'Tumor_hepatocyte','T cells':'T_NK','B cells':'B','CAFs':'CAF','TECs':'Endothelial','unclassified':'Unknown_other'})
else:
 o=pd.read_csv(P/'GSE149614_HCC.metadata.updated.txt.gz',sep='\t');assert o.Cell.is_unique
 sel=o['sample'].str.endswith('T');o=o.loc[sel].copy().reset_index(drop=True);o['donor']=o.patient.astype(str);o['cell_id']=o.Cell
 coarse=o.celltype.map({'Myeloid':'Other_myeloid','Hepatocyte':'Tumor_hepatocyte','T/NK':'T_NK','B':'B','Fibroblast':'CAF','Endothelial':'Endothelial'})
 ann=pd.read_csv(P/'myeloid_cluster_annotation.csv');clusters=ann.loc[ann.included_macrophage_cluster,'cluster'].tolist();coarse.loc[o.celltype.eq('Myeloid')&o['res.3'].isin(clusters)]='TAM'
 rows=[];genes=[]
 with gzip.open(P/'GSE149614_HCC.scRNAseq.S71915.count.txt.gz','rb') as f:
  header=f.readline().strip().decode().split('\t');ind=pd.Index(header).get_indexer(o.Cell);assert (ind>=0).all() and len(set(header))==len(header)
  for k,line in enumerate(f):
   gene,sep,values=line.rstrip(b'\r\n').partition(b'\t');assert sep
   v=np.fromstring(values,dtype=np.int32,sep='\t');assert len(v)==len(header) and (v>=0).all()
   rows.append(sparse.csr_matrix(v[ind].reshape(1,-1)));genes.append(gene.decode())
   if (k+1)%3000==0:print(name,'genes parsed',k+1,flush=True)
 x=sparse.vstack(rows,format='csr').T.tocsr();g=pd.DataFrame({'gene':genes});del rows
assert coarse.notna().all() and (x.data>0).all() and o.cell_id.is_unique
# Merge duplicated symbols by summing counts, preserving a source-to-symbol map.
g.to_csv(R/(name+'_source_gene_map.csv'),index=False);unique=pd.Index(pd.unique(g.gene));col=unique.get_indexer(g.gene);merge=sparse.csr_matrix((np.ones(len(g)),(np.arange(len(g)),col)),shape=(len(g),len(unique)));x=(x@merge).tocsr().astype(np.int32)
tot=np.asarray(x.sum(1)).ravel();ng=np.diff(x.indptr);mt=np.asarray(x[:,unique.str.startswith('MT-')].sum(1)).ravel()/tot*100;qc=(ng>=500)&(mt<20)
pg=unique.get_loc('PGAM5');rawpg=x[:,pg].toarray().ravel();o['coarse_group']=coarse.to_numpy();o['PGAM5_counts']=rawpg;o['total_counts']=tot;o['n_genes']=ng;o['pct_mt']=mt;o['passes_QC']=qc;o.to_csv(R/(name+'_all_tumor_cell_metadata.csv.gz'),index=False)
x=x[qc];o=o.loc[qc].copy().reset_index(drop=True)
# Descriptive competitor subdivision, fixed before training/validation; no PGAM5 feature used.
cyclegenes=['MKI67','TOP2A','UBE2C','CENPF','CDK1','BIRC5'];cyclecount=np.asarray((x[:,unique.get_indexer(cyclegenes)]>0).sum(1)).ravel();o['cycling_flag']=cyclecount>=2
base=o.coarse_group.astype(str).to_numpy();labels=base.copy();ismac=base=='TAM';pos=o.PGAM5_counts.to_numpy()>0
labels[ismac&pos]='TAM_PGAM5_detected';labels[ismac&~pos]='TAM_PGAM5_undetected'
# Cycle competitors are separate, except Unknown/Other-myeloid; target class retains original definition.
for cls in ['Tumor_hepatocyte','T_NK','B','CAF','Endothelial','TAM_PGAM5_undetected']:
 mask=(labels==cls)&o.cycling_flag.to_numpy();labels[mask]='Cycling_'+cls
# object dtype preserves full labels and avoids fixed string truncation.
o['fine_group']=labels;o.index=o.cell_id;o.index.name='cell_ID';v=pd.DataFrame({'gene':unique.to_numpy()},index=unique.to_numpy());a=ad.AnnData(x,obs=o,var=v);a.write_h5ad(R/(name+'_all_tumor_counts_QC.h5ad'),compression='gzip')
o.to_csv(R/(name+'_QC_cell_metadata.csv.gz'));o.groupby(['donor','coarse_group','fine_group'],observed=True).size().rename('cells').reset_index().to_csv(R/(name+'_group_cells_by_donor.csv'),index=False)
if name=='GSE151530':assert int((ismac&pos).sum())==180 and int((ismac&~pos).sum())==4314
else:
 expected=ad.read_h5ad(P/'primary_T_macrophage_counts_QC.h5ad',backed='r');assert set(o.loc[ismac,'cell_id'])==set(expected.obs_names);expected.file.close();assert int((ismac&pos).sum())==167 and int((ismac&~pos).sum())==5984
summary={'dataset':name,'QC_cells':len(o),'genes_after_symbol_merge':len(unique),'duplicate_gene_rows_merged':len(g)-len(unique),'donors':o.donor.nunique(),'PGAM5_positive_TAM':int((ismac&pos).sum()),'coarse_groups':o.coarse_group.value_counts().to_dict(),'fine_groups':o.fine_group.value_counts().to_dict(),'cycling_definition':'At least 2 of MKI67/TOP2A/UBE2C/CENPF/CDK1/BIRC5 detected; competitor subdivision only, exploratory and not independently validated labels','scope':'HCC primary tumor; GSE151530 original author HCC/TAM labels; GSE149614 T-ending primary tumors, prior marker-inferred macrophage clusters retained','limits':'GSE151530 coarse labels lack separately annotated NK/DC/monocyte/hepatocyte classes; Unknown preserved. GSE149614 Hepatocyte includes tumor/hepatocyte-like cells, Other_myeloid is not identified in discovery reference. Fine groups are exploratory marker-defined, not author subtypes.'}
(R/(name+'_preparation_summary.json')).write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
