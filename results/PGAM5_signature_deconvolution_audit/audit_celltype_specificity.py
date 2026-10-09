from pathlib import Path
import sys,gzip,json,os
import pandas as pd,numpy as np
from scipy.io import mmread
from threadpoolctl import threadpool_limits
R=Path(__file__).resolve().parent;name=sys.argv[1];S=Path(os.environ.get('PGAM5_SOURCE_ROOT','/workspace/scratch'))/name
cons=pd.read_csv(R/'primary_consensus_at_least3.csv').gene.tolist()
anchors=['C1QA','C1QB','C1QC','CSF1R','CD68','TYROBP','FCER1G','LST1','AIF1','CD163','MSR1','MRC1','SPP1']
genes=cons+['PGAM5']+anchors
if name=='GSE151530':
 fp=S/'GSE151530_genes.tsv.gz';g=pd.read_csv(fp,sep='\t',header=None,names=['ensembl','gene']);o=pd.read_csv(S/'all_cell_metadata.csv.gz');sel=o.diagnosis.eq('Hepatocellular carcinoma');tam='TAMs';patient='patient_proxy'
else:
 g=pd.read_csv(S/'GSE189903_genes.tsv.gz',sep='\t',header=None,names=['ensembl','gene']);o=pd.read_csv(S/'all_cell_metadata.csv.gz');sel=o.diagnosis.eq('Hepatocellular carcinoma')&o.tissue.isin(['Tumor core','Tumor border']);tam='TAM';patient='patient'
assert o.Cell.is_unique and g.ensembl.is_unique
with threadpool_limits(limits=4):
 with gzip.open(S/(name+'_matrix.mtx.gz'),'rb') as f:x=mmread(f).tocsr()
assert x.shape==(len(g),len(o))
x=x[:,sel.to_numpy()].T.tocsr().astype(np.float64);o=o.loc[sel].copy().reset_index(drop=True)
tot=np.asarray(x.sum(1)).ravel();ng=np.diff(x.indptr);mt=np.asarray(x[:,g.gene.str.startswith('MT-').to_numpy()].sum(1)).ravel()/tot*100
qc=(ng>=500)&(mt<20);x=x[qc];o=o.loc[qc].copy().reset_index(drop=True);tot=tot[qc]
ix=[np.flatnonzero(g.gene.eq(k))[0] for k in genes];assert all(g.gene.eq(k).sum()==1 for k in genes)
z=x[:,ix].toarray()/tot[:,None]*1e4
pos=z[:,genes.index('PGAM5')]>0
labels=o.Type.to_numpy().copy();ismac=o.Type.eq(tam).to_numpy();labels[ismac&pos]='TAM_PGAM5_detected';labels[ismac&~pos]='TAM_PGAM5_undetected';o['reference_group']=labels;o['donor']=o[patient].astype(str).to_numpy();o['total_counts']=tot;o['n_genes']=ng[qc];o['pct_mt']=mt[qc]
np.savez_compressed(R/(name+'_selected_gene_cells.npz'),expression_cp10k=z,genes=np.array(genes));o.to_csv(R/(name+'_QC_cell_metadata.csv.gz'),index=False)
rows=[]
for label in sorted(set(labels)):
 sub=z[labels==label]
 for j,k in enumerate(genes):rows.append({'dataset':name,'group':label,'cells':len(sub),'gene':k,'mean_cp10k':float(sub[:,j].mean()),'detection_fraction':float((sub[:,j]>0).mean()),'role':'consensus_state_candidate' if k in cons else 'definition' if k=='PGAM5' else 'macrophage_lineage_anchor'})
prof=pd.DataFrame(rows);prof.to_csv(R/(name+'_celltype_gene_profiles.csv'),index=False)
piv=prof.pivot(index='gene',columns='group',values='mean_cp10k');nonmac=[c for c in piv if not c.startswith('TAM_PGAM5')];maxother=piv[nonmac].max(1);ratio=piv.TAM_PGAM5_detected/(maxother+1e-9)
a=pd.DataFrame({'gene':piv.index,'PGAM5_positive_TAM_mean_cp10k':piv.TAM_PGAM5_detected,'PGAM5_undetected_TAM_mean_cp10k':piv.TAM_PGAM5_undetected,'max_nonmac_mean_cp10k':maxother,'max_nonmac_group':piv[nonmac].idxmax(1),'positive_TAM_over_max_nonmac_ratio':ratio,'positive_TAM_over_negative_TAM_ratio':piv.TAM_PGAM5_detected/(piv.TAM_PGAM5_undetected+1e-9)})
a['passes_descriptive_2fold_specificity']=a.positive_TAM_over_max_nonmac_ratio.ge(2);a.to_csv(R/(name+'_specificity_audit.csv'),index=False)
print(name,'QC groups',o.reference_group.value_counts().to_dict(),flush=True);print(a[a.gene.isin(cons)].round(3).to_string(index=False),flush=True)
