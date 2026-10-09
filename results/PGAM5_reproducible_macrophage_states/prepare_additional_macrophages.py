from pathlib import Path
import json,os
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
from scipy import sparse
ad.settings.allow_write_nullable_strings=True
R=Path(__file__).resolve().parent
rows=[]
for name in ['GSE242889','GSE202642']:
 a=ad.read_h5ad(Path('/workspace/scratch')/name/'macrophage_counts_QC.h5ad');gene=a.var.gene.astype(str);gene=gene.where(gene.ne('nan')&gene.ne(''),a.var_names);unique=pd.Index(pd.unique(gene));mer=sparse.csr_matrix((np.ones(len(gene)),(np.arange(len(gene)),unique.get_indexer(gene))),shape=(len(gene),len(unique)));a=ad.AnnData((a.X@mer).tocsr().astype(np.int32),obs=a.obs.copy(),var=pd.DataFrame(index=unique));a.obs['cell_id']=a.obs_names.astype(str);a.obs.index.name='barcode_index'
 a.obs['donor_id']=name+':'+(a.obs.patient.astype(str) if name=='GSE242889' else a.obs.sample_name.astype(str));a.obs['donor_kind']='patient' if name=='GSE242889' else 'sample_proxy_not_verified_patient';a.obs['dataset']=name
 c=a.copy();sc.pp.normalize_total(c,target_sum=1e4);sc.pp.log1p(c);sc.pp.highly_variable_genes(c,n_top_genes=2000);c.var.loc[c.var_names=='PGAM5','highly_variable']=False;sc.pp.pca(c,n_comps=30,random_state=0);sc.pp.neighbors(c,n_neighbors=15,n_pcs=30,random_state=0);sc.tl.leiden(c,resolution=.8,flavor='igraph',directed=False,n_iterations=2,random_state=0)
 accepted={};markers=['C1QA','C1QB','C1QC','CD68','CSF1R','TYROBP','CD1C','FCER1A','FCN1','S100A8','S100A9','ALB','APOA1','CD3D','CD79A']
 for k,t in c.obs.groupby('leiden',observed=True):
  sel=c.obs.leiden.eq(k).to_numpy();fr={g:float((a.X[sel,a.var_names.get_loc(g)]>0).mean()) if g in a.var_names else np.nan for g in markers};ok=all(fr[g]>=.4 for g in ['C1QA','C1QB','C1QC']) and fr['CD68']>=.3 and fr['CSF1R']>=.2 and fr['TYROBP']>=.5 and not(fr['CD1C']>=.3 and fr['FCER1A']>=.2);accepted[str(k)]=ok;rows.append({'dataset':name,'QC_cluster':str(k),'cells':len(t),'donor_units':t.donor_id.nunique(),'accepted_macrophage':ok,**{'fraction_'+g:v for g,v in fr.items()}})
 a.obs['macrophage_QC_cluster']=c.obs.leiden.astype(str).to_numpy();a.obs['accepted_macrophage']=a.obs.macrophage_QC_cluster.map(accepted).to_numpy();a.obs.to_csv(R/(name+'_candidate_macrophage_metadata.csv.gz'));clean=a[a.obs.accepted_macrophage].copy();clean.write_h5ad(R/(name+'_marker_consistent_macrophages.h5ad'),compression='gzip');print(name,'original',len(a),'accepted',len(clean),flush=True)
pd.DataFrame(rows).to_csv(R/'additional_cohort_macrophage_marker_QA.csv',index=False)
