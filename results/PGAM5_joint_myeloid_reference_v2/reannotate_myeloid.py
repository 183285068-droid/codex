from pathlib import Path
import sys,json,gzip,os
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
from scipy.io import mmread
from scipy import sparse
from threadpoolctl import threadpool_limits
R=Path(__file__).resolve().parent;name=sys.argv[1];S=Path(os.environ.get('PGAM5_SOURCE_ROOT','/workspace/scratch'));old=Path(os.environ.get('PGAM5_PREVIOUS_REFERENCE_ROOT',str(S/'PGAM5_reference_validation')))
if name!='GSE189903':
 a=ad.read_h5ad(old/(name+'_all_tumor_counts_QC.h5ad'));a.obs['original_group']=a.obs.coarse_group.astype(str);candidate=a.obs.coarse_group.isin(['TAM','Other_myeloid','Unknown_other']).to_numpy()
else:
 P=S/name;g=pd.read_csv(P/'GSE189903_genes.tsv.gz',sep='\t',header=None,names=['ensembl','gene']);o=pd.read_csv(P/'all_cell_metadata.csv.gz');sel=o.diagnosis.eq('Hepatocellular carcinoma')&o.tissue.isin(['Tumor core','Tumor border'])
 with threadpool_limits(limits=4):
  with gzip.open(P/'GSE189903_matrix.mtx.gz','rb') as f:x=mmread(f).tocsr()
 x=x[:,sel.to_numpy()].T.tocsr().astype(np.int32);o=o.loc[sel].copy().reset_index(drop=True);tot=np.asarray(x.sum(1)).ravel();ng=np.diff(x.indptr);mt=np.asarray(x[:,g.gene.str.startswith('MT-').to_numpy()].sum(1)).ravel()/tot*100;qc=(ng>=500)&(mt<20);x=x[qc];o=o.loc[qc].copy().reset_index(drop=True)
 unique=pd.Index(pd.unique(g.gene));merge=sparse.csr_matrix((np.ones(len(g)),(np.arange(len(g)),unique.get_indexer(g.gene))),shape=(len(g),len(unique)));x=(x@merge).tocsr().astype(np.int32)
 o['cell_id']=o.Cell;o['donor']=o.patient.astype(str);o['total_counts']=tot[qc];o['n_genes']=ng[qc];o['pct_mt']=mt[qc];o['original_group']=o.Type.map({'TAM':'TAM','unclassified':'Unknown_other','T-cell':'T_NK','B-cell':'B','CAF':'CAF','TEC':'Endothelial','Malignant cell':'Tumor_hepatocyte'});assert o.original_group.notna().all();o.index=o.cell_id;a=ad.AnnData(x,obs=o,var=pd.DataFrame({'gene':unique},index=unique));candidate=o.original_group.isin(['TAM','Unknown_other']).to_numpy()
pg=a.var_names.get_loc('PGAM5');a.obs['PGAM5_counts']=a.X[:,pg].toarray().ravel();a.obs['dataset']=name;a.obs['donor_id']=name+':'+a.obs.donor.astype(str);a.obs['harmonized_group']=a.obs.original_group.astype(str).to_numpy();a.obs['myeloid_cluster']='not_reclustered'
c=a[candidate].copy();sc.pp.normalize_total(c,target_sum=1e4);sc.pp.log1p(c);sc.pp.highly_variable_genes(c,n_top_genes=2000,flavor='seurat');c.var.loc[c.var_names=='PGAM5','highly_variable']=False;sc.pp.pca(c,n_comps=30,random_state=0);sc.pp.neighbors(c,n_neighbors=15,n_pcs=30,random_state=0);sc.tl.leiden(c,resolution=.8,flavor='igraph',directed=False,n_iterations=2,random_state=0)
c.write_h5ad(R/(name+'_candidate_clustered.h5ad'),compression='gzip')
markers=['C1QA','C1QB','C1QC','CD68','CSF1R','TYROBP','LST1','FCER1G','FCN1','VCAN','S100A8','S100A9','FCGR3A','CD1C','FCER1A','CLEC10A','CLEC9A','XCR1','GZMB','JCHAIN','TCF4','IL3RA','CLEC4C','TPSAB1','TPSB2','CPA3','KIT','CD3D','CD3E','TRAC','GNLY','NKG7','MS4A1','CD79A','COL1A1','DCN','PECAM1','VWF','ALB','APOA1','EPCAM','KRT19']
rows=[];decisions={}
for cluster in c.obs.leiden.cat.categories:
 sel=c.obs.leiden.eq(cluster).to_numpy();r={'dataset':name,'cluster':str(cluster),'cells':int(sel.sum()),'donors':int(c.obs.loc[sel,'donor'].nunique())};f={}
 for g in markers:
  if g in c.var_names:
   v=c.X[sel,c.var_names.get_loc(g)].toarray().ravel();f[g]=float((v>0).mean());r['fraction_'+g]=f[g];r['mean_log1p_'+g]=float(v.mean())
  else:f[g]=np.nan;r['fraction_'+g]=np.nan;r['mean_log1p_'+g]=np.nan
 def allge(gs,t):return all(pd.notna(f[g]) and f[g]>=t for g in gs)
 def lt(g,t):return pd.notna(f[g]) and f[g]<t
 label='Unknown_other'
 if sum(pd.notna(f[g]) and f[g]>=.2 for g in ['TPSAB1','TPSB2','KIT'])>=2:label='Mast'
 elif allge(['GZMB'],.3) and sum(pd.notna(f[g]) and f[g]>=.1 for g in ['TCF4','IL3RA','CLEC4C'])>=2 and lt('CD3D',.2):label='pDC'
 elif allge(['CD1C'],.3) and allge(['FCER1A'],.2) and lt('C1QC',.4):label='DC_enriched'
 elif allge(['FCN1'],.35) and any(pd.notna(f[g]) and f[g]>=.4 for g in ['S100A8','S100A9']) and allge(['TYROBP'],.5) and lt('C1QC',.4):label='Monocyte_enriched'
 elif allge(['CD1C'],.3) and allge(['FCER1A'],.2) and allge(['C1QC'],.4):label='Mixed_APC'
 elif allge(['C1QA','C1QB','C1QC'],.4) and allge(['CSF1R'],.2) and allge(['CD68'],.3) and allge(['TYROBP'],.5):label='Macrophage'
 elif (sum(pd.notna(f[g]) and f[g]>=.3 for g in ['CD3D','CD3E','TRAC'])>=2 or (allge(['GNLY'],.3) and allge(['NKG7'],.4))) and lt('C1QC',.2):label='T_NK'
 elif allge(['MS4A1','CD79A'],.3) and lt('C1QC',.2):label='B'
 elif allge(['COL1A1','DCN'],.3) and lt('TYROBP',.2):label='CAF'
 elif allge(['PECAM1','VWF'],.3) and lt('TYROBP',.2):label='Endothelial'
 elif (allge(['EPCAM','KRT19'],.2) or allge(['ALB','APOA1'],.4)) and lt('TYROBP',.2) and lt('C1QC',.2):label='Tumor_hepatocyte'
 r['assigned_group']=label;r['source_groups_json']=json.dumps(c.obs.loc[sel,'original_group'].value_counts().to_dict());rows.append(r);decisions[str(cluster)]=label
q=pd.DataFrame(rows);q.to_csv(R/(name+'_cluster_annotation_QA.csv'),index=False);mapped=c.obs.leiden.astype(str).map(decisions);a.obs.loc[c.obs_names,'harmonized_group']=mapped.to_numpy();a.obs.loc[c.obs_names,'myeloid_cluster']=c.obs.leiden.astype(str).to_numpy()
# Preserve all original non-candidate identities; define target only within marker-consistent macrophage clusters.
base=a.obs.harmonized_group.astype(str).to_numpy().copy();base[base=='TAM']='Macrophage';a.obs['harmonized_group']=base;labels=base.copy();mac=base=='Macrophage';pos=a.obs.PGAM5_counts.to_numpy()>0;labels[mac&pos]='TAM_PGAM5_detected';labels[mac&~pos]='TAM_PGAM5_undetected'
cyc=['MKI67','TOP2A','UBE2C','CENPF','CDK1','BIRC5'];cycling=np.asarray((a.X[:,a.var_names.get_indexer(cyc)]>0).sum(1)).ravel()>=2;a.obs['cycling_flag']=cycling
for group in ['Tumor_hepatocyte','T_NK','B','CAF','Endothelial','TAM_PGAM5_undetected']:
 labels[(labels==group)&cycling]='Cycling_'+group
assert not np.any(labels=='TAM');a.obs['fine_group']=labels;a.obs.to_csv(R/(name+'_harmonized_cell_metadata.csv.gz'));a.write_h5ad(R/(name+'_harmonized_counts_QC.h5ad'),compression='gzip');a.obs.groupby(['donor_id','original_group','harmonized_group','fine_group'],observed=True).size().rename('cells').reset_index().to_csv(R/(name+'_annotation_transitions.csv'),index=False)
summary={'dataset':name,'cells':len(a),'reclustered_candidates':len(c),'clusters':len(q),'groups':a.obs.harmonized_group.value_counts().to_dict(),'target_positive_cells':int((mac&pos).sum()),'target_undetected_cells':int((mac&~pos).sum()),'annotation':'Frozen marker rules, PGAM5 excluded from graph/HVG; descriptive inferred cluster labels, unknown retained; no cDC subtype claim'};(R/(name+'_annotation_summary.json')).write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True);print(q[['cluster','cells','assigned_group']].to_string(index=False),flush=True)
