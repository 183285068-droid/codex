from pathlib import Path
import json,hashlib,importlib.metadata as md
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.optimize import nnls
from scipy.stats import pearsonr
import legacy_reference_core as core
R=Path(__file__).resolve().parent;S=Path('/workspace/scratch');T='TAM_PGAM5_detected';N='TAM_PGAM5_undetected'
protocol={'target':'Marker-consistent macrophage AND raw symbol-summed PGAM5 >0; zero is RNA undetected, not protein negative','training':['GSE151530','GSE149614'],'external_primary':['GSE189903','GSE242889'],'external_supplementary':'GSE202642 seven HCC sample proxies, not verified seven independent patients','configuration':'Reuse discovery-selected donor_equal_50 reference, no external feature/config tuning. Fixed NNLS and nonnegative ridge alpha=.01, no solver selected from external performance. Internal LODO fixed config is diagnostic, previous nested NNLS analysis remains authoritative for tuning-adjusted internal performance','annotation':'Existing marker-consistent macrophages retained for all five cohorts. For242/202 other cells annotate existing broad clusters using PGAM5-free lineage markers; unresolved cells retained as Unknown_other. Previously used full-data graphs may contain PGAM5; they do not define target membership. No new ambient/doublet correction','mixtures':'2000 cells sampled with replacement; total macrophage fraction .25; target cell fractions0/.005/.01/.02/.05, two repeats per donor. 50% tumor/epithelial proxy,25% remaining lineages; if unavailable donor excluded','units':'equal_RNA_cell_fraction is equalized cell CP10k mixture, an idealized cell fraction, not native bulk. pooled_count_RNA_fraction truth is observed library-count contribution, not calibrated tissue cell fraction','gates':{'standard_MAE_pp_max':1,'standard_Pearson_min':.7,'zero_target_p95_percent_max':.5,'primary_external_donors_min_each':3,'primary_external_Monocyte_DC_testable_donors_min_each':2},'additional_limits':'All cohorts explored previously; retrospective tests. Repeats not biological independent replicates. Passing computational gates still requires bulk platform/abundance-unit calibration and orthogonal biological confirmation before confirmed TCGA cellular abundance. No sequencing-depth matching/covariates, no patient-paired DE'}
protocol['feature_availability']='Frozen training features are intersected with external measured gene symbols; no missing gene imputed as expression zero in fitting. Require >=95% frozen feature coverage. This availability adaptation uses no expression/outcome tuning.'
(R/'protocol.json').write_text(json.dumps(protocol,indent=2));genes=core.genes;models={};pred=[];qa=[];inputs=[];coverage=[];anno=[];grpsummary=[];availability={};audit_cases=[]
for n in ['GSE151530','GSE149614','GSE189903']:
 o=core.D[n]['obs'].copy();o['cell_type']=o.harmonized_group;o['PGAM5_RNA_status']=np.where(o.harmonized_group.eq('Macrophage'),np.where(o.PGAM5_counts.gt(0),'detected','undetected'),'not_target_lineage');o['PGAM5_macrophage_annotation']=np.where(o.PGAM5_RNA_status.eq('detected'),'Macrophage_PGAM5_detected',np.where(o.PGAM5_RNA_status.eq('undetected'),'Macrophage_PGAM5_undetected',o.cell_type));o.to_csv(R/(n+'_all_cell_annotation.csv.gz'),index=False);anno.append(o)
# Fit reference using training data only, same configuration already selected internally in v2.
locked=core.model(core.donors,'donor_equal_50');pd.DataFrame(locked['ref'][locked['ix']],index=genes[locked['ix']],columns=locked['groups']).to_csv(R/'EXPLORATORY_REFERENCE_NOT_VALIDATED.tsv',sep='\t');pd.DataFrame({'gene':genes[locked['ix']],'row_scale':locked['scale']}).to_csv(R/'reference_row_scales.csv',index=False)
(R/'locked_reference_manifest.json').write_text(json.dumps({'groups':locked['groups'],'genes':genes[locked['ix']].tolist(),'training_donors':core.donors,'status':'LOCKED; TCGA_USE_REQUIRES_VALIDATION','solvers':['NNLS','nonnegative_ridge_alpha0.01']},indent=2))
marker=['C1QA','C1QB','C1QC','CD68','CSF1R','TYROBP','FCN1','S100A8','S100A9','CD1C','FCER1A','GZMB','TCF4','IL3RA','CLEC4C','TPSAB1','TPSB2','KIT','CD3D','CD3E','TRAC','GNLY','NKG7','MS4A1','CD79A','MZB1','JCHAIN','COL1A1','DCN','PECAM1','VWF','ALB','APOA1','EPCAM','KRT19']
for n in ['GSE242889','GSE202642']:
 p=S/n/('all_tumor_counts_QC.h5ad' if n=='GSE242889' else 'all_library_counts_QC.h5ad');a=ad.read_h5ad(p);inputs.append(p)
 normp=S/n/('all_tumor_normalized_clustered.h5ad' if n=='GSE242889' else 'all_library_normalized_clustered.h5ad');b=ad.read_h5ad(normp,backed='r');a.obs['broad_cluster']=b.obs.leiden.reindex(a.obs_names).astype(str).to_numpy();b.file.close();inputs.append(normp)
 if n=='GSE202642':a=a[a.obs.library_suffix.astype(int).between(5,11)].copy()
 labels=a.var.gene.astype(str);labels=labels.where(labels.ne('nan')&labels.ne(''),a.var_names);unique=pd.Index(pd.unique(labels));mer=sparse.csr_matrix((np.ones(len(labels)),(np.arange(len(labels)),unique.get_indexer(labels))),shape=(len(labels),len(unique)));x=(a.X@mer).tocsr().astype(np.int32);o=a.obs.copy();o['cell_id']=o.index.astype(str);o['dataset']=n;o['donor']=n+':'+(o.patient.astype(str) if n=='GSE242889' else o.library_suffix.astype(str));o['donor_kind']='patient' if n=='GSE242889' else 'sample_proxy';o['PGAM5_counts']=x[:,unique.get_loc('PGAM5')].toarray().ravel();dec={}
 for k,t in o.groupby('broad_cluster',observed=True):
  ix=o.index.isin(t.index);f={g:float((x[ix,unique.get_loc(g)]>0).mean()) if g in unique else 0. for g in marker}; ge=lambda gs,thr:all(f[g]>=thr for g in gs);lab='Unknown_other'
  if sum(f[g]>=.2 for g in ['TPSAB1','TPSB2','KIT'])>=2:lab='Mast'
  elif f['GZMB']>=.3 and sum(f[g]>=.1 for g in ['TCF4','IL3RA','CLEC4C'])>=2 and f['CD3D']<.2:lab='pDC'
  elif ge(['CD1C'],.3) and ge(['FCER1A'],.2):lab='Mixed_APC' if f['C1QC']>=.4 else 'DC_enriched'
  elif f['FCN1']>=.35 and max(f['S100A8'],f['S100A9'])>=.4 and f['TYROBP']>=.5 and f['C1QC']<.4:lab='Monocyte_enriched'
  elif (sum(f[g]>=.3 for g in ['CD3D','CD3E','TRAC'])>=2 or (f['GNLY']>=.3 and f['NKG7']>=.4)) and f['C1QC']<.2:lab='T_NK'
  elif ge(['MZB1','JCHAIN'],.3) and f['C1QC']<.2:lab='Plasma'
  elif ge(['MS4A1','CD79A'],.3) and f['C1QC']<.2:lab='B'
  elif ge(['COL1A1','DCN'],.3) and f['TYROBP']<.2:lab='CAF'
  elif ge(['PECAM1','VWF'],.3) and f['TYROBP']<.2:lab='Endothelial'
  elif (ge(['EPCAM','KRT19'],.2) or ge(['ALB','APOA1'],.4)) and f['TYROBP']<.2 and f['C1QC']<.2:lab='Tumor_hepatocyte'
  dec[k]=lab;qa.append({'dataset':n,'cluster':k,'cells':len(t),'assigned_nonmacrophage_group':lab,**{'fraction_'+g:v for g,v in f.items()}})
 o['harmonized_group']=o.broad_cluster.map(dec)
 # Target lineage strictly uses the PGAM5-free macrophage re-evaluation from prior state study, never broad-cluster target inference.
 prior=pd.read_csv(S/'PGAM5_reproducible_macrophage_states'/(n+'_macrophage_state_metadata.csv.gz'));macids=set(prior.cell_id.astype(str));assert macids<=set(o.cell_id);o.loc[o.cell_id.isin(macids),'harmonized_group']='Macrophage'
 cycgenes=['MKI67','TOP2A','UBE2C','CENPF','CDK1','BIRC5'];o['cycling_flag']=np.asarray((x[:,unique.get_indexer(cycgenes)]>0).sum(1)).ravel()>=2;fine=o.harmonized_group.to_numpy(dtype=object).copy();mac=o.harmonized_group.eq('Macrophage').to_numpy();pos=o.PGAM5_counts.gt(0).to_numpy();fine[mac&pos]=T;fine[mac&~pos]=N
 for g in ['Tumor_hepatocyte','T_NK','B','CAF','Endothelial',N]:fine[(fine==g)&o.cycling_flag.to_numpy()]='Cycling_'+g
 o['fine_group']=fine;o['cell_type']=o.harmonized_group;o['PGAM5_RNA_status']=np.where(mac,np.where(pos,'detected','undetected'),'not_target_lineage');o['PGAM5_macrophage_annotation']=np.where(mac,np.where(pos,'Macrophage_PGAM5_detected','Macrophage_PGAM5_undetected'),o.cell_type);o.to_csv(R/(n+'_all_cell_annotation.csv.gz'),index=False);anno.append(o)
 common_ix=unique.get_indexer(genes);present=common_ix>=0;mapping=sparse.csr_matrix((np.ones(present.sum()),(common_ix[present],np.flatnonzero(present))),shape=(len(unique),len(genes)));raw=(x@mapping).tocsr().astype(float);total=o.total_counts.to_numpy(dtype=float);core.D[n]={'raw':raw,'norm':raw.multiply(1e4/total[:,None]).tocsr(),'total':total,'obs':o.reset_index(drop=True)}
 coverage.append({'dataset':n,'all_HCC_QC_cells':len(o),'macrophages':int(mac.sum()),'target':int((mac&pos).sum()),'reference_features_available':int(present[locked['ix']].sum()),'reference_features_total':len(locked['ix'])});availability[n]=present;assert present[locked['ix']].mean()>=.95;pd.DataFrame({'missing_reference_gene':genes[locked['ix']][~present[locked['ix']]]}).to_csv(R/(n+'_missing_reference_genes.csv'),index=False);print('Prepared',n,len(o),o.harmonized_group.value_counts().to_dict(),flush=True);del a,x,raw
pd.DataFrame(qa).to_csv(R/'external_broad_lineage_marker_QA.csv',index=False);pd.DataFrame(coverage).to_csv(R/'external_feature_coverage.csv',index=False)
for o in anno:
 for keys,t in o.groupby(['dataset','donor','cell_type','PGAM5_RNA_status'],observed=True):grpsummary.append(dict(zip(['dataset','donor','cell_type','PGAM5_RNA_status'],keys))|{'cells':len(t)})
pd.DataFrame(grpsummary).to_csv(R/'annotation_counts_by_donor.csv',index=False)
def fit(m,y,solver):
 A=m['matrix'];b=y[m['ix']]/m['scale']
 if solver=='NNLS':coef,res=nnls(A,b,maxiter=2000)
 else:
  A2=np.vstack([A,np.sqrt(.01)*np.eye(A.shape[1])]);b2=np.r_[b,np.zeros(A.shape[1])];coef,_=nnls(A2,b2,maxiter=2000);res=np.linalg.norm(A@coef-b)
 return coef/max(coef.sum(),1e-12),res
# Shared controlled mixtures, new target-absent/high-PGAM5 nonmacrophage challenge.
for n in ['joint_discovery','GSE189903','GSE242889','GSE202642']:
 eligible=core.eligible(n);print('Eligible',n,eligible,flush=True)
 for d in eligible:
  m=core.model([k for k in core.donors if k!=d],'donor_equal_50') if n=='joint_discovery' else locked
  if n in availability:
   take=availability[n][m['ix']];m={**m,'ix':m['ix'][take],'scale':m['scale'][take],'matrix':m['matrix'][take]}
  mix=core.mixtures(n,d,replicates=2);data=core.D[n];o=data['obs'];rng=np.random.default_rng(470+sum((n+d).encode()));mask=o.donor.eq(d).to_numpy();pg=genes.tolist().index('PGAM5');nonmac=mask&~o.harmonized_group.eq('Macrophage').to_numpy()&(np.asarray(data['raw'][:,pg].todense()).ravel()>0);hp=np.flatnonzero(nonmac);neg=np.flatnonzero(mask&o.harmonized_group.eq('Macrophage').to_numpy()&o.fine_group.ne(T).to_numpy())
  if len(hp)>=20 and len(neg)>=20:
   for rep in range(2):
    chosen=np.r_[rng.choice(neg,500,replace=True),rng.choice(hp,1500,replace=True)];ye=np.asarray(data['norm'][chosen].mean(0)).ravel();yr=np.asarray(data['raw'][chosen].sum(0)).ravel()/data['total'][chosen].sum()*1e4
    for unit,y in [('equal_RNA_cell_fraction',ye),('pooled_count_RNA_fraction',yr)]:mix.append({'scenario':'zero_PGAM5_detected_nonmacrophages','nominal_fraction':0,'replicate':rep,'unit':unit,'truth':0.,'truth_total_TAM':.25 if unit.startswith('equal') else float(data['total'][chosen[:500]].sum()/data['total'][chosen].sum()),'y':y})
  for z in mix:
   for solver in ['NNLS','nonnegative_ridge_alpha0.01']:
    coef,res=fit(m,z['y'],solver)
    if d==eligible[0] and z['replicate']==0 and ((z['scenario']=='standard' and z['nominal_fraction'] in [0,.05]) or z['scenario']=='zero_PGAM5_detected_nonmacrophages'):
     case_id=len(audit_cases);np.savez_compressed(R/('solver_audit_case_'+str(case_id)+'.npz'),A=m['matrix'],b=z['y'][m['ix']]/m['scale'],weights=coef,target_index=m['groups'].index(T),solver=solver);audit_cases.append({'case':case_id,'dataset':n,'solver':solver,'scenario':z['scenario'],'unit':z['unit']})
    pred.append({k:v for k,v in z.items() if k!='y'}|{'dataset':n,'donor':d,'phase':'discovery_LODO_diagnostic' if n=='joint_discovery' else 'supplementary_external' if n=='GSE202642' else 'primary_external','solver':solver,'predicted':float(coef[m['groups'].index(T)]),'predicted_total_macrophage':float(sum(coef[j] for j,g in enumerate(m['groups']) if g in [T,N,'Cycling_'+N])),'residual':res,'features':len(m['ix'])})
  print('Tested',n,d,len(mix),flush=True)
pd.DataFrame(audit_cases).to_csv(R/'solver_audit_case_manifest.csv',index=False)
p=pd.DataFrame(pred);p.to_csv(R/'validation_predictions.csv',index=False);summary=[]
for keys,t in p.groupby(['dataset','phase','solver','unit','scenario']):
 zero=t[t.truth.eq(0)];corr=float(pearsonr(t.truth,t.predicted).statistic) if t.truth.nunique()>1 and t.predicted.nunique()>1 else np.nan
 summary.append(dict(zip(['dataset','phase','solver','unit','scenario'],keys))|{'donors':t.donor.nunique(),'mixtures':len(t),'MAE_pp':float(abs(t.predicted-t.truth).mean()*100),'Pearson_r':corr,'zero_target_p95_percent':float(zero.predicted.quantile(.95)*100) if len(zero) else np.nan})
su=pd.DataFrame(summary);su.to_csv(R/'validation_summary.csv',index=False);checks=[]
for solver in p.solver.unique():
 for dataset in ['GSE189903','GSE242889']:
  for unit in p.unit.unique():
   st=su[(su.dataset==dataset)&(su.solver==solver)&(su.unit==unit)&(su.scenario=='standard')]
   for crit,col,fn in [('MAE_pp<=1','MAE_pp',lambda v:v<=1),('Pearson>=0.7','Pearson_r',lambda v:v>=.7),('zero_p95_percent<=0.5','zero_target_p95_percent',lambda v:v<=.5),('independent_donors>=3','donors',lambda v:v>=3)]:
    val=float(st.iloc[0][col]) if len(st) else None;checks.append({'solver':solver,'dataset':dataset,'unit':unit,'criterion':crit,'value':val,'passed':bool(val is not None and np.isfinite(val) and fn(val))})
   for scenario in ['zero_Monocyte_enriched','zero_DC_enriched','cycling_zero','zero_PGAM5_detected_nonmacrophages']:
    st=su[(su.dataset==dataset)&(su.solver==solver)&(su.unit==unit)&(su.scenario==scenario)];val=float(st.iloc[0].zero_target_p95_percent) if len(st) else None;count=int(st.iloc[0].donors) if len(st) else 0;checks.append({'solver':solver,'dataset':dataset,'unit':unit,'criterion':scenario+'_p95<=0.5_and_donors>=2','value':val,'donors':count,'passed':bool(val is not None and val<=.5 and count>=2)})
pd.DataFrame(checks).to_csv(R/'predeclared_gate_checks.csv',index=False);status={solver:all(c['passed'] for c in checks if c['solver']==solver) for solver in p.solver.unique()};result={'operational_RNA_annotation_created':True,'solver_external_gates':status,'reference_validated_for_TCGA_cellular_abundance':False,'TCGA_NOT_RUN':True,'limits':protocol['additional_limits'],'interpretation':'Annotation labels raw RNA detection, not a validated protein-positive subtype. No actual TCGA patient abundance generated.'};(R/'result.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
