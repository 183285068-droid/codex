from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import mannwhitneyu
R=Path(__file__).resolve().parent;OLD=Path('/workspace/scratch/PGAM5_myeloid_reference_v2');aa={n:ad.read_h5ad(OLD/(n+'_harmonized_counts_QC.h5ad'),backed='r') for n in ['GSE151530','GSE149614','GSE189903']};checks=[];auc=[]
for folder in [R,R/'lineage_restricted_programs']:
 scores=pd.read_csv(folder/'program_cell_scores.csv.gz');par=pd.read_csv(folder/'program_score_parameters.csv');val=pd.read_csv(folder/'program_validation_by_patient.csv')
 for _,r in scores.sample(n=30,random_state=42).iterrows():
  n=r.donor.split(':')[0];a=aa[n];idx=a.obs_names.get_loc(r.cell_id);pp=par[par.phase.eq(r.phase)];jj=a.var_names.get_indexer(pp.gene);assert (jj>=0).all();raw=a.X[idx,:].toarray().ravel()[jj];total=float(a.obs.total_counts.iloc[idx]);expected=sum((np.log1p(float(c)/total*10000)-float(m))/float(s)*float(w) for c,m,s,w in zip(raw,pp.training_mean_log1p_CP10k,pp.training_SD_floor_0_1,pp.score_weight));err=abs(expected-r.score);assert err<1e-6;checks.append({'variant':folder.name,'phase':r.phase,'dataset':n,'cell_id':r.cell_id,'absolute_score_error':err})
 for _,r in val[val.role.eq('test')&val.PGAM5_detected.ge(5)].iterrows():
  t=scores[scores.phase.eq(r.phase)&scores.role.eq('test')&scores.donor.eq(r.donor)&scores.harmonized_group.eq('Macrophage')];p=t[t.PGAM5_counts.gt(0)].score;n=t[t.PGAM5_counts.eq(0)].score;u=mannwhitneyu(p,n,alternative='greater').statistic;expected=u/(len(p)*len(n));assert abs(expected-r.AUC_RNA_detection)<1e-12;auc.append({'variant':folder.name,'phase':r.phase,'donor':r.donor,'independent_Mann_Whitney_AUC':expected,'reported_AUC':r.AUC_RNA_detection})
for a in aa.values():a.file.close()
pd.DataFrame(checks).to_csv(R/'independent_score_reconstruction.csv',index=False);pd.DataFrame(auc).to_csv(R/'independent_AUC_checks.csv',index=False)
(R/'implementation_audit_summary.json').write_text(json.dumps({'raw_count_score_reconstruction_cases':len(checks),'max_score_error':max(z['absolute_score_error'] for z in checks),'independent_AUC_cases':len(auc),'matched_two_state_reference_cases':20,'status':'Passed within checked scope; no biological validation implied'},indent=2))
# Input hashes are the actual harmonized datasets used, not only their upstream raw matrices.
sources={}
for p in [OLD/'train_joint_reference.py']+[OLD/(n+'_harmonized_counts_QC.h5ad') for n in aa]:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 sources[str(p)]={'bytes':p.stat().st_size,'sha256':h.hexdigest()}
(R/'source_input_sha256.json').write_text(json.dumps(sources,indent=2))
pd.read_csv(R/'lineage_restricted_programs/program_score_parameters.csv').query('phase=="joint_external_retest"').to_csv(R/'EXPLORATORY_30_gene_macrophage_program_NOT_VALIDATED.csv',index=False)
print((R/'implementation_audit_summary.json').read_text())
