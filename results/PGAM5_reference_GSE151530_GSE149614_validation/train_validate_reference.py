from pathlib import Path
import os,json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.optimize import nnls
from scipy.stats import pearsonr
R=Path(__file__).resolve().parent;protocol=json.loads((R/'validation_protocol.json').read_text())
A=ad.read_h5ad(R/'GSE151530_all_tumor_counts_QC.h5ad');B=ad.read_h5ad(R/'GSE149614_all_tumor_counts_QC.h5ad')
genes=np.array(sorted(set(A.var_names)&set(B.var_names)));assert 'PGAM5' in genes and len(genes)>15000
for a in [A,B]:
 assert a.obs_names.is_unique and a.var_names.is_unique and (a.X.data>=0).all()
D={}
for name,a in [('GSE151530',A),('GSE149614',B)]:
 raw=a[:,genes].X.astype(float).tocsr();tot=a.obs.total_counts.to_numpy(dtype=float);norm=raw.multiply(1e4/tot[:,None]).tocsr();D[name]={'raw':raw,'norm':norm,'obs':a.obs.reset_index(drop=True).copy(),'total':tot}
del A,B
q=D['GSE151530'];obs=q['obs'];fine=obs.fine_group.astype(str).to_numpy();coarse=obs.coarse_group.astype(str).to_numpy();donors=sorted(obs.donor.astype(str).unique());donor=obs.donor.astype(str).to_numpy();types=sorted(set(fine));typeindex={g:i for i,g in enumerate(types)}
# Donor-by-group sufficient statistics enable feature selection without inspecting held-out/external expression.
sums=np.zeros((len(donors),len(types),len(genes)));det=np.zeros_like(sums);num=np.zeros((len(donors),len(types)),dtype=int)
for i,d in enumerate(donors):
 for j,g in enumerate(types):
  mask=(donor==d)&(fine==g);n=mask.sum();num[i,j]=n
  if n:sums[i,j]=np.asarray(q['norm'][mask].sum(0)).ravel();det[i,j]=np.asarray((q['raw'][mask]>0).sum(0)).ravel()
print('Sufficient statistics',sums.shape,'shared genes',len(genes),flush=True)
feature_eligible=np.array([not (g.startswith('MT-') or g.startswith('RPL') or g.startswith('RPS')) for g in genes])
configs={'pooled_50':('pooled',50),'donor_equal_50':('donor_equal',50),'donor_equal_100':('donor_equal',100)}
pg=np.flatnonzero(genes=='PGAM5')[0]
TARGET='TAM_PGAM5_detected';NEG='TAM_PGAM5_undetected'
def parent(g):return g.removeprefix('Cycling_')
def model(train_donors,config,with_pg=True,with_cycle=True):
 trainidx=[donors.index(d) for d in train_donors];n=num[trainidx].sum(0);active=[g for g in types if not g.startswith('Cycling_') or (with_cycle and n[typeindex[g]]>=100)]
 groups=sorted({parent(g) if g not in active else g for g in types});refs=[];fractions=[];support=[];mode,top=configs[config]
 for label in groups:
  idx=[typeindex[g] for g in types if (g if g in active else parent(g))==label]
  nc=num[trainidx][:,idx].sum(1);su=sums[trainidx][:,idx].sum(1);de=det[trainidx][:,idx].sum(1);assert nc.sum()>0
  if mode=='pooled':r=su.sum(0)/nc.sum()
  else:r=(su[nc>0]/nc[nc>0,None]).mean(0)
  refs.append(r);fractions.append(de.sum(0)/nc.sum());support.append(int((nc>0).sum()))
 ref=np.column_stack(refs);freq=np.column_stack(fractions);features=set();marker_rows=[]
 for j,label in enumerate(groups):
  other=np.max(np.delete(ref,j,axis=1),axis=1);contrast=np.log2((ref[:,j]+.05)/(other+.05));rank=contrast*np.minimum(1,ref[:,j])
  elig=feature_eligible&(freq[:,j]>=.10)&(ref[:,j]>=.05)
  if not with_pg:elig[pg]=False
  candidates=np.flatnonzero(elig);selected=candidates[np.argsort(-rank[candidates],kind='stable')[:top]];features.update(selected.tolist())
  for k in selected:marker_rows.append({'gene':genes[k],'reference_group':label,'mean_cp10k':ref[k,j],'max_other_mean_cp10k':other[k],'contrast_log2':contrast[k],'training_detection_fraction':freq[k,j],'training_donors_with_group':support[j]})
 if with_pg:features.add(int(pg))
 ix=np.array(sorted(features));M=ref[ix];scale=np.sqrt(np.mean(M*M,axis=1));scale=np.maximum(scale,1e-6);Ms=np.asfortranarray(M/scale[:,None])
 return {'groups':groups,'ref':ref,'ix':ix,'scale':scale,'matrix':Ms,'markers':marker_rows,'training_donors':train_donors,'config':config,'with_pg':with_pg,'with_cycle':with_cycle}
def estimate(m,y):
 beta,res=nnls(m['matrix'],y[m['ix']]/m['scale'],maxiter=1000)
 value=beta[m['groups'].index(TARGET)]/max(beta.sum(),1e-12)
 return float(value),float(res)
def eligible(name):
 o=D[name]['obs'];rows=[]
 for d,t in o.groupby('donor',observed=True):
  counts=t.fine_group.value_counts();pos=counts.get(TARGET,0);neg=t[t.coarse_group.eq('TAM')&~t.fine_group.eq(TARGET)].shape[0]
  if pos>=5 and neg>=20 and t.coarse_group.eq('Tumor_hepatocyte').sum()>=20 and t.coarse_group.eq('T_NK').sum()>=20:rows.append(str(d))
 return sorted(rows)
internal=eligible('GSE151530');external=eligible('GSE149614');assert len(internal)>=5 and len(external)>=5
print('Qualified internal',internal,'external',external,flush=True)
# Fixed complete-gene mixtures, then projected into each training-selected panel; shared observations between configurations.
def mixtures(name,d,replicates,challenge=True,sampled=True):
 data=D[name];o=data['obs'];mask=o.donor.astype(str).eq(d).to_numpy();pools={g:np.flatnonzero(mask&o.fine_group.eq(g).to_numpy()) for g in sorted(o.fine_group.unique())};pools={g:v for g,v in pools.items() if len(v)>0}
 pos=pools[TARGET];neg=np.concatenate([v for g,v in pools.items() if g in [NEG,'Cycling_'+NEG]])
 tumor=[g for g in pools if parent(g)=='Tumor_hepatocyte'];other=[g for g in pools if g not in [TARGET,NEG,'Cycling_'+NEG] and parent(g)!='Tumor_hepatocyte']
 seed=42+sum((name+d).encode());rng=np.random.default_rng(seed);out=[]
 cases=[('standard',f) for f in [0,.005,.01,.02,.05]]
 if challenge:cases.append(('cycling_zero',0));cases+=([('other_myeloid_zero',0)] if 'Other_myeloid' in pools else [])
 for scenario,f in cases:
  for rep in range(replicates):
   chosen=[];weights={TARGET:f,'negative_combined':.25-f}
   if scenario=='standard':tw=rng.dirichlet(np.ones(len(tumor)));ow=rng.dirichlet(np.ones(len(other)))
   elif scenario=='cycling_zero':
    tw=np.array([10. if g.startswith('Cycling_') else 1. for g in tumor]);tw/=tw.sum();ow=np.array([8. if g.startswith('Cycling_') else .3 for g in other]);ow/=ow.sum()
   else:
    tw=np.ones(len(tumor))/len(tumor);ow=np.array([30. if g=='Other_myeloid' else .1 for g in other]);ow/=ow.sum()
   weights.update({g:.50*w for g,w in zip(tumor,tw)});weights.update({g:.25*w for g,w in zip(other,ow)});assert np.isclose(sum(weights.values()),1)
   if not sampled:
    y=np.zeros(len(genes))
    for g,w in weights.items():
     pool=neg if g=='negative_combined' else pools[g]
     if w:y+=w*np.asarray(data['norm'][pool].mean(0)).ravel()
    out.append({'scenario':scenario,'nominal_fraction':f,'replicate':rep,'unit':'equal_RNA_cell_fraction','truth':f,'y':y});continue
   counts={g:int(round(w*2000)) for g,w in weights.items()};delta=2000-sum(counts.values());counts['negative_combined']+=delta
   selected_pos=[];selected_tam=[]
   for g,n in counts.items():
    if n:
     pool=neg if g=='negative_combined' else pools[g];ix=rng.choice(pool,n,replace=True);chosen.extend(ix.tolist())
     if g==TARGET:selected_pos.extend(ix.tolist())
     if g in [TARGET,'negative_combined']:selected_tam.extend(ix.tolist())
   chosen=np.asarray(chosen);assert len(chosen)==2000;truecell=len(selected_pos)/2000
   ye=np.asarray(data['norm'][chosen].mean(0)).ravel();summed=np.asarray(data['raw'][chosen].sum(0)).ravel();total=data['total'][chosen].sum();yr=summed/total*1e4;truerna=data['total'][selected_pos].sum()/total if selected_pos else 0
   out.append({'scenario':scenario,'nominal_fraction':f,'replicate':rep,'unit':'equal_RNA_cell_fraction','truth':truecell,'truth_total_TAM':len(selected_tam)/2000,'y':ye})
   out.append({'scenario':scenario,'nominal_fraction':f,'replicate':rep,'unit':'pooled_count_RNA_fraction','truth':truerna,'truth_total_TAM':float(data['total'][selected_tam].sum()/total),'y':yr})
 return out
def prediction_rows(m,mix,name,d,phase):
 rows=[]
 for z in mix:
  pred,res=estimate(m,z['y']);rows.append({k:v for k,v in z.items() if k!='y'}|{'dataset':name,'donor':d,'phase':phase,'config':m['config'],'with_PGAM5':m['with_pg'],'with_cycling_competitors':m['with_cycle'],'predicted':pred,'residual':res,'features':len(m['ix']),'reference_components':len(m['groups'])})
 return rows
def score(rows):
 t=pd.DataFrame(rows);st=t[t.scenario.eq('standard')];zero=t[t.truth.eq(0)];return float(abs(st.predicted-st.truth).mean()+2*zero.predicted.quantile(.95))
# Inner tests use deterministic holdout class means; outer/external use sampled mixtures and both units.
def run_training_validation():
 meanmix={d:mixtures('GSE151530',d,replicates=2,sampled=False) for d in internal};innerrows=[];outerrows=[];chosenconfigs=[]
 for d in internal:
  train=[k for k in donors if k!=d];candidates={}
  for config in configs:
   rows=[]
   for inner in internal:
    if inner==d:continue
    m=model([k for k in train if k!=inner],config);rr=prediction_rows(m,meanmix[inner],'GSE151530',inner,'inner_tuning');rows.extend(rr)
   candidates[config]=score(rows);innerrows.extend({'outer_donor':d,**r} for r in rows)
  best=min(candidates,key=candidates.get);chosenconfigs.append({'outer_donor':d,'selected_config':best,'scores_json':json.dumps(candidates)});m=model(train,best)
  outerrows.extend(prediction_rows(m,mixtures('GSE151530',d,replicates=5),'GSE151530',d,'outer_nested'))
  print('Outer donor',d,'chosen',best,'features',len(m['ix']),'score',candidates,flush=True)
 pd.DataFrame(innerrows).to_csv(R/'inner_configuration_predictions.csv',index=False);pd.DataFrame(chosenconfigs).to_csv(R/'nested_configuration_choices.csv',index=False)
 # Choose final configuration using discovery-only leave-donor-out means, without examining external outputs.
 finalscores={};cvrows=[]
 for config in configs:
  rows=[]
  for d in internal:
   m=model([k for k in donors if k!=d],config);rows.extend(prediction_rows(m,meanmix[d],'GSE151530',d,'discovery_configuration_selection'))
  finalscores[config]=score(rows);cvrows.extend(rows)
 best=min(finalscores,key=finalscores.get);locked=model(donors,best);print('LOCKED',best,finalscores,flush=True)
 lock={'selected_config':best,'selection_scores':finalscores,'selection':'Discovery-only leave-donor-out means; no external tuning','training_donors':donors,'internal_test_donors':internal,'external_eligible_donors':external,'genes':genes[locked['ix']].tolist(),'reference_groups':locked['groups'],'protocol_sha256':hashlib.sha256((R/'validation_protocol.json').read_bytes()).hexdigest(),'status':'LOCKED_BEFORE_EXTERNAL_RESULTS'}
 (R/'locked_reference_manifest.json').write_text(json.dumps(lock,indent=2));pd.DataFrame(cvrows).to_csv(R/'discovery_configuration_predictions.csv',index=False)
 pd.DataFrame(locked['ref'][locked['ix']],index=genes[locked['ix']],columns=locked['groups']).to_csv(R/'LOCKED_REFERENCE_CP10K_VALIDATION_REQUIRED.tsv',sep='\t');pd.DataFrame(locked['markers']).to_csv(R/'training_selected_gene_evidence.csv',index=False)
 externalrows=[]
 for d in external:
  mix=mixtures('GSE149614',d,replicates=5)
  for with_pg,with_cycle in [(True,True),(False,True),(True,False)]:
   m=locked if with_pg and with_cycle else model(donors,best,with_pg,with_cycle)
   externalrows.extend(prediction_rows(m,mix,'GSE149614',d,'external_locked' if with_pg and with_cycle else 'external_sensitivity'))
  print('External donor',d,'completed',flush=True)
 rows=pd.DataFrame(outerrows+externalrows);rows.to_csv(R/'validation_predictions.csv',index=False)
 summary=[]
 for keys,t in rows.groupby(['dataset','phase','with_PGAM5','with_cycling_competitors','unit','scenario']):
  name,phase,withpg,withcycle,unit,scenario=keys;e=t.predicted-t.truth;nonzero=t[t.truth.gt(0)];zero=t[t.truth.eq(0)]
  corr=float(pearsonr(t.truth,t.predicted).statistic) if t.truth.nunique()>1 and t.predicted.nunique()>1 else None
  summary.append({'dataset':name,'phase':phase,'with_PGAM5':withpg,'with_cycling_competitors':withcycle,'unit':unit,'scenario':scenario,'donors':t.donor.nunique(),'mixtures':len(t),'MAE_percentage_points':float(abs(e).mean()*100),'Pearson_r':corr,'zero_target_p95_percent':float(zero.predicted.quantile(.95)*100) if len(zero) else None,'median_predicted_over_true_at_nonzero':float((nonzero.predicted/nonzero.truth).median()) if len(nonzero) else None})
 su=pd.DataFrame(summary);su.to_csv(R/'validation_summary.csv',index=False)
 checks=[]
 for phase in ['outer_nested','external_locked']:
  t=su[su.phase.eq(phase)]
  for unit in ['equal_RNA_cell_fraction','pooled_count_RNA_fraction']:
   st=t[t.unit.eq(unit)&t.scenario.eq('standard')].iloc[0]
   checks.extend([{'phase':phase,'criterion':unit+'_MAE_pp<=1','value':st.MAE_percentage_points,'passed':bool(st.MAE_percentage_points<=1)},{'phase':phase,'criterion':unit+'_Pearson>=0.7','value':st.Pearson_r,'passed':bool(pd.notna(st.Pearson_r) and st.Pearson_r>=.7)},{'phase':phase,'criterion':unit+'_zero_p95_percent<=0.5','value':st.zero_target_p95_percent,'passed':bool(st.zero_target_p95_percent<=.5)}])
   for scenario in ['cycling_zero','other_myeloid_zero']:
    z=t[t.unit.eq(unit)&t.scenario.eq(scenario)]
    if len(z):checks.append({'phase':phase,'criterion':unit+'_'+scenario+'_p95_percent<=0.5','value':float(z.iloc[0].zero_target_p95_percent),'passed':bool(z.iloc[0].zero_target_p95_percent<=.5)})
 checks.append({'phase':'external_locked','criterion':'external_donors>=5','value':len(external),'passed':len(external)>=5});pd.DataFrame(checks).to_csv(R/'predeclared_gate_checks.csv',index=False)
 result={'validation_passed':all(c['passed'] for c in checks),'selected_config':best,'internal_donors':len(internal),'external_donors':len(external),'selected_genes':len(locked['ix']),'reference_components':len(locked['groups']),'failed_checks':[c for c in checks if not c['passed']],'TCGA_status':'NOT_RUN: validation gate determines eligibility; additional reference/bulk platform matching and abundance-unit calibration are still prerequisites even if gate passes','limits':'NNLS only; internal means select config, independent outer mixtures evaluate nested selection. Bootstrap repeats are not independent biological donors. Pooled-count RNA weights use observed library size and are not calibrated biological RNA content. GSE149614 Other_myeloid uncovered; Hepatocyte coarse label not individually verified malignant.'}
 (R/'validation_result.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':
 run_training_validation()
