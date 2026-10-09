from pathlib import Path
import json,importlib.util
import numpy as np,pandas as pd,anndata as ad
from scipy.optimize import nnls,lsq_linear
R=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('reference_workflow',R/'train_joint_reference.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
lock=json.loads((R/'locked_reference_manifest.json').read_text());m=mod.model(mod.donors,lock['selected_config']);ref=pd.read_csv(R/'LOCKED_REFERENCE_CP10K_VALIDATION_REQUIRED.tsv',sep='\t',index_col=0)
assert ref.index.tolist()==mod.genes[m['ix']].tolist() and ref.columns.tolist()==m['groups'] and np.allclose(ref,m['ref'][m['ix']])
# Exactly matched reference mixtures must be recovered, before interpreting biological holdouts.
rng=np.random.default_rng(20261009);impl=[]
for i in range(20):
 w=rng.dirichlet(np.ones(len(m['groups'])));y=m['matrix']@w;beta,res=nnls(m['matrix'],y);err=float(np.abs(beta-w).max());assert err<1e-6
 impl.append({'test':'exact_reference_mixture','case':i,'max_coefficient_error':err,'residual':res})
coefrows=[];checks=[];pred=pd.read_csv(R/'validation_predictions.csv');baseline=pred[pred.phase.eq('external_locked')]
for donor in mod.external:
 for z in mod.mixtures('GSE189903',donor,replicates=5):
  y=z['y'][m['ix']]/m['scale'];beta,res=nnls(m['matrix'],y);weights=beta/max(beta.sum(),1e-12);value=weights[m['groups'].index(mod.TARGET)]
  existing=baseline[(baseline.donor==donor)&(baseline.unit==z['unit'])&(baseline.scenario==z['scenario'])&(baseline.replicate==z['replicate'])&np.isclose(baseline.nominal_fraction,z['nominal_fraction'])]
  assert len(existing)==1 and np.isclose(value,existing.predicted.iloc[0],atol=1e-12)
  macidx=[i for i,g in enumerate(m['groups']) if g in [mod.TARGET,mod.NEG,'Cycling_'+mod.NEG]]
  row={k:v for k,v in z.items() if k!='y'}|{'donor':donor,'predicted_target':value,'predicted_total_TAM':float(weights[macidx].sum())};row.update({'coefficient_'+g:float(w) for g,w in zip(m['groups'],weights)});coefrows.append(row)
  if z['scenario']=='standard' and z['replicate']==0:
   # Independent bounded least-squares solver checks the NNLS objective and target coefficient.
   alt=lsq_linear(m['matrix'],y,bounds=(0,np.inf),method='bvls',tol=1e-10,max_iter=500);assert alt.success
   ob1=np.linalg.norm(m['matrix']@beta-y);ob2=np.linalg.norm(m['matrix']@alt.x-y);assert np.isclose(ob1,ob2,atol=1e-6)
   ap=alt.x[m['groups'].index(mod.TARGET)]/max(alt.x.sum(),1e-12);assert np.isclose(ap,value,atol=1e-5)
   checks.append({'donor':donor,'unit':z['unit'],'true_fraction':z['truth'],'NNLS_prediction':value,'bounded_solver_prediction':float(ap),'objective_NNLS':float(ob1),'objective_bounded_solver':float(ob2)})
pd.DataFrame(impl).to_csv(R/'exact_reference_recovery_checks.csv',index=False);pd.DataFrame(checks).to_csv(R/'independent_bounded_solver_checks.csv',index=False);co=pd.DataFrame(coefrows);co.to_csv(R/'external_full_component_predictions.csv',index=False)
sumrows=[]
for (unit,scenario),t in co.groupby(['unit','scenario']):
 sumrows.append({'unit':unit,'scenario':scenario,'total_TAM_MAE_percentage_points':float(abs(t.predicted_total_TAM-t.truth_total_TAM).mean()*100),'total_TAM_bias_percentage_points':float((t.predicted_total_TAM-t.truth_total_TAM).mean()*100),'target_MAE_percentage_points':float(abs(t.predicted_target-t.truth).mean()*100)})
pd.DataFrame(sumrows).to_csv(R/'external_total_TAM_control_summary.csv',index=False)

result={'exact_reference_mixtures_checked':len(impl),'max_exact_coefficient_error':max(r['max_coefficient_error'] for r in impl),'independent_bounded_solver_cases':len(checks),'external_predictions_exactly_reproduced':len(co),'status':'Implementation checks passed; scientific validation interpreted separately'}; (R/'implementation_audit_summary.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

# Pure competitor-only pools can be tested even where target-positive counts are insufficient for dose-response validation.
purerows=[]
for name in ['joint_discovery','GSE189903']:
 o=mod.D[name]['obs']
 for d in sorted(o.donor.unique()):
  for g in ['Monocyte_enriched','DC_enriched','pDC','Mixed_APC']:
   sel=o.donor.eq(d)&o.fine_group.eq(g)
   if sel.sum()<20:continue
   mm=m if name=='GSE189903' else mod.model([k for k in mod.donors if k!=d],lock['selected_config'])
   y=np.asarray(mod.D[name]['norm'][sel.to_numpy()].mean(0)).ravel();b,_=nnls(mm['matrix'],y[mm['ix']]/mm['scale']);w=b/max(b.sum(),1e-12)
   purerows.append({'dataset':name,'donor':d,'pure_competitor':g,'cells':int(sel.sum()),'true_target_fraction':0,'estimated_target_fraction':float(w[mm['groups'].index(mod.TARGET)]),'true_total_macrophage_fraction':0,'estimated_total_macrophage_fraction':float(sum(w[i] for i,g in enumerate(mm['groups']) if g in [mod.TARGET,mod.NEG,'Cycling_'+mod.NEG])),'interpretation':'Supplementary pure-cell mean stress test; not a replacement for predeclared mixed-tumor gates'})
pd.DataFrame(purerows).to_csv(R/'pure_competitor_stress_tests.csv',index=False)
