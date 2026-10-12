import json,hashlib,itertools
import numpy as np,pandas as pd,anndata as ad
from scipy.spatial import cKDTree
from scipy.stats import false_discovery_control,ttest_rel
from analyze import ROOT,SAMPLES,PANELS,COMPS,PATIENTS,dense,nuisance
checks=[]
def check(n,v):
 checks.append(dict(check=n,passed=bool(v)))
 if not v:raise AssertionError(n)
def run():
 sec=pd.read_csv(ROOT/'section_effects_descriptive.csv');pat=pd.read_csv(ROOT/'patient_paired_effects.csv');tests=pd.read_csv(ROOT/'patient_level_statistics.csv')
 for row in json.loads((ROOT/'input_manifest.json').read_text()):check(row['path']+' SHA',hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'])
 for item in SAMPLES:
  s=item['sample'];a=ad.read_h5ad(ROOT/'inputs'/f'{s}_plasma_marker_inputs.h5ad');q=pd.read_csv(ROOT/f'{s}_target_profiles.csv.gz',index_col=0);pi=a.obs_names.get_indexer(q.index);raw=dense(a.layers['counts']);norm=np.log1p(raw/a.obs.total_counts.to_numpy()[:,None]*1e4);check(s+' full-library normalization',np.allclose(norm,dense(a.X),rtol=3e-6,atol=3e-6))
  xy=a.obs[['array_col','array_row']].to_numpy(float).copy();xy[:,1]*=np.sqrt(3);pairs=cKDTree(xy).query_pairs(2.000001);neighbors={i:set() for i in range(len(a))}
  for i,j in pairs:neighbors[i].add(j);neighbors[j].add(i)
  high=set(np.flatnonzero(a.obs.PGAM5_mac_group.eq('PGAM5_high')));other=set(np.flatnonzero(a.obs.PGAM5_mac_group.eq('Other_detected')));nh=np.array([len(neighbors[i]&high) for i in pi]);no=np.array([len(neighbors[i]&other) for i in pi]);check(s+' direct target neighbors',np.array_equal(nh,q.n_high_neighbors) and np.array_equal(no,q.n_other_neighbors))
  strict=(raw[pi][:,a.var_names.get_indexer(['MZB1','JCHAIN'])]>0).all(1);check(s+' strict state',np.array_equal(strict,q.strict_plasma));v=np.load(ROOT/f'{s}_residual_audit.npz');C,Q,r=nuisance(a,pi) if len(pi) else (None,None,0);group=np.where(nh>0,np.where(no>0,'Shared','High_only'),np.where(no>0,'Other_only','No_detected_mac_neighbor'));check(s+' target group labels',np.array_equal(group,q.group));reconstructed={}
  for j,(k,gs) in enumerate(PANELS.items()):
   y=dense(a[pi,gs].X).mean(1);check(s+k+' program',np.allclose(y,q[k]))
   res=y-C@np.linalg.lstsq(C,y,rcond=1e-10)[0] if len(pi)-r>=15 else np.full_like(y,np.nan)
   if len(pi)-r>=15:check(s+k+' independent OLS residual',np.allclose(res,v['residuals'][:,j],atol=1e-7))
   reconstructed[(k,'raw')]=y;reconstructed[(k,'adjusted')]=res
  for row in sec[sec['sample']==s].itertuples():
   if row.comparison=='high_any_vs_no_high':left=nh>0;right=nh==0
   else:left=(nh>0)&(no==0);right=(nh==0)&(no>0)
   if row.comparison.startswith('strict'):left&=strict;right&=strict
   z=reconstructed[(row.program,row.method)];check(s+row.comparison+row.method+row.program+' independent section counts',row.n_high==left.sum() and row.n_other==right.sum())
   if row.status=='eligible':check(s+row.comparison+row.method+row.program+' independent section effect',np.isclose(row.delta,z[left].mean()-z[right].mean()))
 for row in pat.itertuples():
  subset=sec[sec.patient.eq(row.patient)&sec.program.eq(row.program)&sec.comparison.eq(row.comparison)&sec.method.eq(row.method)&sec.status.eq('eligible')]
  if row.sensitivity=='primary_only':subset=subset[subset.tumor_status.eq('primary')]
  check(row.patient+row.program+row.comparison+row.method+row.sensitivity+' equal-section aggregation',len(subset)==row.n_eligible_sections and (np.isclose(subset.delta.mean(),row.patient_delta) if len(subset) else pd.isna(row.patient_delta)))
 for row in tests.itertuples():
  p=pat[pat.program.eq(row.program)&pat.comparison.eq(row.comparison)&pat.method.eq(row.method)&pat.sensitivity.eq(row.sensitivity)];p=p[p.patient_delta.notna()];v=p.patient_delta.to_numpy();check(str(row.Index)+' patient count',len(v)==row.n_patients)
  if len(v)==4:
   null=np.array(list(itertools.product([-1,1],repeat=4)))@v/4;pv=(np.abs(null)>=abs(v.mean())-1e-12).mean();check(str(row.Index)+' 16 sign assignments',np.isclose(pv,row.exact_signflip_p));tp=ttest_rel(p.high_section_equal_mean,p.other_section_equal_mean).pvalue;check(str(row.Index)+' paired t sensitivity',np.isclose(tp,row.parametric_paired_t_p))
 for key,ix in tests.groupby(['sensitivity','comparison','method']).groups.items():check(str(key)+' primary-family BH',np.allclose(tests.loc[ix,'family_BH_q_8'],false_discovery_control(tests.loc[ix,'exact_signflip_p'].fillna(1).to_numpy())))
 check('96 planned donor tests',len(tests)==96);check('4 independent patient identifiers',set(pat.patient)==set(PATIENTS))
 (ROOT/'audit_results.json').write_text(json.dumps(dict(checks=len(checks),passed=sum(c['passed'] for c in checks),details=checks),indent=2)+'\n');print('PASS',len(checks))
if __name__=='__main__':run()
