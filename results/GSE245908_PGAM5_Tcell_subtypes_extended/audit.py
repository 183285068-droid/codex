"""Recompute classifications, adjacency, LR summaries and permutation p/q."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
from analyze import R,OLD,PANELS,TYPES,STATES,GROUPS,f
checks=[]
def check(name,ok):
 checks.append(dict(check=name,passed=bool(ok)));assert ok,name
stats=pd.read_csv(R/'spatial_subtype_tests.csv');summary=pd.read_csv(R/'subtype_neighbor_summary.csv');lr=pd.read_csv(R/'LR_by_subtype.csv');thresholds=pd.read_csv(R/'state_thresholds.csv')
for record in json.loads((R/'input_manifest.json').read_text()):
 p=R.parent/record['path'];check('SHA256 '+record['path'],hashlib.sha256(p.read_bytes()).hexdigest()==record['sha256'])
for s in f.SAMPLES:
 a=f.load_sample(s);o=pd.read_csv(R/f'{s}_spot_subtypes.csv',index_col=0);old=pd.read_csv(OLD/f'{s}_spot_states.csv',index_col=0);W=f.graph(o);n=len(o);degree=W.sum(1).A1
 check(s+' unchanged sources',o.source_group.equals(old.source_group));check(s+' unchanged T gate',o.T_enriched.equals(old.T_enriched))
 expected=set();lookup={(int(row),int(col)):i for i,(row,col) in enumerate(o[['array_row','array_col']].to_numpy())}
 for (row,col),i in lookup.items():
  for dr,dc in [(0,2),(0,-2),(1,1),(1,-1),(-1,1),(-1,-1)]:
   j=lookup.get((row+dr,col+dc))
   if j is not None:expected.add((i,j))
 check(s+' hexagonal geometry',set(zip(*W.nonzero()))==expected);check(s+' no self edges',W.diagonal().sum()==0)
 raw=a.layers['counts'];check(s+' full-library denominator',np.allclose(np.asarray(raw.sum(1)).ravel(),o.total_counts));check(s+' nonnegative raw counts',(raw.data>=0).all())
 for k,gs in PANELS.items():
  x=a[:,gs].X.toarray();check(s+' score '+k,np.allclose(x.mean(1),o['score_'+k],atol=1e-7));check(s+' detection '+k,np.array_equal((x>0).sum(1),o['detected_'+k]))
 t=o.T_enriched.to_numpy();cd8=t&(a[:,['CD8A']].X.toarray().ravel()>0)&(o.score_CD8_identity.to_numpy()>o.score_CD4_identity.to_numpy());cd4=t&(o.score_CD4_identity.to_numpy()>0)&(o.score_CD4_identity.to_numpy()>o.score_CD8_identity.to_numpy());reg=cd4&(a[:,['FOXP3']].X.toarray().ravel()>0)&(o.detected_Treg.to_numpy()>=2)
 check(s+' CD8 gate',np.array_equal(cd8,o.gate_All_CD8));check(s+' Treg gate',np.array_equal(reg,o.gate_Treg));check(s+' identity partition',np.array_equal(o[['gate_All_CD8','gate_CD4_other','gate_Treg','gate_T_unresolved']].sum(1).to_numpy(),t.astype(int)))
 for k,p in STATES.items():
  cut=o.loc[cd8,'score_'+p].quantile(.75);mask=cd8&(o['score_'+p]>=cut)&(o['score_'+p]>0)&(o['detected_'+p]>=2);check(s+' subtype '+k,np.array_equal(mask,o['gate_'+k]))
 for k in TYPES:
  nt=W@o['gate_'+k].to_numpy(float);check(s+' neighbor fraction '+k,np.allclose(np.divide(nt,degree,out=np.zeros(n),where=degree>0),o['neighbor_fraction_'+k]))
  for g in GROUPS:
   keep=o.source_group.eq(g).to_numpy();r=summary[(summary['sample']==s)&(summary.subtype==k)&(summary.source_group==g)].iloc[0];check(s+g+k+' edges',int(nt[keep].sum())==r.directed_edges);check(s+g+k+' mean',np.isclose((nt[keep]/degree[keep]).mean(),r.mean_neighbor_fraction))
 edge=pd.read_csv(R/f'{s}_subtype_LR_edges.csv.gz');ix={b:i for i,b in enumerate(o.index)};ri=np.array([ix[b] for b in edge.source_barcode]);ci=np.array([ix[b] for b in edge.target_barcode]);check(s+' expected edge selection',set(zip(ri,ci))=={(i,j) for i,j in expected if o.Mac_enriched.iloc[i] and o.T_enriched.iloc[j]})
 for l,receptor,_ in f.PAIRS:
  ligand=a[:,[l]].layers['counts'].toarray().ravel()/o.total_counts.to_numpy()*1e4;rec=np.min(a[:,receptor.split('_')].layers['counts'].toarray()/o.total_counts.to_numpy()[:,None]*1e4,axis=1);val=np.sqrt(ligand[ri]*rec[ci]);check(s+l+receptor+' LR raw expression',np.allclose(val,edge[l+'--'+receptor]))
  for k in TYPES:
   for g in GROUPS:
    keep=(edge.source_group==g)&o['gate_'+k].to_numpy()[ci];v=val[keep];r=lr[(lr['sample']==s)&(lr.subtype==k)&(lr.source_group==g)&(lr.ligand==l)&(lr.receptor==receptor)].iloc[0];check(s+g+k+l+receptor+' LR aggregation',len(v)==r.n_edges and int((v>0).sum())==r.n_codetected and (np.isnan(r.mean_edge_compatibility) if not len(v) else np.isclose(v.mean(),r.mean_edge_compatibility)))
 for kind in ['subtype_coverage','CD8_neighbor_programs']:
  path=R/f'{s}_{kind}_permutation_audit.npz'
  if not path.exists():continue
  with np.load(path) as p:
   ix=p['source_indices'];h=p['high'];Y=p['outcomes'];obs=h@Y/h.sum()-(1-h)@Y/(len(h)-h.sum());check(s+kind+' observed delta',np.allclose(obs,p['observed_deltas']))
   for j,label in enumerate(p['labels']):
    row=stats[(stats['sample']==s)&(stats.kind==kind)&(stats.feature==label)].iloc[0];ex=(np.abs(p['null_deltas'][:,j]-p['conditional_null_center'][j])>=abs(obs[j]-p['conditional_null_center'][j])-1e-14).sum();pp=ex/len(p['null_deltas']) if str(p['mode'])=='exact_conditional_enumeration' else (1+ex)/(len(p['null_deltas'])+1);check(s+str(label)+' conditional p',np.isclose(pp,row.spatial_conditional_p))
check('24 prespecified rows',len(stats)==24);ok=stats.spatial_conditional_p.notna();check('joint BH',np.allclose(multipletests(stats.loc[ok,'spatial_conditional_p'],method='fdr_bh')[1],stats.loc[ok,'BH_q']))
pd.DataFrame(checks).to_csv(R/'audit_checks.csv',index=False);print('PASS',len(checks),'checks')
