"""Independently reconstruct radius-2 neighborhoods and all proxy/count statistics."""
import json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy.spatial import cKDTree
from scipy.stats import false_discovery_control
from analyze import ROOT,SAMPLES,TIERS,MARKERS,dense
checks=[]
def check(name,v):
 checks.append(dict(check=name,passed=bool(v)))
 if not v:raise AssertionError(name)
def run():
 stat=pd.read_csv(ROOT/'spatial_statistics.csv');summ=pd.read_csv(ROOT/'section_candidate_summary.csv');reviews=[]
 for row in json.loads((ROOT/'input_manifest.json').read_text()):check(row['path']+' SHA',hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'])
 for item in SAMPLES:
  s=item['sample'];a=ad.read_h5ad(ROOT/'inputs'/f'{s}_spatial_marker_subset.h5ad');o=pd.read_csv(ROOT/f'{s}_TLS_proxy_annotations.csv.gz',index_col=0).loc[a.obs_names];raw=dense(a.layers['counts']);norm=np.log1p(raw/a.obs.total_counts.to_numpy()[:,None]*1e4);check(s+' normalized inputs',np.allclose(norm,dense(a.X),rtol=3e-6,atol=3e-6))
  xy=a.obs[['array_col','array_row']].to_numpy(float).copy();xy[:,1]*=np.sqrt(3);edges=cKDTree(xy).query_pairs(2.000001);neighbors={i:set() for i in range(len(a))}
  for i,j in edges:neighbors[i].add(j);neighbors[j].add(i)
  b=o.gate_B.to_numpy(bool);t=o.gate_T.to_numpy(bool);bo=b&~t;to=t&~b;cx=raw[:,a.var_names.get_loc('CXCL13')]>0;mh=(raw[:,a.var_names.get_indexer(MARKERS['MHCII'])]>0).sum(1)>=3;high=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy()
  pg=raw[:,a.var_names.get_loc('PGAM5')];check(s+' original PGAM5 raw counts',np.array_equal(pg,o.PGAM5_raw));pos=o.Mac_RNA_enriched.to_numpy(bool)&(pg>0);expected=pos&o.PGAM5_log1p_CP10k.ge(o.loc[pos,'PGAM5_log1p_CP10k'].quantile(.75)).to_numpy();check(s+' fixed PGAM5 Q75 group',np.array_equal(expected,high))
  bs=norm[:,a.var_names.get_indexer(MARKERS['B'])].mean(1);bd=(raw[:,a.var_names.get_indexer(MARKERS['B'])]>0).sum(1);bg=(bs>=np.quantile(bs,.75))&(bs>0)&(bd>=2)&((raw[:,a.var_names.get_loc('MS4A1')]>0)|(raw[:,a.var_names.get_loc('CD19')]>0));check(s+' B gate',np.array_equal(bg,b))
  cand={k:np.zeros(len(a),bool) for k in TIERS};regions={k:np.zeros(len(a),bool) for k in TIERS}
  for i in np.flatnonzero(b):
   patch={i}|neighbors[i]
   for j in neighbors[i]:patch|=neighbors[j]
   nodes=np.array(sorted(patch));contact=any(bo[j] and any(to[z] and z in patch for z in neighbors[j]) for j in patch);basic=b[nodes].sum()>=2 and t[nodes].sum()>=3 and (b[nodes]|t[nodes]).sum()>=5 and bo[nodes].sum()>=1 and to[nodes].sum()>=1 and contact;support=cx[nodes].sum()>=1 and mh[nodes].sum()>=1;strict=bo[nodes].sum()>=2 and to[nodes].sum()>=3 and cx[nodes].sum()>=2 and mh[nodes].sum()>=1
   for k,v in zip(TIERS,[basic,basic and support,basic and strict]):
    cand[k][i]=v
    if v:regions[k][nodes]=True
  for k in TIERS:
   adj=np.array([not regions[k][i] and any(regions[k][j] for j in neighbors[i]) for i in range(len(a))]);check(s+k+' center reconstruction',np.array_equal(cand[k],o['center_'+k]));check(s+k+' union reconstruction',np.array_equal(regions[k],o['region_'+k]));check(s+k+' exterior adjacency',np.array_equal(adj,o['adjacent_'+k]));row=summ[(summ['sample']==s)&summ.tier.eq(k)].iloc[0];check(s+k+' summary',row.centers==cand[k].sum() and row.high_within==(high&regions[k]).sum() and row.high_adjacent_outside==(high&adj).sum())
  for i in np.flatnonzero(high&regions[TIERS[1]]):reviews.append(dict(sample=s,barcode=a.obs_names[i],PGAM5_raw=int(pg[i]),T_gate=bool(t[i]),B_gate=bool(b[i])))
  for row in stat[stat['sample']==s].itertuples():
   use=o.Mac_RNA_enriched.to_numpy(bool)&o.spatial_degree.gt(0).to_numpy()
   if row.comparison!='high_vs_all_other_mac':use&=o.PGAM5_raw.gt(0).to_numpy()
   if row.comparison.startswith('nonTB'):use&=~b&~t
   if row.endpoint=='adjacent_outside_patch':use&=~regions[row.tier]
   outcome=regions[row.tier] if row.endpoint=='within_patch' else o['adjacent_'+row.tier].to_numpy(bool);check(s+row.tier+row.endpoint+row.comparison+' counts',row.n_high==(use&high).sum() and row.n_comparator==(use&~high).sum() and row.high_supported==(use&high&outcome).sum() and row.comparator_supported==(use&~high&outcome).sum())
   p=ROOT/f'{s}_{row.tier}_{row.endpoint}_{row.comparison}_audit.npz'
   if not p.exists():continue
   v=np.load(p);h=v['high'];y=v['outcome'];obs=100*(h@y/h.sum()-(1-h)@y/(1-h).sum());check(p.name+' observed',np.isclose(obs,v['observed']))
   for g in np.unique(v['groups']):check(p.name+str(g)+' preserved high counts',np.all(v['assignments_first10'][:,v['groups']==g].sum(1)==h[v['groups']==g].sum()))
   for j in range(3):
    hp=v['assignments_first10'][j];val=100*(hp@y/hp.sum()-(1-hp)@y/(1-hp).sum());check(p.name+str(j)+' permutation statistic',np.isclose(val,v['null'][j]))
   if row.status=='tested':
    null=v['null'];pv=min(1,2*min((1+(null>=obs-1e-10).sum())/5000,(1+(null<=obs+1e-10).sum())/5000));check(p.name+' p',np.isclose(pv,row.p))
 check('288 planned tests',len(stat)==288);check('joint BH',np.allclose(stat.joint_q_288,false_discovery_control(stat.p.fillna(1).to_numpy())));ix=stat.tier.eq(TIERS[1])&stat.comparison.eq('high_vs_other_detected');check('primary 32 BH',ix.sum()==32 and np.allclose(stat.loc[ix,'primary_q_32'],false_discovery_control(stat.loc[ix,'p'].fillna(1).to_numpy())))
 check('three supported-patch high sources',len(reviews)==3);check('all three have T/B RNA gate',all(r['T_gate'] or r['B_gate'] for r in reviews))
 (ROOT/'mixed_source_review.json').write_text(json.dumps(dict(status='passed',raw_PGAM5_checked_against_input_layers=True,within_supported_patch_high_sources=len(reviews),all_sources_T_or_B_gate_positive=True,nonTB_high_within_supported_patch=0,sources=reviews),indent=2)+'\n')
 (ROOT/'audit_results.json').write_text(json.dumps(dict(checks=len(checks),passed=sum(c['passed'] for c in checks),details=checks),indent=2)+'\n');print('PASS',len(checks))
if __name__=='__main__':run()
