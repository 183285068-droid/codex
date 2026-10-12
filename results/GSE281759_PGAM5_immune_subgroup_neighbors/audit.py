"""Independent graph, source gate, normalization, denominator and permutation checks."""
import json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy.spatial import cKDTree
from scipy.stats import false_discovery_control
from analyze import ROOT,REF,SAMPLES,TARGETS,PANELS,dense
checks=[]
def check(name,value):
 checks.append({'check':name,'passed':bool(value)})
 if not value:raise AssertionError(name)
def run():
 d=pd.read_csv(ROOT/'neighborhood_counts.csv');s=pd.read_csv(ROOT/'spatial_statistics.csv');lr=pd.read_csv(ROOT/'candidate_LR_codetection.csv');pairs=pd.read_csv(ROOT/'candidate_LR_resource.csv')
 for f in json.loads((ROOT/'reference_manifest.json').read_text()):check(f['path']+' SHA256',hashlib.sha256((ROOT.parent/f['path']).read_bytes()).hexdigest()==f['sha256'])
 for item in SAMPLES:
  sample=item['sample'];a=ad.read_h5ad(REF/f'{sample}_spatial_marker_subset.h5ad');o=pd.read_csv(ROOT/f'{sample}_annotations.csv.gz',index_col=0).loc[a.obs_names];raw=dense(a.layers['counts']);norm=np.log1p(raw/a.obs.total_counts.to_numpy()[:,None]*1e4)
  check(sample+' full-library normalization',np.allclose(norm,dense(a.X),rtol=3e-6,atol=3e-6));check(sample+' unique barcode',a.obs_names.is_unique)
  xy=a.obs[['array_row','array_col']].to_numpy(float).copy();xy[:,0]*=np.sqrt(3);ij=cKDTree(xy).query_pairs(2.000001,output_type='ndarray');expected=set(map(tuple,ij))|set(map(tuple,ij[:,::-1]));coo=a.obsp['spatial_connectivities'].tocoo();check(sample+' independent lattice',set(zip(coo.row,coo.col))==expected)
  neighbors={i:[] for i in range(len(a))}
  for i,j in expected:neighbors[i].append(j)
  m=(o.Macrophage_core>=o.Macrophage_core.quantile(.75))&(o.Macrophage_core_n_detected>=3)&(o.Pan_myeloid>=o.Pan_myeloid.median());pos=m&o.PGAM5_raw.gt(0);high=pos&o.PGAM5_log1p_CP10k.ge(o.loc[pos,'PGAM5_log1p_CP10k'].quantile(.75));check(sample+' source gate',np.array_equal(high,o.PGAM5_mac_group.eq('PGAM5_high')))
  for k,gs in PANELS.items():
   ix=a.var_names.get_indexer(gs);check(sample+' '+k+' score',np.allclose(norm[:,ix].mean(1),o['score_'+k],rtol=3e-6,atol=3e-6))
   target=o['gate_'+k].to_numpy(bool);y=np.array([any(target[j] for j in neighbors[i]) for i in range(len(a))]);check(sample+' '+k+' explicit neighbor',np.array_equal(y,o['neighbor_'+k]))
   for group in ['PGAM5_high','Other_detected','Undetected']:
    mask=m&o.PGAM5_mac_group.eq(group)&o.spatial_degree.gt(0);row=d[(d['sample']==sample)&d.target.eq(k)&d.source_group.eq(group)].iloc[0];check(sample+k+group+' denominator',row.n_sources==mask.sum() and row.with_neighbor==y[mask].sum())
   for pair in pairs.itertuples():
    edges=[(i,j) for i,j in expected if high.iloc[i] and target[j]];lig=a.var_names.get_loc(pair.ligand);ri=a.var_names.get_indexer(pair.receptor_complex.split('_'));supported=sum(raw[i,lig]>0 and (raw[j,ri]>0).all() for i,j in edges);row=lr[(lr['sample']==sample)&lr.target.eq(k)&lr.pair.eq(pair.ligand+'--'+pair.receptor_complex)].iloc[0];check(sample+k+pair.ligand+pair.receptor_complex+' LR denominator',row.eligible_edges==len(edges) and row.supported_edges==supported)
 for path in ROOT.glob('*_audit.npz'):
  v=np.load(path);h=v['high'];y=v['outcomes'];null=v['null'];obs=100*(h@y/h.sum()-(1-h)@y/(1-h).sum());check(path.name+' observed',np.allclose(obs,v['observed']))
  for j,k in enumerate(TARGETS):
   sample=path.name.split('_')[0];kind=path.name[len(sample)+1:-len('_audit.npz')];row=s[(s['sample']==sample)&s.comparison.eq(kind)&s.target.eq(k)].iloc[0]
   if row.status=='tested':
    p=min(1,2*min((1+(null[:,j]>=obs[j]-1e-10).sum())/5000,(1+(null[:,j]<=obs[j]+1e-10).sum())/5000));check(path.name+k+' equal-tail p',np.isclose(p,row.p))
 check('432 planned rows',len(s)==432);check('joint BH',np.allclose(false_discovery_control(s.p.fillna(1).to_numpy()),s.joint_q_432));ix=s.comparison.eq('high_vs_other_detected');check('primary BH',np.allclose(false_discovery_control(s.loc[ix,'p'].fillna(1).to_numpy()),s.loc[ix,'primary_q_144']))
 (ROOT/'audit_results.json').write_text(json.dumps({'checks':len(checks),'passed':sum(x['passed'] for x in checks),'details':checks},indent=2)+'\n');print('PASS',len(checks))
if __name__=='__main__':run()
