"""Original HDF5 gene/library checks and independent graph/permutation recomputation."""
import json,hashlib
import numpy as np,pandas as pd,anndata as ad,h5py
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import false_discovery_control
from analyze import ROOT,RAW,SAMPLES,PANELS,PROGRAMS,COMPS,ISOTYPES,dense,nuisance
checks=[]
def check(name,v):
 checks.append(dict(check=name,passed=bool(v)))
 if not v:raise AssertionError(name)
def run():
 stat=pd.read_csv(ROOT/'program_statistics.csv');iso=pd.read_csv(ROOT/'immunoglobulin_isotype_counts.csv');desc=pd.read_csv(ROOT/'program_group_profiles.csv')
 for r in json.loads((ROOT/'input_manifest.json').read_text()):check(r['path']+' SHA',hashlib.sha256((ROOT/r['path']).read_bytes()).hexdigest()==r['sha256'])
 for sample,gsm in SAMPLES.items():
  a=ad.read_h5ad(ROOT/'inputs'/f'{sample}_functional_marker_inputs.h5ad');o=a.obs;q=pd.read_csv(ROOT/f'{sample}_plasma_target_scores.csv.gz',index_col=0);pi=a.obs_names.get_indexer(q.index);sm=o.source_group.isin(['PGAM5_high','Other_detected']).to_numpy();si=np.flatnonzero(sm);h=o.source_group.eq('PGAM5_high').to_numpy()[si].astype(float);G=a.obsp['spatial_connectivities'];E=G[pi][:,si]
  path=RAW/'inputs'/sample/f'{gsm}_{sample}_filtered_feature_bc_matrix.h5'
  if path.exists():
   with h5py.File(path) as f:
    m=f['matrix'];raw=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:]));ids=[x.decode() for x in m['features/id'][:]];bc=[x.decode() for x in m['barcodes'][:]]
   ci=np.array([bc.index(b) for b in a.obs_names]);gi=np.array([ids.index(g) for g in a.var.gene_ids]);counts=raw[gi][:,ci].T.toarray();total=np.asarray(raw[:,ci].sum(0)).ravel();check(sample+' original 36601 feature library',raw.shape[0]==36601 and np.array_equal(total,o.total_counts));check(sample+' expanded raw gene counts',np.array_equal(counts,dense(a.layers['counts'])));check(sample+' full-library normalized values',np.allclose(np.log1p(counts*1e4/total[:,None]),dense(a.X)))
  xy=o[['array_col','array_row']].to_numpy(float).copy();xy[:,1]*=np.sqrt(3);pairs=cKDTree(xy).query_pairs(2.000001);edges=pairs|{(j,i) for i,j in pairs};coo=G.tocoo();check(sample+' independent hex graph',set(zip(coo.row,coo.col))==edges)
  neighbors={i:set() for i in range(len(a))}
  for i,j in edges:neighbors[i].add(j)
  highset=set(si[h>0]);otherset=set(si[h==0]);nh=np.array([len(neighbors[i]&highset) for i in pi]);no=np.array([len(neighbors[i]&otherset) for i in pi]);check(sample+' explicit high neighbors',np.array_equal(nh,q.n_high_neighbors));check(sample+' explicit other neighbors',np.array_equal(no,q.n_other_detected_neighbors));strict=(dense(a[pi,['MZB1','JCHAIN','DERL3']].layers['counts'])>0).all(1);check(sample+' strict plasma',np.array_equal(strict,q.strict_plasma))
  expected=np.where(nh>0,np.where(no>0,'Shared_high_other','High_only'),np.where(no>0,'Other_detected_only','No_detected_mac_neighbor'));check(sample+' four target groups',np.array_equal(expected,q.neighbor_group))
  for k,gs in PANELS.items():check(sample+k+' program score',np.allclose(dense(a[pi,gs].X).mean(1),q[k]))
  for gp in ['High_only','Shared_high_other','Other_detected_only','No_detected_mac_neighbor']:
   use=expected==gp;den=0;countsiso={}
   for k,gs in ISOTYPES.items():countsiso[k]=int(dense(a[pi,gs].layers['counts'])[use].sum());den+=countsiso[k]
   for k,v in countsiso.items():
    row=iso[(iso['sample']==sample)&iso.group.eq(gp)&iso.isotype.eq(k)].iloc[0];check(sample+gp+k+' constant RNA counts',row.raw_constant_counts==v and row.total_heavy_constant_counts==den)
  for comp in COMPS:
   v=np.load(ROOT/f'{sample}_{comp}_audit.npz');H=v['assignments_first20'];NH=np.asarray(E@H.T).T;nt=np.asarray(E.sum(1)).ravel();NO=nt[None,:]-NH;check(sample+comp+' graph reassignment',np.array_equal(NH,v['target_high_neighbor_counts_first20']))
   for g in np.unique(v['source_strata']):check(sample+comp+str(g)+' source group conservation',np.all(H[:,v['source_strata']==g].sum(1)==h[v['source_strata']==g].sum()))
   def choose(hh,oo):
    if comp==COMPS[1]:return hh>0,hh==0
    l=(hh>0)&(oo==0);r=(hh==0)&(oo>0)
    if comp==COMPS[2]:l&=strict;r&=strict
    return l,r
   left,right=choose(nh,no)
   for j,k in enumerate(PROGRAMS):
    y=q[k].to_numpy();C,Q,rank,names=nuisance(a,pi,set(PANELS[k]));beta=np.linalg.lstsq(C,y,rcond=1e-10)[0];res=y-C@beta;check(sample+comp+k+' full OLS residual',np.allclose(res,v['residuals'][:,j],atol=1e-7));row=stat[(stat['sample']==sample)&stat.comparison.eq(comp)&stat.program.eq(k)].iloc[0];check(sample+comp+k+' observed group sizes',row.n_high_targets==left.sum() and row.n_comparator_targets==right.sum())
    for b in range(3):
     l,r=choose(NH[b],NO[b]);ok=l.sum()>=5 and r.sum()>=10
     if ok:check(sample+comp+k+str(b)+' source-graph null',np.isclose(res[l].mean()-res[r].mean(),v['null'][b,j]))
     else:check(sample+comp+k+str(b)+' invalid permutation missing',np.isnan(v['null'][b,j]))
    nn=v['null'][:,j];vv=nn[np.isfinite(nn)];check(sample+comp+k+' valid null denominator',len(vv)==row.valid_permutations)
    if row.status=='tested':
     obs=res[left].mean()-res[right].mean();pv=min(1,2*min((1+(vv>=obs-1e-12).sum())/(len(vv)+1),(1+(vv<=obs+1e-12).sum())/(len(vv)+1)));check(sample+comp+k+' adjusted effect/p',np.isclose(obs,row.adjusted_residual_delta) and np.isclose(pv,row.p))
 check('66 planned tests',len(stat)==66);check('joint BH66',np.allclose(stat.joint_q_66,false_discovery_control(stat.p.fillna(1).to_numpy())));ix=stat.comparison.eq(COMPS[0]);check('primary BH22',ix.sum()==22 and np.allclose(stat.loc[ix,'primary_q_22'],false_discovery_control(stat.loc[ix,'p'].fillna(1).to_numpy())))
 (ROOT/'audit_results.json').write_text(json.dumps(dict(checks=len(checks),passed=sum(c['passed'] for c in checks),details=checks),indent=2)+'\n');print('PASS',len(checks))
if __name__=='__main__':run()
