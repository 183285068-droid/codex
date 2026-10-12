import json,hashlib
import numpy as np,pandas as pd,h5py
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import false_discovery_control
from analyze import R,BASE,RAW,SAMPLES,PANELS,ANCHORS,PAIRS,GROUPS
checks=[]
def check(name,ok):checks.append(dict(check=name,passed=bool(ok)));assert ok,name
stats=pd.read_csv(R/'spatial_conditional_tests.csv');summary=pd.read_csv(R/'source_immune_neighborhoods.csv');lr=pd.read_csv(R/'LR_by_immune_subtype.csv')
for rec in json.loads((R/'input_manifest.json').read_text()):check('SHA '+rec['path'],hashlib.sha256(open(rec['path'],'rb').read()).hexdigest()==rec['sha256'])
for sample in SAMPLES:
 o=pd.read_csv(R/f'{sample}_immune_spots.csv.gz',index_col=0);baseline=pd.read_csv(BASE/f'{sample}_spot_annotations.csv.gz',index_col=0);check(sample+' frozen Q75 sources',o.source_group.equals(baseline.source_group));coords=np.column_stack([o.array_col*np.sqrt(3)/2,o.array_row*1.5]);pairs=cKDTree(coords).query_pairs(np.sqrt(3)+1e-8);edges={(i,j) for i,j in pairs}|{(j,i) for i,j in pairs};ij=np.array(list(edges));G=sparse.csr_matrix((np.ones(len(ij)),(ij[:,0],ij[:,1])),shape=(len(o),len(o)));deg=G.sum(1).A1
 with h5py.File(RAW/sample/'filtered_feature_bc_matrix.h5') as h:
  m=h['matrix'];raw=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=m['shape'][:]);names=np.array([v.decode() for v in m['features']['name'][:]]);bc=[v.decode() for v in m['barcodes'][:]]
 bi={b:i for i,b in enumerate(bc)};raw=raw[:,[bi[b] for b in o.index]];check(sample+' whole-library normalization',np.allclose(raw.sum(0).A1,o.total_counts));maskdict={}
 def rawgene(g):return raw[names==g].sum(0).A1
 for k,gs in PANELS.items():
  count=np.column_stack([rawgene(g) for g in gs]);score=np.log1p(count/o.total_counts.to_numpy()[:,None]*1e4).mean(1);mask=(score>=np.quantile(score,.75))&(score>0)&((count>0).sum(1)>=2)&(deg>0)
  if k=='T':mask=baseline.T_RNA.to_numpy()
  else:
   anchor=np.column_stack([rawgene(g)>0 for g in ANCHORS[k]]);mask&=anchor.all(1) if k=='Plasma' else anchor.any(1)
  check(sample+k+' score',np.allclose(score,o['score_'+k]));check(sample+k+' gate',np.array_equal(mask,o['gate_'+k]));maskdict[k]=mask;nt=G@mask.astype(float);fraction=np.divide(nt,deg,out=np.zeros(len(o)),where=deg>0);check(sample+k+' geometry/coverage',np.allclose(fraction,o['neighbor_fraction_'+k]))
  for group in GROUPS:
   use=o.source_group.eq(group).to_numpy();row=summary[(summary['sample']==sample)&(summary.subtype==k)&(summary.source_group==group)].iloc[0];check(sample+k+group+' summary',row.n_sources==use.sum() and row.n_edges==nt[use].sum() and np.isclose(row.mean_neighbor_fraction,fraction[use].mean()))
 e=pd.read_csv(R/f'{sample}_directed_LR_edges.csv.gz');index={b:i for i,b in enumerate(o.index)};ri=np.array([index[b] for b in e.source_barcode]);ci=np.array([index[b] for b in e.target_barcode]);union=np.column_stack(list(maskdict.values())).any(1);check(sample+' eligible edge set',set(zip(ri,ci))=={(i,j) for i,j in edges if o.Mac_RNA.iloc[i] and union[j]})
 for l,r in PAIRS:
  ligand=rawgene(l)/o.total_counts.to_numpy()*1e4;receptor=np.column_stack([rawgene(g)/o.total_counts.to_numpy()*1e4 for g in r.split('_')]).min(1);score=np.sqrt(ligand[ri]*receptor[ci]);check(sample+l+r+' raw LR',np.allclose(score,e[l+'--'+r]))
  for k,mask in maskdict.items():
   for group in GROUPS:
    use=e.source_group.eq(group).to_numpy()&mask[ci];v=score[use];row=lr[(lr['sample']==sample)&(lr.subtype==k)&(lr.source_group==group)&(lr.ligand==l)&(lr.receptor==r)].iloc[0];check(sample+l+r+k+group+' LR count',len(v)==row.n_edges and int((v>0).sum())==row.codetected_edges)
 path=R/f'{sample}_immune_coverage_permutation_audit.npz'
 if path.exists():
  with np.load(path) as p:
   ix=p['source_indices'];hh=p['high'];Y=p['outcomes'];check(sample+' test outcomes',np.allclose(Y,o[['neighbor_fraction_'+k for k in PANELS]].iloc[ix]));obs=hh@Y/hh.sum()-(1-hh)@Y/(len(hh)-hh.sum());check(sample+' test delta',np.allclose(obs,p['observed_deltas']))
   for j,label in enumerate(p['labels']):
    row=stats[(stats['sample']==sample)&(stats.feature==label)].iloc[0]
    if pd.isna(row.spatial_conditional_p):continue
    count=(np.abs(p['null_deltas'][:,j]-p['conditional_null_center'][j])>=abs(obs[j]-p['conditional_null_center'][j])-1e-14).sum();pvalue=count/len(p['null_deltas']) if str(p['mode'])=='exact_conditional_enumeration' else (count+1)/(len(p['null_deltas'])+1);check(sample+label+' permutation p',np.isclose(pvalue,row.spatial_conditional_p))
check('63 planned tests',len(stats)==63);valid=stats.spatial_conditional_p.notna();q=false_discovery_control(stats.spatial_conditional_p.fillna(1).to_numpy());check('BH planned family',np.allclose(q[valid],stats.loc[valid,'BH_q_63_planned']));pd.DataFrame(checks).to_csv(R/'audit_checks.csv',index=False);print('PASS',len(checks),'checks')
