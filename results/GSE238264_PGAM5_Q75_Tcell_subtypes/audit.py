import json,hashlib
import numpy as np,pandas as pd,h5py
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import false_discovery_control
from analyze import R,BASE,RAW,SAMPLES,PANELS,TYPES,STATES,PAIRS,GROUPS
checks=[]
def check(name,ok):checks.append(dict(check=name,passed=bool(ok)));assert ok,name
stats=pd.read_csv(R/'spatial_T_subtype_tests.csv');summary=pd.read_csv(R/'T_subtype_neighbors.csv');lr=pd.read_csv(R/'LR_by_T_subtype.csv')
for rec in json.loads((R/'input_manifest.json').read_text()):check('SHA '+rec['path'],hashlib.sha256(open(rec['path'],'rb').read()).hexdigest()==rec['sha256'])
for sample in SAMPLES:
 o=pd.read_csv(R/f'{sample}_T_subtypes.csv.gz',index_col=0);base=pd.read_csv(BASE/f'{sample}_spot_annotations.csv.gz',index_col=0);check(sample+' frozen sources',o.source_group.equals(base.source_group));xy=np.column_stack([o.array_col*np.sqrt(3)/2,o.array_row*1.5]);pairs=cKDTree(xy).query_pairs(np.sqrt(3)+1e-8);edges={(i,j) for i,j in pairs}|{(j,i) for i,j in pairs};ij=np.array(list(edges));G=sparse.csr_matrix((np.ones(len(ij)),(ij[:,0],ij[:,1])),shape=(len(o),len(o)));deg=G.sum(1).A1
 with h5py.File(RAW/sample/'filtered_feature_bc_matrix.h5') as h:
  m=h['matrix'];raw=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=m['shape'][:]);names=np.array([g.decode() for g in m['features']['name'][:]]);bc=[b.decode() for b in m['barcodes'][:]]
 bi={b:i for i,b in enumerate(bc)};raw=raw[:,[bi[b] for b in o.index]];check(sample+' whole-library UMI',np.allclose(raw.sum(0).A1,o.total_counts))
 def count(g):return raw[names==g].sum(0).A1
 for k,gs in PANELS.items():
  c=np.column_stack([count(g) for g in gs]);score=np.log1p(c/o.total_counts.to_numpy()[:,None]*1e4).mean(1);check(sample+k+' score',np.allclose(score,o['score_'+k]));check(sample+k+' detection',np.array_equal((c>0).sum(1),o['detected_'+k]))
 t=base.T_RNA.to_numpy();cd8=t&(count('CD8A')>0)&(o.score_CD8_identity.to_numpy()>o.score_CD4_identity.to_numpy());cd4=t&(o.score_CD4_identity.to_numpy()>0)&(o.score_CD4_identity.to_numpy()>o.score_CD8_identity.to_numpy());reg=cd4&(count('FOXP3')>0)&(o.detected_Treg_program.to_numpy()>=2);masks={'All_T':t,'All_CD8':cd8,'Treg':reg,'CD4_other':cd4&~reg,'T_unresolved':t&~(cd8|cd4)}
 for k,p in STATES.items():masks[k]=cd8&(o['score_'+p]>=o.loc[cd8,'score_'+p].quantile(.75))&(o['score_'+p]>0)&(o['detected_'+p]>=2)
 check(sample+' identity partition',np.array_equal(o[['gate_All_CD8','gate_CD4_other','gate_Treg','gate_T_unresolved']].sum(1).to_numpy(),t.astype(int)))
 for k in TYPES:
  mask=np.asarray(masks[k]);check(sample+k+' gate',np.array_equal(mask,o['gate_'+k]));nt=G@mask.astype(float);fraction=np.divide(nt,deg,out=np.zeros(len(o)),where=deg>0);check(sample+k+' adjacency',np.allclose(fraction,o['neighbor_fraction_'+k]))
  for group in GROUPS:
   use=o.source_group.eq(group).to_numpy();row=summary[(summary['sample']==sample)&(summary.subtype==k)&(summary.source_group==group)].iloc[0];check(sample+k+group+' summary',row.n_edges==nt[use].sum() and row.n_targets==mask.sum() and np.isclose(row.mean_neighbor_fraction,fraction[use].mean()))
 n8=G@cd8.astype(float)
 for p in STATES.values():
  v=np.divide(G@(o['score_'+p].to_numpy()*cd8),n8,out=np.full(len(o),np.nan),where=n8>0);check(sample+p+' continuous CD8 neighbor score',np.allclose(v,o['CD8_neighbor_'+p],equal_nan=True))
 e=pd.read_csv(R/f'{sample}_T_LR_edges.csv.gz');ix={b:i for i,b in enumerate(o.index)};ri=np.array([ix[b] for b in e.source_barcode]);ci=np.array([ix[b] for b in e.target_barcode]);check(sample+' T eligible edge set',set(zip(ri,ci))=={(i,j) for i,j in edges if o.Mac_RNA.iloc[i] and t[j]})
 for l,r in PAIRS:
  lig=count(l)/o.total_counts.to_numpy()*1e4;rec=np.column_stack([count(g)/o.total_counts.to_numpy()*1e4 for g in r.split('_')]).min(1);score=np.sqrt(lig[ri]*rec[ci]);check(sample+l+r+' LR counts',np.allclose(score,e[l+'--'+r]))
  for k in TYPES:
   mask=np.asarray(masks[k])
   for group in GROUPS:
    use=e.source_group.eq(group).to_numpy()&mask[ci];v=score[use];row=lr[(lr['sample']==sample)&(lr.subtype==k)&(lr.source_group==group)&(lr.ligand==l)&(lr.receptor==r)].iloc[0];check(sample+l+r+k+group+' LR aggregation',len(v)==row.n_edges and int((v>0).sum())==row.codetected_edges)
 for kind in ['T_subtype_coverage','CD8_neighbor_programs']:
  path=R/f'{sample}_{kind}_permutation_audit.npz'
  if not path.exists():continue
  with np.load(path) as p:
   h=p['high'];Y=p['outcomes'];obs=h@Y/h.sum()-(1-h)@Y/(len(h)-h.sum());check(sample+kind+' observed delta',np.allclose(obs,p['observed_deltas']))
   for j,label in enumerate(p['labels']):
    row=stats[(stats['sample']==sample)&(stats.kind==kind)&(stats.feature==label)].iloc[0]
    if pd.isna(row.spatial_conditional_p):continue
    n=(np.abs(p['null_deltas'][:,j]-p['conditional_null_center'][j])>=abs(obs[j]-p['conditional_null_center'][j])-1e-14).sum();pv=n/len(p['null_deltas']) if str(p['mode'])=='exact_conditional_enumeration' else (n+1)/(len(p['null_deltas'])+1);check(sample+kind+label+' conditional p',np.isclose(pv,row.spatial_conditional_p))
check('91 fixed tests',len(stats)==91);valid=stats.spatial_conditional_p.notna();q=false_discovery_control(stats.spatial_conditional_p.fillna(1).to_numpy());check('joint BH',np.allclose(q[valid],stats.loc[valid,'BH_q_91_planned']));pd.DataFrame(checks).to_csv(R/'audit_checks.csv',index=False);print('PASS',len(checks),'checks')
