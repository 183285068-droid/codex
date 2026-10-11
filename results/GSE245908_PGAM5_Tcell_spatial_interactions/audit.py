"""Independent HDF5 reconstruction, edge aggregation, null arithmetic and BH audit."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE245908_Tcell_mpl')
import numpy as np,pandas as pd,h5py,scipy.sparse as sp,json,hashlib
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent/'GSE245908_spatial_PGAM5';checks=[]
def check(k,ok,detail=''):
 checks.append({'check':k,'pass':bool(ok),'detail':str(detail)})
 if not ok:raise AssertionError((k,detail))
def close(k,x,y,tol=1e-6):check(k,np.allclose(x,y,rtol=1e-6,atol=tol,equal_nan=True),np.nanmax(np.abs(np.asarray(x)-np.asarray(y))))
if __name__=='__main__':
 programs=json.loads((ROOT/'programs.json').read_text());pairs=pd.read_csv(ROOT/'selected_LR_pairs.csv');tests=pd.read_csv(ROOT/'spatial_conditional_tests.csv');summary=pd.read_csv(ROOT/'sample_state_summary.csv');support=pd.read_csv(ROOT/'permutation_support.csv');resource=pd.read_csv(ROOT/'source_inputs/omni_resource.csv')
 for rec in json.loads((ROOT/'input_manifest.json').read_text()):
  f=ROOT.parent/rec['path'];check('hash_'+rec['path'],hashlib.sha256(f.read_bytes()).hexdigest()==rec['sha256'])
 for p in pairs.itertuples():check('database_'+p.ligand+'_'+p.receptor_complex,((resource.source_genesymbol==p.ligand)&(resource.target_genesymbol==p.receptor_complex)).any())
 for s,gsm in [('CHC20','GSM7850822'),('CHC23','GSM7850823')]:
  o=pd.read_csv(ROOT/f'{s}_spot_states.csv',index_col=0);old=pd.read_csv(BASE/f'{s}_spot_annotations.csv',index_col=0)
  check(s+'_same_QC_IDs',list(o.index)==list(old.index));check(s+'_PGAM5_zeros_preserved',np.array_equal(o.PGAM5_raw,old.PGAM5_raw));close(s+'_PGAM5_expression',o.PGAM5_log1p_CP10k,old.PGAM5_log1p_CP10k)
  with h5py.File(BASE/'inputs'/s/f'{gsm}_{s}_filtered_feature_bc_matrix.h5','r') as f:
   m=f['matrix'];X=sp.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=m['shape'][:]).T.tocsr();b=m['barcodes'].asstr()[:];names=m['features']['name'].asstr()[:]
  ix={v:i for i,v in enumerate(b)};X=X[[ix[v] for v in o.index]].copy();total=np.asarray(X.sum(1)).ravel();close(s+'_full_library_totals',total,o.total_counts,tol=0)
  lookup={g:int(np.flatnonzero(names==g)[0]) for g in set(g for v in programs.values() for g in v)|set(pairs.ligand)|set(g for r in pairs.receptor_complex for g in r.split('_'))}
  for key,gs in programs.items():
   raw=X[:,[lookup[g] for g in gs]].toarray();score=np.log1p(raw/total[:,None]*1e4).mean(1);close(s+'_raw_program_'+key,score,o[key]);check(s+'_program_detection_'+key,np.array_equal((raw>0).sum(1),o[key+'_detected_genes']))
  mac=(o.Macrophage_core>=o.Macrophage_core.quantile(.75))&(o.Pan_myeloid>=o.Pan_myeloid.median())&(o.Macrophage_core_detected_genes>=3)&(o.spatial_degree>0)
  t=(o.T_lineage>=o.T_lineage.quantile(.75))&(o.T_lineage_detected_genes>=2)&(o.spatial_degree>0);check(s+'_mac_gate',np.array_equal(mac,o.Mac_enriched));check(s+'_T_gate',np.array_equal(t,o.T_enriched))
  pos=mac&(o.PGAM5_raw>0);cut=o.loc[pos,'PGAM5_log1p_CP10k'].quantile(.75);high=pos&(o.PGAM5_log1p_CP10k>=cut);check(s+'_high_Q75',np.array_equal(high,o.source_group.eq('PGAM5_high')))
  for key in ['Activation_early','Effector_cytotoxic','Exhaustion_core']:
   q=o.loc[t,key].quantile(.75);expected=t&(o[key]>=q)&(o[key]>0)&(o[key+'_detected_genes']>=2);check(s+'_state_'+key,np.array_equal(expected,o[key+'_supported_high']))
  # Neighbors independently reconstructed from actual hex array geometry.
  xy=np.column_stack([o.array_col,np.sqrt(3)*o.array_row]);n=len(o);geopairs=cKDTree(xy).query_pairs(2.001);ri=np.array([i for i,j in geopairs]+[j for i,j in geopairs]);ci=np.array([j for i,j in geopairs]+[i for i,j in geopairs]);W=sp.csr_matrix((np.ones(len(ri)),(ri,ci)),shape=(n,n));deg=np.asarray(W.sum(1)).ravel();nt=W@t.to_numpy().astype(float);close(s+'_T_neighbor_count',nt,o.T_neighbor_count);close(s+'_T_coverage',np.divide(nt,deg,out=np.zeros(n),where=deg>0),o.T_neighbor_fraction)
  for key in ['Activation_early','Effector_cytotoxic','Exhaustion_core','Exhaustion_expanded']:
   y=np.divide(W@(o[key].to_numpy()*t),nt,out=np.full(n,np.nan),where=nt>0);close(s+'_neighbor_'+key,y,o['T_neighbor_'+key])
  edges=pd.read_csv(ROOT/f'{s}_directed_LR_edges.csv.gz');idx={b:i for i,b in enumerate(o.index)};ridx=np.array([idx[v] for v in edges.source_barcode]);cidx=np.array([idx[v] for v in edges.target_barcode]);check(s+'_edge_no_self',np.all(ridx!=cidx));check(s+'_edge_lineage_gates',mac.iloc[ridx].all() and t.iloc[cidx].all());expected={(i,j) for i,j in zip(ri,ci) if mac.iloc[i] and t.iloc[j]};check(s+'_edge_geometry',set(zip(ridx,cidx))==expected)
  per=pd.read_csv(ROOT/f'{s}_per_source_LR_scores.csv',index_col=0);take=[idx[b] for b in per.index]
  for p in pairs.itertuples():
   l=X[:,lookup[p.ligand]].toarray().ravel()/total*1e4;r=np.min(X[:,[lookup[g] for g in p.receptor_complex.split('_')]].toarray()/total[:,None]*1e4,axis=1);v=np.sqrt(l[ridx]*r[cidx]);key=p.ligand+'--'+p.receptor_complex;close(s+'_edge_score_'+key,v,edges[key]);sv=np.bincount(ridx,weights=v,minlength=n);mean=np.divide(sv,nt,out=np.zeros(n),where=nt>0);close(s+'_source_score_'+key,mean[take],per[key]);check(s+'_complex_zero_'+key,np.all(v[(l[ridx]==0)|(r[cidx]==0)]==0))
  for kind in ['T_coverage','T_programs','LR']:
   z=np.load(ROOT/f'{s}_{kind}_permutation_audit.npz');h=z['high'];Y=z['outcomes'];observed=Y[h==1].mean(0)-Y[h==0].mean(0);close(s+'_'+kind+'_delta',observed,z['observed_deltas'],tol=1e-12)
   for j,key in enumerate(z['labels']):
    null=z['null_deltas'][:,j];center=z['conditional_null_center'][j];count=(abs(null-center)>=abs(observed[j]-center)-1e-14).sum();p=count/len(null) if str(z['mode'])=='exact_conditional_enumeration' else (count+1)/(len(null)+1);row=tests[(tests['sample']==s)&(tests.kind==kind)&(tests.feature==key)].iloc[0];close(s+'_'+kind+'_p_'+key,p,row.spatial_conditional_p,tol=1e-12);close(s+'_'+kind+'_mean_high_'+key,Y[h==1,j].mean(),row.mean_high);close(s+'_'+kind+'_mean_other_'+key,Y[h==0,j].mean(),row.mean_other_detected)
   if str(z['mode'])=='exact_conditional_enumeration':close(s+'_'+kind+'_exact_null_center',z['null_deltas'].mean(0),z['conditional_null_center'],tol=1e-12);check(s+'_'+kind+'_exact_configs',len(z['null_deltas'])==30)
 valid=tests.spatial_conditional_p.notna();p=tests.loc[valid,'spatial_conditional_p'].to_numpy();order=np.argsort(p);adj=np.minimum.accumulate((len(p)*p[order]/np.arange(1,len(p)+1))[::-1])[::-1];q=np.empty(len(p));q[order]=np.minimum(adj,1);close('BH_50_tests',q,tests.loc[valid,'BH_q'],tol=1e-12);check('50_prespecified_tests',len(tests)==50)
 pd.DataFrame(checks).to_csv(ROOT/'independent_audit_checks.csv',index=False);print('Passed',len(checks),'independent checks')
