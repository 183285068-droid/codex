from pathlib import Path
import json,time,sys,hashlib
import numpy as np,pandas as pd
from scTenifold.core._networks import make_networks,manifold_alignment,d_regulation
from scTenifold.core._decomposition import cp_decomposition
from threadpoolctl import threadpool_limits

R=Path(__file__).resolve().parent
x=pd.read_csv(R/'inputs/network_raw_counts.csv.gz',index_col=0)
meta=pd.read_csv(R/'inputs/target_cell_metadata.csv.gz',index_col=0)
assert x.columns.equals(meta.index) and len(meta)==171 and (x.loc['PGAM5']>0).sum()==22
cp=x.div(meta.total_counts,axis='columns')*1e6
protocol=dict(dataset='GSE151530 HCC',target='M2 MKI67/TOP2A cycling macrophages, current 16-group version',cells=171,PGAM5_detected=22,genes=len(x),seeds=[151530,42,101,202,303],normalization='CPM, original full-transcriptome raw library totals; all zero values retained',gene_selection='>=5% detection, exclude MT/RPL/RPS; top1500 normalized dispersions plus PGAM5',network=dict(n_nets=30,n_samp_cells=171,n_comp=3,q=.95,replace=True),tensor=dict(rank=5,max_iter=1000,tol=1e-5,n_decimal=6),alignment_dimensions=30,stability_criteria=dict(minimum_pairwise_distance_rho=.8,minimum_pairwise_top50_Jaccard=.5,all_tensors_converged=True,all_WT_PGAM5_rows_nonzero=True,all_no_KO_controls_zero_significant=True),candidate_criterion='model BH q<0.05 in at least4/5 seeds; descriptive frequency, not a combined p-value')
(R/'protocol.json').write_text(json.dumps(protocol,indent=2))
seeds=[int(v) for v in sys.argv[1:]] or protocol['seeds']
for seed in seeds:
 assert seed in protocol['seeds']
 folder=R/'models'/str(seed);folder.mkdir(parents=True,exist_ok=True)
 if (folder/'metrics.json').exists():continue
 started=time.time();np.random.seed(seed);print('START',seed,flush=True)
 with threadpool_limits(limits=4):
  nets=make_networks(cp,random_state=seed,**protocol['network'])
  print('NETWORKS DONE',seed,flush=True)
  fit=cp_decomposition([n.toarray() for n in nets],K=protocol['tensor']['rank'],max_iter=protocol['tensor']['max_iter'],tol=protocol['tensor']['tol'],random_state=seed)
  print('TENSOR DONE',seed,fit['conv'],len(fit['all_resids']),flush=True)
  w=fit['slice_sum']/len(nets);assert np.isfinite(w).all() and np.abs(w).max()>0
  w=np.round(w/np.abs(w).max(),6).T;np.fill_diagonal(w,0)
  wt=pd.DataFrame(w,index=x.index,columns=x.index);ko=wt.copy();ko.loc['PGAM5']=0
  align=manifold_alignment(wt,ko,d=30);result=d_regulation(align,ko_genes=['PGAM5'])
  np.random.seed(seed);null_align=manifold_alignment(wt,wt,d=30);null=d_regulation(null_align,ko_genes=['PGAM5'])
 np.savez_compressed(folder/'WT_KO_networks.npz',WT=w,KO_delta=ko.to_numpy()-w,KO_SHA256=hashlib.sha256(ko.to_numpy().tobytes()).hexdigest(),genes=x.index.to_numpy(dtype=str))
 np.save(folder/'tensor_residuals.npy',np.asarray(fit['all_resids']))
 align.to_csv(folder/'manifold_alignment.csv.gz');result.to_csv(folder/'differential_regulation.csv',index=False);null.to_csv(folder/'no_KO_control.csv',index=False)
 downstream=result[result.Gene.ne('PGAM5')]
 metrics=dict(seed=seed,significant_downstream_genes=int(downstream['adjusted p-value'].lt(.05).sum()),WT_PGAM5_outgoing_edges=int((wt.loc['PGAM5']!=0).sum()),WT_PGAM5_max_abs_weight=float(wt.loc['PGAM5'].abs().max()),tensor_converged=bool(fit['conv']),tensor_iterations=len(fit['all_resids']),tensor_norm_percent=float(fit['norm_percent']),numerical_noise_threshold=float(np.sqrt(np.finfo(float).eps)*np.abs(align.to_numpy()).max()),max_downstream_distance=float(downstream.Distance.max()),downstream_above_noise=int(downstream.Distance.gt(np.sqrt(np.finfo(float).eps)*np.abs(align.to_numpy()).max()).sum()),no_KO_significant=int(null['adjusted p-value'].lt(.05).sum()),max_null_distance=float(null.Distance.max()),seconds=time.time()-started)
 (folder/'metrics.json').write_text(json.dumps(metrics,indent=2));print('DONE',metrics,flush=True)
