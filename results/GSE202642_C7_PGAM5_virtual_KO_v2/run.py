from pathlib import Path
import json,time,numpy as np,pandas as pd
from scipy.stats import spearmanr
from scTenifold.core._networks import make_networks,manifold_alignment,d_regulation
from scTenifold.core._decomposition import cp_decomposition
from threadpoolctl import threadpool_limits
R=Path(__file__).resolve().parent;src=R/'inputs';src.mkdir(exist_ok=True)
import shutil
old=Path('/workspace/codex/results/GSE202642_C7_PGAM5_virtual_KO')
for f in ['network_raw_counts.csv.gz','target_cell_metadata.csv.gz','program_definitions.json']:
 if not (src/f).exists():shutil.copy2(old/f,src/f)
x=pd.read_csv(src/'network_raw_counts.csv.gz',index_col=0);meta=pd.read_csv(src/'target_cell_metadata.csv.gz',index_col=0);assert x.columns.equals(meta.index);cp=x.div(meta.total_counts,axis='columns')*1e6
seeds=[202642,42,101,202,303]
protocol=dict(cells=403,PGAM5_detected=79,genes=len(x),seeds=seeds,criteria={'minimum_pairwise_distance_rho':.8,'minimum_pairwise_top50_Jaccard':.5,'candidate_FDR_frequency':'at least 4 of 5 seeds','no_KO_significant':0},precision_only=dict(n_nets=10,rank=5,q=.95,n_decimal=6,tol=1e-5),ensemble=dict(n_nets=30,rank=10,q=.95,n_decimal=6,tol=1e-7,max_iter=2000),normalization='CPM original full-library totals, same fixed1501 genes as previous run',network=dict(n_samp_cells=300,replace=True,n_comp=3),alignment_dimensions=30)
(R/'protocol.json').write_text(json.dumps(protocol,indent=2))
rows=[]
for scenario in ['precision_only','ensemble']:
 cfg=protocol[scenario]
 for seed in seeds:
  folder=R/scenario/str(seed);folder.mkdir(parents=True,exist_ok=True)
  if (folder/'metrics.json').exists():rows.append(json.loads((folder/'metrics.json').read_text()));continue
  started=time.time();np.random.seed(seed);print('START',scenario,seed,flush=True)
  with threadpool_limits(limits=4):
   nets=make_networks(cp,n_nets=cfg['n_nets'],n_samp_cells=300,n_comp=3,q=cfg['q'],random_state=seed,replace=True)
   fit=cp_decomposition([n.toarray() for n in nets],K=cfg['rank'],max_iter=cfg.get('max_iter',1000),tol=cfg['tol'],random_state=seed)
   w=fit['slice_sum']/len(nets);w=np.round(w/np.abs(w).max(),cfg['n_decimal']);w=w.T;np.fill_diagonal(w,0);wt=pd.DataFrame(w,index=x.index,columns=x.index);ko=wt.copy();ko.loc['PGAM5']=0
   result=d_regulation(manifold_alignment(wt,ko,d=30),ko_genes=['PGAM5'])
  result.to_csv(folder/'differential_regulation.csv',index=False);np.savez_compressed(folder/'WT_network.npz',weights=w,genes=x.index.to_numpy(dtype=str));np.save(folder/'tensor_residuals.npy',np.asarray(fit['all_resids']));result=result[result.Gene.ne('PGAM5')];metrics=dict(scenario=scenario,seed=seed,significant_genes=int(result['adjusted p-value'].lt(.05).sum()),outgoing_edges=int((wt.loc['PGAM5']!=0).sum()),outgoing_max_weight=float(wt.loc['PGAM5'].abs().max()),tensor_converged=bool(fit['conv']),tensor_iterations=len(fit['all_resids']),seconds=time.time()-started);(folder/'metrics.json').write_text(json.dumps(metrics,indent=2));rows.append(metrics);print(metrics,flush=True)
pd.DataFrame(rows).to_csv(R/'run_metrics.csv',index=False);print('ALL DONE',flush=True)
