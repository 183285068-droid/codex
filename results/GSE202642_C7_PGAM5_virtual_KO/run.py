from pathlib import Path
import json,numpy as np,pandas as pd,anndata as ad,scanpy as sc
from scTenifold.core._base import scTenifoldKnk
from threadpoolctl import threadpool_limits
R=Path(__file__).resolve().parent
if not (R/'C7_raw.h5ad').exists():
 meta=pd.read_csv('/workspace/codex/results/GSE202642_PGAM5_subpopulation_assessment/cell_annotations.csv.gz',index_col=0);sel=meta.index[meta.primary_cluster.eq(7)];b=ad.read_h5ad('/workspace/scratch/GSE202642_HCC_allcell_annotation/HCC_raw_counts_QC.h5ad',backed='r');a=b[sel,:].to_memory();b.file.close();a.obs=meta.loc[sel].copy();assert len(a)==403;assert np.array_equal(a.layers['counts'][:,a.var_names.get_loc('PGAM5')].toarray().ravel(),a.obs.PGAM5_counts);a.X=a.layers['counts'].copy();a.layers.clear();a.write_h5ad(R/'C7_raw.h5ad',compression='gzip')
else:a=ad.read_h5ad(R/'C7_raw.h5ad')
assert np.array_equal(np.asarray(a.X.sum(1)).ravel(),a.obs.total_counts)
z=a.copy();sc.pp.normalize_total(z,target_sum=1e4);sc.pp.log1p(z);sc.pp.highly_variable_genes(z,n_top_genes=1500,flavor='seurat');detect=np.asarray((a.X>0).mean(0)).ravel();eligible=(detect>=.05)&~a.var_names.str.match(r'^(MT-|RPL|RPS)');hv=z.var.loc[eligible].sort_values('dispersions_norm',ascending=False).head(1500).index.tolist();genes=sorted(set(hv+['PGAM5']));pd.DataFrame({'gene':a.var_names,'detection_fraction':detect,'eligible':eligible,'network_selected':a.var_names.isin(genes)}).to_csv(R/'gene_selection.csv',index=False)
x=pd.DataFrame(a[:,genes].X.toarray().T,index=genes,columns=a.obs_names);x.to_csv(R/'network_raw_counts.csv.gz');a.obs.to_csv(R/'target_cell_metadata.csv.gz');cp=x.div(a.obs.total_counts,axis='columns')*1e6
qc=dict(min_lib_size=0,remove_outlier_cells=False,min_percent=0,max_mito_ratio=1,min_exp_sum=0)
for seed in [202642,42]:
 if (R/f'differential_regulation_seed{seed}.csv').exists():continue
 print('SEED',seed,len(genes),flush=True)
 obj=scTenifoldKnk(x,ko_genes=['PGAM5'],qc_kws=qc,nc_kws=dict(n_nets=10,n_samp_cells=300,n_comp=3,q=.95,random_state=seed),td_kws=dict(K=5,max_iter=1000,tol=1e-5,random_state=seed),ma_kws=dict(d=30))
 with threadpool_limits(limits=4):
  obj.run_step('qc');assert obj.QC_dict['WT'].columns.equals(cp.columns);obj.QC_dict['WT']=cp.copy()
  for step in ['nc','td','ko','ma','dr']:obj.run_step(step)
 obj.d_regulation.to_csv(R/f'differential_regulation_seed{seed}.csv',index=False)
 obj.tensor_dict['WT'].to_csv(R/f'WT_network_seed{seed}.csv.gz');obj.tensor_dict['KO'].to_csv(R/f'KO_network_seed{seed}.csv.gz')
 assert np.all(obj.tensor_dict['KO'].loc['PGAM5'].to_numpy()==0)
 assert np.any(obj.tensor_dict['WT'].loc['PGAM5'].to_numpy()!=0)
 print(obj.d_regulation.head(12).to_string(index=False),flush=True)
(R/'protocol.json').write_text(json.dumps(dict(target='Original 15-cluster C7 MKI67/TOP2A',cells=len(a),PGAM5_detected=int(a.obs.PGAM5_counts.gt(0).sum()),samples=a.obs.sample_name.value_counts().to_dict(),genes=len(genes),normalization='CPM using original full-library totals; QC normalization overridden before network construction',gene_filter='>=5% detection, exclude MT/RPL/RPS, top1500 normalized dispersions plus PGAM5',method='scTenifoldpy0.5.1 scTenifoldKnk default row-zero knockout',networks=10,cells_per_network=300,replace=True,PCs=3,sparsity_quantile=.95,tensor_rank=5,manifold_dimensions=30,seeds=[202642,42]),indent=2));print('DONE',flush=True)
