"""Exploratory within-slice association with a graph-spectral spatial null.
Spectral randomization is model-based, not evidence from independent patients.
"""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE245908_mpl')
import numpy as np,pandas as pd,scanpy as sc,scipy.sparse as sp,json
from scipy.linalg import eigh
from scipy.stats import rankdata,spearmanr
from statsmodels.stats.multitest import multipletests
from preprocess import ROOT,WORK,SAMPLES,PROGRAMS,read,array
SEED=245908;N_NULL=4999

def hex_graph(obs):
 rc=obs[['array_row','array_col']].to_numpy().astype(int);lookup={tuple(x):i for i,x in enumerate(rc)};rows=[];cols=[]
 assert len(lookup)==len(obs)
 for i,(r,c) in enumerate(rc):
  for dr,dc in [(0,-2),(0,2),(-1,-1),(-1,1),(1,-1),(1,1)]:
   j=lookup.get((r+dr,c+dc))
   if j is not None:rows.append(i);cols.append(j)
 w=sp.csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(obs),len(obs)))
 assert (w!=w.T).nnz==0 and w.diagonal().sum()==0
 return w

def residual_rank(y,cov):
 y=rankdata(np.round(np.asarray(y,dtype=float),12));c=np.column_stack([np.ones(len(y)),cov]);r=y-c@np.linalg.lstsq(c,y,rcond=None)[0];return r-r.mean()
def covariance(obs):
 return np.column_stack([rankdata(np.log1p(obs.total_counts)),rankdata(obs.n_genes),rankdata(obs.pct_mt)])
def corr(x,y):return float(x@y/np.sqrt((x@x)*(y@y)))
def moran(x,w):
 x=x-x.mean();return float(len(x)/w.sum()*(x@(w@x))/(x@x))
def spectral_basis(w):
 v=w.toarray();means=v.mean(0);v-=means[:,None];v-=means[None,:];v+=means.mean()
 val,vec=eigh(v,overwrite_a=True,check_finite=False,driver='evd');return val,vec

def analyze(a,s):
 for column in a.obs.select_dtypes(include=['floating']).columns:a.obs[column]=a.obs[column].astype(np.float64)
 w0=hex_graph(a.obs);degree0=np.asarray(w0.sum(1)).ravel();keep=degree0>0
 b=a[keep].copy();w=w0[keep][:,keep];degree=np.asarray(w.sum(1)).ravel();wr=sp.diags(1/degree)@w
 cov=covariance(b.obs);x=residual_rank(b.obs.PGAM5_log1p_CP10k.to_numpy(),cov)
 print(s,'spectral eigendecomposition',len(b),flush=True);eig,U=spectral_basis(w)
 rng=np.random.default_rng(SEED+(s=='CHC23'));cx=U.T@x
 signs=rng.choice([-1.,1.],size=(N_NULL,len(b))).astype(np.float32)
 # Null sign flips preserve exact centered norm and graph quadratic form.
 numerical=[]
 for k in range(8):
  xs=U@(cx*signs[k]);numerical.append({'sample':s,'null':k,'norm_relative_error':float(abs(xs@xs-x@x)/(x@x)),'moran_absolute_error':float(abs(moran(xs,w)-moran(x,w))),'mean_absolute_error':float(abs(xs.mean()))})
 out=[];mrows=[];stats_maps={};null_products=[];null_denominators=[];null_labels=[]
 for k in ['PGAM5_log1p_CP10k',*PROGRAMS]:
  y=b.obs[k].to_numpy();I=moran(y,w);null=np.empty(N_NULL)
  for t in range(N_NULL):null[t]=moran(rng.permutation(y),w)
  center=-1/(len(b)-1)
  p=(1+np.count_nonzero(abs(null-center)>=abs(I-center)))/(N_NULL+1)
  mrows.append({'sample':s,'feature':k,'n_spots':len(b),'Moran_I':I,'null_expected_I':center,'two_sided_label_permutation_p':p,'n_permutations':N_NULL})
 for k in PROGRAMS:
  y=b.obs[k].to_numpy();raw=spearmanr(np.round(b.obs.PGAM5_log1p_CP10k,12),np.round(y,12))
  for kind in ['same_spot','adjacent_spot']:
   yy=y if kind=='same_spot' else wr@y
   yr=residual_rank(yy,cov);cy=U.T@yr;den=np.sqrt((x@x)*(yr@yr));rho=corr(x,yr)
   null=signs@(cx*cy)/den
   null_products.append(cx*cy);null_denominators.append(den);null_labels.append(k+'|'+kind)
   p=(1+np.count_nonzero(abs(null)>=abs(rho)))/(N_NULL+1)
   out.append({'sample':s,'program':k,'relationship':kind,'n_spots':len(b),'raw_Spearman_rho':float(raw.statistic) if kind=='same_spot' else float(spearmanr(np.round(b.obs.PGAM5_log1p_CP10k,12),np.round(yy,12)).statistic),'ordinary_independence_p_not_primary':float(raw.pvalue) if kind=='same_spot' else float(spearmanr(np.round(b.obs.PGAM5_log1p_CP10k,12),np.round(yy,12)).pvalue),'QC_adjusted_partial_rank_rho':rho,'spatial_spectral_two_sided_p':float(p),'n_null':N_NULL,'null_sd':float(null.std())})
  stats_maps[k]=wr@y
 np.savez_compressed(ROOT/f'{s}_spectral_null_coefficients.npz',products=np.asarray(null_products),denominators=np.asarray(null_denominators),labels=np.asarray(null_labels),seed=np.asarray(SEED+(s=='CHC23')),n_null=np.asarray(N_NULL),n_spots=np.asarray(len(b)))
 # Accessible graph lists and paired original values allow independent recomputation.
 row,col=w0.nonzero();pd.DataFrame({'barcode_from':a.obs_names[row],'barcode_to':a.obs_names[col]}).to_csv(ROOT/f'{s}_spatial_edges.csv.gz',index=False)
 a.obs['spatial_degree']=degree0
 for k in PROGRAMS:
  dest=np.full(len(a),np.nan);dest[keep]=stats_maps[k];a.obs['neighbor_'+k]=dest
 a.obs.to_csv(ROOT/f'{s}_spot_annotations.csv',index_label='barcode',float_format='%.17g')
 return out,mrows,numerical,{'sample':s,'n_graph_spots':int(keep.sum()),'isolated_spots_excluded_from_statistics':int((~keep).sum()),'undirected_edges':int(w.nnz/2),'mean_neighbors':float(degree.mean()),'min_neighbors':int(degree.min()),'max_neighbors':int(degree.max()),'spectral_residual_moran':moran(x,w)}

def sensitivity():
 records=[]
 for s,g in SAMPLES.items():
  original=read(s,g)
  for scenario,col in [('primary_mt25','pass_qc'),('relaxed_mt30','pass_qc_relaxed_mt30')]:
   a=original[original.obs[col]].copy();sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
   x=array(a[:,['PGAM5']].X).ravel();cov=covariance(a.obs);xr=residual_rank(x,cov)
   variants=dict(PROGRAMS);variants['Macrophage_C1q_only']=['C1QA','C1QB','C1QC'];variants['T_cell_CD3_only']=['CD3D','CD3E','CD3G']
   for k,genes in variants.items():
    genes=[g for g in genes if g in a.var_names];y=array(a[:,genes].X).mean(1)
    records.append({'sample':s,'QC_scenario':scenario,'program':k,'n_spots':len(a),'raw_Spearman_rho':float(spearmanr(np.round(x,12),np.round(y,12)).statistic),'QC_adjusted_partial_rank_rho':corr(xr,residual_rank(y,cov))})
 return pd.DataFrame(records)
if __name__=='__main__':
 allrows=[];allm=[];checks=[];graphs=[]
 for s in SAMPLES:
  a=sc.read_h5ad(WORK/f'{s}_full_processed.h5ad');rows,m,ck,g=analyze(a,s);allrows+=rows;allm+=m;checks+=ck;graphs.append(g)
 df=pd.DataFrame(allrows);df['spatial_BH_q']=multipletests(df.spatial_spectral_two_sided_p,method='fdr_bh')[1];df.to_csv(ROOT/'PGAM5_spatial_associations.csv',index=False)
 df=pd.DataFrame(allm);df['Moran_BH_q']=multipletests(df.two_sided_label_permutation_p,method='fdr_bh')[1];df.to_csv(ROOT/'spatial_autocorrelation.csv',index=False)
 pd.DataFrame(checks).to_csv(ROOT/'spectral_null_numerical_checks.csv',index=False);pd.DataFrame(graphs).to_csv(ROOT/'spatial_graph_summary.csv',index=False)
 sensitivity().to_csv(ROOT/'QC_and_marker_panel_sensitivity.csv',index=False)
 (ROOT/'statistics_settings.json').write_text(json.dumps({'seed':SEED,'null_replicates':N_NULL,'numerical_ties':'Expression and program values rounded to 12 decimal places before ranking; float64 CSV serialization preserves source values.','multiple_testing_family':'all 48 associations (2 slices x 12 programs x 2 relationships)','covariates':'rank(log1p(total UMI)), rank(n_genes), rank(percent mitochondrial RNA), intercept','null_model':'independent random sign flips of coefficients in the complete eigenbasis of H W H, W symmetric binary hex first-neighbor adjacency. Preserves target residual graph spectral power, centered norm and Moran I. This is a model-based exploratory null, not an exact distribution-free patient-level test.','expression':'PGAM5 log1p(CP10k) including zeros; programs average gene log1p(CP10k)','neighbor':'arithmetic mean over observed QC-passing first-ring hex neighbors; zero-neighbor spots excluded','not_estimated':'cell-type fractions, validated PGAM5-positive macrophage abundance, cell-cell ligand-receptor interactions'},indent=2))
 print(df.to_string(index=False),flush=True)
