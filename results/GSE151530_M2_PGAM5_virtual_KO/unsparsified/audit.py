from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parent;checks=[]
def check(name,value):
 checks.append(dict(check=name,passed=bool(value)));assert value,name
a=ad.read_h5ad(R/'inputs/M2_raw.h5ad');x=pd.read_csv(R/'inputs/network_raw_counts.csv.gz',index_col=0);m=pd.read_csv(R/'inputs/target_cell_metadata.csv.gz',index_col=0);c=pd.read_csv(R/'inputs/macrophage_context.csv.gz',index_col=0);protocol=json.loads((R/'protocol.json').read_text())
check('M2_exact_context_cell_ids',a.obs_names.equals(pd.Index(c.index[c.mac_cluster.eq('M2')])))
check('M2_all_171_not_only_detected',len(a)==171 and a.obs.mac_cluster.astype(str).eq('M2').all())
check('full_library_totals_exact',np.array_equal(np.asarray(a.X.sum(1)).ravel(),m.total_counts))
check('input_cell_order_exact',x.columns.equals(a.obs_names) and m.index.equals(a.obs_names))
check('all_selected_raw_counts_exact',np.array_equal(a[:,x.index].X.toarray().T,x.to_numpy()))
check('PGAM5_22_detected_and_149_zeros',int(x.loc['PGAM5'].gt(0).sum())==22 and int(x.loc['PGAM5'].eq(0).sum())==149)
check('PGAM5_in_model_input',('PGAM5' in x.index) and len(x)==1501)
folders=[(str(seed),R/'models'/str(seed)) for seed in protocol['seeds']]+([('high_rank_151530',R/'high_rank_pilot/151530')] if (R/'high_rank_pilot/151530/metrics.json').exists() else [])
for seed,f in folders:
 n=np.load(f/'WT_KO_networks.npz');wt=n['WT'];genes=pd.Index(n['genes']);j=genes.get_loc('PGAM5');rest=np.arange(len(genes))!=j
 if 'KO' in n:ko=n['KO']
 else:
  ko=wt.copy();ko[j]=0;check(str(seed)+'_KO_delta_exact',np.array_equal(n['KO_delta'],ko-wt))
 if 'KO_SHA256' in n:check(str(seed)+'_stored_KO_exact_hash',hashlib.sha256(ko.tobytes()).hexdigest()==str(n['KO_SHA256']))
 check(str(seed)+'_network_genes_order',genes.equals(x.index));check(str(seed)+'_network_finite',np.isfinite(wt).all() and np.isfinite(ko).all());check(str(seed)+'_KO_only_PGAM5_row_zero',np.all(ko[j]==0) and np.array_equal(wt[rest],ko[rest]));check(str(seed)+'_WT_diagonal_zero',np.all(np.diag(wt)==0))
 alignment=pd.read_csv(f/'manifold_alignment.csv.gz',index_col=0);distance=np.linalg.norm(alignment.iloc[:len(genes)].to_numpy()-alignment.iloc[len(genes):].to_numpy(),axis=1);d=pd.read_csv(f/'differential_regulation.csv').set_index('Gene').reindex(genes)
 check(str(seed)+'_saved_alignment_matches_distance',np.allclose(d.Distance,distance,rtol=1e-8,atol=1e-14));check(str(seed)+'_statistics_finite',np.isfinite(d[['Distance','p-value','adjusted p-value']]).all().all());check(str(seed)+'_p_and_q_valid_range',d[['p-value','adjusted p-value']].ge(0).all().all() and d[['p-value','adjusted p-value']].le(1).all().all())
 expected=(d.Distance[rest]**2).mean()
 with np.errstate(divide='ignore',invalid='ignore'):fc=d.Distance**2/expected;p=chi2.sf(fc,df=1)
 noise=np.sqrt(np.finfo(float).eps)*np.abs(alignment.to_numpy()).max();p[d.Distance.le(noise)]=1
 check(str(seed)+'_chi_squared_p_independent',np.allclose(d['p-value'],p,rtol=1e-7,atol=1e-14));check(str(seed)+'_BH_q_independent',np.allclose(d['adjusted p-value'],multipletests(p,method='fdr_bh')[1],rtol=1e-7,atol=1e-14))
 null=pd.read_csv(f/'no_KO_control.csv');check(str(seed)+'_null_all_genes_present',len(null)==len(x) and set(null.Gene)==set(x.index));check(str(seed)+'_null_no_significant_genes',null['adjusted p-value'].ge(.05).all());check(str(seed)+'_null_all_p_equal_one',null['p-value'].eq(1).all());resid=np.load(f/'tensor_residuals.npy');check(str(seed)+'_tensor_residuals_finite',len(resid)>0 and np.isfinite(resid).all())
allres={s:pd.read_csv(R/'models'/str(s)/'differential_regulation.csv').set_index('Gene') for s in protocol['seeds']};dist=pd.DataFrame({s:v.Distance for s,v in allres.items()});qs=pd.DataFrame({s:v['adjusted p-value'] for s,v in allres.items()});out=pd.read_csv(R/'all_genes_perturbation_summary.csv').set_index('gene')
check('consensus_medians_exact',np.allclose(out.median_distance.reindex(dist.index),dist.median(1)))
check('consensus_significance_frequency_exact',out.significant_seed_count.reindex(dist.index).equals(qs.lt(.05).sum(1)))
top=pd.read_csv(R/'TOP20_downstream_genes.csv');thresholds={}
for seed in protocol['seeds']:
 alignment=pd.read_csv(R/'models'/str(seed)/'manifold_alignment.csv.gz',index_col=0);thresholds[seed]=np.sqrt(np.finfo(float).eps)*np.abs(alignment.to_numpy()).max()
above=dist.gt(pd.Series(thresholds),axis='columns').sum(1);valid=above.drop(index='PGAM5').gt(0);expected=dist.drop(index='PGAM5').median(1)[valid].rename_axis('gene').rename('distance').reset_index().sort_values(['distance','gene'],ascending=[False,True]).head(20).gene.tolist();check('top20_excludes_numerical_noise_and_target',top.gene.tolist()==expected and 'PGAM5' not in set(top.gene) and 'PGAM5' in set(out.index));check('consensus_noise_filter_exact',np.array_equal(out.above_numerical_noise_seed_count.reindex(above.index),above))
pd.DataFrame(checks).to_csv(R/'independent_audit_checks.csv',index=False);summary=dict(checks=len(checks),passed=sum(r['passed'] for r in checks),input_cells=171,detected_PGAM5=22,main_models=len(protocol['seeds']),high_rank_pilots=len(folders)-len(protocol['seeds']),no_KO_controls=len(folders));(R/'audit_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
