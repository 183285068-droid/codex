from pathlib import Path
import json,shutil,numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr,fisher_exact
from statsmodels.stats.multitest import multipletests
from scTenifold.core._networks import manifold_alignment,d_regulation
from threadpoolctl import threadpool_limits
R=Path(__file__).resolve().parent
x=pd.read_csv(R/'differential_regulation_seed202642.csv').set_index('Gene');y=pd.read_csv(R/'differential_regulation_seed42.csv').set_index('Gene');merged=x[['Distance','p-value','adjusted p-value']].join(y[['Distance','p-value','adjusted p-value']],rsuffix='_seed42');merged['significant_primary']=merged['adjusted p-value']<.05;merged['significant_seed42']=merged['adjusted p-value_seed42']<.05;merged['robust_both_seeds']=merged.significant_primary&merged.significant_seed42;merged.to_csv(R/'seed_comparison.csv')
sets=json.loads((R/'program_definitions.json').read_text());sets['Heat_shock_proteostasis']=['PGAM5','HSP90AA1','HSP90AB1','HSPA1A','HSPA1B','HSPA5','HSPA8','DNAJB1','HSPH1','HSPB1','HSPD1','BAG3','PPP1R15A'];sets['Immediate_response']=['PGAM5','FOS','FOSB','JUN','JUNB','DUSP1','NR4A1','NR4A2','ATF3','NFKBIA','EGR1','IER2','IER3'];(R/'program_definitions.json').write_text(json.dumps(sets,indent=2))
# Keep PGAM5 in all original programs, but automatic target perturbation is not downstream evidence.
bg=set(x.index)-{'PGAM5'};hit=set(x.index[(x['adjusted p-value']<.05)])-{'PGAM5'};hit2=set(y.index[(y['adjusted p-value']<.05)])-{'PGAM5'};rows=[]
for term,gs in sets.items():
 g=set(gs)&bg;h=g&hit
 if not g:continue
 _,p=fisher_exact([[len(h),len(g-h)],[len(hit-g),len(bg-g-hit)]],alternative='greater');rows.append(dict(program=term,tested_non_target_genes=len(g),primary_hits=len(h),seed42_hits=len(g&hit2),p_value=p,genes=';'.join(sorted(h)),PGAM5_in_program='PGAM5' in gs))
t=pd.DataFrame(rows);t['BH_q']=multipletests(t.p_value,method='fdr_bh')[1];t=t.sort_values('BH_q');t.to_csv(R/'exploratory_program_enrichment.csv',index=False)
wt=pd.read_csv(R/'WT_network_seed202642.csv.gz',index_col=0);ko=pd.read_csv(R/'KO_network_seed202642.csv.gz',index_col=0);assert (ko.loc['PGAM5']==0).all();other=wt.index!='PGAM5';np.testing.assert_array_equal(wt.loc[other],ko.loc[other]);np.random.seed(100)
with threadpool_limits(limits=4):null=d_regulation(manifold_alignment(wt,wt,d=30),ko_genes=['PGAM5'])
null.to_csv(R/'no_KO_control.csv',index=False);nullhits=int((null['adjusted p-value']<.05).sum());assert nullhits==0
meta=pd.read_csv(R/'target_cell_metadata.csv.gz',index_col=0);meta.groupby('sample_name').agg(cells=('PGAM5_counts','size'),PGAM5_detected=('PGAM5_counts',lambda v:(v>0).sum())).to_csv(R/'sample_summary.csv')
z=merged.drop(index='PGAM5');rho,p=spearmanr(z.Distance,z.Distance_seed42);top=x.drop(index='PGAM5').nlargest(20,'Distance');fig,ax=plt.subplots(figsize=(8,7));ax.barh(top.index[::-1],top.Distance.to_numpy()[::-1],color='#397bab');ax.set_xlabel('Network manifold displacement (unsigned)');ax.set_title('PGAM5 virtual KO: primary-run candidates\nNot reproduced as significant in seed 42');fig.tight_layout();fig.savefig(R/'primary_candidate_genes.png',dpi=180);fig.savefig(R/'primary_candidate_genes.pdf');plt.close(fig)
fig,ax=plt.subplots(figsize=(7,6));ax.scatter(z.Distance,z.Distance_seed42,s=12,alpha=.5);ax.set_xscale('symlog',linthresh=1e-12);ax.set_yscale('symlog',linthresh=1e-12);ax.set_xlabel('Distance: seed 202642');ax.set_ylabel('Distance: seed 42');ax.set_title(f'Seed sensitivity | downstream rho={rho:.3f}\nPrimary FDR<0.05: {len(hit)}; seed 42: {len(hit2)}');fig.tight_layout();fig.savefig(R/'seed_sensitivity.png',dpi=180);fig.savefig(R/'seed_sensitivity.pdf');plt.close(fig)
summary=dict(primary_non_target_significant=len(hit),seed42_non_target_significant=len(hit2),both_seeds_non_target_significant=len(hit&hit2),non_target_distance_Spearman=float(rho),no_KO_control_significant=nullhits,no_KO_max_distance=float(null.Distance.max()),PGAM5_WT_outgoing_edges_primary=int((wt.loc['PGAM5']!=0).sum()),PGAM5_WT_outgoing_edges_seed42=int((pd.read_csv(R/'WT_network_seed42.csv.gz',index_col=0).loc['PGAM5']!=0).sum()),WT_KO_other_rows_exactly_equal=True)
(R/'validation_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));print(t.head().to_string(index=False))
