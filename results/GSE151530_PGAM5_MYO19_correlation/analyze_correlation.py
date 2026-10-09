"""Pooled cell correlation, fixed sensitivity analyses, and plots from saved counts."""
from pathlib import Path
import json
import platform
import importlib.metadata as md
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr, rankdata
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R = Path(__file__).resolve().parent
B = 200000
protocol = {'primary':'All3493existing HCC macrophages, Spearman on log1p(CP10k), zeros retained, two-sided scipy asymptotic p plus fixed200000permutation test using absolute rho','secondary_family':['Pearson_all_cells_log1p_CP10k','Spearman_PGAM5_detected_log1p_CP10k','Spearman_both_detected_log1p_CP10k','Spearman_both_detected_raw_counts','Spearman_donor_proxy_mean_CP10k'],'multiple_testing':'One primary gene-pair Spearman test. BH across the five specified secondary asymptotic p-values; co-detected normalized/raw comparisons are diagnostic, not independent validation','permutation':'Fixed200000shuffles of MYO19 ranks across all cells; MonteCarlo p=(absolute-statistic extreme count+1)/(B+1); seed20261009. This cell-level null assumes exchangeable cells, not independent patients','group_aggregation':'Average linear cell CP10k per pre-existing donor_id proxy; all20labels included. Existing suffix a/b/c libraries share donor proxy. No treating27libraries as27independent patients; no expression-selected size threshold','normalization':'Same full-library CP10k and natural log1p as previous analysis; no depth matching, covariate regression or partial correlation','limitations':'Gene expression association is not causal or evidence of protein interaction. Many zeros, low counts, common normalization denominator, cell dependence, unequal sample cell counts and donor proxies affect interpretation. The 22co-detected cells are a selected subset.'}
(R/'protocol.json').write_text(json.dumps(protocol,indent=2))
d = pd.read_csv(R/'macrophage_expression.csv.gz')
assert len(d) == 3493 and d.cell_id.is_unique
x = d.PGAM5_log1p_CP10k.to_numpy()
y = d.MYO19_log1p_CP10k.to_numpy()
primary = spearmanr(x,y)
xc = rankdata(x)-((len(x)+1)/2)
yc = rankdata(y)-((len(y)+1)/2)
den = np.linalg.norm(xc)*np.linalg.norm(yc)
assert np.isclose(np.dot(xc,yc)/den,primary.statistic)
rng = np.random.default_rng(20261009)
extreme = 0
for j in range(B):
    rho = np.dot(xc,rng.permutation(yc))/den
    extreme += abs(rho) >= abs(primary.statistic)-1e-12
perm_p = (extreme+1)/(B+1)
rows = [{'analysis':'Spearman_all_cells_log1p_CP10k','unit':'cell','n':len(d),'correlation':float(primary.statistic),'asymptotic_p':float(primary.pvalue),'permutation_p':perm_p,'permutations':B,'extreme_permutations':int(extreme),'family':'primary','secondary_BH_q':None}]
positive = d.PGAM5_raw_count.gt(0)
both = positive & d.MYO19_raw_count.gt(0)
donors = d.groupby('donor_id',sort=True).agg(macrophages=('cell_id','size'),sample_libraries=('Sample','nunique'),PGAM5_mean_CP10k=('PGAM5_CP10k','mean'),MYO19_mean_CP10k=('MYO19_CP10k','mean'),PGAM5_detected_cells=('PGAM5_detected','sum'),MYO19_detected_cells=('MYO19_detected','sum'))
donors.to_csv(R/'donor_proxy_mean_expression.csv')
specifications = [
    ('Pearson_all_cells_log1p_CP10k','cell',x,y,pearsonr),
    ('Spearman_PGAM5_detected_log1p_CP10k','selected cell',x[positive],y[positive],spearmanr),
    ('Spearman_both_detected_log1p_CP10k','selected cell',x[both],y[both],spearmanr),
    ('Spearman_both_detected_raw_counts','selected cell',d.loc[both,'PGAM5_raw_count'],d.loc[both,'MYO19_raw_count'],spearmanr),
    ('Spearman_donor_proxy_mean_CP10k','donor proxy',donors.PGAM5_mean_CP10k,donors.MYO19_mean_CP10k,spearmanr),
]
secondary = []
for name,unit,xx,yy,method in specifications:
    result = method(xx,yy)
    secondary.append({'analysis':name,'unit':unit,'n':len(xx),'correlation':float(result.statistic),'asymptotic_p':float(result.pvalue),'permutation_p':None,'permutations':None,'extreme_permutations':None,'family':'secondary'})
q = multipletests([v['asymptotic_p'] for v in secondary],method='fdr_bh')[1]
for row,qv in zip(secondary,q): row['secondary_BH_q'] = float(qv)
rows.extend(secondary)
statistics = pd.DataFrame(rows)
statistics.to_csv(R/'correlation_statistics.csv',index=False)
counts = {'cells':len(d),'PGAM5_detected':int(positive.sum()),'MYO19_detected':int(d.MYO19_detected.sum()),'both_detected':int(both.sum()),'neither_detected':int((~positive & ~d.MYO19_detected).sum()),'PGAM5_only':int((positive & ~d.MYO19_detected).sum()),'MYO19_only':int((~positive & d.MYO19_detected).sum())}
(R/'detection_summary.json').write_text(json.dumps(counts,indent=2))
fig,axs = plt.subplots(1,3,figsize=(15,4.8),constrained_layout=True)
ax = axs[0]
ax.scatter(x,y,s=13,alpha=.25,color='#538aa2',edgecolors='none')
ax.scatter(x[both],y[both],s=25,alpha=.8,color='#a36837',edgecolors='none',label='Both detected (n=22)')
ax.legend(loc='upper right',fontsize=8,frameon=False)
ax.set_title(f'All HCC macrophages (n={len(d):,})\nSpearman rho={primary.statistic:.3f}\nasymptotic p={primary.pvalue:.3g}; permutation p={perm_p:.3g}',fontsize=10)
for i,name in [(1,'Spearman_both_detected_log1p_CP10k')]:
    z = statistics[statistics.analysis.eq(name)].iloc[0]
    axs[i].scatter(x[both],y[both],s=35,color='#a36837',alpha=.8,edgecolors='none')
    axs[i].set_title(f'Both genes detected (n={int(both.sum())})\nSpearman rho={z.correlation:.3f}, p={z.asymptotic_p:.3g}\nSelected-subset sensitivity; BH q={z.secondary_BH_q:.3g}',fontsize=10)
z = statistics[statistics.analysis.eq('Spearman_donor_proxy_mean_CP10k')].iloc[0]
dx = np.log1p(donors.PGAM5_mean_CP10k)
dy = np.log1p(donors.MYO19_mean_CP10k)
axs[2].scatter(dx,dy,s=35,color='#4c9273',alpha=.8,edgecolors='none')
axs[2].set_title(f'Mean expression per donor proxy (n={len(donors)})\nSpearman rho={z.correlation:.3f}, p={z.asymptotic_p:.3g}\nBH q={z.secondary_BH_q:.3g}',fontsize=10)
for i,ax in enumerate(axs):
    suffix = 'mean CP10k' if i == 2 else 'CP10k'
    ax.set_xlabel(f'PGAM5 log1p({suffix})')
    ax.set_ylabel(f'MYO19 log1p({suffix})')
    ax.spines[['top','right']].set_visible(False)
    ax.grid(alpha=.12)
fig.suptitle('GSE151530: PGAM5 and MYO19 RNA expression in HCC macrophages',fontsize=13)
fig.savefig(R/'PGAM5_MYO19_correlation.png',dpi=220)
fig.savefig(R/'PGAM5_MYO19_correlation.pdf')
plt.close(fig)
versions = {'python':platform.python_version(),**{k:md.version(k) for k in ['numpy','pandas','scipy','statsmodels','matplotlib','anndata','h5py']}}
(R/'environment_versions.json').write_text(json.dumps(versions,indent=2))
(R/'requirements.txt').write_text('\n'.join(k+'=='+v for k,v in versions.items() if k != 'python')+'\n')
result = {'scope':'Existing GSE151530 HCC macrophages','n_cells':len(d),'n_donor_proxy_labels':len(donors),'Spearman_rho':float(primary.statistic),'Spearman_asymptotic_p':float(primary.pvalue),'Spearman_200000_permutation_p':perm_p,'permutation_extremes':int(extreme),'detection_counts':counts,'interpretation':'Very weak pooled-cell positive association, not strong gene coexpression or a mechanistic/protein conclusion. Sample labels are donor proxies; cell-level p-values do not establish independent patient-level association.'}
(R/'result.json').write_text(json.dumps(result,indent=2))
from plot_main_correlation import render
render(R)
print(statistics.to_string(index=False))
print(json.dumps(result,indent=2))
