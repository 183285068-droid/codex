from pathlib import Path
import nbformat as nb
from nbclient import NotebookClient
R=Path(__file__).resolve().parent
n=nb.v4.new_notebook();n.cells=[nb.v4.new_markdown_cell('# GSE202642: PGAM5 and immune RNA programs within the original 15 groups\nAll program scores retain PGAM5. This executed companion checks portable inputs, reproduces scoring and observed correlations, independently checks BH adjustment, and displays the permutation results. Run `python analyze.py` to rerun all 9,999 permutations; `python plot_report.py` regenerates the figures and report. These are exploratory RNA associations, not causal immune function tests.'),nb.v4.new_code_cell('''from pathlib import Path
import json, numpy as np, pandas as pd
from scipy.stats import spearmanr
from IPython.display import display, Image
R=Path.cwd()
a=pd.read_csv(R/'cell_metadata.csv.gz',index_col=0)
genes=pd.read_csv(R/'program_genes.csv').gene.tolist()
raw=np.load(R/'program_raw_counts.npz')['counts']
panels=json.loads((R/'program_definitions.json').read_text())
assert all(gs.count('PGAM5')==1 for gs in panels.values())
assert np.array_equal(raw[:,genes.index('PGAM5')],a.PGAM5_counts)
assert len(a)==12135 and int((a.PGAM5_counts>0).sum())==868
z=np.log1p(raw/a.total_counts.to_numpy()[:,None]*10000)
x=z[:,genes.index('PGAM5')]
scores=pd.DataFrame({p:z[:,[genes.index(g) for g in gs]].mean(axis=1) for p,gs in panels.items()},index=a.index)
stored=pd.read_csv(R/'cell_program_scores.csv.gz',index_col=0)
assert np.allclose(scores,stored)
display(pd.DataFrame({'program':list(panels),'genes_including_PGAM5':[len(gs) for gs in panels.values()]}))'''),nb.v4.new_markdown_cell('## Self-inclusion-preserving null\nFor each subgroup, permute PGAM5 X, rebuild S=(permuted X + unchanged other member sum)/m, and correlate permuted X with S. This preserves the mechanical shared-gene term. Report observed rho, null median, delta rho, two-sided equal-tail permutation p, and BH q across all 168 valid tests. C1 has no PGAM5 variation; C10 uses exact enumeration of its single detected cell. Null envelopes are not confidence intervals.'),nb.v4.new_code_cell('''results=pd.read_csv(R/'subgroup_program_associations.csv')
for row in results.itertuples():
    ix=np.flatnonzero(a.primary_cluster.to_numpy()==row.cluster)
    if (x[ix]>0).any():
        assert np.isclose(spearmanr(x[ix],scores[row.program].to_numpy()[ix]).statistic,row.rho)
    else:
        assert pd.isna(row.rho) and pd.isna(row.p_calibrated)
valid=results.p_calibrated.notna()
p=results.loc[valid,'p_calibrated'].to_numpy(); order=np.argsort(p)
sorted_q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
q=np.empty(len(p));q[order]=np.minimum(1,sorted_q)
assert np.allclose(q,results.loc[valid,'q_BH'])
assert valid.sum()==168 and len(results)==180
assert np.allclose(results.loc[valid,'delta_rho'],results.loc[valid,'rho']-results.loc[valid,'null_median'])
display(results.loc[(results.q_BH<.05)&(~results.mixed_RNA_flag)&(~results.low_information),['cluster','program','rho','null_median','delta_rho','p_calibrated','q_BH']].sort_values('q_BH').head(20))
print('Input, scoring, 168 observed correlations, BH and delta checks passed.')'''),nb.v4.new_code_cell("display(Image(filename=str(R/'self_inclusion_calibrated_associations.png')))\ndisplay(Image(filename=str(R/'observed_correlations.png')))") ,nb.v4.new_markdown_cell('## Interpretation\nPositive raw correlations can arise from including PGAM5 in the outcome. Interpret calibrated excess associations, with RNA program overlap, cell dependence, sparse detections, mixed RNA groups, and highly uneven patient composition in mind. Patient direction diagnostics are descriptive and subtract only the known shared-term covariance from the full score; they are not patient-level significance tests. No library-depth matching/regression was performed. No conclusion of enhanced or reduced antitumor immunity is established.')]
n.metadata={'kernelspec':{'display_name':'Python (PGAM5)','language':'python','name':'pgam5'},'language_info':{'name':'python','version':'3.12'}}
NotebookClient(n,timeout=180,kernel_name='pgam5',resources={'metadata':{'path':str(R)}}).execute();nb.write(n,R/'analysis.ipynb');print('Notebook executed')
