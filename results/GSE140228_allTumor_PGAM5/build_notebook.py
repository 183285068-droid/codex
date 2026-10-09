from pathlib import Path
import nbformat as nbf
from nbclient import NotebookClient
R=Path(__file__).resolve().parent
nb=nbf.v4.new_notebook();cells=[]
def md(s):cells.append(nbf.v4.new_markdown_cell(s))
def code(s):cells.append(nbf.v4.new_code_cell(s))
md('''# GSE140228: PGAM5-associated macrophage differential expression in all Tumor samples

Scope is literal author `Tissue == Tumor`, including every Tumor sample and all six original `celltype_sub` macrophage labels beginning `Mφ-`. Normal, Blood, Ascites and Lymphnode are excluded. Droplet and Smart-seq2 are analyzed separately because their matrices contain UMI counts and read counts, respectively. **Author Histology labels include CC in Droplet; this is an all-Tumor analysis, not an HCC-only analysis.** Counts of samples must not be interpreted as counts of independent patients; a donor may provide several specimens and occur in both technologies.

## Method and assumptions

PGAM5 raw count >0 defines detected cells; count=0 means RNA nondetection, not protein negativity. QC retains cells with >=500 detected genes and <20% mitochondrial counts. Counts are normalized to 10,000 per cell and log1p transformed. Within each platform, all retained Tumor macrophages are pooled for two-group Wilcoxon tests with tie correction and BH across all genes. No patient pairing, depth matching or sample/depth covariates are used, following the requested analysis.

Candidate thresholds: **FDR <0.05, |Scanpy approximate log2FC| >=1, detection in >=10% of either group**. PGAM5 remains in full DE and threshold-passing tables but is excluded from up/down candidate tables. Ribosomal and mitochondrial genes remain eligible. The approximate log2FC is the ratio of back-transformed mean log expression, not the ratio of arithmetic means. Large fold changes against near-zero controls can reflect pseudocount sensitivity.

These are exploratory associations. Cell-level FDR does not establish reproducibility across patients. Depth, histology and macrophage subtype differences can contribute to the results because the requested tests do not adjust for them. No causal mechanism, validated signature, independent cross-platform replication, ambient-RNA correction or doublet removal is claimed.

Source: [GEO series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140228), [official supplementary files](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140228/suppl/); PMID 31675496. Input filenames and SHA256 hashes are retained in input_sha256.json. Barcode/metadata and count/gene row alignments are asserted in analyze_tumor.py.

The analysis script has been executed on source counts. This notebook independently verifies selected outputs against retained raw macrophage counts, then regenerates the figures.
''')
code('''from pathlib import Path
import json
import numpy as np,pandas as pd,anndata as ad
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from IPython.display import display
R=Path.cwd();summary={k:v for k,v in json.loads((R/'analysis_summary.json').read_text()).items() if k in ['Droplet','Smartseq2']}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':120,'savefig.bbox':'tight'})
rows=[];checks=[]
for tech,s in summary.items():
 p=R/tech;a=ad.read_h5ad(p/'macrophage_counts_QC.h5ad');r=pd.read_csv(p/'full_DE.csv',index_col=0)
 pos=a.obs.PGAM5_counts.to_numpy()>0
 assert a.obs.Tissue.eq('Tumor').all() and a.obs.celltype_sub.str.startswith('Mφ-').all()
 assert a.obs.n_genes.ge(500).all() and a.obs.pct_mt.lt(20).all()
 assert a.obs_names.is_unique and a.var_names.is_unique
 assert int(pos.sum())==s['PGAM5_positive'] and int((~pos).sum())==s['PGAM5_undetected']
 pg=np.flatnonzero(a.var.gene.eq('PGAM5'));assert len(pg)==1 and np.array_equal(a.X[:,pg[0]].toarray().ravel(),a.obs.PGAM5_counts)
 q=multipletests(r.p_value,method='fdr_bh')[1];assert np.allclose(q,r.FDR,atol=1e-12)
 for direction,sign in [('upregulated',1),('downregulated',-1)]:
  c=pd.read_csv(p/(direction+'.csv'),index_col=0)
  expected=r[r.FDR.lt(.05)&(r.log2FC*sign).ge(1)&r[['fraction_PGAM5_positive','fraction_PGAM5_undetected']].max(axis=1).ge(.1)&~r.gene.eq('PGAM5')]
  assert set(c.index)==set(expected.index)
 # Independent SciPy rank test on five strongest nondefining candidates, no continuity correction (Scanpy convention).
 ids=pd.read_csv(p/'upregulated.csv',index_col=0).index[:5]
 totals=np.asarray(a.X.sum(axis=1)).ravel()
 for gene_id in ids:
  j=a.var_names.get_loc(gene_id);counts=a.X[:,j].toarray().ravel()
  v=np.log1p(counts/totals*1e4)
  pv=mannwhitneyu(v[pos],v[~pos],alternative='two-sided',method='asymptotic',use_continuity=False).pvalue
  est=np.log2((np.expm1(v[pos].mean())+1e-9)/(np.expm1(v[~pos].mean())+1e-9))
  assert np.isclose(pv,r.loc[gene_id,'p_value'],rtol=.002,atol=1e-12)
  assert np.isclose(est,r.loc[gene_id,'log2FC'],atol=.002)
  assert np.isclose((counts[pos]>0).mean(),r.loc[gene_id,'fraction_PGAM5_positive'])
  checks.append({'platform':tech,'gene':r.loc[gene_id,'gene'],'SciPy_two_sided_p':pv,'Scanpy_p':r.loc[gene_id,'p_value'],'recomputed_log2FC':est,'reported_log2FC':r.loc[gene_id,'log2FC']})
 rows.append({'platform':tech,'Tumor samples':s['tumor_samples'],'Tumor donors':s['tumor_donors'],'Tumor cells':s['tumor_cells'],'QC macrophages':s['QC_macrophages'],'PGAM5+':s['PGAM5_positive'],'PGAM5 undetected':s['PGAM5_undetected'],'up':s['upregulated_excluding_PGAM5'],'down':s['downregulated_excluding_PGAM5']})
 del a
pd.DataFrame(checks).to_csv(R/'independent_calculation_checks.csv',index=False)
display(pd.DataFrame(rows));display(pd.DataFrame(checks))
print('Scope, identities, PGAM5 counts, BH correction, complete candidate sets and selected rank-test/fold-change calculations verified.')
''')
md('''## Sample composition and PGAM5 detection

Rates use QC macrophages as the denominator. Labels include the author Histology annotation. Multiple specimens may share a donor; both platforms include D20171109. Sample tables retain counts and donor IDs for inspection.
''')
code('''fig,axes=plt.subplots(1,2,figsize=(15,6.2))
for ax,(tech,s) in zip(axes,summary.items()):
 t=pd.read_csv(R/tech/'PGAM5_by_sample.csv').sort_values('rate')
 labels=t.Sample+' ['+t.Histology+']';ax.barh(labels,t.rate*100,color='#355C8A')
 for i,(_,z) in enumerate(t.iterrows()):ax.text(z.rate*100+.4,i,f"{z.PGAM5_positive}/{z.macrophages}",va='center',fontsize=8)
 ax.set_xlim(0,max(5,t.rate.max()*100*1.30));ax.set_xlabel('PGAM5 RNA detection in QC macrophages (%)');ax.set_title(tech);ax.grid(axis='x',alpha=.15)
fig.suptitle('GSE140228 all Tumor: PGAM5 detection by sample',fontsize=13);fig.tight_layout();fig.savefig(R/'PGAM5_detection_by_sample.png');fig.savefig(R/'PGAM5_detection_by_sample.pdf');plt.show()
for tech in summary:
 print(tech);display(pd.read_csv(R/tech/'PGAM5_by_sample.csv'));print('QC histology:',summary[tech]['QC_macrophage_histology'])
''')
md('''## Differential expression

The volcano view includes genes detected in >=10% of either group. Orange marks upregulated candidates, blue marks downregulated candidates, and gray marks the remainder. PGAM5 is not plotted because it defines the comparison. Candidate counts and top genes are shown separately for each technology.
''')
code('''fig,axes=plt.subplots(1,2,figsize=(12,4.7))
for ax,(tech,s) in zip(axes,summary.items()):
 r=pd.read_csv(R/tech/'full_DE.csv',index_col=0);r=r[r[['fraction_PGAM5_positive','fraction_PGAM5_undetected']].max(axis=1).ge(.1)&~r.gene.eq('PGAM5')]
 color=np.where(r.FDR.lt(.05)&r.log2FC.ge(1),'#C68D29',np.where(r.FDR.lt(.05)&r.log2FC.le(-1),'#355C8A','#B8BCC2'))
 ax.scatter(r.log2FC,-np.log10(r.FDR.clip(lower=1e-300)),c=color,s=7,alpha=.6)
 ax.axhline(-np.log10(.05),ls='--',color='#777',lw=.8)
 for x in [-1,1]:ax.axvline(x,ls='--',color='#777',lw=.8)
 ax.set_xlabel('Approximate log2FC (PGAM5 detected / undetected)');ax.set_ylabel('-log10(BH FDR)');ax.set_title(f"{tech}: {s['upregulated_excluding_PGAM5']} up, {s['downregulated_excluding_PGAM5']} down")
fig.tight_layout();fig.savefig(R/'differential_expression.png');fig.savefig(R/'differential_expression.pdf');plt.show()
for tech in summary:
 print(tech,'top up candidates');display(pd.read_csv(R/tech/'top20_upregulated.csv')[['gene','log2FC','FDR','fraction_PGAM5_positive','fraction_PGAM5_undetected']])
''')
md('''## Original macrophage annotation checks

Detection fractions of macrophage and possible nonmyeloid transcripts are reported without post hoc removal of cells. Author macrophage annotations define the selected population; monocyte and dendritic cell labels are excluded. Some macrophages have ALB, CD3D, EPCAM or MS4A1 transcripts; ambient RNA, doublets and phagocytosed material are not distinguished.
''')
code('''fig,axes=plt.subplots(1,2,figsize=(10,5))
markers=['C1QA','C1QB','C1QC','CD68','CSF1R','LST1','TYROBP','SPP1','ALB','EPCAM','CD3D','MS4A1']
for ax,tech in zip(axes,summary):
 q=pd.read_csv(R/tech/'marker_QA.csv').pivot(index='gene',columns='group',values='fraction').reindex(markers)
 im=ax.imshow(q.to_numpy(),vmin=0,vmax=1,cmap='Blues',aspect='auto');ax.set_yticks(range(len(q)),q.index);ax.set_xticks([0,1],['Detected','Undetected']);ax.set_title(tech)
 for i in range(len(q)):
  for j in range(2):
   z=q.iloc[i,j];ax.text(j,i,f'{z:.2f}',ha='center',va='center',color='white' if z>.6 else 'black',fontsize=8)
fig.subplots_adjust(right=.88,wspace=.35);cb=fig.add_axes([.91,.2,.02,.6]);fig.colorbar(im,cax=cb,label='Gene detection fraction');fig.savefig(R/'macrophage_marker_QA.png');fig.savefig(R/'macrophage_marker_QA.pdf');plt.show()
''')
md('''## Interpretation and reproducibility

Upregulated genes are candidate associations with PGAM5 RNA detection under the requested pooled-cell model. The all-Tumor Droplet result mixes HCC and CC histology and should not be labelled HCC-specific. Technology separation avoids mixing UMI and read count distributions but does not make platforms independent patient cohorts. A candidate intersection is a descriptive overlap, not independent validation.

Download the seven official .gz files into the script directory and run `python analyze_tumor.py` using versions listed in environment_versions.txt. The script locates inputs relative to itself. Run `python build_notebook.py` to execute this notebook and regenerate plots from outputs. All filtering, annotation, summary and calculation checks are explicit in the scripts and notebook. Raw source matrices and QC count matrices are not copied into the downloadable result ZIP. Input hashes and reproduction scripts are included; rerunning analyze_tumor.py generates the QC count matrices needed for notebook checks. Executed notebook outputs retain these checks.
''')
nb.cells=cells;nb.metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},'language_info':{'name':'python','version':'3.11'}}
NotebookClient(nb,timeout=600,kernel_name='python3',resources={'metadata':{'path':str(R)}}).execute()
nbf.write(nb,R/'GSE140228_allTumor_PGAM5_analysis.ipynb')
print('Notebook executed successfully',flush=True)
