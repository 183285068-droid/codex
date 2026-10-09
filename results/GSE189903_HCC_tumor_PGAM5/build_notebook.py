from pathlib import Path
import sys
import nbformat as nbf
from nbclient import NotebookClient
R=Path(__file__).resolve().parent
nb=nbf.v4.new_notebook();cells=[]
def md(s):cells.append(nbf.v4.new_markdown_cell(s))
def code(s):cells.append(nbf.v4.new_code_cell(s))
md('''# GSE189903: PGAM5 RNA detection in HCC tumor macrophages

## Scope and methods

User-selected scope: **Hepatocellular carcinoma, Tumor core or Tumor border only**. Adjacent non-tumor tissues and intrahepatic cholangiocarcinoma are excluded. Official GEO family SOFT cancer/tissue characteristics map to S_ID in Info.txt. Author TAM labels define macrophages. There are 16 eligible samples from four HCC patients; samples are separate regional specimens, not independent patients.

PGAM5 raw UMI count >0 defines detected cells, count=0 defines RNA nondetection, not protein negativity. QC retains >=500 detected genes and <20% mitochondrial counts. Counts are normalized to 10,000 per cell and log1p transformed. All selected QC TAMs are pooled for Wilcoxon with tie correction; BH adjustment is across all genes. Following the requested analysis, no patient pairing, depth matching or sample/depth covariates are used.

Thresholds: **FDR <0.05, |Scanpy approximate log2FC| >=1, >=10% detection in either group**. PGAM5 is retained in full and threshold-passing tables, and excluded from up/down candidate lists. Ribosomal and mitochondrial genes remain eligible. Approximate log2FC is the ratio of back-transformed mean log expression, not arithmetic mean expression; folds against near-zero controls may be pseudocount sensitive.

Candidates are exploratory associations. Cell-level FDR does not demonstrate reproducibility across patients. The requested pooled model does not adjust for tumor region, patient, subtype or sequencing depth. No validated signature, causal mechanism, ambient-RNA correction or doublet removal is claimed.

Sources: [GEO GSE189903](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189903), [official supplement directory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189903/suppl/), [family SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189903/soft/GSE189903_family.soft.gz). PMID 36476645.

analyze_HCC_tumor.py has been executed on the official counts. This notebook verifies saved results against raw QC TAM counts and regenerates scientific PNG/PDF figures.
''')
code('''from pathlib import Path
import json
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from IPython.display import display
R=Path.cwd();s=json.loads((R/'analysis_summary.json').read_text())
a=ad.read_h5ad(R/'HCC_tumor_TAM_counts_QC.h5ad')
r=pd.read_csv(R/'GSE189903_HCC_tumor_PGAM5_full_DE.csv',index_col=0)
u=pd.read_csv(R/'GSE189903_HCC_tumor_PGAM5_upregulated.csv',index_col=0)
d=pd.read_csv(R/'GSE189903_HCC_tumor_PGAM5_downregulated.csv',index_col=0)
m=pd.read_csv(R/'sample_diagnosis_tissue_mapping.csv')
assert int(m.included.sum())==16 and m[m.included].diagnosis.eq('Hepatocellular carcinoma').all()
assert m[m.included].tissue.isin(['Tumor core','Tumor border']).all()
assert a.obs.Type.eq('TAM').all() and a.obs.diagnosis.eq('Hepatocellular carcinoma').all()
assert a.obs.tissue.isin(['Tumor core','Tumor border']).all()
assert a.obs.n_genes.ge(500).all() and a.obs.pct_mt.lt(20).all()
assert a.obs_names.is_unique and a.var_names.is_unique
pg=np.flatnonzero(a.var.gene.eq('PGAM5'));assert len(pg)==1
pos=a.X[:,pg[0]].toarray().ravel()>0
assert np.array_equal(pos,a.obs.PGAM5_detected) and int(pos.sum())==s['PGAM5_positive_cells']
assert int((~pos).sum())==s['PGAM5_undetected_cells'] and len(a)==s['QC_TAMs']
assert np.allclose(multipletests(r.p_value,method='fdr_bh')[1],r.FDR,atol=1e-12)
eligible=r.FDR.lt(.05)&r[['fraction_PGAM5_positive','fraction_PGAM5_undetected']].max(axis=1).ge(.1)&~r.gene.eq('PGAM5')
assert set(u.index)==set(r.index[eligible&r.log2FC.ge(1)])
assert set(d.index)==set(r.index[eligible&r.log2FC.le(-1)])
assert len(u)==s['upregulated_excluding_PGAM5'] and len(d)==s['downregulated_excluding_PGAM5']
checks=[];totals=np.asarray(a.X.sum(axis=1)).ravel()
for gene_id in u.index[:10]:
 j=a.var_names.get_loc(gene_id);counts=a.X[:,j].toarray().ravel();v=np.log1p(counts/totals*1e4)
 pv=mannwhitneyu(v[pos],v[~pos],alternative='two-sided',method='asymptotic',use_continuity=False).pvalue
 est=np.log2((np.expm1(v[pos].mean())+1e-9)/(np.expm1(v[~pos].mean())+1e-9))
 assert np.isclose(pv,r.loc[gene_id,'p_value'],rtol=.002,atol=1e-12)
 assert np.isclose(est,r.loc[gene_id,'log2FC'],atol=.002)
 assert np.isclose((counts[pos]>0).mean(),r.loc[gene_id,'fraction_PGAM5_positive'])
 checks.append({'gene':r.loc[gene_id,'gene'],'SciPy_two_sided_p':pv,'Scanpy_p':r.loc[gene_id,'p_value'],'recomputed_log2FC':est,'Scanpy_log2FC':r.loc[gene_id,'log2FC']})
pd.DataFrame(checks).to_csv(R/'independent_calculation_checks.csv',index=False)
print('Scope, raw grouping, QC, all candidate memberships, BH FDR and ten selected independent rank tests/fold changes verified.')
display(pd.DataFrame(checks))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':120,'savefig.bbox':'tight'})
''')
md('''## Cohort and per-sample support

All reported gene tests pool cells across selected samples; descriptive patient/region tables do not introduce patient-paired tests. PGAM5 detection rates have QC TAMs as their denominators.
''')
code('''display(pd.Series(s,name='value').to_frame())
display(m[m.included][['S_ID','Sample','patient','tissue']])
samples=pd.read_csv(R/'PGAM5_by_sample.csv');display(samples)
display(pd.read_csv(R/'PGAM5_by_patient_tissue.csv'))
fig,ax=plt.subplots(figsize=(10,7))
t=samples.sort_values(['patient','Sample']).reset_index(drop=True)
labels=t.Sample+' ['+t.tissue.str.replace('Tumor ','',regex=False)+']'
ax.barh(labels,t.rate*100,color=np.where(t.tissue.eq('Tumor core'),'#355C8A','#C68D29'))
for i,z in t.iterrows():ax.text(z.rate*100+.15,i,f'{z.PGAM5_positive}/{z.TAMs}',va='center',fontsize=8)
ax.set_xlim(0,max(5,t.rate.max()*100*1.25));ax.set_xlabel('PGAM5 RNA detection in QC TAMs (%)')
ax.set_title('GSE189903 HCC: tumor core and border (16 samples, 4 patients)')
ax.grid(axis='x',alpha=.15);fig.tight_layout()
fig.savefig(R/'PGAM5_detection_by_sample.png');fig.savefig(R/'PGAM5_detection_by_sample.pdf');plt.show()
''')
md('''## Differential expression

The volcano includes genes detected in >=10% of either group and excludes PGAM5. Orange marks upregulated candidates, blue marks downregulated candidates, and gray marks the remainder. Dashed lines indicate the FDR and fold thresholds.
''')
code('''q=r[r[['fraction_PGAM5_positive','fraction_PGAM5_undetected']].max(axis=1).ge(.1)&~r.gene.eq('PGAM5')]
color=np.where(q.FDR.lt(.05)&q.log2FC.ge(1),'#C68D29',np.where(q.FDR.lt(.05)&q.log2FC.le(-1),'#355C8A','#B8BCC2'))
fig,ax=plt.subplots(figsize=(8,5))
ax.scatter(q.log2FC,-np.log10(q.FDR.clip(lower=1e-300)),c=color,s=9,alpha=.65)
ax.axhline(-np.log10(.05),ls='--',color='#777',lw=.8)
for x in [-1,1]:ax.axvline(x,ls='--',color='#777',lw=.8)
ax.set_xlabel('Approximate log2FC (PGAM5 detected / undetected)');ax.set_ylabel('-log10(BH FDR)')
ax.set_title(f'GSE189903 HCC tumor TAMs: {len(u)} up, {len(d)} down')
fig.tight_layout();fig.savefig(R/'differential_expression.png');fig.savefig(R/'differential_expression.pdf');plt.show()
display(u[['gene','log2FC','FDR','fraction_PGAM5_positive','fraction_PGAM5_undetected']].head(20))
display(d[['gene','log2FC','FDR','fraction_PGAM5_positive','fraction_PGAM5_undetected']].head(20))
''')
md('''## Macrophage annotation checks

Original TAM labels define the analyzed population. Detection fractions for macrophage and nonmyeloid transcripts are retained without post hoc filtering; these checks do not demonstrate that all cells are singlets or free of ambient RNA.
''')
code('''qa=pd.read_csv(R/'TAM_marker_QA.csv').pivot(index='gene',columns='group',values='fraction')
display(qa.round(3))
q=qa.reindex(['C1QA','C1QB','C1QC','CD68','CSF1R','LST1','TYROBP','SPP1','ALB','EPCAM','CD3D','MS4A1'])
fig,ax=plt.subplots(figsize=(5,5))
im=ax.imshow(q.to_numpy(),vmin=0,vmax=1,cmap='Blues',aspect='auto');ax.set_yticks(range(len(q)),q.index);ax.set_xticks([0,1],['Detected','Undetected'])
for i in range(len(q)):
 for j in range(2):
  z=q.iloc[i,j];ax.text(j,i,f'{z:.2f}',ha='center',va='center',color='white' if z>.6 else 'black',fontsize=9)
ax.set_title('GSE189903: original TAM annotation checks');fig.colorbar(im,ax=ax,label='Gene detection fraction')
fig.tight_layout();fig.savefig(R/'TAM_marker_QA.png');fig.savefig(R/'TAM_marker_QA.pdf');plt.show()
''')
md('''## Interpretation and reproduction

Reported genes are PGAM5 RNA detection-associated candidates under the requested model. They are not a validated macrophage-specific signature. Patient, region, subtype and depth confounding remain possible because the requested test does not adjust for them.

Download matrix, genes, barcodes, Info.txt and family SOFT from the official GEO URLs into the script directory. Install versions in environment_versions.txt, run python analyze_HCC_tumor.py, then python build_notebook.py. Scripts locate inputs relative to themselves. The downloadable ZIP includes tables, scripts, plots, input hashes and executed notebook outputs; full source matrices and QC count matrices remain in the cloud workspace and are regenerated by the analysis script.
''')
nb.cells=cells;nb.metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},'language_info':{'name':'python','version':sys.version.split()[0]}}
NotebookClient(nb,timeout=600,kernel_name='python3',resources={'metadata':{'path':str(R)}}).execute()
nbf.write(nb,R/'GSE189903_HCC_tumor_PGAM5_analysis.ipynb')
print('Notebook executed successfully',flush=True)
