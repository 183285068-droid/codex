from pathlib import Path
import nbformat as nbf
from nbclient import NotebookClient
R=Path('/workspace/scratch/GSE149614');nb=nbf.v4.new_notebook();cells=[]
def md(s):cells.append(nbf.v4.new_markdown_cell(s))
def code(s):cells.append(nbf.v4.new_code_cell(s))
md('''# GSE149614 primary-tumor macrophages: PGAM5 RNA detection

## tl;dr

Ten samples ending in T contain 34,414 primary-tumor cells. Marker-based selection identifies 6,151 cells in macrophage-enriched original myeloid clusters. Of these, 167 have PGAM5 counts >0 and 5,984 have no PGAM5 counts. RNA nondetection is not protein negativity.

With the requested pooled-cell Wilcoxon analysis, 512 upregulated and zero downregulated candidate genes pass BH FDR <0.05, absolute approximate log2FC >=1 and >=10% detection in either group, excluding PGAM5 itself. There are 76 shared upregulated candidates with the previous GSE151530 analysis. These are exploratory associations, not a validated PGAM5 protein-positive signature.
''')
md('''## Context & Methods

### Key assumptions

- Only samples ending in T with official site=Tumor are selected. N (normal), P (PVTT), and L (lymph node) samples are excluded from differential expression.
- Updated GEO metadata labels cells as Myeloid, not macrophages. Marker expression is summarized in the original global clusters using all sites to avoid unstable annotation of tiny tumor-only clusters.
- Macrophage enrichment rule: C1QA, C1QB and C1QC each detected in >=80% of cells, and CD68 in >=75% of cells. Selected clusters: 5, 6, 21, 23, 38, 39, 41, 44, 46. Additional markers CSF1R/CD163/SPP1 support interpretation. This is our marker-based annotation, not an author-provided refined macrophage label.
- Excluded: cluster 10 (DC-enriched/mixed antigen-presenting cells), 16/26 (FCN1/S100A8/A9 monocyte-enriched/transition cells), 52 (GZMB/JCHAIN pDC-enriched), 53 (TPSAB1/TPSB2 mast-cell-enriched). This conservative selection does not cover every possible monocyte-derived inflammatory TAM.
- Basic QC: >=500 detected genes and <20% mitochondrial counts; all 6,151 selected primary-tumor cells meet these cutoffs.
- Grouping: PGAM5 count >0 versus count=0. Total-count normalization to 10,000 followed by log1p. Pooled-cell Wilcoxon with tie correction and BH FDR; no patient pairing, depth matching or sample/depth covariates, as requested.
- PGAM5 is retained in full tables and excluded from candidate lists because it defines the groups. Other gene classes are not removed.
- Scanpy log2FC is an approximate fold change from back-transformed mean log expression. Ranking is by FDR; a top20 shortlist is not a validated signature.
- Technical depth, biological sample/state composition, ambient RNA and doublets remain possible contributors. Macrophage markers are widespread in both groups, but hepatic transcripts are also frequent. Dedicated ambient-RNA correction and doublet detection were not performed.

Preparation and model scripts have been executed and retained. This notebook executes source-backed result checks and figures from their outputs; it does not silently refit models.
''')
code('''from pathlib import Path
import json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display
R=Path.cwd()
s=json.loads((R/'analysis_summary.json').read_text())
cross=json.loads((R/'cross_dataset_summary.json').read_text())
de=pd.read_csv(R/'GSE149614_primaryT_PGAM5_full_DE.csv',index_col=0)
up=pd.read_csv(R/'GSE149614_primaryT_PGAM5_upregulated.csv',index_col=0)
down=pd.read_csv(R/'GSE149614_primaryT_PGAM5_downregulated.csv',index_col=0)
markers=pd.read_csv(R/'original_myeloid_cluster_markers.csv')
annotation=pd.read_csv(R/'myeloid_cluster_annotation.csv')
samples=pd.read_csv(R/'primary_T_macrophages_by_sample.csv')
shared=pd.read_csv(R/'shared_upregulated_candidates_GSE149614_GSE151530.csv')
assert s['PGAM5_positive_cells']+s['PGAM5_undetected_cells']==s['QC_primary_macrophages']
assert len(up)==512 and len(down)==0 and len(shared)==76
assert up.FDR.lt(.05).all() and up.log2FC.ge(1).all() and not up.gene.eq('PGAM5').any()
assert samples['sample'].str.endswith('T').all() and samples['sample'].nunique()==10
assert up[['fraction_PGAM5_positive','fraction_PGAM5_undetected']].max(axis=1).ge(.1).all()
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':120,'savefig.bbox':'tight'})
BLUE='#355C8A';GOLD='#C68D29'
print('Count, scope and differential-threshold checks passed.')
''')
md('''## Data

Sources: [GEO series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614), [updated metadata](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.metadata.updated.txt.gz), [raw processed counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.count.txt.gz). Official GSM4505944 metadata also confirms HCC01T is primary HCC tumor.

The original matrix has 71,915 cells and 25,712 genes. Header IDs match metadata exactly and are unique; all parsed counts are nonnegative integers. Full gzip reading verifies CRC; input SHA256 hashes are recorded in input_sha256.json. Raw EGA sequence access is not needed for these publicly available processed counts.
''')
code('''display(pd.DataFrame({'quantity':['T samples','T cells','T myeloid cells','Selected/QC macrophages','PGAM5 detected','PGAM5 undetected','Upregulated candidates','Downregulated candidates','Shared upregulated with GSE151530'],'value':[10,34414,8209,6151,167,5984,512,0,76]}))
display(samples)
''')
md('''## Results

### 1. Macrophage annotation evidence

This chart uses original all-site myeloid clusters only for annotation. Differential expression uses T samples exclusively. The marker-based rule is independent of PGAM5 detection.
''')
code('''genes=['C1QA','C1QB','C1QC','CD68','CSF1R','CD163','SPP1','FCN1','S100A8','S100A9','CD1C','FCER1A','CLEC9A','GZMB','JCHAIN','TPSAB1']
t=markers.pivot(index='cluster',columns='gene',values='fraction').reindex(columns=genes)
selected=set(s['macrophage_clusters'])
fig,ax=plt.subplots(figsize=(12,7))
im=ax.imshow(t.to_numpy(),cmap='cividis',vmin=0,vmax=1,aspect='auto')
ax.set_xticks(np.arange(len(genes)),genes,rotation=45,ha='right')
ax.set_yticks(np.arange(len(t)),[f"{c} (Mac)" if c in selected else str(c) for c in t.index])
ax.set_ylabel('Original myeloid cluster');ax.set_title('GSE149614 original clusters: marker detection')
fig.colorbar(im,ax=ax,label='Fraction of cells with marker detected',shrink=.8)
fig.tight_layout();fig.savefig(R/'myeloid_annotation_markers.png');fig.savefig(R/'myeloid_annotation_markers.pdf');plt.show()
display(annotation)
''')
md('''### 2. PGAM5 detection in primary-tumor macrophages

All ten T samples contribute macrophages. Percentages are descriptive cell fractions, not estimates of protein-positive macrophages.
''')
code('''t=samples.sort_values('PGAM5_detection_rate',ascending=True).reset_index(drop=True)
fig,ax=plt.subplots(figsize=(9,5))
ax.barh(t['sample'],t.PGAM5_detection_rate*100,color=BLUE)
ax.set_xlabel('PGAM5 RNA detection in selected macrophages (%)')
ax.set_title('GSE149614 primary T samples: PGAM5 detection')
ax.set_xlim(0,8)
for i,row in t.iterrows():ax.text(row.PGAM5_detection_rate*100+.12,i,f"{row.PGAM5_detected}/{row.macrophages}",va='center')
ax.grid(axis='x',alpha=.15)
fig.tight_layout();fig.savefig(R/'PGAM5_detection_primaryT.png');fig.savefig(R/'PGAM5_detection_primaryT.pdf');plt.show()
''')
md('''### 3. Pooled-cell differential expression

The plot shows genes detected in >=10% of cells in either group. Orange points pass the final FDR/effect-size thresholds; blue points do not. PGAM5, the defining gene, is excluded from the plot. This is the requested unpaired cell-level analysis, and no sample-level validation is implied.
''')
code('''d=de[(~de.defining_gene)&de[['fraction_PGAM5_positive','fraction_PGAM5_undetected']].max(axis=1).ge(.1)].dropna(subset=['log2FC','p_value'])
passed=d.FDR.lt(.05)&d.log2FC.abs().ge(1)&d[['fraction_PGAM5_positive','fraction_PGAM5_undetected']].max(axis=1).ge(.1)
fig,ax=plt.subplots(figsize=(9,6))
ax.scatter(d.loc[~passed,'log2FC'],-np.log10(d.loc[~passed,'p_value'].clip(lower=1e-300)),s=7,color=BLUE,alpha=.35)
ax.scatter(d.loc[passed,'log2FC'],-np.log10(d.loc[passed,'p_value'].clip(lower=1e-300)),s=11,color=GOLD,alpha=.65)
ax.axvline(0,color='#666666',lw=.7);ax.axvline(1,color='#666666',lw=.7,ls='--');ax.axvline(-1,color='#666666',lw=.7,ls='--')
ax.set_xlabel('Approximate log2 fold change: PGAM5 detected / undetected')
ax.set_ylabel('-log10(nominal p-value)');ax.set_title('Primary-tumor macrophages: 512 upregulated, 0 downregulated')
ax.grid(alpha=.12)
fig.tight_layout();fig.savefig(R/'PGAM5_primaryT_differential_expression.png');fig.savefig(R/'PGAM5_primaryT_differential_expression.pdf');plt.show()
display(up[['gene','log2FC','FDR','fraction_PGAM5_positive','fraction_PGAM5_undetected']].head(20))
''')
md('''### 4. Comparison with GSE151530

Both cohorts use the same expression/DE thresholds and pooled-cell tests. GSE151530 uses its original TAM annotation; GSE149614 uses marker-defined macrophage-enriched clusters. This population difference must be considered when interpreting the overlap. Shared associations are candidates for further investigation, not proof of a specific PGAM5-driven mechanism.
''')
code('''display(shared[['gene','log2FC_GSE149614','FDR_GSE149614','log2FC_GSE151530','FDR_GSE151530']].head(20))
qa=pd.read_csv(R/'macrophage_and_hepatocyte_marker_QA.csv')
display(qa.pivot(index='gene',columns='group',values='fraction').round(3))
''')
md('''## Takeaways

The requested analysis produces 512 upregulated candidate genes. The 76 shared upregulated genes are a useful next-stage candidate pool; a top-ranked list alone does not establish a gene signature. Hepatic transcripts are frequent despite strong macrophage markers in both groups, so ambient RNA, phagocytosed transcripts and mixed cells remain possible explanations for some associations. They were not resolved in this workflow.

### Reproduction

Download the official metadata and count files, install the scientific dependencies listed in environment_versions.txt, and run prepare_myeloid.py then analyze_primary_macrophages.py. Scripts currently use /workspace/scratch/GSE149614; adjust paths when rerunning elsewhere. All-site counts are used only to annotate original myeloid clusters; DE uses T samples exclusively.

Retained tables support reproducing the displayed figures. The full model can be rerun from official inputs using the scripts. No pathway enrichment or independent signature scoring is claimed. Notebook code cells execute top-to-bottom; PNG/PDF figures are exported for inspection and use.
''')
nb.cells=cells;nb.metadata.kernelspec={'name':'python3','display_name':'Python 3','language':'python'}
NotebookClient(nb,timeout=180,kernel_name='python3',resources={'metadata':{'path':str(R)}}).execute()
nbf.validate(nb);p=R/'GSE149614_primaryT_PGAM5_macrophages.ipynb';nbf.write(nb,p);print('Executed',p)
