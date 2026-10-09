from pathlib import Path
import json,numpy as np,pandas as pd,anndata as ad
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parent;a=ad.read_h5ad(R/'macrophage_counts_QC.h5ad')
pos=a.obs.PGAM5_counts.to_numpy()>0;n1=int(pos.sum());n0=int((~pos).sum())
c1=np.asarray((a.X[pos]>0).sum(0)).ravel();c0=np.asarray((a.X[~pos]>0).sum(0)).ravel()
p=[fisher_exact([[int(x),n1-int(x)],[int(y),n0-int(y)]],alternative='two-sided').pvalue for x,y in zip(c1,c0)]
t=pd.DataFrame({'gene':a.var.gene.to_numpy(),'detected_positive_cells':c1,'detected_undetected_cells':c0,'Fisher_two_sided_p':p,'Fisher_BH_FDR':multipletests(p,method='fdr_bh')[1]},index=a.var_names);t.index.name='gene_ID'
t.to_csv(R/'exact_detection_Fisher_all_genes.csv')
u=pd.read_csv(R/'GSE146115_HCC_tumor_PGAM5_upregulated.csv',index_col=0)
d=pd.read_csv(R/'GSE146115_HCC_tumor_PGAM5_downregulated.csv',index_col=0)
for label,x in [('upregulated',u),('downregulated',d)]:
 z=x.join(t.drop(columns='gene'));z.to_csv(R/(label+'_with_exact_detection_check.csv'))
 z[z.Fisher_BH_FDR.lt(.05)].to_csv(R/(label+'_passing_exact_detection_check.csv'))
s={'positive_cells':n1,'undetected_cells':n0,'tested_genes':len(t),'non_PGAM5_genes_with_Fisher_FDR005':int((~t.gene.eq('PGAM5')&t.Fisher_BH_FDR.lt(.05)).sum()),'up_candidates_passing_Fisher_FDR005':int(t.loc[u.index].Fisher_BH_FDR.lt(.05).sum()),'down_candidates_passing_Fisher_FDR005':int(t.loc[d.index].Fisher_BH_FDR.lt(.05).sum()),'method':'Supplementary two-sided Fisher exact detection-frequency test, BH across all human source rows. Different estimand from expression-rank Wilcoxon; does not replace requested primary results.'}
(R/'exact_detection_summary.json').write_text(json.dumps(s,indent=2));print(json.dumps(s,indent=2))
