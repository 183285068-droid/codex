from pathlib import Path
from functools import lru_cache
import json
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
R=Path('/workspace/scratch/GSE125449');a=ad.read_h5ad(R/'HCC_TAM_counts_QC.h5ad');pos=a.obs.PGAM5_counts.to_numpy()>0;npos=int(pos.sum());nneg=int((~pos).sum())
kpos=np.asarray((a.X[pos]>0).sum(0)).ravel().astype(int);kneg=np.asarray((a.X[~pos]>0).sum(0)).ravel().astype(int)
@lru_cache(None)
def exact(kp,kn):
 if kp+kn==0:return 1.
 return float(fisher_exact([[kp,npos-kp],[kn,nneg-kn]],alternative='two-sided').pvalue)
p=np.array([exact(int(kp),int(kn)) for kp,kn in zip(kpos,kneg)])
r=pd.DataFrame({'gene':a.var.gene,'detected_positive_cells':kpos,'detected_undetected_cells':kneg,'fraction_positive':kpos/npos,'fraction_undetected':kneg/nneg,'Fisher_two_sided_p':p,'Fisher_BH_FDR':multipletests(p,method='fdr_bh')[1]},index=a.var_names)
r.to_csv(R/'exact_detection_Fisher_all_genes.csv')
u=pd.read_csv(R/'GSE125449_HCC_PGAM5_upregulated.csv',index_col=0)
joined=u.join(r[['detected_positive_cells','detected_undetected_cells','Fisher_two_sided_p','Fisher_BH_FDR']]);joined.to_csv(R/'Wilcoxon_candidates_with_exact_detection_check.csv')
credible=joined[joined.Fisher_BH_FDR.lt(.05)];credible.to_csv(R/'Wilcoxon_candidates_passing_exact_detection_check.csv')
s={'test':'Two-sided Fisher exact comparison of gene detection proportions, BH across all 21324 tested genes; supplementary, not a substitute for the requested expression Wilcoxon test','positive_cells':npos,'undetected_cells':nneg,'Wilcoxon_up_candidates':len(u),'Wilcoxon_candidates_passing_Fisher_FDR05':len(credible),'Fisher_FDR05_genes_excluding_PGAM5':int((r.Fisher_BH_FDR.lt(.05)&r.gene.ne('PGAM5')).sum()),'candidate_genes_passing':credible.gene.tolist(),'one_positive_zero_negative_exact_p':exact(1,0),'one_positive_zero_negative_Wilcoxon_p':float(joined.loc[(joined.detected_positive_cells==1)&(joined.detected_undetected_cells==0),'p_value'].iloc[0])}
(R/'exact_detection_summary.json').write_text(json.dumps(s,indent=2));print(json.dumps(s,indent=2))
