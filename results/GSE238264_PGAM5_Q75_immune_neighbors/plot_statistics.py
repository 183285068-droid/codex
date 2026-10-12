from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/gse245908_mpl')
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;d=pd.read_csv(R/'spatial_conditional_tests.csv');types=['T','B','Plasma','NK','cDC1','cDC2','pDC','Monocyte','Neutrophil'];samples=['HCC1R','HCC2R','HCC3R','HCC4R','HCC5NR','HCC6NR','HCC7NR'];fig,axs=plt.subplots(2,4,figsize=(18,9))
for ax,sample in zip(axs.ravel(),samples):
 t=d[d['sample']==sample].set_index('feature').loc[types];ax.axvline(0,c='#777',lw=.8);ax.scatter(t.delta_high_minus_other*100,np.arange(9),c='#2678a8',s=18);ax.set_yticks(range(9),types);ax.set_ylim(8.5,-.5);ax.set_xlim(-70,70);ax.set_title(sample);ax.set_xlabel('High − other detected (pp)')
 for i,row in enumerate(t.itertuples()):
  label='Not testable' if pd.isna(row.spatial_conditional_p) else f'p={row.spatial_conditional_p:.3f}; q={row.BH_q_63_planned:.2f}'
  ax.text(.99,i,label,transform=ax.get_yaxis_transform(),ha='right',va='bottom',fontsize=7)
axs.ravel()[-1].axis('off');fig.suptitle('Q75 immune neighborhood contrasts | 63 planned spatial tests\nDifferences are unadjusted; p values condition on location/depth strata');fig.tight_layout();fig.savefig(R/'figures/spatial_immune_tests.png',dpi=180,bbox_inches='tight');fig.savefig(R/'figures/spatial_immune_tests.pdf',bbox_inches='tight');plt.close(fig)
