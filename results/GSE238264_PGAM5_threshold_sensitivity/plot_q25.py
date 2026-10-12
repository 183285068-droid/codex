from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/gse245908_mpl')
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;QS=[.25];METRICS=['mac_fraction','T_fraction','fibroblast_score'];d=pd.read_csv(R/'slice_metrics_by_threshold.csv');t=pd.read_csv(R/'slice_response_tests.csv')
fig,axs=plt.subplots(1,3,figsize=(14,4.5));axs=np.atleast_2d(axs)
for i,q in enumerate(QS):
 for j,metric in enumerate(METRICS):
  ax=axs[i,j];dd=d[d['quantile']==q];factor=100 if metric.endswith('fraction') else 1
  for k,g in enumerate(['responder','non-responder']):
   ss=dd[dd.response==g];xx=k+np.linspace(-.15,.15,len(ss));ax.scatter(xx,ss[metric]*factor,c='#2678a8' if k==0 else '#b95e27')
   for x,row in zip(xx,ss.itertuples()):ax.annotate(f'{row.sample} (n={row.n_high})',(x,getattr(row,metric)*factor),xytext=(3,14 if k==1 and row.sample=='HCC7NR' and getattr(row,metric)<.01 else 3),textcoords='offset points',fontsize=7)
  rr=t[(t['quantile']==q)&(t.metric==metric)].iloc[0];ax.set_title(f'Q{int(q*100)} | {metric}\np={rr.p_exact:.3f}; joint q={rr.BH_q_12_tests:.3f}',fontsize=9);ax.set_xticks([0,1],['Responder (4)','Non-responder (3)']);ax.set_xlim(-.4,1.7);ax.set_ylim(bottom=0);ax.set_ylabel('%' if factor==100 else 'Mean log1p(CP10k)')
fig.suptitle('Adopted PGAM5 Q25 threshold | each dot is one slice/patient\nSource counts shown; low-count slice estimates remain unstable');fig.tight_layout();fig.savefig(R/'figures/Q25_response_comparisons.png',dpi=180,bbox_inches='tight');fig.savefig(R/'figures/Q25_response_comparisons.pdf',bbox_inches='tight');plt.close(fig)
