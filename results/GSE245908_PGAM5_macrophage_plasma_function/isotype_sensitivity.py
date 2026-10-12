"""Post-hoc descriptive count-concentration review; no added statistical tests."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE245908_plasma_mpl')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/GSE245908_plasma_cache')
import numpy as np,pandas as pd,anndata as ad,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze import ROOT,SAMPLES,ISOTYPES,dense
rows=[];summ=[]
for sample in SAMPLES:
 a=ad.read_h5ad(ROOT/'inputs'/f'{sample}_functional_marker_inputs.h5ad');q=pd.read_csv(ROOT/f'{sample}_plasma_target_scores.csv.gz',index_col=0);a=a[q.index];counts=np.column_stack([dense(a[:,gs].layers['counts']).sum(1) for gs in ISOTYPES.values()]);total=counts.sum(1)
 out=q[['neighbor_group','strict_plasma']].copy()
 for j,k in enumerate(ISOTYPES):out[k+'_raw']=counts[:,j]
 out['total_heavy_constant_raw']=total;out.to_csv(ROOT/f'{sample}_target_isotype_counts.csv.gz')
 for group in ['High_only','Other_detected_only']:
  ix=np.flatnonzero(q.neighbor_group.eq(group));largest=ix[np.argmax(total[ix])];igm=counts[ix,list(ISOTYPES).index('IgM')];summ.append(dict(sample=sample,group=group,n_targets=len(ix),IgM_total=int(igm.sum()),IgM_max=int(igm.max()),IgM_max_fraction=float(igm.max()/igm.sum()) if igm.sum() else np.nan,largest_total_barcode=q.index[largest],largest_total_heavy_counts=int(total[largest]),total_heavy_counts=int(total[ix].sum()),post_hoc_descriptive_only=True))
  for mode,ids in [('all_targets',ix),('exclude_largest_heavy_target',ix[ix!=largest])]:
   aggregate=counts[ids].sum(0);den=aggregate.sum()
   for j,k in enumerate(ISOTYPES):rows.append(dict(sample=sample,group=group,mode=mode,n_targets=len(ids),isotype=k,count=int(aggregate[j]),total_count=int(den),fraction=aggregate[j]/den if den else np.nan,post_hoc_descriptive_only=True))
pd.DataFrame(summ).to_csv(ROOT/'isotype_count_concentration.csv',index=False);d=pd.DataFrame(rows);d.to_csv(ROOT/'isotype_leave_largest_out.csv',index=False)
fig,axes=plt.subplots(1,2,figsize=(14,5));colors=['#285A8E','#B5822C','#A8536A','#607C44','#777777'];groups=[('High_only','all_targets'),('Other_detected_only','all_targets'),('High_only','exclude_largest_heavy_target'),('Other_detected_only','exclude_largest_heavy_target')]
for ax,sample in zip(axes,SAMPLES):
 bottom=np.zeros(4)
 for k,color in zip(ISOTYPES,colors):
  vals=np.array([100*d[(d['sample']==sample)&d.group.eq(g)&d['mode'].eq(m)&d.isotype.eq(k)].iloc[0].fraction for g,m in groups]);ax.bar(range(4),vals,bottom=bottom,color=color,label=k);bottom+=vals
 ax.set_xticks(range(4),['High-only','Other-only','High-only\nminus max','Other-only\nminus max'],fontsize=9);ax.set_ylim(0,100);ax.set_title(sample);ax.set_ylabel('Heavy-constant RNA count share (%)')
axes[1].legend(fontsize=8);fig.suptitle('One target can dominate isotype RNA | Post-hoc leave-largest-count-target-out description');fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(ROOT/'figures'/f'isotype_count_concentration_sensitivity.{ext}',dpi=180,bbox_inches='tight')
plt.close(fig);print(pd.DataFrame(summ).to_string(index=False))
