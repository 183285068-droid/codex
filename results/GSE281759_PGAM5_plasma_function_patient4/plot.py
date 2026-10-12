import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE281759_plasma_patient_mpl')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/GSE281759_plasma_patient_cache')
import numpy as np,pandas as pd,anndata as ad,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze import ROOT,PANELS,PATIENTS,COMPS,SAMPLES
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
def save(fig,name):
 fig.savefig(ROOT/'figures'/f'{name}.png',dpi=180,bbox_inches='tight');fig.savefig(ROOT/'figures'/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
def run():
 p=pd.read_csv(ROOT/'patient_paired_effects.csv');t=pd.read_csv(ROOT/'patient_level_statistics.csv');colors=['#285A8E','#B5822C','#A8536A','#607C44'];base=p[p.comparison.eq(COMPS[0])&p.sensitivity.eq('all_tumor_sections')]
 raw=base[base.method.eq('raw')];fig,axes=plt.subplots(2,4,figsize=(14,7));
 for ax,k in zip(axes.flat,PANELS):
  z=raw[raw.program.eq(k)].set_index('patient').reindex(PATIENTS)
  for j,patient in enumerate(PATIENTS):ax.plot([0,1],[z.loc[patient,'other_section_equal_mean'],z.loc[patient,'high_section_equal_mean']],marker=['o','s','^','D'][j],color=colors[j],lw=1.2,label=patient)
  ax.set_xticks([0,1],['Other-only','High-only']);ax.set_title(k);ax.set_ylabel('Equal-section patient mean RNA');ax.spines[['top','right']].set_visible(False)
 axes.flat[0].legend(fontsize=7);fig.suptitle('Four independent patient pairs | Same eligible sections in each pair; points are patients');fig.tight_layout();save(fig,'four_patient_paired_profiles')
 fig,axes=plt.subplots(1,2,figsize=(12,6),sharey=True)
 for ax,method in zip(axes,['raw','adjusted']):
  mat=base[base.method.eq(method)].pivot(index='program',columns='patient',values='patient_delta').reindex(index=PANELS,columns=PATIENTS);limit=max(.05,float(np.nanmax(np.abs(mat))))
  im=ax.imshow(mat,cmap='RdBu_r',vmin=-limit,vmax=limit,aspect='auto');ax.set_yticks(range(len(PANELS)),PANELS);ax.set_xticks(range(4),PATIENTS);ax.set_title(method+' patient paired effects')
  for i in range(len(PANELS)):
   for j in range(4):ax.text(j,i,f'{mat.iloc[i,j]:+.3f}',ha='center',va='center',fontsize=8,color='white' if abs(mat.iloc[i,j])>.65*limit else '#252525')
  fig.colorbar(im,ax=ax,shrink=.7,label='High minus other; mean log1p(CP10k) score')
 fig.suptitle('One effect per patient | Raw and adjusted panels have separate labelled color scales');fig.tight_layout();save(fig,'patient_raw_adjusted_effect_heatmaps')
 main=t[t.comparison.eq(COMPS[0])&t.method.eq('raw')&t.sensitivity.eq('all_tumor_sections')].set_index('program').reindex(PANELS);fig,ax=plt.subplots(figsize=(11,6));y=np.arange(len(PANELS));lo=main.t_model_mean_CI95_low;hi=main.t_model_mean_CI95_high;mean=main.mean_patient_delta
 ax.errorbar(mean,y,xerr=np.vstack([mean-lo,hi-mean]),fmt='o',color='#285A8E',capsize=4);ax.axvline(0,color='#555555',lw=.8);ax.set_yticks(y,PANELS);ax.invert_yaxis();ax.set_xlabel('Mean of four patient effects; t-model 95% interval (df=3)');ax.set_title('Parametric sensitivity, not a distribution-free interval\nExact patient sign-flip p and BH q reported at right')
 for i,r in enumerate(main.itertuples()):ax.text(1.01,i,f'p={r.exact_signflip_p:.3f}; q={r.family_BH_q_8:.3f}',transform=ax.get_yaxis_transform(),va='center',fontsize=8)
 fig.tight_layout();save(fig,'patient_effect_intervals_and_exact_tests')
 st=pd.read_csv(ROOT/'RNA_state_counts.csv');fig,axes=plt.subplots(1,2,figsize=(12,5));stateorder=['JCHAIN_MZB1_both','JCHAIN_only','MZB1_only','Neither']
 for ax,gp in zip(axes,['High_only','Other_only']):
  sub=st[st.group.eq(gp)].groupby(['patient','state']).n_targets.sum().unstack().reindex(index=PATIENTS,columns=stateorder).fillna(0);den=sub.sum(axis=1);bottom=np.zeros(4)
  for state,color in zip(stateorder,colors):
   vals=100*sub[state].to_numpy()/den.to_numpy();ax.bar(range(4),vals,bottom=bottom,color=color,label=state);bottom+=vals
  ax.set_xticks(range(4),[f'{pp}\nn={int(den[pp])} spots' for pp in PATIENTS]);ax.set_ylim(0,100);ax.set_title(gp);ax.set_ylabel('Within-group RNA-state spot proportion (%)')
 axes[1].legend(fontsize=7,loc='upper right');fig.suptitle('JCHAIN/MZB1 RNA co-detection states | Descriptive only, not IgA/IgG subtypes');fig.tight_layout();save(fig,'patient_RNA_state_composition')
 counts=pd.read_csv(ROOT/'section_target_counts.csv');agg=counts.groupby('patient')[['high_only','other_only','high_any','plasma_targets']].sum().reindex(PATIENTS);fig,ax=plt.subplots(figsize=(10,5))
 for j,k in enumerate(agg.columns):ax.bar(np.arange(4)+(j-1.5)*.2,agg[k],width=.19,color=colors[j],label=k)
 ax.set_xticks(range(4),PATIENTS);ax.set_ylabel('Descriptive target spot counts (not statistical n)');ax.set_title('Patient support | Statistical n=4 patients, not total spots');ax.legend(fontsize=8);fig.tight_layout();save(fig,'patient_target_support')
 for patient in PATIENTS:
  items=[i for i in SAMPLES if i['patient']==patient];fig,axes=plt.subplots(1,len(items),figsize=(3.4*len(items),4),squeeze=False)
  for ax,item in zip(axes[0],items):
   sample=item['sample'];a=ad.read_h5ad(ROOT/'inputs'/f'{sample}_plasma_marker_inputs.h5ad');q=pd.read_csv(ROOT/f'{sample}_target_profiles.csv.gz',index_col=0);im=a.uns['spatial'][sample]['images']['lowres'];sf=a.uns['scalefactors']['tissue_lowres_scalef'];xy=a.obsm['spatial']*sf;ix=a.obs_names.get_indexer(q.index);ax.imshow(im)
   for gp,color in [('High_only',colors[0]),('Other_only',colors[1]),('Shared',colors[2]),('No_detected_mac_neighbor','#777777')]:
    idx=ix[q.group.eq(gp).to_numpy()];ax.scatter(xy[idx,0],xy[idx,1],s=8,color=color,label=gp)
   h=a.obs.PGAM5_mac_group.eq('PGAM5_high').to_numpy();ax.scatter(xy[h,0],xy[h,1],s=16,facecolors='none',edgecolors=colors[2],label='PGAM5-high source');ax.set_title(sample,fontsize=9);ax.axis('off')
  axes[0,0].legend(fontsize=5);fig.suptitle(patient+' | Real H&E / Visium positions; mixed RNA targets');fig.tight_layout();save(fig,patient+'_plasma_target_spatial_maps')
if __name__=='__main__':run()
