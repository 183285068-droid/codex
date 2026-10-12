import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE281759_TLS_mpl')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/GSE281759_TLS_cache')
import numpy as np,pandas as pd,anndata as ad,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze import ROOT,SAMPLES,TIERS
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
def save(fig,name):
 fig.savefig(ROOT/'figures'/f'{name}.png',dpi=180,bbox_inches='tight');fig.savefig(ROOT/'figures'/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
def run():
 summ=pd.read_csv(ROOT/'section_candidate_summary.csv');stats=pd.read_csv(ROOT/'spatial_statistics.csv');cov=pd.read_csv(ROOT/'TLS_marker_coverage.csv');samples=[s['sample'] for s in SAMPLES]
 fig,axes=plt.subplots(1,2,figsize=(14,6));x=np.arange(len(samples));colors=['#777777','#285A8E','#B5822C']
 for i,k in enumerate(TIERS):
  d=summ[summ.tier.eq(k)].set_index('sample').reindex(samples);axes[0].bar(x+(i-1)*.25,d.centers,width=.24,color=colors[i],label=k)
  axes[1].bar(x+(i-1)*.25,100*d.high_within/d.high_all,width=.24,color=colors[i],label=k)
 axes[0].set_ylabel('Qualifying B-centered neighborhoods (overlap allowed)');axes[1].set_ylabel('PGAM5-high sources within proxy patch union (%)')
 for ax in axes:ax.set_xticks(x,samples,rotation=55,ha='right');ax.spines[['top','right']].set_visible(False)
 axes[1].set_ylim(0,100);axes[0].legend(fontsize=7);fig.suptitle('GSE281759 | T/B aggregate RNA patches and CXCL13 support\nExploratory proxies; no pathological TLS identification');fig.tight_layout();save(fig,'candidate_patch_overview')
 primary=stats[stats.tier.eq(TIERS[1])&stats.comparison.eq('high_vs_other_detected')];fig,ax=plt.subplots(figsize=(12,6))
 for i,endpoint in enumerate(['within_patch','adjacent_outside_patch']):
  d=primary[primary.endpoint.eq(endpoint)].set_index('sample').reindex(samples);ax.bar(x+(i-.5)*.3,d.delta_pp,width=.28,color=['#285A8E','#B5822C'][i],label=endpoint)
  for j,row in enumerate(d.itertuples()):
   if row.status!='tested':ax.text(j+(i-.5)*.3,1,'NA',rotation=90,ha='center',fontsize=7)
 ax.set_xticks(x,samples,rotation=50,ha='right');ax.axhline(0,color='#555555',lw=.8);ax.set_ylabel('High minus other detected sources (percentage points)');ax.set_title('CXCL13/MHCII-supported proxy patches | Conditional spatial contrast\nNA: unavailable comparison or constant null; primary BH q in result table');ax.legend(fontsize=9);fig.tight_layout();save(fig,'primary_spatial_contrasts')
 cc=cov[cov['sample'].eq(samples[0])].drop_duplicates(['panel','gene']);fig,ax=plt.subplots(figsize=(10,5));by=cc.groupby('panel').present.agg(['sum','count']);x2=np.arange(len(by));ax.barh(x2,by['count'],color='#d4d7db',label='Requested illustrative marker set');ax.barh(x2,by['sum'],color='#285A8E',label='Available in 179-marker input');ax.set_yticks(x2,by.index);ax.invert_yaxis();ax.set_xlabel('Genes');ax.set_title('TLS-related marker coverage | Missing FDC/HEV and chemokine markers\nPresence in file does not imply RNA detected or confirm TLS identity');ax.legend(fontsize=8);fig.tight_layout();save(fig,'TLS_marker_coverage')
 for patient in ['A009','B013','C004','C008']:
  items=[i for i in SAMPLES if i['patient']==patient];fig,axes=plt.subplots(len(items),3,figsize=(12,3.5*len(items)),squeeze=False)
  for r,item in enumerate(items):
   sample=item['sample'];a=ad.read_h5ad(ROOT/'inputs'/f'{sample}_spatial_marker_subset.h5ad');o=pd.read_csv(ROOT/f'{sample}_TLS_proxy_annotations.csv.gz',index_col=0).loc[a.obs_names];xy=a.obsm['spatial']*a.uns['scalefactors']['tissue_lowres_scalef'];im=a.uns['spatial'][sample]['images']['lowres'];high=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy()
   for c,k in enumerate(TIERS):
    ax=axes[r,c];ax.imshow(im);region=o['region_'+k].to_numpy(bool);center=o['center_'+k].to_numpy(bool);t=o.gate_T.to_numpy(bool);b=o.gate_B.to_numpy(bool)
    ax.scatter(xy[region,0],xy[region,1],s=20,color='#56A2AC',alpha=.2,label='Proxy patch union');ax.scatter(xy[t,0],xy[t,1],s=3,color='#285A8E',alpha=.6,label='T RNA');ax.scatter(xy[b,0],xy[b,1],s=10,color='#B5822C',label='B RNA');ax.scatter(xy[center,0],xy[center,1],s=30,marker='s',facecolors='none',edgecolors='#2a7777',label='Qualifying center');ax.scatter(xy[high,0],xy[high,1],s=18,facecolors='none',edgecolors='#A8536A',linewidths=.8,label='PGAM5-high Mac RNA');ax.set_title(sample+' | '+k,fontsize=8);ax.axis('off')
    if r==0:ax.legend(fontsize=5,loc='upper right',framealpha=.85)
  fig.suptitle(patient+' | Actual H&E / Visium positions; RNA proxies, not TLS boundaries',fontsize=11);fig.tight_layout();save(fig,patient+'_TLS_proxy_spatial_maps')
if __name__=='__main__':run()
