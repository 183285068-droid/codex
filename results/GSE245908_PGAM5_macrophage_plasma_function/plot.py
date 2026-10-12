import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE245908_plasma_mpl')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/GSE245908_plasma_cache')
import numpy as np,pandas as pd,anndata as ad,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from analyze import ROOT,SAMPLES,PANELS,PROGRAMS,COMPS
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
def save(fig,name):
 fig.savefig(ROOT/'figures'/f'{name}.png',dpi=180,bbox_inches='tight');fig.savefig(ROOT/'figures'/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
def run():
 s=pd.read_csv(ROOT/'program_statistics.csv');d=pd.read_csv(ROOT/'program_group_profiles.csv');iso=pd.read_csv(ROOT/'immunoglobulin_isotype_counts.csv');gs=pd.read_csv(ROOT/'target_group_summary.csv');genes=pd.read_csv(ROOT/'marker_gene_descriptive_profiles.csv')
 primary=s[s.comparison.eq(COMPS[0])];mat=primary.pivot(index='program',columns='sample',values='standardized_delta').reindex(PROGRAMS);fig,ax=plt.subplots(figsize=(9,7));im=ax.imshow(mat,cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto');ax.set_yticks(range(len(PROGRAMS)),PROGRAMS);ax.set_xticks([0,1],mat.columns)
 for i in range(len(PROGRAMS)):
  for j in range(2):ax.text(j,i,f'{mat.iloc[i,j]:+.2f}',ha='center',va='center',color='white' if abs(mat.iloc[i,j])>.65 else '#252525')
 ax.set_title('Plasma-related RNA near high-only vs other-detected-only Mac-RNA\nNuisance-residual contrasts; all primary BH q=1.0');fig.colorbar(im,ax=ax,label='Residual mean difference / all plasma-target score SD');fig.tight_layout();save(fig,'adjusted_program_contrasts')
 fig,axes=plt.subplots(1,2,figsize=(14,6),sharey=True,sharex=True);x=np.arange(len(PROGRAMS));groups=['High_only','Other_detected_only','No_detected_mac_neighbor'];colors=['#285A8E','#B5822C','#777777']
 for ax,sample in zip(axes,SAMPLES):
  for j,g in enumerate(groups):
   v=d[(d['sample']==sample)&d.group.eq(g)].set_index('program').reindex(PROGRAMS);ax.barh(x+(j-1)*.23,v['mean'],height=.22,color=colors[j],label=g+f' (n={int(v.n_targets.iloc[0])})')
  ax.set_yticks(x,PROGRAMS);ax.set_title(sample);ax.set_xlabel('Mean log1p(full-library CP10k), zeros retained');ax.legend(fontsize=7);ax.spines[['top','right']].set_visible(False)
 axes[0].invert_yaxis();fig.suptitle('Descriptive RNA profiles | Program means are not measured functional activity');fig.tight_layout();save(fig,'functional_program_profiles')
 fig,axes=plt.subplots(1,2,figsize=(12,5));icolors=['#285A8E','#B5822C','#A8536A','#607C44','#777777'];isotypes=['IgG','IgA','IgM','IgD','IgE']
 for ax,sample in zip(axes,SAMPLES):
  grp=['High_only','Other_detected_only'];bottom=np.zeros(2)
  for k,color in zip(isotypes,icolors):
   v=iso[(iso['sample']==sample)&iso.isotype.eq(k)].set_index('group').reindex(grp);vals=100*v.count_fraction.to_numpy();ax.bar(np.arange(2),vals,bottom=bottom,color=color,label=k);bottom+=vals
  for j,g in enumerate(grp):
   row=iso[(iso['sample']==sample)&iso.group.eq(g)].iloc[0];ax.text(j,102,f'n={int(row.n_targets)} spots\n{int(row.total_heavy_constant_counts)} heavy-constant UMIs',ha='center',fontsize=8)
  ax.set_xticks([0,1],['High-only neighbor','Other-detected-only neighbor']);ax.set_ylim(0,118);ax.set_title(sample);ax.set_ylabel('Share of selected heavy-chain constant-gene counts (%)')
 axes[1].legend(fontsize=8);fig.suptitle('Immunoglobulin constant-region RNA composition | Not antibody secretion or class-switch rates');fig.tight_layout();save(fig,'immunoglobulin_isotype_composition')
 marker=['PRDM1','XBP1','IRF4','TNFRSF17','MZB1','JCHAIN','DERL3','SDC1','IGKC','IL10','TGFB1','CD274','IGHG1','IGHA1','IGHM'];fig,axes=plt.subplots(1,2,figsize=(14,6),sharey=True,sharex=True)
 for ax,sample in zip(axes,SAMPLES):
  z=genes[genes['sample']==sample].set_index('gene').reindex(marker);y=np.arange(len(marker));ax.scatter(z.high_mean,y-.12,s=z.high_detect_pct*1.5,color='#285A8E',label='High-only (area = detected %)');ax.scatter(z.other_mean,y+.12,s=z.other_detect_pct*1.5,color='#B5822C',label='Other-only (area = detected %)');ax.set_yticks(y,marker);ax.set_xlabel('Mean log1p(CP10k)');ax.set_title(sample);ax.legend(fontsize=7);ax.spines[['top','right']].set_visible(False)
 axes[0].invert_yaxis();fig.suptitle('Marker profiles | Descriptive; zero-size dot = no RNA detection, no gene-level significance test');fig.tight_layout();save(fig,'plasma_marker_profiles')
 fig,ax=plt.subplots(figsize=(9,5));x=np.arange(2)
 for j,col in enumerate(['plasma_spots','high_any','high_only','other_only','strict_plasma']):ax.bar(x+(j-2)*.15,gs[col],width=.14,label=col,color=['#777777','#285A8E','#597EA2','#B5822C','#A8536A'][j])
 ax.set_xticks(x,gs['sample']);ax.set_ylabel('Mixed Visium target spots');ax.set_title('Plasma-target support and group sizes\nStrict MZB1+JCHAIN+DERL3 gate yields 0 high-adjacent targets');ax.legend(fontsize=8);fig.tight_layout();save(fig,'target_group_and_strict_support')
 spatial_max={k:max(pd.read_csv(ROOT/f'{sample}_plasma_target_scores.csv.gz')[k].max() for sample in SAMPLES) for k in ['Plasma_differentiation','Antibody_heavy_chain']}
 for sample in SAMPLES:
  a=ad.read_h5ad(ROOT/'inputs'/f'{sample}_functional_marker_inputs.h5ad');o=a.obs;q=pd.read_csv(ROOT/f'{sample}_plasma_target_scores.csv.gz',index_col=0);im=a.uns['spatial'][sample]['images']['lowres'];sf=a.uns['spatial'][sample]['scalefactors'];xy=a.obsm['spatial']*sf['tissue_lowres_scalef'];idx=a.obs_names.get_indexer(q.index);source=o.source_group.eq('PGAM5_high').to_numpy();e=a.obsp['spatial_connectivities'].tocoo();target=o.gate_Plasma_related.to_numpy(bool);use=source[e.row]&target[e.col]
  fig,axes=plt.subplots(1,3,figsize=(13,4.8))
  for ax in axes:ax.imshow(im);ax.axis('off')
  for group,color in [('High_only','#285A8E'),('Shared_high_other','#A8536A'),('Other_detected_only','#B5822C'),('No_detected_mac_neighbor','#777777')]:
   ix=idx[q.neighbor_group.eq(group).to_numpy()];axes[0].scatter(xy[ix,0],xy[ix,1],s=15,color=color,label=group)
  axes[0].scatter(xy[source,0],xy[source,1],s=25,facecolors='none',edgecolors='#A8536A',label='PGAM5-high source');axes[0].add_collection(LineCollection(np.stack([xy[e.row[use]],xy[e.col[use]]],axis=1),colors='#222222',lw=.6,alpha=.6));axes[0].legend(fontsize=5);axes[0].set_title('Actual source-target first-ring edges',fontsize=9)
  for ax,k in zip(axes[1:],['Plasma_differentiation','Antibody_heavy_chain']):
   sc=ax.scatter(xy[idx,0],xy[idx,1],c=q[k],s=18,cmap='viridis',vmin=0,vmax=spatial_max[k]);ax.set_title(k,fontsize=9);fig.colorbar(sc,ax=ax,shrink=.6,label='RNA score (shared scale)')
  fig.suptitle(sample+' | H&E / real Visium coordinates; labels are mixed-spot RNA proxies');fig.tight_layout();save(fig,sample+'_plasma_spatial_profiles')
if __name__=='__main__':run()
