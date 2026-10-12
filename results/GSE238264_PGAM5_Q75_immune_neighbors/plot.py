import os,json
os.environ.setdefault('MPLCONFIGDIR','/tmp/gse245908_mpl')
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from analyze import R,RAW,SAMPLES,PANELS,GROUPS
F=R/'figures';types=list(PANELS);summary=pd.read_csv(R/'source_immune_neighborhoods.csv');lr=pd.read_csv(R/'LR_by_immune_subtype.csv');plt.rcParams.update({'font.size':10})
def save(fig,name):fig.savefig(F/(name+'.png'),dpi=170,bbox_inches='tight');fig.savefig(F/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
fig,axs=plt.subplots(2,4,figsize=(18,9));colors=['#de6027','#2678a8','#888888']
for ax,sample in zip(axs.ravel(),SAMPLES):
 for j,(group,c) in enumerate(zip(GROUPS,colors)):
  d=summary[(summary['sample']==sample)&(summary.source_group==group)].set_index('subtype').loc[types];ax.barh(np.arange(9)+(j-1)*.23,d.mean_neighbor_fraction*100,height=.23,color=c,label=f'{group} (n={int(d.n_sources.iloc[0])})')
 ax.set_yticks(range(9),types);ax.invert_yaxis();ax.set_xlim(0,100);ax.set_title(sample);ax.set_xlabel('Mean neighbor coverage (%)');ax.legend(fontsize=7)
axs.ravel()[-1].axis('off');fig.suptitle('PGAM5 Q75 myeloid sources → immune RNA neighborhoods\nSource-weighted means; overlapping proxy gates, not cell proportions');fig.tight_layout();save(fig,'immune_neighborhood_coverage')
selected=['SPP1--CD44','CXCL12--CXCR4','CXCL13--CXCR5','TNFSF13B--TNFRSF13C','TNFSF13--TNFRSF17','PVR--TIGIT','LGALS9--HAVCR2','TGFB1--TGFBR1_TGFBR2'];fig,axs=plt.subplots(2,4,figsize=(23,12))
for ax,sample in zip(axs.ravel(),SAMPLES):
 d=lr[(lr['sample']==sample)&lr.source_group.eq('PGAM5_high')].copy();d['pair']=d.ligand+'--'+d.receptor;matrix=d.pivot(index='subtype',columns='pair',values='codetection_pct').reindex(index=types,columns=selected).to_numpy();heat=ax.imshow(np.ma.masked_invalid(matrix),vmin=0,vmax=100,cmap='Blues',aspect='auto');ax.set_xticks(range(len(selected)),selected,rotation=45,ha='right',fontsize=8);ax.set_yticks(range(9),[k+' (e='+str(int(d[d.subtype==k].n_edges.iloc[0]))+')' for k in types],fontsize=8);ax.set_title(sample)
 for i in range(9):
  for j in range(len(selected)):ax.text(j,i,'NA' if np.isnan(matrix[i,j]) else f'{matrix[i,j]:.0f}',ha='center',va='center',fontsize=7,color='white' if not np.isnan(matrix[i,j]) and matrix[i,j]>55 else 'black')
axs.ravel()[-1].axis('off');fig.suptitle('Selected candidate LR codetection in High-source edges (%)\nNA = no eligible edges; e = all edges; no interaction probabilities');fig.tight_layout(rect=[0,0,.94,.94]);cax=fig.add_axes([.96,.3,.01,.4]);fig.colorbar(heat,cax=cax,label='Codetected edges (%)');save(fig,'immune_LR_codetection')
checks=[]
for sample in SAMPLES:
 o=pd.read_csv(R/f'{sample}_immune_spots.csv.gz');root=RAW/sample/'spatial';im=np.array(Image.open(root/'tissue_hires_image.png').convert('RGB'));sc=json.loads((root/'scalefactors_json.json').read_text());xy=o[['pixel_col','pixel_row']].to_numpy()*sc['tissue_hires_scalef'];radius=sc['spot_diameter_fullres']*sc['tissue_hires_scalef']/2;assert (xy>=0).all() and (xy[:,0]<im.shape[1]).all() and (xy[:,1]<im.shape[0]).all();high=o.source_group.eq('PGAM5_high').to_numpy();fig,axs=plt.subplots(3,3,figsize=(12,12))
 for ax,k in zip(axs.ravel(),types):
  ax.imshow(im,origin='upper');mask=o['gate_'+k].to_numpy();ax.add_collection(PatchCollection([Circle(p,radius) for p in xy[mask]],facecolors='#007bb6',edgecolors='#007bb6',alpha=.5,lw=.4));ax.add_collection(PatchCollection([Circle(p,radius*1.15) for p in xy[high]],facecolors='none',edgecolors='#f06016',lw=1.3));ax.set_xlim(max(0,xy[:,0].min()-40),min(im.shape[1],xy[:,0].max()+40));ax.set_ylim(min(im.shape[0],xy[:,1].max()+40),max(0,xy[:,1].min()-40));ax.set_axis_off();ax.set_title(f'{k} candidates={mask.sum()}\nQ75 High={high.sum()}',fontsize=10)
 fig.suptitle(sample+' | Q75 source + immune RNA candidates fused with H&E\nOrange rings: High source; translucent blue: target candidate',fontsize=13);fig.subplots_adjust(top=.91,hspace=.14,wspace=.06,bottom=.02);save(fig,sample+'_HE_immune_neighbors');checks.append(dict(sample=sample,official_hires_scale=sc['tissue_hires_scalef'],spot_radius=radius,all_coordinates_inside_image=True))
(R/'HE_registration_checks.json').write_text(json.dumps(checks,indent=2));print('9 PNG/PDF figure pairs saved')
