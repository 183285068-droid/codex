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
checks=[]
for sample in SAMPLES:
 o=pd.read_csv(R/f'{sample}_immune_spots.csv.gz');root=RAW/sample/'spatial';im=np.array(Image.open(root/'tissue_hires_image.png').convert('RGB'));sc=json.loads((root/'scalefactors_json.json').read_text());xy=o[['pixel_col','pixel_row']].to_numpy()*sc['tissue_hires_scalef'];radius=sc['spot_diameter_fullres']*sc['tissue_hires_scalef']/2;assert (xy>=0).all() and (xy[:,0]<im.shape[1]).all() and (xy[:,1]<im.shape[0]).all();high=o.source_group.eq('PGAM5_high').to_numpy();fig,axs=plt.subplots(3,3,figsize=(12,12))
 for ax,k in zip(axs.ravel(),types):
  ax.imshow(im,origin='upper');mask=o['gate_'+k].to_numpy();ax.add_collection(PatchCollection([Circle(p,radius) for p in xy[mask]],facecolors='#007bb6',edgecolors='#007bb6',alpha=.5,lw=.4));ax.add_collection(PatchCollection([Circle(p,radius*1.15) for p in xy[high]],facecolors='none',edgecolors='#f06016',lw=1.3));ax.set_xlim(max(0,xy[:,0].min()-40),min(im.shape[1],xy[:,0].max()+40));ax.set_ylim(min(im.shape[0],xy[:,1].max()+40),max(0,xy[:,1].min()-40));ax.set_axis_off();ax.set_title(f'{k} candidates={mask.sum()}\nQ75 High={high.sum()}',fontsize=10)
 fig.suptitle(sample+' | Q75 source + immune RNA candidates fused with H&E\nOrange rings: High source; translucent blue: target candidate',fontsize=13);fig.subplots_adjust(top=.91,hspace=.14,wspace=.06,bottom=.02);save(fig,sample+'_HE_immune_neighbors');checks.append(dict(sample=sample,official_hires_scale=sc['tissue_hires_scalef'],spot_radius=radius,all_coordinates_inside_image=True))
(R/'HE_registration_checks.json').write_text(json.dumps(checks,indent=2));print('9 PNG/PDF figure pairs saved')
