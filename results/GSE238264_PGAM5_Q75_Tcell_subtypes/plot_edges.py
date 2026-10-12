import os,json
os.environ.setdefault('MPLCONFIGDIR','/tmp/gse245908_mpl')
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from matplotlib.collections import LineCollection,PatchCollection
from matplotlib.patches import Circle
from analyze import R,RAW
s='HCC4R';o=pd.read_csv(R/f'{s}_T_subtypes.csv.gz',index_col=0);e=pd.read_csv(R/f'{s}_T_LR_edges.csv.gz');root=RAW/s/'spatial';im=np.array(Image.open(root/'tissue_hires_image.png').convert('RGB'));sc=json.loads((root/'scalefactors_json.json').read_text());xy=o[['pixel_col','pixel_row']].to_numpy()*sc['tissue_hires_scalef'];radius=sc['spot_diameter_fullres']*sc['tissue_hires_scalef']/2;lookup={b:i for i,b in enumerate(o.index)};fig,axs=plt.subplots(2,2,figsize=(11,11));high=o.source_group.eq('PGAM5_high').to_numpy()
for ax,k,pair in zip(axs.ravel(),['All_CD8','All_CD8','CD8_exhaustion','Treg'],['CXCL12--CXCR4','SPP1--CD44','TGFB1--TGFBR1_TGFBR2','NECTIN2--TIGIT']):
 ax.imshow(im,origin='upper');mask=o['gate_'+k].to_numpy();ax.add_collection(PatchCollection([Circle(p,radius) for p in xy[mask]],facecolors='#007bb6',edgecolors='#007bb6',alpha=.5,lw=.4));ax.add_collection(PatchCollection([Circle(p,radius*1.15) for p in xy[high]],facecolors='none',edgecolors='#f06016',lw=1.1));d=e[e.source_group.eq('PGAM5_high')&e['target_'+k]&(e[pair]>0)];segments=[np.array([xy[lookup[a]],xy[lookup[b]]]) for a,b in zip(d.source_barcode,d.target_barcode)];ax.add_collection(LineCollection(segments,colors='#00c6d4',lw=2,zorder=5));ax.set_xlim(max(0,xy[:,0].min()-40),min(im.shape[1],xy[:,0].max()+40));ax.set_ylim(min(im.shape[0],xy[:,1].max()+40),max(0,xy[:,1].min()-40));ax.set_axis_off();ax.set_title(k.replace('_',' ')+' | '+pair+'\n'+f'Edges={len(d)} | sources={d.source_barcode.nunique()} | targets={d.target_barcode.nunique()}',fontsize=10)
fig.suptitle('HCC4R | Q75 candidate T-subtype LR edges over H&E\nOrange: High source; blue: target proxy; cyan: adjacent codetected RNA',fontsize=13);fig.subplots_adjust(top=.9,hspace=.18,wspace=.07,bottom=.02);fig.savefig(R/'figures/HCC4R_HE_candidate_T_LR_edges.png',dpi=180,bbox_inches='tight');fig.savefig(R/'figures/HCC4R_HE_candidate_T_LR_edges.pdf',bbox_inches='tight');plt.close(fig)
