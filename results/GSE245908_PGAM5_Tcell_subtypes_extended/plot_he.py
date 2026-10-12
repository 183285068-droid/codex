"""Overlay measured Visium coordinates on supplied histology, using official scale."""
from pathlib import Path
import os,gzip,json,hashlib
os.environ.setdefault('MPLCONFIGDIR','/tmp/gse245908_mpl')
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection,LineCollection
from matplotlib.lines import Line2D
R=Path(__file__).resolve().parent;BASE=R.parent/'GSE245908_spatial_PGAM5';F=R/'figures';checks=[];D={}
TYPES=['All_CD8','CD8_exhaustion','CD8_activation','CD8_effector','CD8_naive_memory','CD4_other','Treg','T_unresolved'];LABELS=['All CD8 candidates','CD8 exhaustion','CD8 activation','CD8 effector','CD8 naive/memory','CD4 other','Treg candidates','T unresolved']
plt.rcParams.update({'font.size':10})
def save(fig,name):
 fig.savefig(F/(name+'.png'),dpi=200,bbox_inches='tight');fig.savefig(F/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
for s,gsm in [('CHC20','GSM7850822'),('CHC23','GSM7850823')]:
 root=BASE/'inputs'/s;ip=root/f'{gsm}_{s}_tissue_hires_image.png.gz';sp=root/f'{gsm}_{s}_scalefactors_json.json.gz';pp=root/f'{gsm}_{s}_tissue_positions_list.csv.gz'
 with gzip.open(ip,'rb') as z:im=np.array(Image.open(z).convert('RGB'))
 with gzip.open(sp,'rt') as z:scale=json.load(z)
 o=pd.read_csv(R/f'{s}_spot_subtypes.csv',index_col=0);pos=pd.read_csv(pp,header=None,index_col=0).loc[o.index];assert np.array_equal(pos[[4,5]].to_numpy(),o[['pixel_row','pixel_col']].to_numpy());sf=scale['tissue_hires_scalef'];xy=o[['pixel_col','pixel_row']].to_numpy()*sf;radius=scale['spot_diameter_fullres']*sf/2
 assert np.isfinite(xy).all() and (xy>=0).all() and (xy[:,0]<im.shape[1]).all() and (xy[:,1]<im.shape[0]).all()
 check=dict(sample=s,image_shape=list(im.shape),tissue_hires_scalef=sf,spot_radius_hires_pixels=radius,n_QC_coordinates=len(xy),coordinates_match_original_positions=True,all_coordinates_inside_image=True,transform='x=pixel_col*tissue_hires_scalef; y=pixel_row*tissue_hires_scalef; image origin upper; no rotation, mirroring or manual alignment',inputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ip,sp,pp]});checks.append(check)
 D[s]=(o,im,xy,radius)
def base(ax,s):
 o,im,xy,r=D[s];ax.imshow(im,origin='upper');pad=70;ax.set_xlim(max(0,xy[:,0].min()-pad),min(im.shape[1],xy[:,0].max()+pad));ax.set_ylim(min(im.shape[0],xy[:,1].max()+pad),max(0,xy[:,1].min()-pad));ax.set_aspect('equal');ax.set_axis_off();return o,xy,r
def circles(ax,xy,r,mask,color,fill=True):
 c=PatchCollection([Circle(p,r) for p in xy[np.asarray(mask)]],facecolors=color if fill else 'none',edgecolors=color,linewidths=.65 if fill else 1.1,alpha=.48 if fill else 1,zorder=3 if fill else 5);ax.add_collection(c)
for s in D:
 fig,axs=plt.subplots(2,4,figsize=(17,9))
 for ax,k,label in zip(axs.ravel(),TYPES,LABELS):
  o,xy,r=base(ax,s);t=o['gate_'+k];h=o.source_group.eq('PGAM5_high');circles(ax,xy,r,t,'#006bb6');circles(ax,xy,r*1.15,h,'#ed5c16',False);ax.set_title(f'{label}\nTargets={t.sum()} | High sources={h.sum()}',fontsize=11)
 fig.suptitle(s+' | T subtype RNA candidates fused with supplied H&E histology',fontsize=15);fig.legend(handles=[Line2D([],[],marker='o',ls='',color='#006bb6',label='T subtype candidate (translucent spot)'),Line2D([],[],marker='o',ls='',markerfacecolor='none',markeredgecolor='#ed5c16',color='#ed5c16',label='PGAM5-high myeloid RNA source (ring)')],loc='lower center',ncol=2);fig.subplots_adjust(top=.87,bottom=.07,hspace=.25,wspace=.05);save(fig,s+'_HE_T_subtypes')
fig,axs=plt.subplots(2,2,figsize=(11,11))
for row,s in enumerate(D):
 o,xy,r=base(axs[row,0],s);axs[row,0].set_title(s+' | supplied H&E tissue image');o,xy,r=base(axs[row,1],s)
 for g,c in [('Undetected','#555555'),('Other_detected','#0078ac'),('PGAM5_high','#ed5c16')]:circles(axs[row,1],xy,r,o.source_group.eq(g),c,False)
 axs[row,1].set_title(s+' | PGAM5 groups within myeloid RNA gate')
fig.legend(handles=[Line2D([],[],marker='o',ls='',markerfacecolor='none',color=c,label=g) for g,c in [('PGAM5 high','#ed5c16'),('Other detected','#0078ac'),('Undetected','#555555')]],loc='lower center',ncol=3);fig.suptitle('H&E registration and PGAM5 source distribution',fontsize=15);fig.subplots_adjust(top=.92,bottom=.06,wspace=.05,hspace=.15);save(fig,'HE_PGAM5_source_groups')
fig,axs=plt.subplots(2,2,figsize=(11,11))
for row,s in enumerate(D):
 o,_,_,_=D[s];e=pd.read_csv(R/f'{s}_subtype_LR_edges.csv.gz');lookup={b:i for i,b in enumerate(o.index)}
 for ax,pair in zip(axs[row],['SPP1--CD44','PVR--TIGIT']):
  o,xy,r=base(ax,s);keep=e.source_group.eq('PGAM5_high')&e.target_All_CD8&(e[pair]>0);d=e[keep];segments=[np.array([xy[lookup[x]],xy[lookup[y]]]) for x,y in zip(d.source_barcode,d.target_barcode)];ax.add_collection(LineCollection(segments,colors='#00d4dd',linewidths=1.8,zorder=6));circles(ax,xy,r,o.gate_All_CD8,'#006bb6');circles(ax,xy,r*1.15,o.source_group.eq('PGAM5_high'),'#ed5c16',False);ax.set_title(f'{s} | {pair}\nCodetected CD8 edges={len(d)}; source spots={d.source_barcode.nunique()}')
fig.legend(handles=[Line2D([],[],color='#00a8ad',lw=2,label='Adjacent codetected RNA candidate edge'),Line2D([],[],marker='o',ls='',color='#006bb6',label='CD8 candidate'),Line2D([],[],marker='o',ls='',markerfacecolor='none',color='#ed5c16',label='High source')],loc='lower center',ncol=3);fig.suptitle('Candidate LR edges fused with H&E | not causal communication',fontsize=14);fig.subplots_adjust(top=.9,bottom=.06,wspace=.08,hspace=.22);save(fig,'HE_CD8_candidate_LR_edges')
(R/'HE_registration_checks.json').write_text(json.dumps(checks,indent=2));print('PASS: original coordinates and image bounds checked for both sections; 4 PNG/PDF pairs saved')
