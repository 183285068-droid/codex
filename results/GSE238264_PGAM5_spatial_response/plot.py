from pathlib import Path
import os,json
os.environ.setdefault('MPLCONFIGDIR','/tmp/gse245908_mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt,numpy as np,pandas as pd
from PIL import Image
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from analyze import R,INPUT,SAMPLES,PRIMARY
F=R/'figures';d=pd.read_csv(R/'patient_spatial_metrics.csv');tests=pd.read_csv(R/'response_association_tests.csv');COLORS={'responder':'#2678a8','non-responder':'#b95e27'};plt.rcParams.update({'font.size':11})
def save(fig,name):
 fig.savefig(F/(name+'.png'),dpi=180,bbox_inches='tight');fig.savefig(F/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
fig,axs=plt.subplots(2,2,figsize=(13,9));titles=['PGAM5 detected / macrophage RNA spots','Macrophage RNA neighbor coverage','T RNA neighbor coverage','Neighbor fibroblast RNA program']
for ax,metric,title in zip(axs.ravel(),PRIMARY,titles):
 factor=100 if metric.endswith('fraction') else 1
 for i,response in enumerate(['responder','non-responder']):
  dd=d[(d.response==response)&d[metric].notna()];x=i+np.linspace(-.12,.12,len(dd));ax.scatter(x,dd[metric]*factor,c=COLORS[response],s=60)
  for xx,row in zip(x,dd.itertuples()):ax.annotate(row.sample,(xx,getattr(row,metric)*factor),xytext=(5,5),textcoords='offset points',fontsize=9)
  ax.hlines(dd[metric].mean()*factor,i-.18,i+.18,color=COLORS[response],lw=2)
 row=tests[tests.metric==metric].iloc[0];ax.set_xticks([0,1],[f'Responder ({int(((d.response=="responder")&d[metric].notna()).sum())}/4 available)',f'Non-responder ({int(((d.response=="non-responder")&d[metric].notna()).sum())}/3 available)']);ax.set_xlim(-.45,1.5);ax.set_ylim(bottom=0);ax.set_ylabel('%' if factor==100 else 'Mean log1p(CP10k)');ax.set_title(title+f'\nExact p={row.p_exact:.3f}; BH q={row.BH_q:.3f}' if pd.notna(row.p_exact) else title+'\nInsufficient high-source support for complete 7-patient test')
fig.suptitle('GSE238264 | post-neoadjuvant spatial RNA features vs pathological response\nEach dot is one patient; horizontal lines are group means');fig.tight_layout();save(fig,'patient_response_primary_comparisons')
checks=[]
for sample in SAMPLES:
 o=pd.read_csv(R/f'{sample}_spot_annotations.csv.gz',index_col=0);root=INPUT/sample;image=Image.open(root/'spatial/tissue_hires_image.png');im=np.array(image.convert('RGB'));scale=json.loads((root/'spatial/scalefactors_json.json').read_text());sf=scale['tissue_hires_scalef'];xy=o[['pixel_col','pixel_row']].to_numpy()*sf;rad=scale['spot_diameter_fullres']*sf/2;assert (xy>=0).all() and (xy[:,0]<image.width).all() and (xy[:,1]<image.height).all();checks.append(dict(sample=sample,scale=sf,spot_radius_hires_pixels=rad,image_dimensions=list(image.size),all_coordinates_inside_image=True))
 fig,axs=plt.subplots(2,2,figsize=(11,10));high=o.source_group.eq('PGAM5_high');pad=40
 def base(ax):
  ax.imshow(im,origin='upper');ax.set_xlim(max(0,xy[:,0].min()-pad),min(image.width,xy[:,0].max()+pad));ax.set_ylim(min(image.height,xy[:,1].max()+pad),max(0,xy[:,1].min()-pad));ax.set_aspect('equal');ax.set_axis_off()
 def circle(ax,mask,c,fill=False):ax.add_collection(PatchCollection([Circle(p,rad) for p in xy[np.asarray(mask)]],facecolors=c if fill else 'none',edgecolors=c,alpha=.48 if fill else 1,linewidths=.7 if fill else 1.1,zorder=4))
 for ax in axs.ravel():base(ax)
 axs[0,0].set_title('Supplied H&E tissue image')
 for g,c in [('Undetected','#555'),('Other_detected','#2678a8'),('PGAM5_high','#ed5c16')]:circle(axs[0,1],o.source_group.eq(g),c)
 axs[0,1].set_title(f'Myeloid RNA sources | High={high.sum()}')
 circle(axs[1,0],o.T_RNA,'#2678a8',True);circle(axs[1,0],high,'#ed5c16');axs[1,0].set_title('T RNA candidates (blue) + high source rings')
 fib=o.score_Fibroblast.to_numpy();strong=(fib>=np.quantile(fib,.75))&(o.detected_Fibroblast>=2)&(fib>0);circle(axs[1,1],strong,'#765f97',True);circle(axs[1,1],high,'#ed5c16');axs[1,1].set_title('Fibroblast RNA candidates (purple) + high rings')
 response=d.set_index('sample').loc[sample,'response'];fig.suptitle(sample+' | '+response+' | post-neoadjuvant HCC',fontsize=15);fig.legend(handles=[Line2D([],[],marker='o',ls='',markerfacecolor='none',color=c,label=g) for g,c in [('PGAM5 high','#ed5c16'),('Other detected','#2678a8'),('Undetected','#555')]],loc='lower center',ncol=3);fig.subplots_adjust(top=.91,bottom=.07,hspace=.18,wspace=.08);save(fig,sample+'_HE_spatial_features')
(R/'HE_registration_checks.json').write_text(json.dumps(checks,indent=2));print('8 PNG/PDF pairs saved, all HE coordinates in bounds')
