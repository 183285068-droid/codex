"""Standalone PNG/PDF figures with genuine Space Ranger coordinates."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE245908_mpl')
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Patch
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize,ListedColormap
import numpy as np,pandas as pd,scanpy as sc,json
from scipy.stats import spearmanr
import gzip
from PIL import Image
from io import BytesIO
from preprocess import ROOT,WORK,SAMPLES,PROGRAMS,array
from annotate_regions import LABELS
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':13,'axes.labelsize':11,'pdf.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white'})
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
COLORS=['#2166AC','#E6AB02','#1B9E77','#B358A6','#D95F02','#7570B3']
DATA={s:sc.read_h5ad(ROOT/f'{s}_spatial_marker_subset.h5ad') for s in SAMPLES}
def save(fig,name):
 fig.savefig(OUT/f'{name}.png',dpi=220,bbox_inches='tight');fig.savefig(OUT/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
def image_axis(ax,s):
 a=DATA[s];meta=a.uns['spatial'][s];sf=meta['scalefactors']['tissue_hires_scalef'];xy=a.obsm['spatial']*sf
 if 'hires' in meta['images']:im=meta['images']['hires']
 else:im=np.asarray(Image.open(BytesIO(gzip.decompress((ROOT/'inputs'/s/f'{SAMPLES[s]}_{s}_tissue_hires_image.png.gz').read_bytes()))))
 ax.imshow(im,alpha=.68)
 ax.set_xlim(xy[:,0].min()-35,xy[:,0].max()+35);ax.set_ylim(xy[:,1].max()+35,xy[:,1].min()-35);ax.set_aspect('equal');ax.axis('off')
 return xy,meta['scalefactors']['spot_diameter_fullres']*sf/2

def spatial(ax,s,values=None,title=None,vmin=0,vmax=None,cmap='YlOrRd',categorical=False,zero_grey=False):
 xy,r=image_axis(ax,s)
 if values is not None:
  coll=PatchCollection([Circle(p,r) for p in xy],linewidth=0,alpha=.88,rasterized=True)
  if categorical:coll.set_facecolor([COLORS[int(str(v)[1:])] for v in values])
  else:
   cm=plt.get_cmap(cmap);norm=Normalize(vmin=vmin,vmax=vmax);co=cm(norm(values))
   if zero_grey:co[np.asarray(values)==0]=[.73,.75,.78,.9]
   coll.set_facecolor(co)
  ax.add_collection(coll)
  if not categorical:
   sm=plt.cm.ScalarMappable(norm=norm,cmap=cm);return sm
 if title:ax.set_title(title)

def cluster_figures():
 for s,a in DATA.items():
  fig,axs=plt.subplots(1,3,figsize=(18,6));spatial(axs[0],s,title=f'{s}: H&E tissue image')
  spatial(axs[1],s,a.obs.cluster,title=f'{s}: expression regions',categorical=True);axs[1].set_title(f'{s}: spatial expression regions')
  xy=a.obsm['X_umap'];axs[2].scatter(xy[:,0],xy[:,1],c=[COLORS[int(v[1:])] for v in a.obs.cluster],s=6,alpha=.75,rasterized=True,linewidths=0)
  for c in LABELS[s]:
   med=np.median(xy[a.obs.cluster==c],axis=0);axs[2].text(*med,c,ha='center',va='center',fontweight='bold',bbox={'fc':'white','ec':'none','alpha':.8,'pad':2})
  axs[2].set_xlabel('UMAP 1');axs[2].set_ylabel('UMAP 2');axs[2].set_title(f'{s}: UMAP of mixed spots');axs[2].spines[['top','right']].set_visible(False)
  leg=[Patch(fc=COLORS[int(c[1:])],label=f'{c}: {v[0]} | {v[1]}') for c,v in LABELS[s].items()]
  fig.legend(handles=leg,loc='lower center',ncol=2,bbox_to_anchor=(.5,-.04),frameon=False,fontsize=10)
  fig.suptitle(f'GSE245908 {s}: {len(a):,} QC-passing spots | provisional RNA-region annotations',fontsize=16)
  fig.subplots_adjust(bottom=.22,wspace=.19,top=.87);save(fig,f'{s}_H_E_spatial_clusters_UMAP')

def QCfigure():
 fig,axs=plt.subplots(2,3,figsize=(14,8))
 for row,s in enumerate(SAMPLES):
  d=pd.read_csv(ROOT/f'{s}_all_spot_QC.csv')
  for col,(k,label) in enumerate([('total_counts','Total UMI'),('n_genes','Detected genes'),('pct_mt','Mitochondrial RNA (%)')]):
   axs[row,col].hist(d[k],bins=45,color='#2166AC',alpha=.85);axs[row,col].set_xlabel(label);axs[row,col].set_ylabel('Spots');axs[row,col].set_title(f'{s}: {d.pass_qc.sum():,}/{len(d):,} retained');axs[row,col].spines[['top','right']].set_visible(False)
   cutoff={'total_counts':500,'n_genes':200,'pct_mt':25}[k];axs[row,col].axvline(cutoff,c='#D95F02',ls='--',lw=1.5)
 fig.suptitle('GSE245908: spot QC (UMI >= 500, genes >= 200, MT <= 25%)',fontsize=16);fig.tight_layout(rect=[0,0,1,.95]);save(fig,'QC_distributions')

def PGAM5figure():
 maxraw=max(a.obs.PGAM5_raw.max() for a in DATA.values());maxlog=max(a.obs.PGAM5_log1p_CP10k.max() for a in DATA.values());maxmean=max(pd.read_csv(ROOT/f'{s}_cluster_profiles.csv').PGAM5_mean_log1p_CP10k.max() for s in SAMPLES)
 fig,axs=plt.subplots(2,3,figsize=(15,10))
 for row,(s,a) in enumerate(DATA.items()):
  d=pd.read_csv(ROOT/f'{s}_cluster_profiles.csv').set_index('cluster');means=a.obs.cluster.map(d.PGAM5_mean_log1p_CP10k).astype(float)
  for col,(values,label,mx) in enumerate([(a.obs.PGAM5_raw,'PGAM5 raw UMI',maxraw),(a.obs.PGAM5_log1p_CP10k,'PGAM5 log1p(CP10k)',maxlog),(means,'PGAM5 region mean log1p(CP10k)',maxmean)]):
   sm=spatial(axs[row,col],s,values,vmax=mx,zero_grey=True);axs[row,col].set_title(f'{s}: {label}',fontsize=12);fig.colorbar(sm,ax=axs[row,col],fraction=.038,pad=.02)
   if col==2:
    xy=a.obsm['spatial']*a.uns['spatial'][s]['scalefactors']['tissue_hires_scalef']
    for c in LABELS[s]:
     med=np.median(xy[a.obs.cluster==c],axis=0);axs[row,col].text(*med,c,ha='center',fontsize=10,fontweight='bold',bbox={'fc':'white','ec':'none','alpha':.85,'pad':1})
 fig.suptitle('GSE245908: PGAM5 spatial expression | shared scales across slices',fontsize=16)
 fig.text(.02,.015,'Grey = no PGAM5 detected. Region means include zeros and are assigned to every spot in that region.\nSpot mixtures cannot identify which cell type expresses PGAM5. S0-S5 names are supplied in region_names.csv.',fontsize=11)
 fig.subplots_adjust(left=.03,right=.97,top=.9,bottom=.09,wspace=.18,hspace=.16);save(fig,'PGAM5_spatial_expression')

def program_figures():
 selected=['Macrophage_core','T_cell','Cytotoxic_T_NK','Fibroblast_matrix','Hepatocyte','THBS1_TREM1_myeloid']
 for start in [0,3]:
  subset=selected[start:start+3];fig,axs=plt.subplots(2,3,figsize=(15,10))
  for col,k in enumerate(subset):
   mx=max(a.obs[k].max() for a in DATA.values())
   for row,(s,a) in enumerate(DATA.items()):
    sm=spatial(axs[row,col],s,a.obs[k],vmax=mx,cmap='viridis');axs[row,col].set_title(f'{s}: {k.replace("_"," ")}',fontsize=12);fig.colorbar(sm,ax=axs[row,col],fraction=.038,pad=.02)
  fig.suptitle('GSE245908: spatial RNA programs | common scale for each program',fontsize=16)
  fig.text(.03,.02,'Program = mean gene log1p(CP10k), including zeros. Scores are RNA signals, not cell fractions.\nTHBS1/TREM1 panel is exploratory; its genes also occur outside myeloid cells.',fontsize=11)
  fig.subplots_adjust(left=.03,right=.97,top=.9,bottom=.095,wspace=.18,hspace=.15);save(fig,f'spatial_programs_{start//3+1}')

def abundancebars():
 for s in SAMPLES:
  d=pd.read_csv(ROOT/f'{s}_cluster_profiles.csv');fig,axs=plt.subplots(1,2,figsize=(16,6))
  names=[f'{r.cluster}: {LABELS[s][r.cluster][0]}\n{LABELS[s][r.cluster][1]}' for r in d.itertuples()];y=np.arange(len(d))
  for ax,k,label in [(axs[0],'PGAM5_mean_log1p_CP10k','Mean PGAM5 log1p(CP10k), including zeros'),(axs[1],'PGAM5_detection_pct','PGAM5-detected spots (%)')]:
   ax.barh(y,d[k],color=[COLORS[i] for i in y],height=.66);ax.set_yticks(y,labels=names if ax==axs[0] else d.cluster);ax.invert_yaxis();ax.set_xlabel(label);ax.set_xlim(0,d[k].max()*1.24);ax.spines[['top','right']].set_visible(False)
   for i,r in enumerate(d.itertuples()):ax.text(getattr(r,k)+d[k].max()*.02,i,f'{getattr(r,k):.3f}' if k.endswith('CP10k') else f'{getattr(r,k):.1f}% ({r.PGAM5_detected_spots}/{r.n_spots})',va='center',fontsize=10)
  fig.suptitle(f'{s}: PGAM5 expression by RNA region',fontsize=16);fig.tight_layout(rect=[0,0,1,.95]);save(fig,f'{s}_PGAM5_by_region')

def dotplots():
 values=[]
 for s in SAMPLES:
  a=sc.read_h5ad(WORK/f'{s}_full_processed.h5ad')
  genes=list(dict.fromkeys([g for v in LABELS[s].values() for g in v[0].split('/')]+['PGAM5','C1QA','CD68','CSF1R','LST1','LYZ','CD3D','TRAC','NKG7','COL1A1','ALB','GPC3','EPCAM','MKI67']))
  genes=[g for g in genes if g in a.var_names];x=array(a[:,genes].X);cs=list(LABELS[s]);means=np.array([x[a.obs.cluster==c].mean(0) for c in cs]);rates=np.array([(x[a.obs.cluster==c]>0).mean(0) for c in cs]);scaled=(means-means.min(0))/np.where(np.ptp(means,axis=0)>0,np.ptp(means,axis=0),1)
  for i,c in enumerate(cs):
   for j,g in enumerate(genes):values.append({'sample':s,'cluster':c,'gene':g,'mean_log1p_CP10k':float(means[i,j]),'detection_pct':float(rates[i,j]*100),'within_gene_scaled_mean':float(scaled[i,j])})
  fig,ax=plt.subplots(figsize=(17,5));yy,xx=np.indices(means.shape);dots=ax.scatter(xx.ravel(),yy.ravel(),s=(rates.ravel()*170)+1,c=scaled.ravel(),cmap='YlOrRd',vmin=0,vmax=1,edgecolors='#666666',linewidths=.2)
  ax.set_xticks(np.arange(len(genes)),labels=genes,rotation=60,ha='right');ax.set_yticks(np.arange(len(cs)),labels=[f'{c}: {LABELS[s][c][0]}' for c in cs]);ax.set_xlim(-.7,len(genes)-.3);ax.set_ylim(len(cs)-.4,-.6);ax.set_title(f'{s}: region markers | mixed-spot expression',fontsize=16);ax.spines[['top','right']].set_visible(False)
  fig.colorbar(dots,ax=ax,pad=.02,fraction=.015,label='Within-gene scaled mean (0-1)');handles=[ax.scatter([],[],s=v*170+1,c='white',edgecolors='#666666',label=f'{v*100:.0f}%') for v in [.25,.5,1]];ax.legend(handles=handles,title='Detected spots',loc='upper left',bbox_to_anchor=(1.09,1),frameon=False)
  fig.tight_layout();save(fig,f'{s}_region_marker_dotplot')
 pd.DataFrame(values).to_csv(ROOT/'region_marker_dotplot_values.csv',index=False)

def associations():
 d=pd.read_csv(ROOT/'PGAM5_spatial_associations.csv');keys=list(PROGRAMS);cols=[(s,k) for s in SAMPLES for k in ['same_spot','adjacent_spot']]
 matrix=np.zeros((len(keys),len(cols)));q=matrix.copy()
 for i,k in enumerate(keys):
  for j,(s,rel) in enumerate(cols):
   r=d[(d['sample']==s)&(d.program==k)&(d.relationship==rel)].iloc[0];matrix[i,j]=r.QC_adjusted_partial_rank_rho;q[i,j]=r.spatial_BH_q
 fig,ax=plt.subplots(figsize=(10,9));im=ax.imshow(matrix,cmap='PuOr_r',vmin=-.25,vmax=.25,aspect='auto');ax.set_xticks(range(len(cols)),labels=[f'{s}\n'+('Same spot' if k=='same_spot' else 'Adjacent spots') for s,k in cols]);ax.set_yticks(range(len(keys)),labels=[k.replace('_',' ') for k in keys])
 for i in range(len(keys)):
  for j in range(len(cols)):ax.text(j,i,f'{matrix[i,j]:+.3f}\nq={q[i,j]:.3f}',ha='center',va='center',fontsize=9,color='white' if abs(matrix[i,j])>.19 else '#111111',fontweight='bold' if q[i,j]<.05 else 'normal')
 fig.colorbar(im,ax=ax,fraction=.04,pad=.03,label='QC-adjusted partial rank correlation')
 ax.set_title('PGAM5 vs spatial RNA programs',fontsize=16,pad=18);fig.text(.02,.018,'BH q across all 48 tests; model-based graph-spectral randomization (4,999 nulls).\nRanks adjusted for UMI depth, detected genes and mitochondrial RNA. Two slices; exploratory associations.',fontsize=10);fig.tight_layout(rect=[0,.065,1,1]);save(fig,'PGAM5_spatial_association_heatmap')
 fig,axs=plt.subplots(2,3,figsize=(14,8));selected=['Macrophage_core','T_cell','Cytotoxic_T_NK']
 for row,(s,a) in enumerate(DATA.items()):
  a=a[a.obs.spatial_degree>0]
  for col,k in enumerate(selected):
   ax=axs[row,col];x=a.obs[k];y=a.obs.PGAM5_log1p_CP10k;ax.scatter(x,y,s=5,c='#2166AC',alpha=.23,rasterized=True,linewidths=0);r=d[(d['sample']==s)&(d.program==k)&(d.relationship=='same_spot')].iloc[0];ax.set_title(f'{s}: raw rho={r.raw_Spearman_rho:+.3f}\nQC-adjusted rho={r.QC_adjusted_partial_rank_rho:+.3f}, spatial q={r.spatial_BH_q:.3f}',fontsize=11);ax.set_xlabel(k.replace('_',' ')+' program');ax.set_ylabel('PGAM5 log1p(CP10k)');ax.spines[['top','right']].set_visible(False)
 fig.suptitle('PGAM5 and immune RNA programs | graph-connected QC spots, including zeros',fontsize=16);fig.tight_layout(rect=[0,0,1,.93]);save(fig,'PGAM5_immune_signal_scatter')
if __name__=='__main__':
 QCfigure();cluster_figures();PGAM5figure();program_figures();abundancebars();dotplots();associations();print('Figures complete',flush=True)
