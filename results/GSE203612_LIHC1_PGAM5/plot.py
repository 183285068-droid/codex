from pathlib import Path
import json
import numpy as np,pandas as pd,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Circle,Patch
from matplotlib.collections import PatchCollection
from scipy.stats import rankdata
from association import residualizer
ROOT=Path(__file__).resolve().parent;FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
a=ad.read_h5ad(ROOT/'singlecell_marker_subset.h5ad');b=ad.read_h5ad(ROOT/'NYU_LIHC1_Vis_spatial_marker_subset.h5ad');rules=json.loads((ROOT/'annotation_rules.json').read_text())
cats=sorted(a.obs.expression_cluster.unique(),key=lambda x:int(x[1:]));colors={c:plt.get_cmap('tab20')(i) for i,c in enumerate(cats)}
labels={c:c+': '+v[0]+' | '+v[1] for c,v in rules['singlecell'].items()}
def save(fig,name):
    fig.savefig(FIG/f'{name}.png',dpi=180,bbox_inches='tight');fig.savefig(FIG/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
fig,axs=plt.subplots(1,3,figsize=(17,8));xy=a.obsm['X_umap']
for c in cats:
    mask=a.obs.expression_cluster.eq(c);axs[0].scatter(*xy[mask].T,c=[colors[c]],s=6,rasterized=True);axs[0].text(*np.median(xy[mask],axis=0),c,fontsize=8,weight='bold',bbox=dict(facecolor='white',alpha=.65,edgecolor='none',pad=1))
axs[0].set_title('14 RNA groups (provisional labels)')
v=a.obs.PGAM5_log1p_CP10k.to_numpy();ix=np.argsort(v);p=axs[1].scatter(*xy[ix].T,c=v[ix],s=7,cmap='YlOrRd',rasterized=True);plt.colorbar(p,ax=axs[1],label='PGAM5 log1p(CP10k)');axs[1].set_title('PGAM5 expression; zeros retained')
for lib,c in zip(['lib1','lib2','lib3'],['#487b9b','#c88054','#719164']):
    m=a.obs.library.eq(lib);axs[2].scatter(*xy[m].T,c=c,s=6,label=lib,rasterized=True)
axs[2].legend(frameon=False);axs[2].set_title('Technical library (same patient)')
for ax in axs:ax.set(xlabel='UMAP 1',ylabel='UMAP 2')
fig.suptitle('GSE203612 | NYU_LIHC1 | 1,783 QC single cells',fontsize=15)
fig.legend(handles=[Patch(facecolor=colors[c],label=labels[c]) for c in cats],loc='lower left',bbox_to_anchor=(.02,.055),ncol=3,frameon=False,fontsize=8)
fig.text(.03,.018,'RNA labels do not prove malignancy. C4: cycling lineage unresolved; C5: monocyte/DC continuum; C8: T/NK-related. No cross-patient replication.',fontsize=8)
fig.subplots_adjust(top=.88,bottom=.34,wspace=.35);save(fig,'singlecell_UMAP_PGAM5_libraries')
genes=['PGAM5','MYO19','ALB','CYP3A4','CYP2E1','GPC3','AFP','C1QA','C1QB','CD163','MRC1','CD5L','MARCO','FCN1','VCAN','CD1C','CD3D','TRAC','IL7R','CCL5','FOXP3','IKZF2','TIGIT','CTLA4','PDCD1','TOX','CXCL13','JCHAIN','MZB1','MS4A1','CD79A','VWF','KDR','COL1A2','NOTCH3','FCGR3B','S100A8','KIT','CPA3','MKI67','TOP2A']
genes=[g for g in genes if g in a.var_names];x=a[:,genes].X.toarray();raw=a[:,genes].layers['counts'].toarray();means=np.vstack([x[a.obs.expression_cluster.eq(c)].mean(0) for c in cats]);pct=np.vstack([(raw[a.obs.expression_cluster.eq(c)]>0).mean(0)*100 for c in cats])
xx,yy=np.meshgrid(range(len(genes)),range(len(cats)));fig,ax=plt.subplots(figsize=(19,8));p=ax.scatter(xx.ravel(),yy.ravel(),s=pct.ravel()*1.15,c=means.ravel(),cmap='YlOrRd',vmin=0,linewidths=.25,edgecolors='#999999')
ax.set_xticks(range(len(genes)),genes,rotation=90);ax.set_yticks(range(len(cats)),[labels[c] for c in cats]);ax.set_ylim(len(cats)-.5,-.5);ax.set_xlim(-.5,len(genes)-.5)
plt.colorbar(p,ax=ax,label='Mean log1p(CP10k), including zeros')
for pctv in [25,50,100]:ax.scatter([],[],s=pctv*1.15,c='#aaaaaa',label=f'{pctv}%')
ax.legend(loc='upper center',bbox_to_anchor=(.55,-.20),ncol=3,frameon=False);ax.set_title('NYU_LIHC1 single-cell marker dotplot | area = detected fraction; color = actual mean')
fig.tight_layout(rect=[0,.10,1,1]);save(fig,'singlecell_marker_dotplot')
fig,axs=plt.subplots(1,2,figsize=(12,5));d=pd.read_csv(ROOT/'singlecell_cluster_profiles.csv').set_index('cluster').loc[cats]
axs[0].barh(range(14),d.PGAM5_detection_pct,color=[colors[c] for c in cats]);axs[0].set_yticks(range(14),[c+': '+rules['singlecell'][c][1] for c in cats]);axs[0].invert_yaxis();axs[0].set(xlabel='PGAM5 RNA detected (%)',title='PGAM5 across RNA groups')
for k,r in enumerate(d.itertuples()):axs[0].text(r.PGAM5_detection_pct+.08,k,f'{r.PGAM5_detected}/{r.n_cells}',va='center',fontsize=8)
axs[0].set_xlim(0,d.PGAM5_detection_pct.max()+3)
mac=a.obs.expression_cluster.isin(['C3','C6']);axs[1].scatter(*xy[~mac].T,c='#d3d7dd',s=5,rasterized=True)
for c in ['C3','C6']:
    m=a.obs.expression_cluster.eq(c);axs[1].scatter(*xy[m].T,c=[colors[c]],s=8,label=c+': '+rules['singlecell'][c][0],rasterized=True)
positive=a.obs.Mac_RNA_cluster_confirmed&a.obs.PGAM5_raw.gt(0);axs[1].scatter(*xy[positive].T,facecolors='none',edgecolors='#bf2537',s=60,linewidths=1,label='6 PGAM5-detected lineage-supported candidates')
axs[1].legend(loc='lower right',fontsize=7,frameon=False);axs[1].set(xlabel='Original all-cell UMAP 1',ylabel='Original all-cell UMAP 2',title='Macrophage-like groups; all 6 detections = 1 UMI')
fig.text(.03,.01,'C3/C6: 325 cells, 8 PGAM5 detections; additional lineage gate: 246 candidates, 6 detections. This is not evidence for a stable PGAM5-positive subtype.',fontsize=8)
fig.tight_layout(rect=[0,.04,1,1]);save(fig,'PGAM5_detection_and_macrophage_candidates')
im=b.uns['spatial']['NYU_LIHC1_Vis']['images']['lowres'];sf=b.uns['scalefactors'];sxy=b.obsm['spatial']*sf['tissue_lowres_scalef'];rad=sf['spot_diameter_fullres']*sf['tissue_lowres_scalef']/2

def base(ax):ax.imshow(im,alpha=.7);ax.set_xlim(0,im.shape[1]);ax.set_ylim(im.shape[0],0);ax.axis('off')
def draw(ax,mask,values=None,colors=None,vmax=None):
    ids=np.flatnonzero(np.asarray(mask));p=PatchCollection([Circle(sxy[i],rad) for i in ids],linewidths=0,alpha=.95,rasterized=True)
    if colors is not None:p.set_facecolor(colors)
    else:p.set_array(np.asarray(values)[ids]);p.set_cmap('YlOrRd');p.set_norm(Normalize(0,vmax))
    ax.add_collection(p);return p

def continuous(ax,col,title,mask=None):
    base(ax);m=np.ones(len(b),bool) if mask is None else np.asarray(mask);draw(ax,~m,colors='#cbd0d5')
    vm=max(.01,float(b.obs.loc[m,col].quantile(.99)));p=draw(ax,m,b.obs[col],vmax=vm);plt.colorbar(p,ax=ax,fraction=.032,pad=.015,extend='max');ax.set_title(title,fontsize=10)
fig,axs=plt.subplots(3,3,figsize=(13,13));base(axs[0,0])
spcats=['S0','S1','S2']
for k,c in enumerate(spcats):draw(axs[0,0],b.obs.expression_cluster.eq(c),colors=plt.get_cmap('tab10')(k))
axs[0,0].set_title('Three section-specific mixed RNA regions')
for ax,col,title,mask in zip(list(axs.ravel())[1:],['PGAM5_log1p_CP10k','Macrophage_core','T_lineage','Exhaustion_core','Exhaustion_no_CTLA4_TIGIT','Activation_early','Treg_related','IFN_response'],['PGAM5 RNA (all spots)','Macrophage-core RNA','T-lineage RNA','Exhaustion-core RNA (T gate)','Exhaustion without CTLA4/TIGIT (T gate)','Early activation RNA (T gate)','Treg-related RNA (T gate)','Interferon-response RNA'],[None,None,None,b.obs.T_RNA_enriched,b.obs.T_RNA_enriched,b.obs.T_RNA_enriched,b.obs.T_RNA_enriched,None]):continuous(ax,col,title,mask)
fig.suptitle('GSE203612 | NYU_LIHC1 Visium | 1,634 retained spots | genuine H&E coordinates',fontsize=14)
fig.legend(handles=[Patch(facecolor=plt.get_cmap('tab10')(k),label=c+': '+rules['spatial'][c][0]+' | '+rules['spatial'][c][1]) for k,c in enumerate(spcats)],loc='lower left',bbox_to_anchor=(.04,.05),fontsize=8,frameon=False)
fig.text(.05,.015,'Colors: PGAM5 log1p(CP10k); programs: mean log1p(CP10k), zeros retained. Each panel has its own Q99 display limit. Grey = outside T gate.\nMixed spots: signal source and causal communication are unresolved. Labels do not classify every spot as malignant.',fontsize=8)
fig.subplots_adjust(top=.94,bottom=.15,hspace=.16,wspace=.25);save(fig,'Visium_spatial_RNA_maps')
# Original spatial-cluster expression embedding and marker dotplot use identical dot conventions.
gs=[g for g in ['PGAM5','MYO19','ALB','CYP3A4','GLUL','APOE','CD74','HLA-DRB1','C1QA','CD68','CD3D','TRAC','FOXP3','PDCD1','TIGIT','CTLA4','TOX','JCHAIN','MS4A1','DCN','COL1A1','VWF','MKI67'] if g in b.var_names]
x=b[:,gs].X.toarray();raw=b[:,gs].layers['counts'].toarray();means=np.vstack([x[b.obs.expression_cluster.eq(c)].mean(0) for c in spcats]);pct=np.vstack([(raw[b.obs.expression_cluster.eq(c)]>0).mean(0)*100 for c in spcats])
fig,(ax,dx)=plt.subplots(1,2,figsize=(14,5),gridspec_kw={'width_ratios':[1,2]});z=b.obsm['X_umap']
for k,c in enumerate(spcats):
    m=b.obs.expression_cluster.eq(c);ax.scatter(*z[m].T,c=[plt.get_cmap('tab10')(k)],s=5,label=c,rasterized=True)
ax.legend(frameon=False);ax.set(xlabel='UMAP 1',ylabel='UMAP 2',title='Spot expression embedding (not tissue coordinates)')
xx,yy=np.meshgrid(range(len(gs)),range(3));p=dx.scatter(xx.ravel(),yy.ravel(),s=pct.ravel()*1.2,c=means.ravel(),cmap='YlOrRd',vmin=0,linewidths=.25,edgecolors='#999999');dx.set_xticks(range(len(gs)),gs,rotation=90);dx.set_yticks(range(3),[c+': '+rules['spatial'][c][0] for c in spcats]);dx.set_ylim(2.5,-.5);plt.colorbar(p,ax=dx,label='Mean log1p(CP10k)')
for v in [25,50,100]:dx.scatter([],[],s=v*1.2,c='#aaaaaa',label=f'{v}%')
dx.legend(loc='upper center',bbox_to_anchor=(.5,-.32),ncol=3,frameon=False);dx.set_title('Dot area: detected %; color: actual mean including zeros');fig.subplots_adjust(top=.86,bottom=.30,wspace=.38);save(fig,'Visium_UMAP_marker_dotplot')
tab=pd.read_csv(ROOT/'association_statistics.csv',float_precision='round_trip');data=pd.read_csv(ROOT/'NYU_LIHC1_Vis_source_neighbor_RNA.csv',float_precision='round_trip');data=data[data.T_neighbor_count>0]
fig,(ax,bx)=plt.subplots(1,2,figsize=(13,5));q,_=residualizer(data);xr=rankdata(data.PGAM5_log1p_CP10k);yr=rankdata(data.neighbor_Exhaustion_core);xr-=q@(q.T@xr);yr-=q@(q.T@yr);xr/=xr.std(ddof=1);yr/=yr.std(ddof=1)
for group,col in [('Undetected','#8794a0'),('Other_detected','#df9b43'),('PGAM5_high','#bd3445')]:
    m=data.PGAM5_mac_group.eq(group);ax.scatter(xr[m],yr[m],s=15,c=col,label=group,rasterized=True)
r=tab.iloc[0];xx=np.linspace(xr.min(),xr.max(),50);ax.plot(xx,r.statistic*xx,c='#444444',lw=1);ax.legend(fontsize=8,frameon=False)
ax.set(xlabel='Source PGAM5 adjusted rank / SD',ylabel='Neighbor exhaustion adjusted rank / SD',title=f'n=163; partial rho={r.statistic:.4f}; p={r.permutation_p:.4f}')
t=tab[(tab.cohort=='all_mac_including_zero')&(tab.contrast=='partial_Spearman')]
for i,row in enumerate(t.itertuples()):
    bx.barh(i,row.statistic,color='#bd3445' if row.statistic>0 else '#337ea0');bx.text(.10,i,f'p={row.permutation_p:.3g}; q={row.BH_q_7_planned_tests:.2g}',va='center',fontsize=9)
bx.set_yticks(range(len(t)),['Exhaustion core','Expanded exhaustion','Early activation','Cytotoxic effector','No CTLA4/TIGIT']);bx.invert_yaxis();bx.set_xlim(-.01,.34);bx.set_xlabel('Adjusted partial Spearman rho');bx.set_title('4,999 spatial / UMI conditional permutations')
fig.suptitle('NYU_LIHC1: PGAM5 in Mac-RNA sources vs neighboring T-RNA programs',fontsize=13)
fig.text(.03,.015,'9 quality / composition covariates; PGAM5 zeros retained. q: 7 planned tests (2 untestable). One patient, no causal or macrophage-specific RNA attribution.\nHigh vs other detected test not run: only 19 detected sources with T neighbors (4 high, 15 other); raw co-localization does not establish exhaustion promotion.',fontsize=8)
fig.tight_layout(rect=[0,.10,1,.92]);save(fig,'Visium_PGAM5_T_exhaustion_association')
fig,axs=plt.subplots(1,2,figsize=(10,6));continuous(axs[0],'PGAM5_log1p_CP10k','PGAM5: Mac-RNA source spots with T neighbor',b.obs.Mac_RNA_enriched&b.obs.T_neighbor_count.gt(0));continuous(axs[1],'Exhaustion_core','Exhaustion-core RNA: T-RNA spots',b.obs.T_RNA_enriched)
for i in np.flatnonzero(b.obs.PGAM5_mac_group.eq('PGAM5_high')&b.obs.T_neighbor_count.gt(0)):axs[0].add_patch(Circle(sxy[i],rad*1.25,facecolor='none',edgecolor='#21b8c4',lw=1))
fig.text(.03,.025,'Cyan rings: 4 high-PGAM5 Mac-RNA sources with T neighbors. None has an exhaustion-supported T-RNA neighbor under the fixed gate.\nSpatial RNA scores are not cell percentages or proven cell-to-cell interactions.',fontsize=8);fig.subplots_adjust(bottom=.15,wspace=.2);save(fig,'Visium_PGAM5_Tex_spatial_neighborhoods')
print('Generated 7 PNG/PDF pairs',flush=True)
