"""Inspectably source-backed Visium/H&E maps; expression groups are mixed RNA regions."""
from pathlib import Path
import gzip,io
import numpy as np
import pandas as pd
import anndata as ad
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Circle, Patch
from matplotlib.collections import PatchCollection
from preprocess import ROOT,WORK,SAMPLES,PROGRAMS
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
IMG=ROOT/'tissue_images';IMG.mkdir(exist_ok=True)
DATA={};LABELS=[]
for item in SAMPLES:
    s=item['sample'];a=ad.read_h5ad(ROOT/f'{s}_spatial_marker_subset.h5ad')
    path=WORK/'inputs'/f'{item["gsm"]}_{s}_tissue_lowres_image.png.gz'
    im=np.asarray(Image.open(io.BytesIO(gzip.decompress(path.read_bytes()))).convert('RGB'))
    Image.fromarray(im).save(IMG/f'{s}_tissue_lowres_image.png')
    # Conventional portable spatial object, with original full-resolution coordinates.
    a.uns['spatial']={s:{'images':{'lowres':im},'scalefactors':dict(a.uns['scalefactors']),
                       'metadata':{'source':'GEO GSE281759','sample':s}}}
    a.write_h5ad(ROOT/f'{s}_spatial_marker_subset.h5ad',compression='gzip')
    markers=pd.read_csv(ROOT/f'{s}_top100_markers.csv');names={}
    for c in a.obs.expression_cluster.cat.categories:
        d=markers[(markers.group==c)&markers.logfoldchanges.gt(0)&markers.names.isin(a.var_names)]
        d=d[~d.names.str.match(r'^(MT-|RP[SL]\d)')]
        gs=d.names.head(2).tolist();names[c]=c+': '+('/'.join(gs) if gs else 'Mixed RNA')
        LABELS.append(dict(sample=s,patient=item['patient'],cluster=c,marker_pair='/'.join(gs),annotation='Mixed-spot RNA region; provisional',n_spots=int(a.obs.expression_cluster.eq(c).sum())))
    DATA[s]=(a,im,names)
pd.DataFrame(LABELS).to_csv(ROOT/'expression_group_names.csv',index=False)
VARS=['PGAM5_log1p_CP10k','Macrophage_core','T_lineage','Activation_early','Exhaustion_core']
VMAX={v:max(float(a.obs[v].quantile(.99)) for a,_,_ in DATA.values()) for v in VARS}
# Only display clipping for readability; all numeric results retain full values.
def backdrop(ax,a,im):
    ax.imshow(im,alpha=.70);ax.set_xlim(0,im.shape[1]);ax.set_ylim(im.shape[0],0);ax.axis('off')
def dots(ax,a,mask,colors,cmap=None,norm=None,alpha=.9):
    sf=a.uns['scalefactors'];xy=a.obsm['spatial']*sf['tissue_lowres_scalef']
    radius=sf['spot_diameter_fullres']*sf['tissue_lowres_scalef']/2
    ix=np.flatnonzero(mask)
    coll=PatchCollection([Circle(xy[i],radius) for i in ix],linewidths=0,alpha=alpha,rasterized=True)
    if cmap is None:coll.set_facecolor(colors)
    else:coll.set_array(np.asarray(colors));coll.set_cmap(cmap);coll.set_norm(norm)
    ax.add_collection(coll);return coll
def continuous(ax,a,im,var,title,mask=None):
    backdrop(ax,a,im);m=np.ones(len(a),bool) if mask is None else np.asarray(mask)
    if mask is not None:dots(ax,a,~m,'#cbd0d5',alpha=.65)
    coll=dots(ax,a,m,a.obs[var].to_numpy()[m],cmap='YlOrRd',norm=Normalize(0,VMAX[var]))
    ax.set_title(title,fontsize=10);plt.colorbar(coll,ax=ax,fraction=.035,pad=.015,extend='max').set_label('Mean log1p(CP10k)' if var!='PGAM5_log1p_CP10k' else 'log1p(CP10k)',fontsize=8)
def clusters(ax,a,im,names):
    backdrop(ax,a,im);cats=list(a.obs.expression_cluster.cat.categories)
    palette=plt.get_cmap('tab10')
    for k,c in enumerate(cats):
        m=a.obs.expression_cluster.eq(c).to_numpy();dots(ax,a,m,palette(k))

    ax.set_title('Within-section expression groups')
def proximity(ax,a,im):
    backdrop(ax,a,im);dots(ax,a,np.ones(len(a),bool),'#cbd0d5',alpha=.6)
    for mask,col in [(a.obs.Exhaustion_core_supported,'#7455a5'),
                     (a.obs.Mac_RNA_enriched&a.obs.T_exhaustion_neighbor_count.gt(0),'#2e92ad'),
                     (a.obs.PGAM5_mac_group.eq('PGAM5_high')&a.obs.T_exhaustion_neighbor_count.gt(0),'#bf2537')]:dots(ax,a,np.asarray(mask),col)
    ax.set_title('Mac RNA / exhaustion-T RNA neighbors')
def save(fig,name):
    fig.savefig(FIG/f'{name}.png',dpi=170,bbox_inches='tight');fig.savefig(FIG/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
for item in SAMPLES:
    s=item['sample'];a,im,names=DATA[s]
    fig,axs=plt.subplots(2,3,figsize=(13,9))
    clusters(axs[0,0],a,im,names)
    continuous(axs[0,1],a,im,'PGAM5_log1p_CP10k','PGAM5 RNA: all retained spots')
    continuous(axs[0,2],a,im,'Macrophage_core','Macrophage-core RNA program')
    continuous(axs[1,0],a,im,'T_lineage','T-lineage RNA program')
    continuous(axs[1,1],a,im,'Activation_early','Early activation: T-RNA gate',a.obs.T_RNA_enriched)
    continuous(axs[1,2],a,im,'Exhaustion_core','Exhaustion-core: T-RNA gate',a.obs.T_RNA_enriched)
    fig.suptitle(f'GSE281759 | {s} | {item["patient"]} | {item["tumor_status"]} | n={len(a):,}',fontsize=15)
    fig.legend(handles=[Patch(facecolor=plt.get_cmap('tab10')(k),label=v) for k,v in enumerate(names.values())],loc='lower left',bbox_to_anchor=(.03,.073),ncol=3,frameon=False,fontsize=8)
    fig.text(.05,.015,'Genuine Visium coordinates on original low-resolution H&E. Program means include zeros. Grey: outside T-RNA gate.\nShared scales; values above the displayed maxima are clipped for plotting only. S groups are section-specific mixed RNA regions.',fontsize=8)
    fig.subplots_adjust(top=.91,bottom=.17,hspace=.2,wspace=.22);save(fig,f'{s}_spatial_RNA_maps')
    # Match prior dotplot conventions: size=nonzero fraction; color=mean including zeros.
    genes=['PGAM5','MYO19','C1QA','C1QB','CD68','CSF1R','LST1','TYROBP','CD3D','TRAC','CD69','IL2RA','ICOS','IFNG','GZMB','PDCD1','LAG3','TIGIT','TOX','CXCL13','MS4A1','JCHAIN','PECAM1','COL1A1','ALB','EPCAM','MKI67','SPP1']
    genes=[g for g in genes if g in a.var_names];cats=list(names)
    fig,(ax,dx)=plt.subplots(1,2,figsize=(17,5.6),gridspec_kw={'width_ratios':[1,2.6]})
    xy=a.obsm['X_umap']
    for k,c in enumerate(cats):
        m=a.obs.expression_cluster.eq(c).to_numpy();ax.scatter(*xy[m].T,s=3,c=[plt.get_cmap('tab10')(k)],rasterized=True)
        ax.text(*np.median(xy[m],axis=0),c,weight='bold',bbox=dict(facecolor='white',alpha=.75,edgecolor='none',pad=1))
    ax.set(xlabel='UMAP 1',ylabel='UMAP 2',title='Expression UMAP (not spatial coordinates)')
    xv=a[:,genes].X.toarray();counts=a[:,genes].layers['counts'].toarray()
    means=np.vstack([xv[a.obs.expression_cluster.eq(c)].mean(0) for c in cats]);pct=np.vstack([(counts[a.obs.expression_cluster.eq(c)]>0).mean(0)*100 for c in cats])
    xx,yy=np.meshgrid(np.arange(len(genes)),np.arange(len(cats)))
    sc=dx.scatter(xx.ravel(),yy.ravel(),s=pct.ravel()*1.1,c=means.ravel(),cmap='YlOrRd',vmin=0,vmax=max(1,float(means.max())),linewidths=.25,edgecolors='#999999')
    dx.set_xticks(np.arange(len(genes)),genes,rotation=90);dx.set_yticks(np.arange(len(cats)),[names[c] for c in cats]);dx.set_ylim(len(cats)-.5,-.5);dx.set_xlim(-.6,len(genes)-.4)
    dx.set_title('Dot area: RNA detected % | Color: mean log1p(CP10k), including zeros')
    plt.colorbar(sc,ax=dx,fraction=.025,pad=.02).set_label('Mean log1p(CP10k)')
    for p in [25,50,100]:dx.scatter([],[],s=p*1.1,c='#aaaaaa',label=f'{p}%')
    dx.legend(loc='upper center',bbox_to_anchor=(.5,-.31),ncol=3,frameon=False)
    fig.suptitle(f'{s} | mixed-spot RNA groups and marker expression',fontsize=14);fig.subplots_adjust(top=.85,bottom=.29,wspace=.52);save(fig,f'{s}_UMAP_marker_dotplot')
    print('Figures:',s,flush=True)
for patient in ['A009','B013','C004','C008']:
    items=[i for i in SAMPLES if i['patient']==patient]
    fig,axs=plt.subplots(len(items),3,figsize=(11,3.5*len(items)),squeeze=False)
    for row,it in enumerate(items):
        a,im,names=DATA[it['sample']]
        continuous(axs[row,0],a,im,'PGAM5_log1p_CP10k',it['sample']+' | PGAM5 RNA')
        continuous(axs[row,1],a,im,'Macrophage_core','Macrophage-core RNA')
        proximity(axs[row,2],a,im)
    fig.suptitle(f'{patient} | spatial RNA patterns across tumor sections',fontsize=15)
    fig.text(.03,.012,'Neighbor map: purple = exhaustion-supported T-RNA spot; blue = Mac-RNA spot with such neighbor; red = PGAM5-high subset.\nRNA gates on mixed spots, not cell identities or validated communication. C004T1P: primary; C004T1R: recurrent.',fontsize=8)
    fig.subplots_adjust(top=.96,bottom=.065,hspace=.18,wspace=.2);save(fig,f'{patient}_PGAM5_immune_spatial_overview')
q=pd.read_csv(ROOT/'sample_QC_summary.csv');stats=pd.read_csv(ROOT/'PGAM5_spatial_autocorrelation.csv')
fig,axs=plt.subplots(1,3,figsize=(15,6));ys=np.arange(len(q))
axs[0].barh(ys,q.qc_spots,color='#3589aa',label='Main QC (MT <=25%, UMI >=500)');axs[0].barh(ys,q.qc_mt5_author_spots,color='#dfb151',alpha=.8,label='Author rule (MT <=5%, genes >=200)')
axs[0].set_yticks(ys,q['sample']);axs[0].invert_yaxis();axs[0].set(xlabel='Spots retained',title='QC sensitivity');axs[0].legend(fontsize=7,loc='lower right')
axs[1].barh(ys,q.PGAM5_detection_pct,color='#c95c49');axs[1].set_yticks(ys,q['sample']);axs[1].invert_yaxis();axs[1].set(xlabel='PGAM5 RNA detected (%)',title='All main-QC spots')
d=stats[stats.method.eq('UMI_MT_stratified')].set_index('sample').loc[q['sample']]
axs[2].barh(ys,d.Moran_I,color='#527d66');axs[2].set_yticks(ys,q['sample']);axs[2].invert_yaxis();axs[2].set(xlabel="Moran's I",title='PGAM5 spatial autocorrelation')
for k,r in enumerate(d.itertuples()):axs[2].text(max(0,r.Moran_I)+.005,k,f'q={r.BH_q_32_tests:.3g}',va='center',fontsize=8)
fig.suptitle('GSE281759 | 16 tumor sections, 4 patients',fontsize=15);fig.text(.03,.01,'q: BH across 32 tests, 4,999 permutations; conditional null preserves UMI-quintile / MT-tertile expression assignments. Sections are not independent patients.',fontsize=8)
fig.tight_layout(rect=[0,.045,1,.95]);save(fig,'cohort_QC_PGAM5_spatial_statistics')
lr=pd.read_csv(ROOT/'candidate_LR_codetection_summary.csv');d=lr[lr.source_group.eq('PGAM5_high')&lr.target_state.eq('Exhaustion_supported')]
pairs=pd.read_csv(ROOT/'selected_LR_pairs.csv');order=[r.ligand+'--'+r.receptor_complex for r in pairs.itertuples()]
mat=d.pivot(index='pair',columns='sample',values='edge_codetection_pct').reindex(index=order,columns=q['sample'])
fig,ax=plt.subplots(figsize=(11,8));cmap=plt.get_cmap('YlOrRd').copy();cmap.set_bad('#d2d7dd')
plot=ax.imshow(mat,cmap=cmap,vmin=0,vmax=100,aspect='auto');ax.set_xticks(range(16),q['sample'],rotation=90);ax.set_yticks(range(21),order)
plt.colorbar(plot,ax=ax,label='Eligible spatial edges with ligand + all receptor subunits (%)')
ax.set_title('PGAM5-high Mac-RNA sources -> exhaustion-supported T-RNA neighbors\nCandidate RNA co-detection only; grey = no eligible edges')
fig.text(.04,.012,'Targets and sources are mixed spots. Repeated edges are not independent observations. No communication significance or causality test.\nTIGIT/PDCD1/CTLA4 also belong to the T-state panel: candidate support is not independent validation of exhaustion.',fontsize=8)
fig.tight_layout(rect=[0,.06,1,1]);save(fig,'PGAM5_high_Tex_candidate_LR_codetection')
