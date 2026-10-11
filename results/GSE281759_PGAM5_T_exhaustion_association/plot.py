"""Association effects, adjusted scatterplots and genuine-coordinate spatial maps."""
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Circle,Patch
from matplotlib.collections import PatchCollection
from scipy.stats import rankdata
from analyze import ROOT,REF,SAMPLES,COV,residualizer
FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
table=pd.read_csv(ROOT/'association_statistics.csv',float_precision='round_trip')
primary=table[(table.cohort=='all_mac_including_zero')&(table.program=='Exhaustion_core')].set_index('sample')
NAMES=[x['sample'] for x in SAMPLES]
GROUPS=['PGAM5_high','Other_detected','Undetected'];COLORS=['#bd3445','#df9b43','#758492']
def save(fig,name):
    fig.savefig(FIG/f'{name}.png',dpi=180,bbox_inches='tight');fig.savefig(FIG/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
fig,ax=plt.subplots(figsize=(11,7));ys=np.arange(16);d=primary.loc[NAMES]
for k,r in enumerate(d.itertuples()):
    if r.status=='tested':
        ax.barh(k,r.statistic,color='#bd3445' if r.statistic>0 else '#337ea0',height=.65)
        ax.text(.34,k,f'n={r.n_sources}; p={r.permutation_p:.3g}; q={r.primary_BH_q_16_sections:.3g}',fontsize=9,va='center')
    else:ax.text(.02,k,f'Not tested: n={r.n_sources} (<30)',va='center',color='#777777')
ax.axvline(0,color='#999999',lw=.7);ax.set_yticks(ys,NAMES);ax.invert_yaxis();ax.set_xlim(-.18,.82)
ax.set_xlabel('Adjusted partial Spearman rho (PGAM5 zeros included)');ax.set_title('PGAM5 in Mac-RNA source spots vs neighboring T-RNA exhaustion\n9 covariates; 4,999 local-space / UMI conditional permutations',fontsize=12)
fig.text(.02,.012,'q: BH across 16 planned primary section tests. Red: positive; blue: negative association. Sections are nested within 4 patients.\nMixed Visium spots: neither macrophage-specific PGAM5 expression nor a causal effect is identified.',fontsize=8)
fig.tight_layout(rect=[0,.055,1,1]);save(fig,'primary_PGAM5_neighbor_exhaustion_association')
progs=['Exhaustion_core','Exhaustion_expanded','Activation_early','Effector_cytotoxic']
a=table[(table.cohort=='all_mac_including_zero')&(table.contrast=='partial_Spearman')]
mat=a.pivot(index='sample',columns='program',values='statistic').reindex(index=NAMES,columns=progs)
fig,ax=plt.subplots(figsize=(10,9));cm=plt.get_cmap('RdBu_r').copy();cm.set_bad('#d6d9de');lim=max(.3,float(np.nanmax(np.abs(mat))))
p=ax.imshow(mat,cmap=cm,vmin=-lim,vmax=lim,aspect='auto');ax.set_yticks(range(16),NAMES);ax.set_xticks(range(4),['Exhaustion core','Exhaustion expanded','Early activation','Cytotoxic effector'],rotation=20,ha='right')
for i,s in enumerate(NAMES):
    for j,pn in enumerate(progs):
        r=a[(a['sample']==s)&(a.program==pn)].iloc[0]
        ax.text(j,i,f'rho={r.statistic:.2f}\nq={r.BH_q_96_planned_tests:.2g}' if r.status=='tested' else 'Not tested',ha='center',va='center',fontsize=8)
plt.colorbar(p,ax=ax,label='Adjusted partial Spearman rho');ax.set_title('Neighbor T-RNA program associations with source PGAM5 RNA\nq: BH across all 96 planned primary / secondary / sensitivity tests')
fig.text(.04,.01,'All Mac-RNA sources with >=1 T-RNA neighbor; zeros retained. Grey: insufficient sources. No patient-level pooled p value or causal inference.',fontsize=8)
fig.tight_layout(rect=[0,.055,1,1]);save(fig,'PGAM5_neighbor_T_program_association_heatmap')
fig,axs=plt.subplots(4,4,figsize=(15,12));rng=np.random.default_rng(281759)
fig2,axs2=plt.subplots(4,4,figsize=(15,12))
for k,it in enumerate(SAMPLES):
    s=it['sample'];df=pd.read_csv(ROOT/f'{s}_source_neighbor_RNA.csv',float_precision='round_trip');df=df[df.T_neighbor_count>0]
    ax=axs.ravel()[k];ax2=axs2.ravel()[k]
    r=primary.loc[s]
    if r.status=='tested':
        q,rank=residualizer(df);xr=rankdata(df.PGAM5_log1p_CP10k);yr=rankdata(df.neighbor_Exhaustion_core)
        xr-=q@(q.T@xr);yr-=q@(q.T@yr);xr/=np.std(xr,ddof=1);yr/=np.std(yr,ddof=1)
        for g,c in zip(GROUPS,COLORS):
            choose=df.PGAM5_mac_group.eq(g);ax.scatter(xr[choose],yr[choose],c=c,s=7,alpha=.65,rasterized=True)
        xx=np.linspace(xr.min(),xr.max(),30);ax.plot(xx,r.statistic*xx,c='#555555',lw=1)
        ax.set_title(f'{s}: rho={r.statistic:.2f}, q={r.primary_BH_q_16_sections:.2g}',fontsize=9)
    else:ax.text(.5,.5,f'{s}\nn={len(df)}: not tested',ha='center',va='center',transform=ax.transAxes)
    ax.set_xlabel('PGAM5 adjusted rank / SD',fontsize=8);ax.set_ylabel('Neighbor exhaustion adjusted rank / SD',fontsize=8)
    data=[df.loc[df.PGAM5_mac_group.eq(g),'neighbor_Exhaustion_core'].to_numpy() for g in GROUPS]
    box=ax2.boxplot(data,positions=[0,1,2],widths=.5,patch_artist=True,showfliers=False)
    for patch,c in zip(box['boxes'],COLORS):patch.set_facecolor(c);patch.set_alpha(.35)
    for j,(v,c) in enumerate(zip(data,COLORS)):ax2.scatter(j+rng.uniform(-.17,.17,len(v)),v,c=c,s=4,alpha=.45,rasterized=True)
    rg=table[(table['sample']==s)&(table.contrast=='adjusted_high_minus_other')].iloc[0]
    qstr=f'adjusted H vs D q={rg.BH_q_96_planned_tests:.2g}' if rg.status=='tested' else 'H vs D not tested (small n)'
    ax2.set_title(s+'\n'+qstr,fontsize=9);ax2.set_xticks([0,1,2],[f'H\nn={len(data[0])}',f'D\nn={len(data[1])}',f'0\nn={len(data[2])}']);ax2.set_ylabel('Mean neighbor exhaustion\nlog1p(CP10k)',fontsize=8)
fig.suptitle('Adjusted source-PGAM5 / neighboring T-RNA exhaustion associations',fontsize=14)
fig.text(.025,.012,'Red = PGAM5-high; orange = other detected; grey = undetected. Residual ranks remove 9 quality / RNA-composition covariates.\nLines summarize association; no causal direction, confidence interval, or independent-cell identity is implied.',fontsize=8)
fig.tight_layout(rect=[0,.07,1,.95]);save(fig,'adjusted_PGAM5_exhaustion_scatter_by_section')
fig2.suptitle('Neighbor T-RNA exhaustion scores by source PGAM5 group',fontsize=14)
fig2.text(.025,.012,'H: high among PGAM5-detected Mac-RNA spots; D: other detected; 0: undetected. Raw scores include zeros.\nAnnotated q applies to the adjusted residual-rank H vs D test, not to unadjusted boxes or the 0 group. Shared T targets create dependence handled approximately by local permutation.',fontsize=8)
fig2.tight_layout(rect=[0,.07,1,.95]);save(fig2,'PGAM5_groups_neighbor_exhaustion_boxplots')
# Two-map rows use original tissue and full-resolution pixel coordinates, never UMAP.
def draw(ax,a,im,mask,values,vmax):
    ax.imshow(im,alpha=.7);xy=a.obsm['spatial']*a.uns['scalefactors']['tissue_lowres_scalef'];rad=a.uns['scalefactors']['spot_diameter_fullres']*a.uns['scalefactors']['tissue_lowres_scalef']/2
    background=PatchCollection([Circle(p,rad) for p in xy],facecolors='#bbc2c9',linewidths=0,alpha=.55,rasterized=True);ax.add_collection(background)
    ids=np.flatnonzero(np.asarray(mask));pc=PatchCollection([Circle(xy[i],rad) for i in ids],linewidths=0,alpha=.95,rasterized=True)
    pc.set_array(np.asarray(values)[ids]);pc.set_cmap('YlOrRd');pc.set_norm(Normalize(0,vmax));ax.add_collection(pc)
    ax.set_xlim(0,im.shape[1]);ax.set_ylim(im.shape[0],0);ax.axis('off');plt.colorbar(pc,ax=ax,fraction=.035,pad=.02)
objects={it['sample']:ad.read_h5ad(REF/f'{it["sample"]}_spatial_marker_subset.h5ad') for it in SAMPLES}
pgmax=max(float(a.obs.PGAM5_log1p_CP10k.quantile(.99)) for a in objects.values());exmax=max(float(a.obs.Exhaustion_core.quantile(.99)) for a in objects.values())
for pat in ['A009','B013','C004','C008']:
    its=[it for it in SAMPLES if it['patient']==pat];fig,axs=plt.subplots(len(its),2,figsize=(9,3.8*len(its)),squeeze=False)
    for row,it in enumerate(its):
        a=objects[it['sample']];im=a.uns['spatial'][it['sample']]['images']['lowres'];o=a.obs
        draw(axs[row,0],a,im,o.Mac_RNA_enriched&o.T_neighbor_count.gt(0),o.PGAM5_log1p_CP10k,pgmax)
        draw(axs[row,1],a,im,o.T_RNA_enriched,o.Exhaustion_core,exmax)
        # Ring sources high in PGAM5 and directly adjacent to exhaustion-supported T-RNA.
        xy=a.obsm['spatial']*a.uns['scalefactors']['tissue_lowres_scalef'];rad=a.uns['scalefactors']['spot_diameter_fullres']*a.uns['scalefactors']['tissue_lowres_scalef']/2
        high=o.PGAM5_mac_group.eq('PGAM5_high')&o.T_exhaustion_neighbor_count.gt(0)
        rings=PatchCollection([Circle(p,rad*1.15) for p in xy[high]],facecolors='none',edgecolors='#24b4c3',linewidths=.6,rasterized=True);axs[row,0].add_collection(rings)
        axs[row,0].set_title(it['sample']+' | PGAM5 in Mac-RNA sources with T neighbor',fontsize=9)
        axs[row,1].set_title('T-RNA spots | exhaustion-core RNA',fontsize=10)
    fig.suptitle(pat+' | PGAM5 / T-exhaustion spatial RNA neighborhoods',fontsize=14)
    fig.text(.035,.012,'True Visium positions on original H&E. Grey = outside displayed RNA gate. Cyan rings = high-PGAM5 Mac-RNA source next to exhaustion-supported T-RNA.\nColor: log1p(CP10k) for PGAM5; mean log1p(CP10k) for exhaustion. Shared scales; display clipping at cohort-wide maximum section Q99.\nCo-localization alone is not evidence that macrophage PGAM5 promotes T-cell exhaustion.',fontsize=8)
    fig.subplots_adjust(top=.955,bottom=.075,hspace=.25,wspace=.18);save(fig,pat+'_PGAM5_Tex_spatial_neighborhoods')
print('Generated 8 PNG/PDF figure pairs',flush=True)
