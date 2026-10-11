"""Reviewed standalone scientific figures; coordinates never invented."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE245908_immune_mpl')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/GSE245908_immune_cache')
import numpy as np,pandas as pd,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Patch
from matplotlib.collections import PatchCollection,LineCollection
from matplotlib.colors import Normalize
from analyze import ROOT,TARGETS,PANELS,SAMPLES,PAIRS
FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
NAMES={'All_T':'All T-lineage (baseline)','T_exhaustion':'Exhaustion-supported T','T_activation':'Early activation-supported T','CD8_cytotoxic_T':'CD8 / cytotoxic T-related','Naive_memory_T':'CCR7/TCF7 naive-memory T','Treg_related':'FOXP3/IKZF2 Treg-related','NK_related':'KLRD1/GNLY NK-related','B_related':'MS4A1/CD79A B-related','Plasma_related':'MZB1/JCHAIN plasma-related','cDC2_related':'CD1C/FCER1A cDC2-related','cDC1_related':'CLEC9A/XCR1 cDC1-related','pDC_related':'CLEC4C/IL3RA pDC-related','Inflammatory_monocyte':'FCN1/VCAN inflammatory monocyte','Neutrophil_related':'FCGR3B/CSF3R neutrophil-related','Mast_related':'TPSB2/CPA3 mast-related'}
BLUE='#2b7d9e';GOLD='#d89437';PURPLE='#7b61a5';GREY='#cbd1d7'
DATA={s:ad.read_h5ad(ROOT/f'{s}_immune_marker_subset.h5ad') for s in SAMPLES}
counts=pd.read_csv(ROOT/'subgroup_neighborhood_counts.csv');stats=pd.read_csv(ROOT/'spatial_conditional_tests.csv');lr=pd.read_csv(ROOT/'candidate_LR_by_immune_subgroup.csv')
def save(fig,name):
    fig.savefig(FIG/f'{name}.png',dpi=165,bbox_inches='tight');fig.savefig(FIG/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(15,8.5));y=np.arange(len(TARGETS))
for ax,s in zip(axs,SAMPLES):
    for group,off,col in [('PGAM5_high',-.17,GOLD),('Other_detected',.17,BLUE)]:
        d=counts[(counts['sample']==s)&(counts.source_group==group)].set_index('target').reindex(TARGETS)
        ax.barh(y+off,d.source_adjacency_pct,height=.32,color=col,label=group.replace('_',' '))
        if group=='PGAM5_high':
            for j,r in enumerate(d.itertuples()):ax.text(r.source_adjacency_pct+1,j+off,f'{r.sources_with_neighbor}/{r.n_sources}',va='center',fontsize=8)
    ax.set_yticks(y,[NAMES[k] for k in TARGETS]);ax.invert_yaxis();ax.set(xlim=(0,105),xlabel='Sources with >=1 direct target neighbor (%)',title=s+' | descriptive neighborhood coverage');ax.legend(fontsize=9)
fig.suptitle('PGAM5-high macrophage-associated spots and immune RNA subgroup neighbors',fontsize=15)
fig.text(.025,.015,'Labels are multi-gene RNA gates on mixed Visium spots; subgroups may overlap. Counts are source spots, never cell counts.\nAll-T is a baseline; nominal proximity does not imply selective enrichment. No-target rows describe detection limits.',fontsize=9)
fig.tight_layout(rect=[0,.065,1,.95]);save(fig,'immune_subgroup_neighborhood_coverage')

fig,axs=plt.subplots(1,2,figsize=(12,9));categories=['neighbor_fraction','neighbor_present','source_not_target_fraction'];titles=['Fraction\nΔ (pp)','Presence\nΔ (pp)','Fraction\nstrict (pp)']
for ax,s in zip(axs,SAMPLES):
    matrix=np.full((len(TARGETS),3),np.nan);texts=np.empty(matrix.shape,dtype=object)
    for i,k in enumerate(TARGETS):
        for j,kind in enumerate(categories):
            r=stats[(stats['sample']==s)&(stats.target==k)&(stats.kind==kind)].iloc[0]
            if r.status=='tested':matrix[i,j]=100*r.statistic;texts[i,j]=f'{100*r.statistic:+.1f}\nq={r.BH_q_90_planned_tests:.2f}'
            else:texts[i,j]='Not\ntestable'
    cmap=plt.get_cmap('PuOr_r').copy();cmap.set_bad('#dfe3e7');im=ax.imshow(matrix,cmap=cmap,vmin=-30,vmax=30,aspect='auto')
    ax.set_yticks(range(len(TARGETS)),[NAMES[k] for k in TARGETS]);ax.set_xticks(range(3),titles,fontsize=9);ax.set_title(s)
    for (i,j),v in np.ndenumerate(texts):ax.text(j,i,v,ha='center',va='center',fontsize=8,color='white' if np.isfinite(matrix[i,j]) and abs(matrix[i,j])>22 else '#22282e')
fig.suptitle('Spatially conditional subgroup neighborhood tests',fontsize=15)
fig.text(.03,.012,'Delta: High minus other; strict: exclude sources matching that target gate. q: BH across 90 planned tests.\nConditional deviations tested about exact stratum expectation, not zero. Grey: insufficient/constant/degenerate data. Only 2 patients.',fontsize=9)
fig.subplots_adjust(top=.93,bottom=.12,left=.28,right=.88,wspace=1.15)
fig.colorbar(im,cax=fig.add_axes([.925,.20,.02,.60]),label='Raw High-minus-other difference (percentage points)')
save(fig,'immune_subgroup_spatial_statistics')

def backdrop(ax,a):
    s=a.obs['sample'].iloc[0];im=a.uns['spatial'][s]['images']['lowres'];ax.imshow(im,alpha=.7);ax.set_xlim(0,im.shape[1]);ax.set_ylim(im.shape[0],0);ax.axis('off')
def dots(ax,a,mask,col,edge=False):
    sf=a.uns['spatial'][a.obs['sample'].iloc[0]]['scalefactors'];xy=a.obsm['spatial']*sf['tissue_lowres_scalef'];r=sf['spot_diameter_fullres']*sf['tissue_lowres_scalef']/2
    collection=PatchCollection([Circle(xy[i],r) for i in np.flatnonzero(mask)],linewidths=.85 if edge else 0,rasterized=True,alpha=.95)
    collection.set_facecolor('none' if edge else col);collection.set_edgecolor(col if edge else 'none');ax.add_collection(collection)
for s,a in DATA.items():
    fig,axs=plt.subplots(4,4,figsize=(14,14));o=a.obs;high=o.source_group.eq('PGAM5_high').to_numpy();G=a.obsp['spatial_connectivities'].tocoo();sf=a.uns['spatial'][s]['scalefactors'];xy=a.obsm['spatial']*sf['tissue_lowres_scalef']
    for ax,k in zip(axs.ravel(),TARGETS):
        backdrop(ax,a);dots(ax,a,np.ones(len(a),bool),GREY);target=o['gate_'+k].to_numpy(bool);dots(ax,a,target,BLUE)
        keep=high[G.row]&target[G.col];si,ti=G.row[keep],G.col[keep]
        if len(si):ax.add_collection(LineCollection(np.stack([xy[si],xy[ti]],axis=1),colors=PURPLE,linewidths=.8,alpha=.9,rasterized=True))
        dots(ax,a,high,GOLD,edge=True)
        ax.set_title(f'{NAMES[k]}\n{len(np.unique(si))} High sources; {len(si)} edges',fontsize=9)
    axs.ravel()[-1].axis('off')
    axs.ravel()[-1].text(.03,.7,'Original H&E coordinates\nBlue: target RNA gate\nGold rings: PGAM5-high Mac-RNA\nPurple lines: direct spot neighbors\n\nRNA gates can overlap.\nZero edges != absent cell type.\nNo cell-specific communication\nis established.',transform=axs.ravel()[-1].transAxes,fontsize=11,va='top')
    fig.suptitle(f'GSE245908 | {s} | PGAM5-high Mac-RNA / immune subgroup neighborhoods',fontsize=15)
    fig.text(.03,.012,'Main QC and original PGAM5 Q75 labels are unchanged. Non-self hex-grid neighbors only; no spatial location inferred from UMAP.\nTarget programs are provisional mixed-spot RNA annotations; not pathology ROIs, pure single cells, or interaction probabilities.',fontsize=9)
    fig.subplots_adjust(top=.93,bottom=.06,hspace=.23,wspace=.1);save(fig,s+'_immune_subgroup_spatial_neighborhoods')
    # Multi-gene marker evidence: mean INCLUDING zeros, dot area = nonzero fraction.
    marker=['CD3D','TRAC','CD8A','CD8B','PDCD1','LAG3','TOX','TIGIT','CD69','ICOS','CCR7','TCF7','FOXP3','IKZF2','KLRD1','GNLY','MS4A1','CD79A','MZB1','JCHAIN','CD1C','FCER1A','CLEC9A','XCR1','CLEC4C','IL3RA','FCN1','VCAN','FCGR3B','CSF3R','TPSB2','CPA3']
    x=a[:,marker].X.toarray();c=a[:,marker].layers['counts'].toarray();means=[];pct=[]
    for k in TARGETS:
        m=o['gate_'+k].to_numpy(bool);means.append(x[m].mean(0) if m.any() else np.full(len(marker),np.nan));pct.append((c[m]>0).mean(0)*100 if m.any() else np.zeros(len(marker)))
    means=np.asarray(means);pct=np.asarray(pct);xx,yy=np.meshgrid(np.arange(len(marker)),np.arange(len(TARGETS)))
    fig,ax=plt.subplots(figsize=(16,8));scat=ax.scatter(xx.ravel(),yy.ravel(),s=pct.ravel(),c=means.ravel(),cmap='YlOrRd',vmin=0,vmax=max(1,np.nanmax(means)),edgecolors='#8b9299',linewidths=.3)
    ax.set_xticks(range(len(marker)),marker,rotation=90);ax.set_yticks(range(len(TARGETS)),[NAMES[k]+' (n='+str(int(o['gate_'+k].sum()))+')' for k in TARGETS]);ax.set_xlim(-.7,len(marker)-.3);ax.set_ylim(len(TARGETS)-.5,-.5);ax.set_title(s+' | immune RNA gate marker evidence')
    plt.colorbar(scat,ax=ax,label='Mean log1p(CP10k), including zeros',fraction=.025,pad=.02)
    for n in [25,50,100]:ax.scatter([],[],s=n,c='#999999',label=f'{n}%')
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.16),ncol=3,frameon=False,title='RNA detected %')
    fig.text(.025,.012,'Marker gene means computed from full-library normalization. Gates use these markers: this dotplot documents the rule, not independent validation.\nEmpty pDC gate means no supported spots under the fixed rule, not proof of absent pDCs. CD8A and immune programs can mix across cells.',fontsize=9)
    fig.tight_layout(rect=[0,.11,1,1]);save(fig,s+'_immune_subgroup_marker_dotplot')

fixed=[l+'--'+r for l,r in PAIRS];fig,axs=plt.subplots(1,2,figsize=(17,9.5))
for ax,s in zip(axs,SAMPLES):
    d=lr[(lr['sample']==s)&(lr.source_group=='PGAM5_high')];matrix=d.pivot(index='target',columns='pair',values='codetection_pct').reindex(index=TARGETS,columns=fixed)
    cm=plt.get_cmap('YlOrRd').copy();cm.set_bad('#d4dae0');im=ax.imshow(matrix,cmap=cm,vmin=0,vmax=100,aspect='auto');ax.set_yticks(range(len(TARGETS)),[NAMES[k] for k in TARGETS]);ax.set_xticks(range(len(fixed)),fixed,rotation=90,fontsize=8);ax.set_title(s+' | High sources to immune subgroup targets')
fig.suptitle('Candidate ligand-receptor RNA co-detection (descriptive)',fontsize=15)
fig.text(.025,.012,'15 fixed pairs from LIANA snapshot 1.10.0; no communication significance test. Grey: no eligible edges; yellow zero: edges without RNA support.\nAll-T includes T subgroups; subgroup gates and targets can overlap. A 100% value from 1-2 edges is not replicated evidence. Mixed-spot receptor attribution remains unresolved.',fontsize=9)
fig.subplots_adjust(left=.24,right=.89,top=.92,bottom=.32,wspace=1.35)
fig.colorbar(im,cax=fig.add_axes([.935,.32,.016,.60]),label='Eligible edges with ligand + all receptor subunits (%)')
save(fig,'immune_subgroup_candidate_LR_codetection')
print('Rendered',len(list(FIG.glob('*.png'))),'PNG/PDF pairs',flush=True)
