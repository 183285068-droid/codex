"""Original H&E + genuine Visium coordinates, and fully source-backed static figures."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE281759_SPP1_mpl')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/GSE281759_SPP1_cache')
from pathlib import Path
import json
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Circle,Patch
from matplotlib.collections import PatchCollection,LineCollection
from analyze import ROOT,REF,SAMPLES,augment,PROGRAMS,LR_NAMES

FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
BLUE='#287ea3';ORANGE='#db963c';PURPLE='#765ca2';OLIVE='#6f8744';GREY='#cbd0d5'
DATA={i['sample']:augment(ad.read_h5ad(REF/f'{i["sample"]}_spatial_marker_subset.h5ad')) for i in SAMPLES}
V={k:max(float(a.obs[k].quantile(.99)) for a in DATA.values()) for k in ['SPP1_log1p_CP10k','T_lineage','Exhaustion_core']}
def save(fig,name):
    fig.savefig(FIG/f'{name}.png',dpi=160,bbox_inches='tight');fig.savefig(FIG/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
def image(ax,a):
    im=a.uns['spatial'][a.obs['sample'].iloc[0]]['images']['lowres'];ax.imshow(im,alpha=.65);ax.set_xlim(0,im.shape[1]);ax.set_ylim(im.shape[0],0);ax.axis('off')
def dots(ax,a,mask,col=None,values=None,var=None,edge=False):
    sf=a.uns['scalefactors'];xy=a.obsm['spatial']*sf['tissue_lowres_scalef'];r=sf['spot_diameter_fullres']*sf['tissue_lowres_scalef']/2
    ix=np.flatnonzero(mask);pc=PatchCollection([Circle(xy[i],r) for i in ix],rasterized=True,linewidths=.7 if edge else 0,alpha=.9)
    if values is None:pc.set_facecolor('none' if edge else col);pc.set_edgecolor(col if edge else 'none')
    else:pc.set_array(np.asarray(values)[ix]);pc.set_cmap('YlOrRd');pc.set_norm(Normalize(0,V[var]))
    ax.add_collection(pc);return pc
def continuous(ax,a,k,mask,title):
    image(ax,a);dots(ax,a,np.ones(len(a),bool),GREY)
    pc=dots(ax,a,mask,values=a.obs[k].to_numpy(),var=k)
    ax.set_title(title,fontsize=10);plt.colorbar(pc,ax=ax,fraction=.035,pad=.015,extend='max').set_label('log1p(CP10k)',fontsize=8)
def states(ax,a):
    image(ax,a);o=a.obs;t=o.T_RNA_enriched.to_numpy(bool);ex=o.Exhaustion_core_supported.to_numpy(bool);act=o.Activation_early_supported.to_numpy(bool)
    dots(ax,a,np.ones(len(a),bool),GREY);dots(ax,a,t,BLUE);dots(ax,a,ex,PURPLE);dots(ax,a,act,OLIVE);dots(ax,a,ex&act,'#222222')
    ax.set_title('T-RNA states (can overlap)',fontsize=10)
def neighbors(ax,a):
    image(ax,a);o=a.obs;m=o.Mac_RNA_enriched.to_numpy(bool)&o.SPP1_raw.gt(0).to_numpy();t=o.T_RNA_enriched.to_numpy(bool)
    dots(ax,a,np.ones(len(a),bool),GREY,);dots(ax,a,t,BLUE);dots(ax,a,o.Exhaustion_core_supported.to_numpy(bool),PURPLE);dots(ax,a,o.Activation_early_supported.to_numpy(bool),OLIVE)
    g=a.obsp['spatial_connectivities'].tocoo();keep=m[g.row]&t[g.col];s,j=g.row[keep],g.col[keep]
    xy=a.obsm['spatial']*a.uns['scalefactors']['tissue_lowres_scalef'];cols=np.repeat('#7e8892',len(s)).astype(object)
    cols[o.Exhaustion_core_supported.to_numpy(bool)[j]]=PURPLE;cols[o.Activation_early_supported.to_numpy(bool)[j]]=OLIVE
    ax.add_collection(LineCollection(np.stack([xy[s],xy[j]],axis=1),colors=cols,linewidths=.6,alpha=.7,rasterized=True))
    dots(ax,a,m,ORANGE,edge=True)
    ax.set_title(f'SPP1+ Mac-RNA to T neighbors\n{len(s)} edges; {len(np.unique(s))} sources',fontsize=10)

for pat in ['A009','B013','C004','C008']:
    items=[i for i in SAMPLES if i['patient']==pat]
    fig,axs=plt.subplots(len(items),4,figsize=(15,3*len(items)),squeeze=False)
    for row,it in enumerate(items):
        a=DATA[it['sample']];o=a.obs
        continuous(axs[row,0],a,'SPP1_log1p_CP10k',np.ones(len(a),bool),it['sample']+' | SPP1: all spots')
        continuous(axs[row,1],a,'SPP1_log1p_CP10k',o.Mac_RNA_enriched.to_numpy(bool),'SPP1: Mac-RNA gate')
        states(axs[row,2],a);neighbors(axs[row,3],a)
    fig.suptitle(f'GSE281759 | {pat} | SPP1 macrophage-associated / T-RNA spatial neighborhoods',fontsize=15)
    fig.legend(handles=[Patch(facecolor=BLUE,label='T-RNA'),Patch(facecolor=PURPLE,label='Exhaustion-supported T-RNA'),Patch(facecolor=OLIVE,label='Activation-supported T-RNA'),Patch(facecolor='#222222',label='Both states'),Patch(facecolor='none',edgecolor=ORANGE,label='SPP1-detected Mac-RNA')],loc='lower center',bbox_to_anchor=(.5,.026),ncol=3,frameon=False,fontsize=9)
    fig.text(.035,.007,'Original H&E / true Visium coordinates; shared expression scales with display-only Q99 clipping. Lines: direct non-self spot neighbors.\nMixed spots, not single cells or proven communication. C004T1P: primary; C004T1R: recurrent.',fontsize=9)
    fig.subplots_adjust(top=.92 if len(items)>2 else .84,bottom=.09 if len(items)>2 else .19,hspace=.14,wspace=.15)
    save(fig,pat+'_SPP1_mac_T_spatial_neighborhoods')

summary=pd.read_csv(ROOT/'section_summary.csv');stats=pd.read_csv(ROOT/'interaction_statistics.csv');order=summary['sample'].tolist();y=np.arange(len(order))
fig,axs=plt.subplots(1,3,figsize=(16,8.5))
axs[0].barh(y,100*summary.SPP1_mac_detected/summary.mac_RNA_spots,color=ORANGE);axs[0].set(xlabel='SPP1 detected / Mac-RNA spots (%)',title='SPP1 in macrophage-associated spots',xlim=(0,110))
for k,r in enumerate(summary.itertuples()):axs[0].text(100*r.SPP1_mac_detected/r.mac_RNA_spots+1,k,f'{r.SPP1_mac_detected}/{r.mac_RNA_spots}',va='center',fontsize=8)
rates=pd.read_csv(ROOT/'source_group_neighbor_summary.csv')
for ax,state in zip(axs[1:],['All_T','Exhaustion']):
    for group,delta,col in [('SPP1_detected',-.17,ORANGE),('Undetected',.17,BLUE)]:
        dd=rates[(rates.source_group==group)&(rates.target_state==state)].set_index('sample').reindex(order)
        ax.barh(y+delta,dd.neighbor_pct,height=.32,color=col,label=group.replace('_',' '))
    ax.set(xlabel='Source spots with a direct neighbor (%)',title=state.replace('_',' ')+' RNA adjacency',xlim=(0,105));ax.legend(fontsize=9)
for ax in axs:ax.set_yticks(y,order);ax.invert_yaxis();ax.grid(axis='x',alpha=.15)
fig.suptitle('GSE281759 | SPP1 detection and direct T-RNA neighborhood proportions',fontsize=15)
fig.text(.035,.015,'Adjacency excludes isolated sources; different group denominators; a missing comparator is blank, never zero. See permutation tests.\n16 sections belong to 4 patients; C004T1R is recurrent. Raw count >0 defines SPP1 detection.',fontsize=9)
fig.tight_layout(rect=[0,.07,1,.95]);save(fig,'SPP1_detection_T_neighborhood_overview')

def heatmap(family,tests,name,title,label,limit):
    matrix=np.full((len(order),len(tests)),np.nan);labels=np.empty(matrix.shape,dtype=object)
    for i,s in enumerate(order):
        for j,(endpoint,cohort,short) in enumerate(tests):
            r=stats[(stats['sample']==s)&(stats.family==family)&(stats.endpoint==endpoint)&(stats.cohort==cohort)].iloc[0]
            if r.status=='tested':matrix[i,j]=r.statistic;labels[i,j]=f'{r.statistic:+.1f}\nq={r.BH_q_176_planned_tests:.2f}' if family=='adjacency' else f'{r.statistic:+.2f}\nq={r.BH_q_176_planned_tests:.2f}'
            else:labels[i,j]='Not\ntestable'
    fig,ax=plt.subplots(figsize=(max(11,len(tests)*2),10));cm=plt.get_cmap('PuOr_r').copy();cm.set_bad('#e0e3e7');im=ax.imshow(matrix,cmap=cm,vmin=-limit,vmax=limit,aspect='auto')
    ax.set_xticks(range(len(tests)),[t[2] for t in tests]);ax.set_yticks(range(len(order)),order)
    for (i,j),v in np.ndenumerate(labels):ax.text(j,i,v,ha='center',va='center',fontsize=8,color='white' if np.isfinite(matrix[i,j]) and abs(matrix[i,j])>limit*.68 else '#20262b')
    plt.colorbar(im,ax=ax,label=label,fraction=.025,pad=.02);ax.set_title(title,fontsize=14,pad=16)
    fig.text(.035,.012,'q: BH across 176 planned tests; 4,999 within-tile/UMI/Mac-core permutations. Grey: missing comparison or insufficient/constant data.\nSpot-level conditional association; sections are not independent patients. Positive: greater T adjacency / stronger neighbor program.',fontsize=9)
    fig.tight_layout(rect=[0,.06,1,1]);save(fig,name)
adjtests=[('All_T','detected_vs_undetected','All T\nSPP1+ vs zero'),('Exhaustion','detected_vs_undetected','Exhaustion\nSPP1+ vs zero'),('Activation','detected_vs_undetected','Activation\nSPP1+ vs zero'),('Exhaustion','nonT_detected_vs_undetected','Exhaustion\nexclude T+ sources'),('Exhaustion','high_vs_other_detected','Exhaustion\nhigh vs other detected')]
heatmap('adjacency',adjtests,'SPP1_T_adjacency_statistics','SPP1 macrophage-associated / T-RNA adjacency: conditional permutation tests','Difference in source-neighbor proportion (percentage points)',40)
corrtests=[(p,'all_mac_including_zero',n) for p,n in zip(PROGRAMS,['Exhaustion\ncore','Exhaustion\nexpanded','Early\nactivation','Cytotoxic\neffector','Exhaustion without\nCTLA4/TIGIT'])]+[('Exhaustion_core','nonT_mac_including_zero','Exhaustion core\nexclude T+ sources')]
heatmap('correlation',corrtests,'SPP1_neighbor_T_program_association','SPP1 RNA and direct T-neighbor programs: partial Spearman','Partial Spearman rho (including SPP1 zeros)',.6)

pt=pd.read_csv(ROOT/'patient_descriptive_effect_summary.csv');fig,axs=plt.subplots(1,4,figsize=(16,4.8))
spec=[('adjacency','All_T','detected_vs_undetected','All-T adjacency\nSPP1+ minus zero (pp)'),('adjacency','Exhaustion','detected_vs_undetected','Exhaustion adjacency\nSPP1+ minus zero (pp)'),('adjacency','Activation','detected_vs_undetected','Activation adjacency\nSPP1+ minus zero (pp)'),('correlation','Exhaustion_core','all_mac_including_zero','Exhaustion program\npartial Spearman rho')]
for ax,(fam,ep,coh,label) in zip(axs,spec):
    for scope,offset,marker,col in [('all_sections',-.10,'o',BLUE),('primary_only',.10,'s',ORANGE)]:
        dd=pt[(pt.scope==scope)&(pt.family==fam)&(pt.endpoint==ep)&(pt.cohort==coh)].set_index('patient').reindex(['A009','B013','C004','C008'])
        ax.scatter(np.arange(4)+offset,dd.median_tested_section_statistic,c=col,marker=marker,s=55,label=scope.replace('_',' '))
    ax.axhline(0,color='#7e8892',lw=.8);ax.set_xticks(range(4),['A009','B013','C004','C008']);ax.set_title(label,fontsize=11);ax.set_ylabel('Median of testable section effects')
axs[0].legend(fontsize=9);fig.suptitle('Patient summaries: descriptive medians, without pooling spots as patients',fontsize=14)
fig.text(.035,.02,'4 independent patients. Each section contributes equally; untestable sections omitted, never set to zero. Primary-only excludes C004T1R.\nNo patient-level significance or confidence interval is claimed.',fontsize=9);fig.tight_layout(rect=[0,.14,1,.92]);save(fig,'patient_descriptive_effect_summary')

lr=pd.read_csv(ROOT/'candidate_LR_codetection_summary.csv');fig,axs=plt.subplots(2,1,figsize=(13,10.5))
for ax,state in zip(axs,['All_T','Exhaustion']):
    dd=lr[(lr.source_group=='SPP1_detected')&(lr.target_state==state)]
    mat=dd.pivot(index='pair',columns='sample',values='edge_codetection_pct').reindex(index=LR_NAMES,columns=order)
    cm=plt.get_cmap('YlOrRd').copy();cm.set_bad('#d8dde2');im=ax.imshow(mat,cmap=cm,vmin=0,vmax=100,aspect='auto')
    ax.set_xticks(range(len(order)),order,rotation=90);ax.set_yticks(range(len(LR_NAMES)),LR_NAMES);ax.set_title('SPP1-detected Mac-RNA -> '+state.replace('_',' ')+' neighbors')
    plt.colorbar(im,ax=ax,label='Eligible edges with ligand + all receptor subunits (%)',fraction=.025,pad=.02)
fig.suptitle('Candidate ligand-receptor RNA co-detection on direct spatial edges',fontsize=14)
fig.text(.03,.012,'Fixed 8 pairs from LIANA resource snapshot 1.10.0; descriptive co-detection, no CellChat/LIANA communication significance.\nCD44 is broad-lineage; SPP1 positivity is selected by definition. Mixed spots and repeated targets prevent cell-specific attribution. Grey: no eligible edges.',fontsize=9)
fig.tight_layout(rect=[0,.07,1,.95]);save(fig,'SPP1_mac_T_candidate_LR_codetection')
print('Rendered',len(list(FIG.glob('*.png'))),'PNG/PDF pairs',flush=True)
