"""Standalone research figures from the saved PGAM5-free graph and real cells."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

R=Path(__file__).resolve().parent

def main():
    obs=pd.read_csv(R/'cell_annotations.csv.gz',index_col=0)
    pc=np.load(R/'primary_PCA_model.npz')['PC']
    graph=sparse.load_npz(R/'primary_neighbor_graph.npz')
    model=ad.AnnData(np.zeros((len(obs),1),dtype=np.float32),obs=obs.copy())
    model.obsm['X_pca']=pc;model.obsp['connectivities']=graph
    model.uns['neighbors']={'connectivities_key':'connectivities','distances_key':'distances',
                            'params':{'n_neighbors':20,'method':'umap','metric':'euclidean','n_pcs':30,'use_rep':'X_pca'}}
    sc.tl.umap(model,min_dist=.35,random_state=0,maxiter=400)
    xy=model.obsm['X_umap'];pd.DataFrame({'cell_id':obs.cell_id,'UMAP1':xy[:,0],'UMAP2':xy[:,1]}).to_csv(R/'UMAP_coordinates.csv.gz',index=False)
    association=pd.read_csv(R/'cluster_PGAM5_association.csv',dtype={'cluster':str})
    association=association[association.geometry.eq('primary')].sort_values('cluster',key=lambda z:z.astype(int))
    labels=obs.primary_cluster.astype(str).to_numpy()
    candidate=str(association[(association.pooled_detection_rate_ratio>1)&(association.evaluable_samples>=3)&
                             (association.fraction_evaluable_samples_detection_enriched>=.7)].sort_values('pooled_detection_rate_ratio',ascending=False).iloc[0].cluster)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
    cmap=plt.get_cmap('tab20');colors={str(i):cmap(i) for i in range(15)}
    colors[candidate]='#ad4e22'
    fig,axes=plt.subplots(2,2,figsize=(13.6,10.2),gridspec_kw={'width_ratios':[1.1,1]})
    ax=axes[0,0]
    for cluster in association.cluster:
        sel=labels==cluster
        ax.scatter(xy[sel,0],xy[sel,1],s=3,color=colors[cluster],rasterized=True,alpha=.65)
        center=np.median(xy[sel],axis=0)
        ax.text(center[0],center[1],cluster,ha='center',va='center',fontsize=9,
                bbox={'facecolor':'white','alpha':.8,'edgecolor':'none','pad':1})
    ax.set(title='A  Macrophage clusters (PGAM5 excluded from features)',xlabel='UMAP1',ylabel='UMAP2')
    ax=axes[0,1];pos=obs.PGAM5_detected.to_numpy(bool)
    ax.scatter(xy[~pos,0],xy[~pos,1],s=3,color='#ccd1d5',rasterized=True,alpha=.4,label=f'Undetected (n={int((~pos).sum())})')
    ax.scatter(xy[pos,0],xy[pos,1],s=7,color='#b9552c',rasterized=True,alpha=.8,label=f'Detected (n={int(pos.sum())})')
    ax.set(title='B  PGAM5 RNA detection spans multiple clusters',xlabel='UMAP1',ylabel='UMAP2')
    ax.legend(loc='lower left',frameon=False,fontsize=9)
    ax=axes[1,0];y=np.arange(len(association))
    ax.barh(y,association.PGAM5_detection_fraction*100,
            color=[colors[c] if c==candidate else '#7c98b1' for c in association.cluster],height=.7)
    ax.axvline(pos.mean()*100,color='#444444',ls='--',lw=1,label=f'All macrophages: {pos.mean():.2%}')
    ax.set(yticks=y,yticklabels=['C'+c+f' (n={n})' for c,n in zip(association.cluster,association.cells)],
           xlabel='PGAM5 RNA-detected cells within cluster (%)',title='C  Within-cluster detection rate')
    ax.invert_yaxis();ax.legend(loc='upper right',frameon=False,fontsize=8.5)
    ax=axes[1,1]
    ax.barh(y,association.capture_of_all_PGAM5_detected*100,
            color=[colors[c] if c==candidate else '#7c98b1' for c in association.cluster],height=.7)
    ax.set(yticks=y,yticklabels=['C'+c for c in association.cluster],xlabel='Share of all 868 PGAM5 RNA-detected cells (%)',
           title='D  No single cluster contains most detected cells')
    ax.invert_yaxis()
    for ax in axes.ravel():ax.spines[['top','right']].set_visible(False)
    fig.suptitle('GSE202642 HCC: PGAM5 RNA detection and macrophage clustering\n7 patients (user-confirmed); 12,135 macrophages; cycle genes/cells retained',fontsize=12,y=.982)
    fig.text(.5,.024,'UMAP is descriptive. Cluster evidence also uses non-PGAM5 markers, patient consistency and resampling.\n'
             'Highlighted C'+candidate+': PGAM5-enriched cycling macrophage candidate; it is not the complete PGAM5-detected category.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.885,bottom=.12,hspace=.36,wspace=.38,left=.11,right=.97)
    for ext in ['png','pdf']:fig.savefig(R/f'PGAM5_subpopulation_overview.{ext}',dpi=200)
    plt.close(fig)
    # Patient-level evidence for the selected candidate, keeping sparse patients visible.
    table=pd.read_csv(R/'cluster_PGAM5_by_sample.csv',dtype={'cluster':str})
    t=table[table.geometry.eq('primary')&table.cluster.eq(candidate)].sort_values('sample_name')
    fig,axes=plt.subplots(1,2,figsize=(12.5,5.8))
    x=np.arange(len(t));width=.34
    ax=axes[0]
    ax.bar(x-width/2,t.cluster_detection_rate*100,width,color='#ad4e22',label='C'+candidate+' candidate')
    ax.bar(x+width/2,t.other_detection_rate*100,width,color='#7c98b1',label='Other macrophages')
    ax.set(xticks=x,xticklabels=t.sample_name,xlabel='Patient',ylabel='PGAM5 RNA detection (%)',
           title='PGAM5 enrichment in all seven patients')
    ax.tick_params(axis='x',rotation=35)
    for j,row in enumerate(t.itertuples()):
        ax.text(j,max(row.cluster_detection_rate,row.other_detection_rate)*100+1,f'n={row.cluster_cells}',ha='center',fontsize=8)
    ax.legend(frameon=False,fontsize=9,loc='upper right')
    ax=axes[1]
    ax.bar(x-width/2,t.PGAM5_CP10k_cluster_mean,width,color='#ad4e22',label='C'+candidate+' candidate')
    ax.bar(x+width/2,t.PGAM5_CP10k_other_mean,width,color='#7c98b1',label='Other macrophages')
    ax.set(xticks=x,xticklabels=t.sample_name,xlabel='Patient',ylabel='Mean PGAM5 CP10k (including zeros)',
           title='Normalized PGAM5 expression, including zero counts')
    ax.tick_params(axis='x',rotation=35);ax.legend(frameon=False,fontsize=9,loc='upper right')
    for ax in axes:ax.spines[['top','right']].set_visible(False)
    fig.suptitle('GSE202642: patient evidence for C'+candidate+' cycling macrophage candidate',fontsize=12,y=.98)
    fig.text(.5,.035,'All seven patients shown; five have >=20 candidate cells and >=20 other macrophages.\n'
             'Two small patient subsets (n=7 and n=12) are descriptive. No depth matching or regression performed.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.84,bottom=.25,left=.08,right=.98,wspace=.32)
    for ext in ['png','pdf']:fig.savefig(R/f'candidate_patient_consistency.{ext}',dpi=200)
    plt.close(fig)
    # Source-backed dot plot: color is mean log1p CP10k; size is raw detection fraction.
    genes=pd.read_csv(R/'portable_feature_genes.csv').gene.to_numpy()
    raw=sparse.load_npz(R/'portable_raw_feature_counts.npz')
    norm=raw.multiply((1e4/obs.total_counts.to_numpy())[:,None]).tocsr();log=norm.copy();log.data=np.log1p(log.data)
    selected=['C1QA','C1QB','C1QC','CD68','CSF1R','MSR1','PGAM5','MKI67','TOP2A','CENPF','PCLAF','TYMS','CDK1','TPX2','PRC1','GTSE1']
    points=[]
    for yi,cluster in enumerate(association.cluster):
        sel=labels==cluster
        for xi,gene in enumerate(selected):
            j=int(np.flatnonzero(genes==gene)[0]);points.append({'cluster':cluster,'gene':gene,'x':xi,'y':yi,
                'detection_fraction':float((raw[sel,j]>0).mean()),'mean_log1p_CP10k':float(log[sel,j].mean())})
    dots=pd.DataFrame(points);dots.to_csv(R/'marker_dotplot_values.csv',index=False)
    fig,ax=plt.subplots(figsize=(12.6,7))
    v=ax.scatter(dots.x,dots.y,s=4+100*dots.detection_fraction,c=dots.mean_log1p_CP10k,cmap='viridis',vmin=0,
                  vmax=dots.mean_log1p_CP10k.max(),edgecolors='none')
    ax.set(xticks=range(len(selected)),xticklabels=selected,yticks=range(len(association)),
           yticklabels=['C'+c for c in association.cluster],ylim=(len(association)-.4,-.6),
           title='Macrophage identity and cycling program (PGAM5 not used in clustering)')
    ax.tick_params(axis='x',rotation=45);ax.axhspan(int(candidate)-.45,int(candidate)+.45,color='#ad4e22',alpha=.06,zorder=-1)
    fig.colorbar(v,ax=ax,label='Mean log1p(CP10k)',pad=.02)
    for fraction in [.1,.5,1.]:ax.scatter([],[],s=4+100*fraction,c='#777777',label=f'{fraction:.0%} detected')
    ax.legend(frameon=False,bbox_to_anchor=(1.19,1),loc='upper left',fontsize=9)
    ax.spines[['top','right']].set_visible(False)
    fig.text(.5,.025,'C'+candidate+' is a cycling RNA cluster supported by macrophage markers; PGAM5 RNA is detected in only a subset.\n'
             'Dot size: raw RNA detection fraction. Color: mean log1p full-library CP10k. Values are exported for audit.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.91,bottom=.24,left=.07,right=.80)
    for ext in ['png','pdf']:fig.savefig(R/f'macrophage_cycling_marker_dotplot.{ext}',dpi=200)
    plt.close(fig)
    print('Rendered overview, patient evidence and marker dot plot. Candidate C'+candidate,flush=True)

if __name__=='__main__':main()
