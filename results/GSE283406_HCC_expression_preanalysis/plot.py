from pathlib import Path
import gzip
from io import BytesIO
import json
import os

os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE283406_mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import anndata as ad
import numpy as np
import pandas as pd
from PIL import Image

ROOT=Path(__file__).resolve().parent
WORK=Path(os.environ.get('GSE283406_WORK',ROOT/'work'))
OUT=ROOT/'figures'
OUT.mkdir(exist_ok=True)
PROGRAMS=json.loads((ROOT/'programs.json').read_text())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42})
DATA={s:ad.read_h5ad(ROOT/f'{s}_AVITI_expression_marker_subset.h5ad') for s in ['J259B','J259X']}
LABELS={'Macrophage_core':'Macrophage-associated', 'Pan_myeloid':'Myeloid',
        'T_lineage':'T-lineage','B_cell':'B-lineage','Plasma':'Plasma-associated',
        'DC':'DC-associated','Monocyte_inflammatory':'Inflammatory myeloid',
        'Fibroblast_matrix':'Matrix-associated','Endothelial':'Endothelial-associated',
        'Hepatocyte':'Hepatocytic','Epithelial':'Epithelial','Cycling':'Cycling',
        'Hypoxia':'Hypoxia-associated',
        'Acute_phase':'Acute phase'}
DESCRIPTIONS={
    ('J259B','E0'):'Metabolic / matrix-mixed',
    ('J259B','E1'):'Immunoglobulin-associated',
    ('J259B','E2'):'Hepatic secretory / acute-phase',
    ('J259B','E3'):'Matrix-associated',
    ('J259X','E0'):'Hepatic metabolic / T-marker-mixed',
    ('J259X','E1'):'Immunoglobulin / matrix-mixed',
    ('J259X','E2'):'Hepatic metabolic / secretory',
}
PALETTE=['#4477AA','#EE6677','#228833','#CCBB44','#66CCEE','#AA3377','#BBBBBB',
         '#332288','#44AA99','#999933','#882255','#6699CC','#DDCC77','#CC6677','#117733','#AA4499']

def save(fig,name):
    fig.savefig(OUT/f'{name}.png',dpi=210,bbox_inches='tight')
    fig.savefig(OUT/f'{name}.pdf',bbox_inches='tight')
    plt.close(fig)

def names():
    records=[]
    for section,a in DATA.items():
        markers=pd.read_csv(ROOT/f'{section}_AVITI_top100_markers.csv')
        for cluster in sorted(a.obs.expression_cluster.unique(),key=lambda x:int(x[1:])):
            mask=a.obs.expression_cluster.eq(cluster).to_numpy()
            m=markers[markers.group.eq(cluster) & markers.logfoldchanges.gt(0) & markers.pct_nz_group.ge(.10)]
            candidates=[g for g in m.names if g in a.var_names and not g.startswith(('MT-','RPL','RPS'))]
            selected=[]
            for gene in candidates:
                x=a[:,[gene]].X.toarray().ravel()
                if x[mask].mean() > x[~mask].mean():
                    selected.append(gene)
                if len(selected)==2:break
            assert len(selected)==2
            standardized={k:(a.obs.loc[mask,k].mean()-a.obs[k].mean())/(a.obs[k].std()+1e-8) for k in LABELS}
            best=max(standardized,key=standardized.get)
            description=DESCRIPTIONS.get((section,cluster),'Marker-associated mixed')
            label=f'{cluster}: {"/".join(selected)} | {description} RNA'
            for gene in selected:
                x=a[:,[gene]].X.toarray().ravel()
                records.append(dict(section=section,expression_cluster=cluster,cluster_name=label,
                    marker=gene,cluster_mean_log1p_CP10k=x[mask].mean(),other_mean_log1p_CP10k=x[~mask].mean(),
                    RNA_program=best,program_mean=a.obs.loc[mask,best].mean(),
                    program_all_spot_mean=a.obs[best].mean(),n_spots=int(mask.sum())))
    evidence=pd.DataFrame(records)
    evidence.to_csv(ROOT/'expression_cluster_naming_evidence.csv',index=False)
    annotation=evidence.drop_duplicates(['section','expression_cluster']).drop(columns=['marker','cluster_mean_log1p_CP10k','other_mean_log1p_CP10k'])
    annotation.to_csv(ROOT/'expression_cluster_names.csv',index=False)
    return annotation

ANNOTATION=names()

def clusters():
    fig,axs=plt.subplots(1,2,figsize=(15,7.8))
    for ax,(section,a) in zip(axs,DATA.items()):
        xy=a.obsm['X_umap'];cl=sorted(a.obs.expression_cluster.unique(),key=lambda x:int(x[1:]))
        ann=ANNOTATION[ANNOTATION.section.eq(section)].set_index('expression_cluster')
        for i,c in enumerate(cl):
            mask=a.obs.expression_cluster.eq(c).to_numpy()
            ax.scatter(xy[mask,0],xy[mask,1],s=4,c=PALETTE[i%len(PALETTE)],rasterized=True)
            center=np.median(xy[mask],axis=0)
            ax.text(*center,c,ha='center',va='center',fontsize=10,bbox=dict(facecolor='white',alpha=.75,edgecolor='none',pad=1))
        ax.set_title(f'{section}: {len(a):,} QC spots, {len(cl)} expression groups\n'+('HCC tumor section' if section=='J259X' else 'HCC + adjacent section; no spot-level pathology labels'))
        ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2');ax.spines[['top','right']].set_visible(False)
        text='\n'.join(ann.loc[c,'cluster_name'] for c in cl)
        ax.text(0,-.12,text,transform=ax.transAxes,va='top',fontsize=10,linespacing=1.55)
    fig.suptitle('GSE283406 AVITI: expression UMAP (not tissue coordinates)',fontsize=17)
    fig.text(.04,.012,'Each point is a mixed-cell Visium spot. Names are provisional RNA annotations; group IDs are section-specific.',fontsize=10)
    fig.subplots_adjust(top=.85,bottom=.29,wspace=.22)
    save(fig,'AVITI_expression_clusters_UMAP')

def features():
    keys=['PGAM5_log1p_CP10k','Macrophage_core','T_lineage','Exhaustion_core']
    fig,axs=plt.subplots(2,4,figsize=(18,9))
    for j,key in enumerate(keys):
        maximum=max(a.obs.loc[a.obs.T_RNA_enriched,key].max() if key=='Exhaustion_core' else a.obs[key].max() for a in DATA.values())
        norm=Normalize(0,maximum)
        for i,(section,a) in enumerate(DATA.items()):
            ax=axs[i,j];xy=a.obsm['X_umap'];v=a.obs[key].to_numpy()
            mask=a.obs.T_RNA_enriched.to_numpy() if key=='Exhaustion_core' else np.ones(len(a),dtype=bool)
            ax.scatter(xy[:,0],xy[:,1],s=3,c='#D4D8DC',rasterized=True)
            ids=np.flatnonzero(mask);ids=ids[np.argsort(v[ids])]
            scatter=ax.scatter(xy[ids,0],xy[ids,1],c=v[ids],s=4,cmap='YlOrRd',norm=norm,rasterized=True)
            ax.set_title(f'{section}: {key.replace("_"," ")}'+ ('\nT-RNA gate only' if key=='Exhaustion_core' else ''),fontsize=11)
            ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2');ax.spines[['top','right']].set_visible(False)
            fig.colorbar(scatter,ax=ax,fraction=.045,pad=.025)
    fig.suptitle('GSE283406: PGAM5 and immune RNA programs in expression UMAP',fontsize=17)
    fig.text(.025,.014,'Not tissue coordinates. PGAM5: log1p(CP10k); programs: marker mean log1p(CP10k), including zeros.\nGrey in exhaustion panels = outside T-RNA gate; exhaustion-associated transcripts do not confirm terminal exhaustion.',fontsize=10)
    fig.tight_layout(rect=[0,.07,1,.94]);save(fig,'PGAM5_immune_expression_UMAP')

def dotplot():
    genes=['PGAM5','C1QA','C1QB','CD68','CSF1R','LST1','LYZ','CD3D','TRAC','CD79A','JCHAIN','COL1A1','VWF','ALB','KRT8','MKI67','CA9','CD74','PDCD1','TIGIT']
    rows=[]
    fig,axs=plt.subplots(2,1,figsize=(15,8))
    vmax=0
    for section,a in DATA.items():
        for cluster in sorted(a.obs.expression_cluster.unique(),key=lambda x:int(x[1:])):
            mask=a.obs.expression_cluster.eq(cluster).to_numpy()
            for gene in genes:
                x=a[:,[gene]].X.toarray().ravel();raw=a[:,[gene]].layers['counts'].toarray().ravel()
                mean=float(x[mask].mean());vmax=max(vmax,mean)
                rows.append(dict(section=section,expression_cluster=cluster,gene=gene,
                    mean_log1p_CP10k=mean,detection_percent=float(100*(raw[mask]>0).mean()),n_spots=int(mask.sum())))
    table=pd.DataFrame(rows);table.to_csv(ROOT/'marker_dotplot_values.csv',index=False)
    for ax,(section,a) in zip(axs,DATA.items()):
        cl=sorted(a.obs.expression_cluster.unique(),key=lambda x:int(x[1:]));d=table[table.section.eq(section)]
        ix=d.gene.map({g:i for i,g in enumerate(genes)});iy=d.expression_cluster.map({c:i for i,c in enumerate(cl)})
        points=ax.scatter(ix,iy,s=d.detection_percent*1.8,c=d.mean_log1p_CP10k,cmap='YlOrRd',norm=Normalize(0,vmax),edgecolors='#707070',linewidths=.2)
        ax.set_xticks(range(len(genes)),genes,rotation=45,ha='right');ax.set_yticks(range(len(cl)),cl);ax.invert_yaxis()
        ax.set_title(section+' AVITI expression groups');ax.set_xlim(-.6,len(genes)-.4);ax.spines[['top','right']].set_visible(False)
        fig.colorbar(points,ax=ax,label='Mean log1p(CP10k)',fraction=.023,pad=.03)
    for pct in [10,50,100]:axs[0].scatter([],[],s=pct*1.8,c='none',edgecolor='#707070',label=f'{pct}%')
    axs[0].legend(title='Detected spots',loc='upper left',bbox_to_anchor=(1.08,1),frameon=False)
    fig.suptitle('GSE283406: expression marker dotplot (mixed spots)',fontsize=17)
    fig.text(.03,.008,'Dot area = detected spot percentage; colour = mean log1p(CP10k), including zeros. Gene-expression groups are not pure cell types.',fontsize=10)
    fig.tight_layout(rect=[0,.045,.91,.95]);save(fig,'AVITI_marker_dotplot')

def clusterPGAM5():
    fig,axs=plt.subplots(2,2,figsize=(13,9))
    for i,section in enumerate(DATA):
        d=pd.read_csv(ROOT/f'{section}_AVITI_cluster_profiles.csv')
        d=d.sort_values('expression_cluster',key=lambda x:x.str[1:].astype(int))
        for j,(key,label) in enumerate([('PGAM5_mean_log1p_CP10k','Mean log1p(CP10k)'),('PGAM5_detection_pct','Detected spots (%)')]):
            ax=axs[i,j];bars=ax.bar(d.expression_cluster,d[key],color='#C65D6D')
            for bar,value in zip(bars,d[key]):ax.text(bar.get_x()+bar.get_width()/2,bar.get_height(),f'{value:.2f}',ha='center',va='bottom',fontsize=9)
            ax.set_title(section);ax.set_ylabel(label);ax.set_ylim(0,1.18*max(pd.read_csv(ROOT/f'{s}_AVITI_cluster_profiles.csv')[key].max() for s in DATA))
            ax.spines[['top','right']].set_visible(False)
    fig.suptitle('GSE283406: PGAM5 expression by AVITI expression group',fontsize=17)
    fig.text(.03,.01,'Zeros included. Descriptive comparisons only; these groups are not spatial domains or single-cell macrophage subgroups.',fontsize=10)
    fig.tight_layout(rect=[0,.05,1,.95]);save(fig,'PGAM5_expression_group_comparison')

def tissue():
    fig,axs=plt.subplots(1,2,figsize=(11,6))
    for ax,section,gsm in zip(axs,['J259B','J259X'],['GSM8661791','GSM8661792']):
        path=WORK/'inputs'/f'{gsm}_{section}.tissue_hires_image.png.gz'
        image=np.asarray(Image.open(BytesIO(gzip.decompress(path.read_bytes()))))
        ax.imshow(image);ax.axis('off');ax.set_title(section+ (' (HCC + adjacent)' if section=='J259B' else ' (HCC tumor)'))
    fig.suptitle('GSE283406: publicly deposited tissue images',fontsize=16)
    fig.text(.035,.022,'Original images only. No barcode-to-pixel coordinates or pathology masks were provided; no expression overlay reconstructed.',fontsize=10)
    fig.tight_layout(rect=[0,.055,1,.93]);save(fig,'original_tissue_images_no_expression_overlay')

def technical():
    fig,axs=plt.subplots(2,3,figsize=(15,9))
    summaries=pd.read_csv(ROOT/'paired_platform_summary.csv').set_index('section')
    samples={'J259B':('GSM8661789','GSM8661791'),'J259X':('GSM8661790','GSM8661792')}
    for i,(section,gsms) in enumerate(samples.items()):
        for gsm,platform,color in zip(gsms,['NextSeq','AVITI'],['#4477AA','#CC6677']):
            d=pd.read_csv(ROOT/f'{gsm}_QC_expression_annotations.csv')
            axs[i,0].hist(np.log10(d.total_counts),bins=np.linspace(2.6,4.7,35),density=True,alpha=.45,color=color,label=f'{platform} n={len(d)}')
            axs[i,1].hist(d.n_genes,bins=np.linspace(0,10000,35),density=True,alpha=.45,color=color)
        axs[i,0].set_title(section+': QC UMI distribution');axs[i,0].set_xlabel('log10(total UMI)');axs[i,0].set_ylabel('Density');axs[i,0].legend(frameon=False)
        axs[i,1].set_title(section+': detected genes');axs[i,1].set_xlabel('Number of detected features');axs[i,1].set_ylabel('Density')
        d=pd.read_csv(ROOT/f'{section}_paired_platform_expression.csv');r=summaries.loc[section]
        axs[i,2].scatter(d.NextSeq_PGAM5_log1p_CP10k,d.AVITI_PGAM5_log1p_CP10k,s=6,alpha=.4,c='#4477AA',rasterized=True)
        maximum=max(d.NextSeq_PGAM5_log1p_CP10k.max(),d.AVITI_PGAM5_log1p_CP10k.max())*1.04
        axs[i,2].plot([0,maximum],[0,maximum],color='#777777',linewidth=.8)
        axs[i,2].set(xlim=(-.03,maximum),ylim=(-.03,maximum),xlabel='NextSeq PGAM5 log1p(CP10k)',ylabel='AVITI PGAM5 log1p(CP10k)')
        axs[i,2].set_title(f'{section}: same-spot PGAM5\nSpearman rho={r.PGAM5_normalized_spearman:.3f}')
        axs[i,2].text(.03,.97,f'Shared QC spots: {len(d)}\nDetected only in AVITI: {int(r.PGAM5_only_AVITI)}\nDetected only in NextSeq: {int(r.PGAM5_only_NextSeq)}',transform=axs[i,2].transAxes,va='top',fontsize=10)
        for ax in axs[i]:ax.spines[['top','right']].set_visible(False)
    fig.suptitle('GSE283406: QC and paired sequencing measurements',fontsize=17)
    fig.text(.03,.012,'The two platforms measure the same section/library barcodes; they are technical comparisons, not independent patients. Zeros retained.',fontsize=10)
    fig.tight_layout(rect=[0,.05,1,.94]);save(fig,'QC_and_paired_platform_comparison')

if __name__=='__main__':
    clusters();features();dotplot();clusterPGAM5();tissue();technical()
    print(ANNOTATION.to_string(index=False))
