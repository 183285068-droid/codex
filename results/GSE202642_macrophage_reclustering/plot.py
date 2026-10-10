from pathlib import Path
import scanpy as sc,pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt,json
R=Path(__file__).resolve().parent
a=sc.read_h5ad(R/'macrophages_full.h5ad');names=json.loads((R/'names.json').read_text());a.obs['subtype']=a.obs.mac_cluster.astype(str).map(names).astype('category');xy=a.obsm['X_umap'];groups=sorted(a.obs.mac_cluster.unique(),key=lambda c:int(c[1:]));colors=plt.get_cmap('tab20')
def save(fig,n):fig.tight_layout();fig.savefig(R/(n+'.png'),dpi=180);fig.savefig(R/(n+'.pdf'));plt.close(fig)
def axis(ax):ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2')
fig,ax=plt.subplots(figsize=(11,7))
for i,c in enumerate(groups):
 ix=a.obs.mac_cluster.eq(c).to_numpy();ax.scatter(*xy[ix].T,s=3,color=colors(i%20),label=f'{c} {names[c]} (n={ix.sum()})',linewidths=0,rasterized=True);ax.text(*np.median(xy[ix],axis=0),c,ha='center',bbox=dict(facecolor='white',alpha=.8,edgecolor='none',pad=1),fontsize=9)
axis(ax);ax.set_title('GSE202642 | re-clustered macrophages');ax.legend(loc='center left',bbox_to_anchor=(1,0.5),fontsize=8,markerscale=3);save(fig,'UMAP_named_macrophages')
fig,axs=plt.subplots(1,2,figsize=(14,6))
for i,sample in enumerate(a.obs.sample_name.cat.categories):
 ix=a.obs.sample_name.eq(sample).to_numpy();axs[0].scatter(*xy[ix].T,s=3,color=colors(i),label=sample,linewidths=0,rasterized=True)
axs[0].legend(fontsize=8,markerscale=3);axs[0].set_title('Harmony UMAP: samples')
for i,c in enumerate(groups):
 ix=a.obs.mac_cluster.eq(c).to_numpy();axs[1].scatter(*a.obsm['X_umap_unintegrated'][ix].T,s=3,color=colors(i%20),linewidths=0,rasterized=True)
axs[1].set_title('Unintegrated UMAP: same cluster labels')
for ax in axs:axis(ax)
save(fig,'UMAP_sample_and_integration')
j=a.var_names.get_loc('PGAM5');raw=a.layers['counts'][:,j].toarray().ravel();cp=raw/a.obs.total_counts.to_numpy()*10000;log=np.log1p(cp);a.obs['PGAM5_CP10k']=cp;a.obs['PGAM5_log1p_CP10k']=log
fig,axs=plt.subplots(1,2,figsize=(14,6));pos=np.flatnonzero(raw>0);pos=pos[np.argsort(log[pos])]
for ax in axs:ax.scatter(*xy.T,s=2,color='#d5d5d5',linewidths=0,rasterized=True);axis(ax)
p=axs[0].scatter(*xy[pos].T,s=4,c=log[pos],vmin=0,vmax=log.max(),cmap='viridis',linewidths=0,rasterized=True);fig.colorbar(p,ax=axs[0],label='PGAM5 log1p(CP10k)');axs[0].set_title('PGAM5 normalized RNA expression');axs[1].scatter(*xy[pos].T,s=4,color='#c6363d',linewidths=0,rasterized=True);axs[1].set_title(f'PGAM5 RNA detected: {len(pos)}/{len(a)} cells')
for ax in axs:
 for c in groups:ax.text(*np.median(xy[a.obs.mac_cluster.eq(c)],axis=0),c,fontsize=8,bbox=dict(facecolor='white',alpha=.8,edgecolor='none',pad=1))
save(fig,'UMAP_PGAM5')
genes=['C1QA','CSF1R','FOLR2','CD163','LYVE1','MRC1','TREM2','SPP1','APOE','MARCO','CD5L','HLA-DRA','CD74','IL1B','CXCL8','FCN1','S100A8','MKI67','TOP2A','ISG15','IFIT1','MT1G','MT2A','PGAM5']
f=sc.pl.dotplot(a,genes,groupby='mac_cluster',standard_scale='var',show=False,return_fig=True);f.savefig(R/'marker_dotplot.png',dpi=180);f.savefig(R/'marker_dotplot.pdf')
rows=[]
for c in groups:
 ix=a.obs.mac_cluster.eq(c).to_numpy();samples=a.obs.loc[ix,'sample_name'].value_counts(normalize=True);rows.append(dict(cluster=c,name=names[c],cells=int(ix.sum()),patients=int((samples>0).sum()),largest_patient_fraction=float(samples.max()),PGAM5_detected_cells=int((raw[ix]>0).sum()),PGAM5_detection_percent=float((raw[ix]>0).mean()*100),PGAM5_mean_CP10k=float(cp[ix].mean())))
pd.DataFrame(rows).to_csv(R/'cluster_summary.csv',index=False);a.obs.to_csv(R/'cell_annotations.csv.gz');a.obs.groupby(['mac_cluster','sample_name'],observed=True).size().rename('cells').reset_index().to_csv(R/'cluster_sample_counts.csv',index=False)
portable=list(dict.fromkeys(genes+pd.read_csv(R/'top50_markers.csv.gz').names.tolist()));portable=[g for g in portable if g in a.var_names];a[:,portable].copy().write_h5ad(R/'macrophages_marker_subset.h5ad',compression='gzip')
print(pd.DataFrame(rows).to_string(index=False))
