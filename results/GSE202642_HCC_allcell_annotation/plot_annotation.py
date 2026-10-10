from pathlib import Path
import numpy as np,pandas as pd,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R=Path(__file__).resolve().parent

def main():
 a=pd.read_csv(R/'cell_annotations.csv.gz',index_col=0,dtype={'cluster':str});u=pd.read_csv(R/'UMAP_coordinates.csv.gz',index_col=0);assert a.index.equals(u.index);defs=pd.read_csv(R/'cluster_annotations.csv',dtype={'cluster':str});types=defs.cell_type_EN.drop_duplicates().tolist();colors={t:plt.get_cmap('tab20')(i%20) for i,t in enumerate(types)};rng=np.random.default_rng(202642);order=rng.permutation(len(a))
 for key,fn in [('cell_type_EN','UMAP_cell_types'),('cluster','UMAP_clusters'),('cluster_label','UMAP_named_clusters'),('sample_name','UMAP_samples')]:
  fig,ax=plt.subplots(figsize=((16,9) if key=='cluster_label' else (12,8)));cats=(types if key=='cell_type_EN' else defs.cluster.tolist() if key=='cluster' else defs.cluster_label.tolist() if key=='cluster_label' else sorted(a[key].unique()));pal=colors if key=='cell_type_EN' else {t:plt.get_cmap('turbo')(i/max(1,len(cats)-1)) for i,t in enumerate(cats)}
  for t in cats:
   ix=order[a.iloc[order][key].eq(t).to_numpy()];ax.scatter(u.iloc[ix].UMAP1,u.iloc[ix].UMAP2,s=1.6,c=[pal[t]],rasterized=True,alpha=.7,linewidths=0)
  if key=='cell_type_EN':
   for t in cats:
    v=u[a[key].eq(t)];ax.text(v.UMAP1.median(),v.UMAP2.median(),t,fontsize=9,ha='center',bbox={'facecolor':'white','alpha':.8,'edgecolor':'none','pad':1})
  elif key in ['cluster','cluster_label']:
   for t in cats:
    v=u[a[key].eq(t)];ax.text(v.UMAP1.median(),v.UMAP2.median(),('C'+t if key=='cluster' else t.split(' ')[0]),fontsize=9,ha='center',bbox={'facecolor':'white','alpha':.8,'edgecolor':'none','pad':1})
  handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=pal[t],label=('C'+t if key=='cluster' else t),markersize=5) for t in cats];ax.legend(handles=handles,bbox_to_anchor=(1.02,1),loc='upper left',fontsize=8,frameon=False,ncol=2 if len(cats)>25 else 1);ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2');ax.set_title('GSE202642 | 7 HCC tumor samples | Harmony + Leiden');ax.set_xticks([]);ax.set_yticks([]);fig.tight_layout();fig.savefig(R/(fn+'.png'),dpi=180);fig.savefig(R/(fn+'.pdf'));plt.close(fig)
 fig,axs=plt.subplots(1,2,figsize=(17,7));
 for ax,cols,title in zip(axs,[['unintegrated_UMAP1','unintegrated_UMAP2'],['UMAP1','UMAP2']],['Unintegrated PCA neighbors','Sample-adjusted Harmony neighbors']):
  for t in types:
   ix=order[a.iloc[order].cell_type_EN.eq(t).to_numpy()];ax.scatter(u.iloc[ix][cols[0]],u.iloc[ix][cols[1]],s=1.3,c=[colors[t]],alpha=.7,rasterized=True,linewidths=0)
  ax.set_title(title);ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2')
 fig.legend(handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=colors[t],label=t,markersize=5) for t in types],loc='lower center',ncol=5,fontsize=8,frameon=False);fig.tight_layout(rect=[0,.12,1,1]);fig.savefig(R/'UMAP_integration_comparison.png',dpi=180);fig.savefig(R/'UMAP_integration_comparison.pdf');plt.close(fig)
 markers=['CD3D','CD3E','TRAC','IL7R','CD8A','FOXP3','NKG7','GNLY','MS4A1','CD79A','MZB1','JCHAIN','C1QA','C1QB','CSF1R','CD68','LST1','FCN1','S100A8','CD14','CD1C','FCER1A','CLEC9A','LAMP3','IL3RA','TCF4','PECAM1','VWF','COL1A1','DCN','RGS5','ACTA2','ALB','APOA1','TTR','EPCAM','KRT19','TPSAB1','KIT','FCGR3B','CSF3R','MKI67','TOP2A']
 e=pd.read_csv(R/'cluster_marker_evidence.csv',dtype={'cluster':str});joined=e.merge(defs[['cluster','cell_type_EN']],on='cluster');rows=[]
 for (t,g),d in joined.groupby(['cell_type_EN','gene']):rows.append({'cell_type':t,'gene':g,'fraction':np.average(d.fraction,weights=d.cells),'mean_log1p_CP10k':np.average(d.mean_log1p_CP10k,weights=d.cells)})
 pd.DataFrame(rows).to_csv(R/'celltype_marker_evidence.csv',index=False)
 for level,evidence,key,groups,fn in [('type',pd.DataFrame(rows),'cell_type',types,'celltype_marker_dotplot'),('cluster',e,'cluster',defs.cluster.tolist(),'cluster_marker_dotplot')]:
  mean=evidence.pivot(index=key,columns='gene',values='mean_log1p_CP10k').reindex(index=groups,columns=markers);frac=evidence.pivot(index=key,columns='gene',values='fraction').reindex(index=groups,columns=markers);scaled=mean/mean.max(axis=0).replace(0,1);fig,ax=plt.subplots(figsize=(19,max(7,len(groups)*.34)));xx,yy=np.meshgrid(np.arange(len(markers)),np.arange(len(groups)));pts=ax.scatter(xx.ravel(),yy.ravel(),s=frac.to_numpy().ravel()*100,c=scaled.to_numpy().ravel(),cmap='Reds',vmin=0,vmax=1,linewidths=.2,edgecolors='#999999');ax.set_xticks(range(len(markers)),markers,rotation=90,fontsize=8);labs=groups if level=='type' else ['C'+c+' '+defs.set_index('cluster').loc[c,'cell_type_EN'] for c in groups];ax.set_yticks(range(len(groups)),labs,fontsize=9);ax.invert_yaxis();ax.set_xlim(-.8,len(markers)-.2);ax.set_title('Marker RNA evidence | dot size: detected fraction | color: per-gene scaled mean');fig.colorbar(pts,ax=ax,shrink=.5,label='Scaled mean log1p(CP10k)');fig.tight_layout();fig.savefig(R/(fn+'.png'),dpi=180);fig.savefig(R/(fn+'.pdf'));plt.close(fig)
 counts=a.groupby(['sample_name','cell_type_EN']).size().unstack(fill_value=0).reindex(columns=types);counts.to_csv(R/'sample_celltype_counts.csv');frac=counts.div(counts.sum(axis=1),axis=0);frac.to_csv(R/'sample_celltype_fractions.csv');fig,ax=plt.subplots(figsize=(12,6));frac.plot.bar(stacked=True,color=[colors[t] for t in types],ax=ax);ax.set_ylabel('Fraction among QC-passing tumor cells');ax.set_xlabel('HCC sample');ax.legend(bbox_to_anchor=(1.02,1),fontsize=8,frameon=False);ax.set_title('Sample composition (RNA-inferred labels)');fig.tight_layout();fig.savefig(R/'sample_celltype_composition.png',dpi=180);fig.savefig(R/'sample_celltype_composition.pdf');plt.close(fig)
 print('Figures complete',flush=True)
if __name__=='__main__':main()
