import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE281759_distribution_mpl');os.environ.setdefault('XDG_CACHE_HOME','/tmp/GSE281759_distribution_cache')
import numpy as np,pandas as pd,anndata as ad,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from analyze import R,BASE,SAMPLES,METRICS
NAMES=['High-high neighbor edges','Median nearest-high spacing','Largest high component fraction','3x3 regional heterogeneity','Median boundary proxy: high - other'];COL=['#285A8E','#B5822C','#A8536A','#607C44'];P=['A009','B013','C004','C008']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
def save(f,n):
 for ext in ['png','pdf']:f.savefig(R/'figures'/f'{n}.{ext}',dpi=160,bbox_inches='tight')
 plt.close(f)
def component_plot():
 components=[]
 for item in SAMPLES:
  o=pd.read_csv(R/f'{item["sample"]}_annotations.csv.gz',index_col=0);sizes=o.loc[o.PGAM5_mac_group.eq('PGAM5_high'),'high_component_id'].value_counts()
  components.extend(dict(**item,component_id=int(cid),component_size=int(size)) for cid,size in sizes.items())
 cc=pd.DataFrame(components);cc.to_csv(R/'high_connected_components.csv',index=False);freq=cc.component_size.value_counts().sort_index();f,ax=plt.subplots(figsize=(8,5));ax.bar(freq.index.astype(str),freq.to_numpy(),color=COL[0]);ax.set_xlabel('Number of high RNA spots in a direct-neighbor component');ax.set_ylabel('Components across 16 sections');ax.set_title('First-ring component size distribution\nA spot is not an individual macrophage')
 for i,v in enumerate(freq):ax.text(i,v+3,str(v),ha='center')
 f.tight_layout();save(f,'high_component_size_distribution')
def run():
 d=pd.read_csv(R/'section_summary.csv');s=pd.read_csv(R/'spatial_distribution_statistics.csv');p=pd.read_csv(R/'patient_descriptive.csv');bins=pd.read_csv(R/'spatial_bin_counts.csv');main=s[s.null_kind.eq('quality_stratified')]
 f,axs=plt.subplots(1,3,figsize=(15,7),sharey=True);y=np.arange(len(d));colors=[COL[P.index(v)] for v in d.patient]
 for ax,col,title in zip(axs,['n_high','high_pct_mac','high_with_high_neighbor'],['High source spots','High / Mac-RNA spots (%)','High spots with another high neighbor']):ax.barh(y,d[col],color=colors);ax.set_title(title);ax.set_xlabel('Count' if '%' not in title else 'Percent');ax.spines[['top','right']].set_visible(False)
 axs[0].set_yticks(y,d['sample']);axs[0].invert_yaxis();f.suptitle('Fixed PGAM5-high macrophage-RNA distribution by section\nSection-specific Q75 among detected sources; proportions are not resolved cell fractions');f.tight_layout();save(f,'section_high_distribution_summary')
 f,axs=plt.subplots(1,5,figsize=(22,9),sharey=True)
 for j,k in enumerate(METRICS):
  z=main[main.metric==k].set_index('sample').reindex(d['sample']);ax=axs[j]
  for i,row in enumerate(z.itertuples()):
   ax.plot([row.null_q025,row.null_q975],[i,i],color='#BBBBBB',linewidth=3);ax.scatter(row.observed,i,color=COL[P.index(row.patient)],s=32,zorder=3)
   if row.primary_q_80<.05:ax.text(row.observed,i,' *',va='center',fontsize=13)
  ax.set_title(NAMES[j],fontsize=10);ax.set_xlabel(['Edges','Grid-neighbor spacing','Fraction','Statistic','Lattice hops'][j]);ax.spines[['top','right']].set_visible(False)
 axs[0].set_yticks(y,d['sample']);axs[0].invert_yaxis();f.suptitle('Observed spatial metrics vs quality-stratified random labels\nGrey: central 95% conditional null range (not a confidence interval); *: primary 80-test BH q<0.05');f.tight_layout();save(f,'spatial_metrics_vs_conditional_null')
 mat=main.pivot(index='sample',columns='metric',values='primary_q_80').reindex(index=d['sample'],columns=METRICS)
 for row in main.itertuples():
  if row.status!='tested':mat.loc[row.sample,row.metric]=np.nan
 f,ax=plt.subplots(figsize=(12,8));im=ax.imshow(mat,aspect='auto',vmin=0,vmax=1,cmap='Blues');ax.set_xticks(range(5),NAMES,rotation=30,ha='right');ax.set_yticks(range(len(d)),d['sample'])
 for i in range(len(d)):
  for j in range(5):ax.text(j,i,f'{mat.iloc[i,j]:.3f}' if np.isfinite(mat.iloc[i,j]) else 'NA',ha='center',va='center',fontsize=8)
 ax.set_title('Primary BH q: 80 preplanned section-metric tests\nNA: constant null or insufficient groups');f.colorbar(im,ax=ax,shrink=.6);f.tight_layout();save(f,'spatial_distribution_FDR')
 component_plot()
 for patient in P:
  items=[v for v in SAMPLES if v['patient']==patient];f,axs=plt.subplots(len(items),3,figsize=(13,3.5*len(items)),squeeze=False)
  for row,item in enumerate(items):
   name=item['sample'];a=ad.read_h5ad(BASE/f'{name}_spatial_marker_subset.h5ad');o=pd.read_csv(R/f'{name}_annotations.csv.gz',index_col=0).loc[a.obs_names];xy=a.obsm['spatial']*a.uns['scalefactors']['tissue_lowres_scalef'];image=a.uns['spatial'][name]['images']['lowres'];high=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy();mac=o.Mac_RNA_enriched.to_numpy(bool);pos=o.source_candidate_detected.to_numpy(bool);e=a.obsp['spatial_connectivities'].tocoo();edges=high[e.row]&high[e.col];component=o.loc[high,'high_component_id'].value_counts();cluster=o.loc[high,'high_component_id'].map(component).to_numpy()>1
   for ax in axs[row]:ax.imshow(image);ax.axis('off')
   ax=axs[row,0];ax.scatter(xy[mac,0],xy[mac,1],s=6,color='#BBBBBB',alpha=.65,label='Mac RNA');ax.scatter(xy[pos&~high,0],xy[pos&~high,1],s=9,color=COL[0],alpha=.75,label='Other PGAM5 detected');ax.scatter(xy[high,0],xy[high,1],s=22,facecolors='none',edgecolors='#B33232',linewidths=1,label='PGAM5 high');ax.set_title(name+' | high / other detected sources',fontsize=10)
   if row==0:ax.legend(loc='upper right',fontsize=7,framealpha=.8)
   ax=axs[row,1];hx=xy[high];ax.scatter(hx[~cluster,0],hx[~cluster,1],s=18,color='#999999',label='Singleton high');ax.scatter(hx[cluster,0],hx[cluster,1],s=22,color=COL[1],label='High in adjacent component')
   if edges.any():ax.add_collection(LineCollection(np.stack([xy[e.row[edges]],xy[e.col[edges]]],axis=1),colors='#B33232',linewidths=1,alpha=.7))
   ax.set_title(name+' | direct high-high components',fontsize=10)
   if row==0:ax.legend(loc='upper right',fontsize=7,framealpha=.8)
   ax=axs[row,2];z=bins[bins['sample']==name].set_index('tile');rate={t:100*v.n_high/v.n_detected if v.n_detected else np.nan for t,v in z.iterrows()};v=np.array([rate[t] for t in o.loc[pos,'spatial_bin']]);ax.scatter(xy[pos,0],xy[pos,1],s=10,c=v,cmap='viridis',vmin=0,vmax=100);ax.scatter(xy[high,0],xy[high,1],s=19,facecolors='none',edgecolors='#222222',linewidths=.6)
   for t,b in z.iterrows():
    mask=o.spatial_bin.eq(t).to_numpy()
    if mask.any() and b.n_detected:
     center=xy[mask].mean(0);ax.text(center[0],center[1],f'{int(b.n_high)}/{int(b.n_detected)}',fontsize=7,ha='center',color='black',bbox=dict(facecolor='white',alpha=.7,edgecolor='none',pad=1))
   ax.set_title(name+' | bin high / detected Mac RNA (%)',fontsize=10)
  f.suptitle(patient+' | True Visium H&E coordinates; all points are real spots; 3x3 bins are not pathological ROIs',fontsize=12);f.tight_layout(rect=(0,0,.9,.98));cax=f.add_axes([.92,.22,.014,.55]);f.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,100),cmap='viridis'),cax=cax,label='Bin high / detected Mac RNA (%)');save(f,patient+'_high_distribution_spatial_maps')
if __name__=='__main__':run()
