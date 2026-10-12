import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE283406_candidates_mpl')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/GSE283406_candidates_cache')
import numpy as np,pandas as pd,anndata as ad,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze import ROOT,PANELS,PAIRS
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
def save(fig,name):
 fig.savefig(ROOT/'figures'/f'{name}.png',dpi=180,bbox_inches='tight');fig.savefig(ROOT/'figures'/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
def run():
 d=pd.read_csv(ROOT/'same_spot_RNA_overlap.csv');lr=pd.read_csv(ROOT/'candidate_LR_independent_expression.csv');marker=pd.read_csv(ROOT/'target_marker_evidence.csv');order=list(PANELS);fig,axes=plt.subplots(1,2,figsize=(13,7),sharey=True)
 for ax,section in zip(axes,['J259X','J259B']):
  for j,(group,color) in enumerate([('PGAM5_high','#285A8E'),('Other_detected','#B5822C')]):
   v=d[d.section.eq(section)&d.source_group.eq(group)].set_index('target').reindex(order);ax.barh(np.arange(len(order))+(j-.5)*.3,v.same_spot_overlap_pct,height=.28,color=color,label=group+f' (n={int(v.n_sources.iloc[0])})')
   for y,row in enumerate(v.itertuples()):ax.text(row.same_spot_overlap_pct+1,y+(j-.5)*.3,f'{int(row.same_spot_overlap)}/{int(row.n_sources)}',va='center',fontsize=8)
  ax.set_yticks(range(len(order)),order);ax.set_xlim(0,105);ax.set_xlabel('Source spots also meeting immune RNA gate (%)');ax.set_title(section+(' (tumor primary)' if section=='J259X' else ' (mixed tissue context only)'));ax.legend(fontsize=8);ax.spines[['top','right']].set_visible(False)
 axes[0].invert_yaxis();fig.suptitle('Same-spot RNA co-occurrence | NOT neighboring cells or interaction probability');fig.tight_layout();save(fig,'same_spot_immune_RNA_overlap')
 primary=lr[lr.section.eq('J259X')];fig,axes=plt.subplots(1,2,figsize=(13,6),gridspec_kw={'width_ratios':[1,2.8]});pairs=[l+'--'+r for l,r in PAIRS];src=primary.drop_duplicates('pair').set_index('pair').reindex(pairs)
 axes[0].barh(range(len(pairs)),src.high_ligand_detect_pct,color='#285A8E');axes[0].set_yticks(range(len(pairs)),pairs);axes[0].invert_yaxis();axes[0].set_xlim(0,100);axes[0].set_xlabel('Ligand detection in 21 high sources (%)')
 for j,row in enumerate(src.itertuples()):axes[0].text(row.high_ligand_detect_pct+2,j,f'{int(row.high_ligand_detected)}/21',va='center',fontsize=8)
 mat=primary.pivot(index='target',columns='pair',values='target_receptor_detect_pct').reindex(index=order,columns=pairs);cm=plt.get_cmap('Blues').copy();cm.set_bad('#d9dde0');im=axes[1].imshow(mat,cmap=cm,vmin=0,vmax=100,aspect='auto');axes[1].set_xticks(range(len(pairs)),pairs,rotation=50,ha='right');nt=primary.drop_duplicates('target').set_index('target').n_target_spots;axes[1].set_yticks(range(len(order)),[f'{k} (n={nt[k]})' for k in order]);axes[1].set_title('Receptor detection within overlapping RNA target gates');fig.colorbar(im,ax=axes[1],label='Target RNA spots with receptor (%)');fig.suptitle('J259X candidate ligand/receptor expression | No spatial pairing, score or significance');fig.tight_layout();save(fig,'candidate_LR_independent_expression')
 for section in ['J259X','J259B']:
  a=ad.read_h5ad(ROOT/'inputs'/f'{section}_AVITI_expression_marker_subset.h5ad');o=pd.read_csv(ROOT/f'{section}_immune_candidate_annotations.csv.gz',index_col=0).loc[a.obs_names];xy=a.obsm['X_umap'];fig,axes=plt.subplots(3,3,figsize=(12,11));high=o.PGAM5_mac_RNA_group.eq('PGAM5_high').to_numpy()
  for ax,k in zip(axes.flat,order):
   ax.scatter(xy[:,0],xy[:,1],s=2,color='#d6d9de',alpha=.6,rasterized=True);gate=o['gate_'+k].to_numpy(bool);ax.scatter(xy[gate,0],xy[gate,1],s=5,color='#285A8E',alpha=.7,rasterized=True);ax.scatter(xy[high,0],xy[high,1],s=28,facecolors='none',edgecolors='#B5822C',linewidths=1.0);ax.set_title(k+f' | {int(gate.sum())} RNA targets',fontsize=9);ax.set_xlabel('Expression UMAP 1');ax.set_ylabel('Expression UMAP 2')
  fig.suptitle(section+' | Expression embedding only; UMAP proximity is NOT tissue proximity\nBlue = immune RNA gate; gold rings = PGAM5-high Mac-RNA sources');fig.tight_layout();save(fig,section+'_immune_expression_UMAP')
 genes=['CD3D','TRAC','PDCD1','LAG3','CD69','ICOS','IFNG','PRF1','MS4A1','CD79A','MZB1','JCHAIN','CD1C','FCER1A','CLEC9A','XCR1','S100A8','VCAN'];z=marker[marker.section.eq('J259X')];fig,ax=plt.subplots(figsize=(12,6));means=z.pivot(index='target',columns='gene',values='mean_expression').reindex(index=order,columns=genes);pct=z.pivot(index='target',columns='gene',values='detected_pct').reindex(index=order,columns=genes);xx,yy=np.meshgrid(range(len(genes)),range(len(order)));sc=ax.scatter(xx.ravel(),yy.ravel(),s=np.nan_to_num(pct.to_numpy()).ravel(),c=means.to_numpy().ravel(),cmap='viridis');ax.set_xticks(range(len(genes)),genes,rotation=50,ha='right');ax.set_yticks(range(len(order)),order);ax.invert_yaxis();fig.colorbar(sc,ax=ax,label='Mean log1p(CP10k), zeros retained');ax.set_title('J259X gate-marker documentation (only each gate defining genes shown)\nDot area = RNA detection%; gate markers are not independent identity validation');fig.tight_layout();save(fig,'J259X_marker_gate_evidence')
if __name__=='__main__':run()
