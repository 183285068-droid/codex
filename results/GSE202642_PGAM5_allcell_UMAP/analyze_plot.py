from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;src=R/'inputs' if (R/'inputs').exists() else Path('/workspace/codex/results/GSE202642_HCC_allcell_annotation')
def main():
 a=pd.read_csv(src/'cell_annotations.csv.gz',index_col=0,dtype={'cluster':str});u=pd.read_csv(src/'UMAP_coordinates.csv.gz',index_col=0);d=pd.read_csv(src/'cluster_annotations.csv',dtype={'cluster':str});assert a.index.equals(u.index) and len(a)==56719
 genes=pd.read_csv(src/'marker_genes.csv').gene.tolist();raw=sparse.load_npz(src/'marker_raw_counts.npz');pg=raw[:,genes.index('PGAM5')].toarray().ravel();assert np.array_equal(pg,a.PGAM5_counts);assert a.total_counts.gt(0).all();a['PGAM5_CP10k']=pg/a.total_counts.to_numpy()*10000;a['PGAM5_log1p_CP10k']=np.log1p(a.PGAM5_CP10k);a['PGAM5_detected']=pg>0
 rows=[]
 for c in sorted(a.cluster.unique(),key=int):
  v=a[a.cluster.eq(c)];rows.append(dict(cluster=c,cells=len(v),PGAM5_detected_cells=int(v.PGAM5_detected.sum()),PGAM5_detection_percent=float(v.PGAM5_detected.mean()*100),mean_PGAM5_CP10k_all_cells=float(v.PGAM5_CP10k.mean()),median_PGAM5_CP10k_all_cells=float(v.PGAM5_CP10k.median()),mean_PGAM5_log1p_CP10k_all_cells=float(v.PGAM5_log1p_CP10k.mean()),mean_PGAM5_CP10k_detected_only=float(v.loc[v.PGAM5_detected,'PGAM5_CP10k'].mean()),mean_raw_PGAM5_counts=float(v.PGAM5_counts.mean()),median_total_counts=float(v.total_counts.median()),patients_with_cells=int(v.sample_name.nunique()),largest_patient_fraction=float(v.sample_name.value_counts(normalize=True).max())))
 s=pd.DataFrame(rows).merge(d[['cluster','subtype_EN','subtype_CN','cell_type_EN','boundary_or_mixed_flag','annotation_confidence']],on='cluster',validate='one_to_one');s['mean_expression_rank']=s.mean_PGAM5_CP10k_all_cells.rank(ascending=False,method='min').astype(int);s['detection_rank']=s.PGAM5_detection_percent.rank(ascending=False,method='min').astype(int);s.to_csv(R/'PGAM5_by_cluster.csv',index=False)
 a[['sample_name','cluster','subtype_EN','cell_type_EN','PGAM5_counts','total_counts','PGAM5_CP10k','PGAM5_log1p_CP10k','PGAM5_detected']].join(u).to_csv(R/'PGAM5_cell_expression_and_UMAP.csv.gz');p=a.groupby(['cluster','sample_name'],observed=True).agg(cells=('PGAM5_counts','size'),detected=('PGAM5_detected','sum'),mean_CP10k=('PGAM5_CP10k','mean'),mean_log1p_CP10k=('PGAM5_log1p_CP10k','mean'),median_total_counts=('total_counts','median')).reset_index();p['detection_percent']=p.detected/p.cells*100;p.to_csv(R/'PGAM5_by_cluster_patient.csv',index=False)
 base=s.set_index('cluster');order=np.random.default_rng(202642).permutation(len(a));coords=u[['UMAP1','UMAP2']].to_numpy();pos=np.flatnonzero(pg>0);pos=pos[np.argsort(a.PGAM5_log1p_CP10k.to_numpy()[pos])]
 def labels(ax):
  for c in sorted(a.cluster.unique(),key=int):
   xy=np.median(coords[a.cluster.eq(c).to_numpy()],axis=0);ax.text(*xy,'C'+c,ha='center',va='center',fontsize=8,bbox={'facecolor':'white','alpha':.8,'edgecolor':'none','pad':1})
 def axes(ax):ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2')
 def save(fig,fn):fig.tight_layout();fig.savefig(R/(fn+'.png'),dpi=180);fig.savefig(R/(fn+'.pdf'));plt.close(fig)
 fig,axs=plt.subplots(1,2,figsize=(17,7));axs[0].scatter(coords[order,0],coords[order,1],s=1.5,c='#d5d5d5',linewidths=0,rasterized=True);pt=axs[0].scatter(coords[pos,0],coords[pos,1],s=3,c=a.PGAM5_log1p_CP10k.to_numpy()[pos],cmap='viridis',vmin=0,vmax=a.PGAM5_log1p_CP10k.max(),linewidths=0,rasterized=True);fig.colorbar(pt,ax=axs[0],shrink=.7,label='PGAM5 log1p(CP10k)');axs[0].set_title('Per-cell normalized PGAM5 RNA\nGray: raw PGAM5 count = 0; full color range');axs[1].scatter(coords[order,0],coords[order,1],s=1.5,c='#d5d5d5',linewidths=0,rasterized=True);axs[1].scatter(coords[pos,0],coords[pos,1],s=3,c='#c6363d',linewidths=0,rasterized=True);axs[1].set_title(f'PGAM5 RNA detected: {len(pos):,}/{len(a):,} cells\nRed: raw count > 0; gray: undetected')
 for ax in axs:labels(ax);axes(ax)
 save(fig,'PGAM5_expression_UMAP')
 fig,axs=plt.subplots(1,2,figsize=(17,7))
 for ax,field,title,barlabel in [(axs[0],'mean_PGAM5_CP10k_all_cells','Cluster mean PGAM5 expression (zeros included)','Mean PGAM5 CP10k'),(axs[1],'PGAM5_detection_percent','Cluster PGAM5 RNA detection rate','Detected cells (%)')]:
  values=a.cluster.map(base[field]).to_numpy();pt=ax.scatter(coords[order,0],coords[order,1],s=1.5,c=values[order],cmap='viridis',vmin=0,vmax=float(s[field].max()),linewidths=0,rasterized=True);fig.colorbar(pt,ax=ax,shrink=.7,label=barlabel);ax.set_title(title);labels(ax);axes(ax)
 save(fig,'PGAM5_cluster_mean_and_detection_UMAP')
 ranked=s.sort_values('mean_PGAM5_CP10k_all_cells',ascending=False);ys=np.arange(len(ranked));labs=['C'+r.cluster+' '+r.subtype_EN+(' [qualified/mixed]' if r.boundary_or_mixed_flag else '') for r in ranked.itertuples()];fig,axs=plt.subplots(1,2,figsize=(17,10),sharey=True);axs[0].barh(ys,ranked.mean_PGAM5_CP10k_all_cells,color='#3977aa');axs[1].barh(ys,ranked.PGAM5_detection_percent,color='#b97028');axs[0].set_yticks(ys,labs,fontsize=9);axs[0].invert_yaxis();axs[0].set_xlabel('Mean PGAM5 CP10k (all cells, including zeros)');axs[1].set_xlabel('PGAM5 RNA detected cells (%)');axs[0].set_title('Normalized expression');axs[1].set_title('Detection rate')
 for i,r in enumerate(ranked.itertuples()):axs[0].text(r.mean_PGAM5_CP10k_all_cells+.002,i,f'{r.mean_PGAM5_CP10k_all_cells:.3f}',va='center',fontsize=8);axs[1].text(r.PGAM5_detection_percent+.2,i,f'{r.PGAM5_detection_percent:.2f}% ({r.PGAM5_detected_cells}/{r.cells})',va='center',fontsize=8)
 axs[0].set_xlim(0,s.mean_PGAM5_CP10k_all_cells.max()*1.22);axs[1].set_xlim(0,s.PGAM5_detection_percent.max()*1.4);fig.suptitle('GSE202642 | C0-C20 | same cluster definitions and original UMAP coordinates');save(fig,'PGAM5_cluster_comparison')
 hashes={f:hashlib.sha256((src/f).read_bytes()).hexdigest() for f in ['cell_annotations.csv.gz','cluster_annotations.csv','UMAP_coordinates.csv.gz','marker_genes.csv','marker_raw_counts.npz']};(R/'validation.json').write_text(json.dumps(dict(cells=len(a),clusters=len(s),PGAM5_detected_cells=len(pos),source_PGAM5_counts_exact_match=True,UMAP_coordinates_unchanged=True,zeros_included_in_means=True,normalization_denominator='Original full-library counts, not marker subset',color_clipping=False,source_file_hashes=hashes),indent=2));print(s.sort_values('mean_expression_rank').head()[['cluster','mean_PGAM5_CP10k_all_cells','PGAM5_detection_percent']].to_string(index=False),flush=True)
if __name__=='__main__':main()
