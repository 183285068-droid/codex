import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE283406_T_mpl');os.environ.setdefault('XDG_CACHE_HOME','/tmp/GSE283406_T_cache')
import numpy as np,pandas as pd,anndata as ad,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze import R,SOURCE,K,PANELS,PAIRS,dense,input_path
NAMES=['T lineage reference','Exhaustion supported','Early activation supported','Cytotoxic supported'];COL=['#285A8E','#B5822C','#999999']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
def save(f,n):
 for ext in ['png','pdf']:f.savefig(R/'figures'/f'{n}.{ext}',dpi=180,bbox_inches='tight')
 plt.close(f)
def run():
 d=pd.read_csv(R/'same_barcode_T_program_overlap.csv');c=pd.read_csv(R/'conditional_T_program_descriptive.csv');lr=pd.read_csv(R/'ligand_receptor_independent_expression.csv');a=ad.read_h5ad(input_path());o=pd.read_csv(R/'J259X_T_annotations.csv.gz',index_col=0).loc[a.obs_names]
 f,axs=plt.subplots(1,2,figsize=(14,6));y=np.arange(4)
 for i,group in enumerate(['PGAM5_high','Other_detected','Undetected']):
  z=d[d.source_group.eq(group)].set_index('target').reindex(K);axs[0].barh(y+(i-1)*.22,z.same_barcode_pct,height=.2,color=COL[i],label=f'{group} (n={int(z.n_sources.iloc[0])})')
  zc=c[c.source_group.eq(group)].set_index('target').reindex(K);axs[1].barh(y+(i-1)*.22,zc.state_pct,height=.2,color=COL[i],label=f'{group} / T coexpressing (n={int(zc.n_T_coexpressing_sources.iloc[0])})')
 for ax in axs:ax.set_yticks(y,NAMES);ax.invert_yaxis();ax.set_xlim(0,100);ax.legend(loc='lower right',fontsize=8);ax.spines[['top','right']].set_visible(False)
 axs[0].set_title('Same-barcode RNA overlap / all Mac-RNA sources');axs[1].set_title('Program support / T-coexpressing Mac-RNA sources');axs[0].set_xlabel('Source spots (%)');axs[1].set_xlabel('T-coexpressing source spots (%)');f.suptitle('J259X tumor: mixed RNA cooccurrence only; no cell pairing, spatial adjacency or p values');f.tight_layout();save(f,'T_same_barcode_coexpression')
 signatures=pd.read_csv(R/'T_program_signatures.csv');f,axs=plt.subplots(1,2,figsize=(15,7),sharey=True);labels=[s.replace('_T','').replace('_',' ').replace('+',' + ') for s in signatures.signature];y=np.arange(len(labels));axs[0].barh(y,signatures.n_T_targets,color=COL[0]);axs[1].barh(y,signatures.high_same_barcode,color=COL[1]);axs[0].set_yticks(y,labels);axs[0].invert_yaxis();axs[0].set_title('All T-lineage RNA targets (n=575)');axs[1].set_title('High Mac-RNA / T SAME-barcode overlap (n=7)');axs[0].set_xlabel('Spots');axs[1].set_xlabel('Spots');f.suptitle('T-state combinations: mutually exclusive descriptive signatures\nNo-supported-state does not establish naive/memory identity');f.tight_layout();save(f,'T_state_combinations')
 f,axs=plt.subplots(1,2,figsize=(14,6),gridspec_kw={'width_ratios':[1,2]});source=lr[lr.source_group.eq('PGAM5_high')].drop_duplicates('pair').set_index('pair');pairs=[l+'--'+r for l,r in PAIRS];src=source.reindex(pairs).ligand_pct.to_numpy()[:,None];target=source.index
 im=axs[0].imshow(src,vmin=0,vmax=100,cmap='Blues',aspect='auto');axs[0].set_yticks(range(5),pairs);axs[0].set_xticks([0],['High-source ligand (n=21)']);axs[0].set_title('Source expression')
 m=lr[lr.source_group.eq('PGAM5_high')].pivot(index='pair',columns='target',values='receptor_pct').reindex(index=pairs,columns=K);axs[1].imshow(m,vmin=0,vmax=100,cmap='Blues',aspect='auto');axs[1].set_yticks(range(5),pairs);axs[1].set_xticks(range(4),[NAMES[j]+f' (n={int(o["gate_"+k].sum())})' for j,k in enumerate(K)],rotation=35,ha='right');axs[1].set_title('Receptor detection within overlapping T RNA gates')
 for i in range(5):
  axs[0].text(0,i,f'{src[i,0]:.1f}',ha='center',va='center')
  for j in range(4):axs[1].text(j,i,f'{m.iloc[i,j]:.1f}',ha='center',va='center',fontsize=9)
 f.suptitle('Separate ligand and receptor detection (%)\nNo source-target pairing, interaction score or communication significance');f.tight_layout();save(f,'T_ligand_receptor_independent_expression')
 xy=a.obsm['X_umap'];high=o.PGAM5_mac_RNA_group.eq('PGAM5_high').to_numpy();f,axs=plt.subplots(1,3,figsize=(15,5))
 for j,k in enumerate(K[1:]):
  mask=o['gate_'+k].to_numpy(bool);ax=axs[j];ax.scatter(xy[:,0],xy[:,1],s=3,color='#DDDDDD',alpha=.6);ax.scatter(xy[mask,0],xy[mask,1],s=10,color=COL[0],label=f'{NAMES[j+1]} (n={int(mask.sum())})');ax.scatter(xy[high,0],xy[high,1],s=24,facecolors='none',edgecolors='#B33232',linewidths=.9,label='PGAM5-high Mac RNA (n=21)');ax.set_xlabel('Expression UMAP 1');ax.set_ylabel('Expression UMAP 2');ax.legend(fontsize=7);ax.set_title(NAMES[j+1])
 f.suptitle('J259X expression UMAP: NOT tissue coordinates or spatial neighbors; no interaction edges');f.tight_layout();save(f,'T_expression_UMAP_not_spatial')
 gs=list(dict.fromkeys(g for panel in PANELS.values() for g in panel));f,ax=plt.subplots(figsize=(16,5));raw=dense(a.layers['counts']);x=dense(a.X)
 for row,k in enumerate(K):
  mask=o['gate_'+k].to_numpy(bool)
  for j,g in enumerate(gs):
   ix=a.var_names.get_loc(g);pct=100*(raw[mask,ix]>0).mean();mean=x[mask,ix].mean();ax.scatter(j,row,s=pct*1.6,c=mean,cmap='viridis',vmin=0,vmax=2,edgecolors='#444444',linewidths=.2)
 ax.set_xticks(range(len(gs)),gs,rotation=60,ha='right');ax.set_yticks(range(4),[NAMES[i]+f' (n={int(o["gate_"+k].sum())})' for i,k in enumerate(K)]);ax.invert_yaxis();ax.set_title('Defining T-program RNA markers; not independent cell identity validation\nDot area: detection 0-100%; color: mean log1p(CP10k), capped at 2');f.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,2),cmap='viridis'),ax=ax,shrink=.6);f.tight_layout();save(f,'T_marker_evidence')
if __name__=='__main__':run()
