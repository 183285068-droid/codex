import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE281759_Tfocused_mpl');os.environ.setdefault('XDG_CACHE_HOME','/tmp/GSE281759_Tfocused_cache')
import numpy as np,pandas as pd,anndata as ad,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from analyze import R,BASE,K,SAMPLES,dense
NAMES=['CD8 exhaustion supported','CD8 early activation supported','CD8 memory related','CD8 effector related','Treg supported (separate T gate)'];COL=['#285A8E','#B5822C','#A8536A','#607C44'];P=['A009','B013','C004','C008']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
def save(f,n):
 for ext in ['png','pdf']:f.savefig(R/'figures'/f'{n}.{ext}',dpi=160,bbox_inches='tight')
 plt.close(f)
def heat(m,title,n,lo=0,hi=100,cmap='Blues'):
 f,ax=plt.subplots(figsize=(13,8));im=ax.imshow(m,aspect='auto',vmin=lo,vmax=hi,cmap=cmap);ax.set_xticks(range(len(m.columns)),m.columns,rotation=35,ha='right');ax.set_yticks(range(len(m)),m.index)
 for i in range(len(m)):
  for j in range(len(m.columns)):
   v=m.iloc[i,j];ax.text(j,i,f'{v:.1f}' if np.isfinite(v) else 'NA',ha='center',va='center',fontsize=9)
 ax.set_title(title,pad=15);f.colorbar(im,ax=ax,shrink=.65);f.tight_layout();save(f,n)
def run():
 p=pd.read_csv(R/'patient_descriptive.csv');s=pd.read_csv(R/'spatial_statistics.csv');lr=pd.read_csv(R/'candidate_LR_codetection.csv')
 f,ax=plt.subplots(figsize=(12,6));y=np.arange(5)
 for i,patient in enumerate(P):
  z=p[(p.patient==patient)&p.source_group.eq('PGAM5_high')].set_index('target').reindex(K);ax.barh(y+(i-1.5)*.18,z.neighbor_pct,height=.17,color=COL[i],label=f'{patient}: {int(z.n_sources.iloc[0])} high sources')
 ax.set_yticks(y,NAMES);ax.invert_yaxis();ax.set_xlim(0,20);ax.set_xlabel('High Mac-RNA sources with >=1 direct target neighbor (%)');ax.legend(loc='lower right',fontsize=9);ax.set_title('Focused T-RNA neighborhood coverage by patient\nFour overlapping CD8 programs; Treg uses a separate T-lineage gate');f.tight_layout();save(f,'patient_focused_T_coverage')
 main=s[s.comparison.eq('high_vs_other_detected')];m=main.pivot(index='sample',columns='target',values='delta_pp').reindex(index=[v['sample'] for v in SAMPLES],columns=K);m.columns=NAMES;heat(m,'High minus other PGAM5-detected sources (percentage points)\nNA: not testable; 80 planned primary tests','section_focused_T_differences',-20,20,'RdBu_r')
 pooled=lr.groupby(['target','pair'])[['eligible_edges','supported_edges']].sum();pooled['pct']=100*pooled.supported_edges/pooled.eligible_edges.replace(0,np.nan);pooled.to_csv(R/'pooled_LR_descriptive.csv');m=pooled.pct.unstack().reindex(K);m.index=NAMES;heat(m,'Candidate ligand-source / receptor-target RNA co-detection (%)\nActual adjacent edges; denominators differ; no communication significance','focused_T_LR_codetection')
 gs=['CD8A','CD8B','PDCD1','TOX','TCF7','SLAMF6','HAVCR2','ENTPD1','CD69','IL2RA','TNFRSF9','ICOS','CCR7','IL7R','GZMK','ITGAE','ZNF683','PRF1','GZMB','FGFBP2','CX3CR1','FOXP3','CTLA4','IKZF2','TIGIT'];records={k:[] for k in K};ov=np.zeros((5,5),int);counts=np.zeros(5,int)
 for item in SAMPLES:
  a=ad.read_h5ad(R/f'{item["sample"]}_focused_T_markers.h5ad');g=a.obs[['gate_'+k for k in K]].to_numpy(int);ov+=g.T@g;counts+=g.sum(0);ix=a.var_names.get_indexer(gs)
  for k in K:
   mask=a.obs['gate_'+k].to_numpy(bool)&a.obs.near_high_source.to_numpy(bool);records[k].append((dense(a.layers['counts'])[mask][:,ix],dense(a.X)[mask][:,ix]))
 f,ax=plt.subplots(figsize=(16,6));dots=[]
 for row,k in enumerate(K):
  raw=np.concatenate([v for v,n in records[k]]);x=np.concatenate([n for v,n in records[k]])
  for j,g in enumerate(gs):
   if len(raw):
    pct=100*(raw[:,j]>0).mean();mean=x[:,j].mean();ax.scatter(j,row,s=pct*1.7,c=mean,cmap='viridis',vmin=0,vmax=2,edgecolors='#444444',linewidths=.2);dots.append(dict(target=k,gene=g,n_unique_high_near_targets=len(raw),detection_pct=pct,mean_log1p_CP10k=mean))
 ax.set_yticks(range(5),[NAMES[j]+f' (n={sum(len(v) for v,n in records[k])})' for j,k in enumerate(K)]);ax.set_xticks(range(len(gs)),gs,rotation=60,ha='right');ax.invert_yaxis();ax.set_title('Defining RNA markers in unique targets adjacent to high sources\nDot area: detection 0-100%; color: mean log1p(CP10k), capped at 2; n=0 rows empty');f.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,2),cmap='viridis'),ax=ax,shrink=.6);f.tight_layout();save(f,'focused_T_marker_evidence');pd.DataFrame(dots).to_csv(R/'marker_dotplot_values.csv',index=False)
 pd.DataFrame(ov,index=K,columns=K).to_csv(R/'focused_program_overlap.csv');m=pd.DataFrame(100*ov/counts[:,None],index=NAMES,columns=NAMES);heat(m,'Fraction of row-category targets also meeting column gate (%)\nOverlapping mixed RNA; Treg-CD8 overlap does not identify CD8 regulatory cells','focused_T_overlap')
 c=pd.read_csv(R/'CD8_four_program_signatures.csv').groupby('signature')[['n_CD8_targets','n_high_near_targets']].sum();c.to_csv(R/'CD8_four_program_signature_pooled.csv');f,ax=plt.subplots(figsize=(13,7));labels=[v.replace('CD8_','').replace('+',' + ') for v in c.index];ax.barh(labels,c.n_high_near_targets,color=COL[0]);ax.invert_yaxis();ax.set_xlabel('Unique high-near CD8 targets');ax.set_title('Observed combinations of the four CD8 programs\nNo-four-program label is not a naive/resting identity');f.tight_layout();save(f,'CD8_four_program_combinations')
 for patient in P:
  items=[v for v in SAMPLES if v['patient']==patient];f,axs=plt.subplots(len(items),5,figsize=(21,3.2*len(items)),squeeze=False)
  for row,item in enumerate(items):
   name=item['sample'];a=ad.read_h5ad(BASE/f'{name}_spatial_marker_subset.h5ad');o=pd.read_csv(R/f'{name}_annotations.csv.gz',index_col=0).loc[a.obs_names];xy=a.obsm['spatial']*a.uns['scalefactors']['tissue_lowres_scalef'];high=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy();e=a.obsp['spatial_connectivities'].tocoo()
   for j,k in enumerate(K):
    ax=axs[row,j];target=o['gate_'+k].to_numpy(bool);ax.imshow(a.uns['spatial'][name]['images']['lowres']);ax.scatter(xy[target,0],xy[target,1],s=11,color=COL[0],alpha=.8);use=high[e.row]&target[e.col]
    if use.any():ax.add_collection(LineCollection(np.stack([xy[e.row[use]],xy[e.col[use]]],axis=1),colors='#333333',linewidths=.7,alpha=.6))
    ax.scatter(xy[high,0],xy[high,1],s=20,facecolors='none',edgecolors='#B33232',linewidths=.8);ax.set_title(name+' | '+k.replace('CD8_','CD8 ').replace('_',' '),fontsize=9);ax.axis('off')
  f.suptitle(patient+' | True Visium H&E coordinates; red rings: high Mac RNA; blue: targets; lines: direct neighbors',fontsize=12);f.tight_layout();save(f,patient+'_focused_T_spatial_maps')
if __name__=='__main__':run()
