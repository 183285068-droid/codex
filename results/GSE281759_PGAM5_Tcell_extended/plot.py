from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE281759_extended_mpl');os.environ.setdefault('XDG_CACHE_HOME','/tmp/GSE281759_extended_cache')
import numpy as np,pandas as pd,anndata as ad,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from analyze import R,BASE,K,SAMPLES
P=['A009','B013','C004','C008'];COL=['#285A8E','#B5822C','#A8536A','#607C44'];LABEL=['CD4 related','CD8 supported','Treg supported','Naive like','Central memory like','Effector memory like']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
def save(f,n):
 for ext in ['png','pdf']:f.savefig(R/'figures'/f'{n}.{ext}',dpi=160,bbox_inches='tight')
 plt.close(f)
def heat(m,title,n,cmap='Blues',lo=0,hi=100):
 f,ax=plt.subplots(figsize=(12,8));im=ax.imshow(m.to_numpy(float),aspect='auto',cmap=cmap,vmin=lo,vmax=hi);ax.set_xticks(range(len(m.columns)),m.columns,rotation=35,ha='right');ax.set_yticks(range(len(m)),m.index)
 for i in range(len(m)):
  for j in range(len(m.columns)):
   v=m.iloc[i,j];ax.text(j,i,f'{v:.1f}' if np.isfinite(v) else 'NA',ha='center',va='center',fontsize=9)
 ax.set_title(title,pad=15);f.colorbar(im,ax=ax,shrink=.6);f.tight_layout();save(f,n)
def run():
 d=pd.read_csv(R/'neighborhood_counts.csv');p=pd.read_csv(R/'patient_descriptive.csv');s=pd.read_csv(R/'spatial_statistics.csv');lr=pd.read_csv(R/'candidate_LR_codetection.csv')
 f,ax=plt.subplots(figsize=(12,6));x=np.arange(len(K))
 for i,patient in enumerate(P):
  z=p[(p.patient==patient)&p.source_group.eq('PGAM5_high')].set_index('target').reindex(K);ax.barh(x+(i-1.5)*.18,z.neighbor_pct,height=.17,color=COL[i],label=f'{patient} (n={int(z.n_sources.iloc[0])} high sources)')
 ax.set_yticks(x,LABEL);ax.invert_yaxis();ax.set_xlim(0,100);ax.set_xlabel('High source spots with >=1 adjacent target (%)');ax.legend(loc='lower right',fontsize=9);ax.set_title('Extended T-RNA neighborhood coverage by patient\nOverlapping mixed-spot programs; no resolved cell identity');f.tight_layout();save(f,'patient_extended_T_coverage')
 main=s[s.comparison.eq('high_vs_other_detected')];m=main.pivot(index='sample',columns='target',values='delta_pp').reindex(index=[i['sample'] for i in SAMPLES],columns=K);m.columns=LABEL
 heat(m,'High minus other PGAM5-detected source coverage (percentage points)\nNA: not testable; FDR from 96 preplanned primary tests','section_extended_T_differences','RdBu_r',-50,50)
 q=main.pivot(index='sample',columns='target',values='primary_q_96').reindex(index=m.index,columns=K);q.columns=LABEL
 heat(q,'BH q values: 96 planned primary tests\nNot-testable p included as 1 only for family correction','section_extended_T_FDR','Blues',0,1)
 pool=lr.groupby(['target','pair'])[['eligible_edges','supported_edges']].sum();pool['pct']=100*pool.supported_edges/pool.eligible_edges.replace(0,np.nan);pool.to_csv(R/'pooled_LR_descriptive.csv');heat(pool.pct.unstack().reindex(K),'Ligand-source / receptor-target RNA co-detection (%)\nAdjacent high source edges; descriptive, not communication significance','extended_T_LR_codetection')
 ov=[]
 for item in SAMPLES:
  o=pd.read_csv(R/f'{item["sample"]}_annotations.csv.gz',index_col=0)
  for k in K:
   for k2 in K:ov.append(dict(patient=item['patient'],sample=item['sample'],first=k,second=k2,both=int((o['gate_'+k]&o['gate_'+k2]).sum())))
 ov=pd.DataFrame(ov);ov.to_csv(R/'target_overlap_counts.csv',index=False);heat(ov.groupby(['first','second']).both.sum().unstack().reindex(index=K,columns=K),'Overlapping T-program target counts (diagonal = program total)\nCD4/CD8 and naive/memory labels are not exclusive','extended_T_target_overlap','Blues',0,float(ov.groupby(['first','second']).both.sum().max()))
 for patient in P:
  items=[i for i in SAMPLES if i['patient']==patient];f,axs=plt.subplots(len(items),6,figsize=(22,3*len(items)),squeeze=False)
  for row,item in enumerate(items):
   sample=item['sample'];a=ad.read_h5ad(BASE/f'{sample}_spatial_marker_subset.h5ad');o=pd.read_csv(R/f'{sample}_annotations.csv.gz',index_col=0).loc[a.obs_names];high=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy();xy=a.obsm['spatial']*a.uns['scalefactors']['tissue_lowres_scalef'];image=a.uns['spatial'][sample]['images']['lowres'];e=a.obsp['spatial_connectivities'].tocoo()
   for j,k in enumerate(K):
    target=o['gate_'+k].to_numpy(bool);use=high[e.row]&target[e.col];ax=axs[row,j];ax.imshow(image);ax.scatter(xy[target,0],xy[target,1],s=9,color=COL[0],alpha=.7)
    if use.any():ax.add_collection(LineCollection(np.stack([xy[e.row[use]],xy[e.col[use]]],axis=1),colors='#333333',linewidths=.6,alpha=.5))
    ax.scatter(xy[high,0],xy[high,1],s=20,facecolors='none',edgecolors='#B33232',linewidths=.8);ax.set_title(sample+' | '+LABEL[j],fontsize=10);ax.axis('off')
  f.suptitle(patient+' | True Visium H&E coordinates; red rings: high Mac RNA; blue: T programs; lines: direct neighbors',fontsize=13);f.tight_layout();save(f,patient+'_extended_T_spatial_maps')
if __name__=='__main__':run()
