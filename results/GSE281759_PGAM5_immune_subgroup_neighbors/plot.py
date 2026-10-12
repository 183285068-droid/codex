from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE281759_immune_mpl')
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import anndata as ad
from analyze import ROOT,REF,SAMPLES,TARGETS
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
def save(fig,name):
 fig.savefig(ROOT/'figures'/f'{name}.png',dpi=180,bbox_inches='tight');fig.savefig(ROOT/'figures'/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
def heat(data,title,name,cmap='Blues',vmin=0,vmax=100,suffix=''):
 fig,ax=plt.subplots(figsize=(12,7));im=ax.imshow(data.to_numpy(float),aspect='auto',cmap=cmap,vmin=vmin,vmax=vmax)
 ax.set_xticks(range(len(data.columns)),data.columns,rotation=50,ha='right');ax.set_yticks(range(len(data)),data.index)
 for i in range(len(data)):
  for j in range(len(data.columns)):
   v=data.iloc[i,j]
   if np.isfinite(v):ax.text(j,i,f'{v:.0f}{suffix}',ha='center',va='center',fontsize=8,color='white' if abs(v)>vmax*.65 else '#252525')
 fig.colorbar(im,ax=ax,shrink=.7);ax.set_title(title,pad=16);fig.tight_layout();save(fig,name)
def run():
 d=pd.read_csv(ROOT/'neighborhood_counts.csv');s=pd.read_csv(ROOT/'spatial_statistics.csv');p=pd.read_csv(ROOT/'patient_neighborhood_descriptive.csv');lr=pd.read_csv(ROOT/'candidate_LR_codetection.csv')
 h=d[d.source_group.eq('PGAM5_high')];matrix=h.pivot(index='sample',columns='target',values='neighbor_pct').reindex(index=[i['sample'] for i in SAMPLES],columns=TARGETS)
 heat(matrix,'PGAM5-high macrophage-RNA sources with >=1 adjacent immune-RNA target (%)\nOverlapping mixed-spot labels; denominator: connected high sources','section_neighborhood_coverage')
 q=s[s.comparison.eq('high_vs_other_detected')];mat=q.pivot(index='sample',columns='target',values='delta_pp').reindex(index=matrix.index,columns=TARGETS)
 heat(mat,'High minus other PGAM5-detected macrophage-RNA sources (percentage points)\nNA: insufficient groups; FDR results in spatial_statistics.csv','section_neighbor_differences','RdBu_r',-50,50)
 fig,ax=plt.subplots(figsize=(10,6));sub=p[p.source_group.eq('PGAM5_high')];x=np.arange(len(TARGETS));colors=['#285A8E','#B5822C','#A8536A','#607C44']
 for i,patient in enumerate(['A009','B013','C004','C008']):
  z=sub[sub.patient.eq(patient)].set_index('target').reindex(TARGETS);ax.barh(x+(i-1.5)*.19,z.neighbor_pct,height=.18,color=colors[i],label=f'{patient} (n={int(z.n_sources.iloc[0])} sources)')
 ax.set_yticks(x,TARGETS);ax.invert_yaxis();ax.set_xlim(0,100);ax.set_xlabel('High source spots with target neighbor (%)');ax.set_title('Patient-level descriptive spatial coverage');ax.legend(fontsize=9);ax.spines[['top','right']].set_visible(False);fig.tight_layout();save(fig,'patient_neighborhood_coverage')
 pooled=lr.groupby(['target','pair'])[['eligible_edges','supported_edges']].sum();pooled['pct']=100*pooled.supported_edges/pooled.eligible_edges.replace(0,np.nan)
 lm=pooled.pct.unstack().reindex(TARGETS)
 totals=lr.groupby('target').eligible_edges.sum()/lr.pair.nunique()
 lm.index=[f'{k} (n={int(totals[k])} edges)' for k in TARGETS]
 heat(lm,'Candidate ligand-source / receptor-target RNA co-detection (%)\nAdjacent PGAM5-high sources only; descriptive, no communication significance','candidate_LR_codetection')
 for patient in ['A009','B013','C004','C008']:
  items=[i for i in SAMPLES if i['patient']==patient];fig,axes=plt.subplots(len(items),4,figsize=(14,3.2*len(items)),squeeze=False)
  for r,item in enumerate(items):
   sample=item['sample'];a=ad.read_h5ad(REF/f'{sample}_spatial_marker_subset.h5ad');o=pd.read_csv(ROOT/f'{sample}_annotations.csv.gz',index_col=0).loc[a.obs_names];xy=a.obsm['spatial']*a.uns['scalefactors']['tissue_lowres_scalef'];image=a.uns['spatial'][sample]['images']['lowres'];high=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy()
   for c,keys in enumerate([['Exhausted_T','Cytotoxic_T','Activated_T'],['B_related','Plasma_related'],['cDC2_related','cDC1_related'],['Inflammatory_myeloid']]):
    ax=axes[r,c];ax.imshow(image);ax.set_title(sample+' | '+', '.join(keys),fontsize=8)
    for j,k in enumerate(keys):
     mask=o['gate_'+k].to_numpy(bool);ax.scatter(xy[mask,0],xy[mask,1],s=7,color=colors[j],alpha=.65,label=k)
    edge=a.obsp['spatial_connectivities'].tocoo();target=o[['gate_'+k for k in keys]].any(axis=1).to_numpy();use=high[edge.row]&target[edge.col]
    segments=np.stack([xy[edge.row[use]],xy[edge.col[use]]],axis=1)
    if len(segments):ax.add_collection(LineCollection(segments,colors='#292929',linewidths=.5,alpha=.4))
    ax.scatter(xy[high,0],xy[high,1],s=19,facecolors='none',edgecolors='#B33232',linewidths=.9,label='PGAM5-high Mac RNA');ax.axis('off')
    if r==0:ax.legend(fontsize=6,loc='upper right',framealpha=.85)
  fig.suptitle(patient+' | Actual Visium coordinates on H&E; lines = direct neighbor edges; RNA labels may overlap',fontsize=11);fig.tight_layout();save(fig,patient+'_immune_spatial_maps')
if __name__=='__main__':run()
