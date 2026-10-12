import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE281759_fullimmune_mpl');os.environ.setdefault('XDG_CACHE_HOME','/tmp/GSE281759_fullimmune_cache')
import numpy as np,pandas as pd,anndata as ad,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from analyze import R,BASE,K,SAMPLES,dense
P=['A009','B013','C004','C008'];COL=['#285A8E','#B5822C','#A8536A','#607C44']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
def label(k):return k.replace('_',' ')
def save(f,n):
 for ext in ['png','pdf']:f.savefig(R/'figures'/f'{n}.{ext}',dpi=160,bbox_inches='tight')
 plt.close(f)
def heat(m,title,name,cmap='Blues',lo=0,hi=100):
 f,ax=plt.subplots(figsize=(max(14,len(m.columns)*.7),max(9,len(m)*.4)));im=ax.imshow(m.to_numpy(float),aspect='auto',cmap=cmap,vmin=lo,vmax=hi);ax.set_xticks(range(len(m.columns)),m.columns,rotation=65,ha='right');ax.set_yticks(range(len(m)),m.index)
 for i in range(len(m)):
  for j in range(len(m.columns)):
   v=m.iloc[i,j];ax.text(j,i,f'{v:.1f}' if np.isfinite(v) else 'NA',ha='center',va='center',fontsize=7)
 ax.set_title(title,pad=18);f.colorbar(im,ax=ax,shrink=.65);f.tight_layout();save(f,name)
def run():
 p=pd.read_csv(R/'patient_descriptive.csv');s=pd.read_csv(R/'spatial_statistics.csv');lr=pd.read_csv(R/'candidate_LR_codetection.csv');zones=pd.read_csv(R/'unique_target_neighborhood_zones.csv');raw_objects=[]
 f,ax=plt.subplots(figsize=(14,18));y=np.arange(len(K))
 for i,patient in enumerate(P):
  z=p[(p.patient==patient)&p.source_group.eq('PGAM5_high')].set_index('target').reindex(K);ax.barh(y+(i-1.5)*.18,z.neighbor_pct,height=.17,color=COL[i],label=f'{patient}: {int(z.n_sources.iloc[0])} high source spots')
 ax.set_yticks(y,[label(k) for k in K]);ax.invert_yaxis();ax.set_xlim(0,100);ax.set_xlabel('High Mac-RNA sources with >=1 direct target neighbor (%)');ax.legend(loc='lower right',fontsize=10);ax.set_title('Full-gene immune RNA neighborhood coverage by patient\nOverlapping mixed-spot programs; spatial candidate evidence, not resolved cell interaction');f.tight_layout();save(f,'patient_immune_coverage')
 main=s[s.comparison.eq('high_vs_other_detected')];m=main.pivot(index='sample',columns='target',values='delta_pp').reindex(index=[i['sample'] for i in SAMPLES],columns=K);m.columns=[label(k) for k in K];heat(m,'High minus other PGAM5-detected sources (percentage points)\nNA: insufficient groups or invariant conditional null; FDR in spatial_statistics.csv','section_immune_differences','RdBu_r',-50,50)
 pool=lr.groupby(['target','pair'])[['eligible_edges','supported_edges']].sum();pool['pct']=100*pool.supported_edges/pool.eligible_edges.replace(0,np.nan);pool.to_csv(R/'pooled_LR_descriptive.csv');lm=pool.pct.unstack().reindex(K);lm.index=[label(k) for k in K];heat(lm,'Actual adjacent ligand-source / receptor-target RNA co-detection (%)\nNA: no eligible edges; category denominators differ; no communication significance','immune_LR_codetection')
 ag=zones.groupby(['target','zone']).n_targets_in_zone.sum().unstack().reindex(K);ag.to_csv(R/'unique_target_zone_pooled.csv');cols=['High_near_including_shared','Other_only_near','Neither_high_nor_other_near'];pct=100*ag[cols].div(ag[cols].sum(1),axis=0)
 f,ax=plt.subplots(figsize=(14,16));left=np.zeros(len(K));names=['Near high (may also near other)','Near other detected only','Near neither high nor other detected'];colors=[COL[0],COL[1],'#CCCCCC']
 for col,name,color in zip(cols,names,colors):ax.barh(np.arange(len(K)),pct[col],left=left,color=color,label=name);left+=pct[col].to_numpy()
 ax.set_yticks(range(len(K)),[label(k)+f' (n={int(ag.loc[k,cols].sum())})' for k in K]);ax.invert_yaxis();ax.set_xlim(0,100);ax.set_xlabel('Unique target spots within each overlapping RNA category (%)');ax.set_title('Target neighborhood zones; zones partition each category\nNeither zone can still be near PGAM5-undetected sources');ax.legend(loc='lower right',fontsize=9);f.tight_layout();save(f,'immune_unique_target_zones')
 overlap=np.zeros((len(K),len(K)),int);rows=np.zeros(len(K),int)
 for item in SAMPLES:
  a=ad.read_h5ad(R/f'{item["sample"]}_immune_markers.h5ad');g=a.obs[['gate_'+k for k in K]].to_numpy(int);overlap+=g.T@g;rows+=g.sum(0);raw_objects.append(a)
 ov=pd.DataFrame(overlap,index=K,columns=K);ov.to_csv(R/'immune_program_overlap_counts.csv');op=100*ov.div(pd.Series(rows,index=K).replace(0,np.nan),axis=0);op.index=[label(k) for k in K];op.columns=[label(k) for k in K];heat(op,'Fraction of row-category target spots also meeting column-category gate (%)\nAsymmetric row denominators; overlaps do not establish same-cell identity','immune_program_overlap')
 selected=['All_T','CD4_related','CD8_supported','Treg_supported','NK_related','B_related','Plasma_related','cDC1_related','cDC2_related','pDC_related','Classical_monocyte_related','FCGR3A_myeloid_related','Neutrophil_related','Mast_related','SPP1_mac_related','TREM2_mac_related']
 gs=['CD3D','CD4','CD8A','CD8B','FOXP3','NCR1','NCR3','NCAM1','MS4A1','CD19','MZB1','JCHAIN','CLEC9A','XCR1','CD1C','CLEC4C','LILRA4','FCN1','CD14','FCGR3A','MS4A7','FCGR3B','CSF3R','TPSAB1','TPSB2','CPA3','SPP1','TREM2'];dot=[];f,ax=plt.subplots(figsize=(18,10))
 for row,k in enumerate(selected):
  data=[];norm=[]
  for a in raw_objects:
   mask=a.obs['gate_'+k].to_numpy(bool)&a.obs.near_high_source.to_numpy(bool);ix=a.var_names.get_indexer(gs);data.append(dense(a.layers['counts'])[mask][:,ix]);norm.append(dense(a.X)[mask][:,ix])
  raw=np.concatenate(data);x=np.concatenate(norm);n=len(raw)
  for j,g in enumerate(gs):
   if n:
    pct=100*(raw[:,j]>0).mean();mean=x[:,j].mean();ax.scatter(j,row,s=pct*1.4,c=mean,cmap='viridis',vmin=0,vmax=2,edgecolors='#444444',linewidths=.2);dot.append(dict(target=k,gene=g,n_high_near_targets=n,detection_pct=pct,mean_log1p_CP10k=mean))
  ax.text(len(gs)+.1,row,f'n={n}',va='center',fontsize=8)
 ax.set_xticks(range(len(gs)),gs,rotation=60,ha='right');ax.set_yticks(range(len(selected)),[label(k) for k in selected]);ax.invert_yaxis();ax.set_xlim(-1,len(gs)+2);ax.set_title('Defining marker RNA in unique target spots adjacent to high sources\nDot area: detection 0-100%; color: mean log1p(CP10k), capped at 2; not independent validation');f.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,2),cmap='viridis'),ax=ax,shrink=.6);f.tight_layout();save(f,'high_near_immune_marker_evidence');pd.DataFrame(dot).to_csv(R/'marker_dotplot_values.csv',index=False)
 groups=[['CD8_supported','Treg_supported','NK_related'],['B_related','Plasma_related'],['cDC1_related','cDC2_related','pDC_related'],['Classical_monocyte_related','FCGR3A_myeloid_related'],['Neutrophil_related','Mast_related'],['SPP1_mac_related','TREM2_mac_related']]
 for patient in P:
  items=[i for i in SAMPLES if i['patient']==patient];f,axs=plt.subplots(len(items),6,figsize=(23,3.2*len(items)),squeeze=False)
  for row,item in enumerate(items):
   sample=item['sample'];a=ad.read_h5ad(BASE/f'{sample}_spatial_marker_subset.h5ad');o=pd.read_csv(R/f'{sample}_annotations.csv.gz',index_col=0).loc[a.obs_names];xy=a.obsm['spatial']*a.uns['scalefactors']['tissue_lowres_scalef'];image=a.uns['spatial'][sample]['images']['lowres'];high=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy();e=a.obsp['spatial_connectivities'].tocoo()
   for j,keys in enumerate(groups):
    ax=axs[row,j];ax.imshow(image);target=np.zeros(len(o),bool)
    for n,k in enumerate(keys):
     mask=o['gate_'+k].to_numpy(bool);target|=mask;ax.scatter(xy[mask,0],xy[mask,1],s=10,color=COL[n],alpha=.7,label=label(k))
    use=high[e.row]&target[e.col]
    if use.any():ax.add_collection(LineCollection(np.stack([xy[e.row[use]],xy[e.col[use]]],axis=1),colors='#333333',linewidths=.6,alpha=.5))
    ax.scatter(xy[high,0],xy[high,1],s=20,facecolors='none',edgecolors='#B33232',linewidths=.8);ax.set_title(sample,fontsize=10);ax.axis('off')
    if row==0:ax.legend(loc='upper right',fontsize=6,framealpha=.8)
  f.suptitle(patient+' | True Visium H&E coordinates; red rings: high Mac RNA; colored overlapping targets; lines: direct neighbors',fontsize=13);f.tight_layout();save(f,patient+'_immune_spatial_maps')
if __name__=='__main__':run()
