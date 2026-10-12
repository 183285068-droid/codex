import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE281759_CD8_mpl');os.environ.setdefault('XDG_CACHE_HOME','/tmp/GSE281759_CD8_cache')
import numpy as np,pandas as pd,anndata as ad,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from analyze import R,BASE,K,ALL,SAMPLES,dense
P=['A009','B013','C004','C008'];COL=['#285A8E','#B5822C','#A8536A','#607C44']
LABEL=['Naive like','Central memory like','GZMK memory like','Cytotoxic effector','Terminal effector like','Exhaustion supported','Progenitor exhaustion like','Terminal exhaustion like','Resident memory like','Early activation','Cycling','IFN response']
NAMES=dict(zip(K,LABEL));NAMES['All_CD8_reference']='All CD8 reference'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
def save(f,n):
 for ext in ['png','pdf']:f.savefig(R/'figures'/f'{n}.{ext}',dpi=160,bbox_inches='tight')
 plt.close(f)
def heat(m,title,name,cmap='Blues',lo=0,hi=100):
 f,ax=plt.subplots(figsize=(14,9));im=ax.imshow(m.to_numpy(float),aspect='auto',cmap=cmap,vmin=lo,vmax=hi);ax.set_xticks(range(len(m.columns)),m.columns,rotation=45,ha='right');ax.set_yticks(range(len(m)),m.index)
 for i in range(len(m)):
  for j in range(len(m.columns)):
   v=m.iloc[i,j];ax.text(j,i,f'{v:.1f}' if np.isfinite(v) else 'NA',ha='center',va='center',fontsize=8)
 ax.set_title(title,pad=18);f.colorbar(im,ax=ax,shrink=.65);f.tight_layout();save(f,name)
def run():
 p=pd.read_csv(R/'patient_descriptive.csv');s=pd.read_csv(R/'spatial_statistics.csv');lr=pd.read_csv(R/'candidate_LR_codetection.csv')
 f,ax=plt.subplots(figsize=(13,9));y=np.arange(len(ALL))
 for i,patient in enumerate(P):
  z=p[(p.patient==patient)&p.source_group.eq('PGAM5_high')].set_index('target').reindex(ALL);ax.barh(y+(i-1.5)*.18,z.neighbor_pct,height=.17,color=COL[i],label=f'{patient}: {int(z.n_sources.iloc[0])} high sources')
 ax.set_yticks(y,[NAMES[k] for k in ALL]);ax.invert_yaxis();ax.set_xlim(0,30);ax.set_xlabel('High Mac-RNA sources with >=1 direct CD8-program neighbor (%)');ax.legend(loc='lower right',fontsize=9);ax.set_title('Fine CD8-supported RNA programs: patient coverage\nPrograms overlap; cell identity/function unresolved in mixed Visium spots');f.tight_layout();save(f,'patient_CD8_fine_coverage')
 main=s[s.comparison.eq('high_vs_other_detected')];m=main.pivot(index='sample',columns='target',values='delta_pp').reindex(index=[i['sample'] for i in SAMPLES],columns=K);m.columns=LABEL
 heat(m,'High minus other PGAM5-detected sources (percentage points)\nNA: insufficient groups or invariant conditional null','section_CD8_fine_differences','RdBu_r',-20,20)
 pool=lr.groupby(['target','pair'])[['eligible_edges','supported_edges']].sum();pool['pct']=100*pool.supported_edges/pool.eligible_edges.replace(0,np.nan);pool.to_csv(R/'pooled_LR_descriptive.csv');lm=pool.pct.unstack().reindex(ALL);lm.index=[NAMES[k] for k in ALL]
 heat(lm,'Actual adjacent ligand-source / receptor-target RNA co-detection (%)\nNA: no eligible edges; different denominators; no communication significance','CD8_fine_LR_codetection')
 z=pd.read_csv(R/'CD8_target_neighborhood_programs.csv');agg=z.groupby(['zone','target'])[['n_CD8','n_program']].sum();agg['pct']=100*agg.n_program/agg.n_CD8.replace(0,np.nan);agg.to_csv(R/'CD8_target_zone_pooled.csv');zm=agg.pct.unstack('zone').reindex(K);zm.index=LABEL
 heat(zm,'T-program frequency within unique CD8 target neighborhood zones (%)\nHigh-near includes targets shared with other detected sources; descriptive only','CD8_unique_target_zones')
 ov=pd.read_csv(R/'program_overlap.csv').groupby(['first','second']).both.sum().unstack().reindex(index=K,columns=K);ov.index=LABEL;ov.columns=LABEL
 heat(ov,'Overlapping fine CD8-program target counts\nDiagonal = total target count; overlaps do not establish single-cell coexpression','CD8_program_overlap','Blues',0,max(1,float(ov.to_numpy().max())))
 edges=pd.read_csv(R/'high_source_CD8_directed_edges.csv.gz');unique=edges[['sample','target_barcode','signature']].drop_duplicates();signatures=unique.signature.value_counts();table=pd.DataFrame({'signature_id':[f'S{i+1:02d}' for i in range(len(signatures))],'signature':signatures.index,'unique_high_near_CD8_targets':signatures.values});table.to_csv(R/'high_near_exact_signature_summary.csv',index=False)
 membership=np.array([[int(k in sig.split('+')) for k in K] for sig in signatures.index]);f,(left,right)=plt.subplots(1,2,figsize=(15,10),gridspec_kw={'width_ratios':[2,5]},sharey=True);pos=np.arange(len(table));left.barh(pos,table.unique_high_near_CD8_targets,color=COL[0]);left.set_yticks(pos,table.signature_id);left.set_xlabel('Unique CD8 targets');left.invert_yaxis();left.set_title('Observed signature counts')
 for row in range(len(table)):
  left.text(table.unique_high_near_CD8_targets.iloc[row]+.15,row,str(table.unique_high_near_CD8_targets.iloc[row]),va='center',fontsize=9)
  right.scatter(np.arange(len(K)),np.full(len(K),row),s=45,c=['#285A8E' if v else '#DDDDDD' for v in membership[row]],edgecolors='none')
 right.set_xticks(range(len(K)),LABEL,rotation=60,ha='right');right.set_title('Blue: supported program; grey: gate not met');right.tick_params(axis='y',labelleft=False);right.set_xlim(-.6,len(K)-.4);f.suptitle('Finest descriptive combinations in unique CD8 targets adjacent to high sources\nSignature IDs map to high_near_exact_signature_summary.csv; combinations are not validated cell subtypes',fontsize=12);f.tight_layout();save(f,'CD8_high_near_exact_signatures')
 # Dot plot: anchors and distinguishing markers, not an independent identity validation.
 gs=['CD8A','CD8B','CCR7','SELL','TCF7','LEF1','IL7R','GZMK','PRF1','GZMB','FGFBP2','CX3CR1','PDCD1','TOX','SLAMF6','HAVCR2','ENTPD1','CXCL13','ITGAE','ZNF683','CXCR6','CD69','IL2RA','TNFRSF9','MKI67','TOP2A','IFIT1','MX1'];records={k:[] for k in ALL}
 for item in SAMPLES:
  a=ad.read_h5ad(R/f'{item["sample"]}_CD8_fine_markers.h5ad');ix=a.var_names.get_indexer(gs)
  for k in ALL:
   mask=a.obs['gate_'+k].to_numpy(bool)&a.obs.CD8_near_high.to_numpy(bool)
   v=dense(a.layers['counts'])[mask][:,ix];n=dense(a.X)[mask][:,ix];records[k].append((v,n))
 f,ax=plt.subplots(figsize=(17,9));dotrows=[]
 for row,k in enumerate(ALL):
  raw=np.concatenate([v for v,n in records[k]],axis=0);norm=np.concatenate([n for v,n in records[k]],axis=0)
  for j,g in enumerate(gs):
   if len(raw):
    pct=100*(raw[:,j]>0).mean();mean=norm[:,j].mean();ax.scatter(j,row,s=pct*1.8,c=mean,cmap='viridis',vmin=0,vmax=2,edgecolors='#444444',linewidths=.2);dotrows.append(dict(target=k,gene=g,n_targets=len(raw),detection_pct=pct,mean_log1p_CP10k=mean))
 ax.set_xticks(range(len(gs)),gs,rotation=60,ha='right');ax.set_yticks(range(len(ALL)),[NAMES[k]+f' (n={sum(len(v) for v,n in records[k])})' for k in ALL]);ax.invert_yaxis();ax.set_xlim(-1,len(gs));ax.set_title('Markers in unique CD8 targets adjacent to high sources (defining genes, not validation)\nDot area: detection 0-100%; color: mean log1p(CP10k), capped at 2; n=0 rows are empty');f.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,2),cmap='viridis'),ax=ax,shrink=.6);f.tight_layout();save(f,'CD8_fine_marker_evidence');pd.DataFrame(dotrows).to_csv(R/'marker_dotplot_values.csv',index=False)
 for patient in P:
  items=[i for i in SAMPLES if i['patient']==patient]
  for part,keys in enumerate([K[:6],K[6:]],start=1):
   f,axs=plt.subplots(len(items),6,figsize=(23,3*len(items)),squeeze=False)
   for row,item in enumerate(items):
    sample=item['sample'];a=ad.read_h5ad(BASE/f'{sample}_spatial_marker_subset.h5ad');o=pd.read_csv(R/f'{sample}_annotations.csv.gz',index_col=0).loc[a.obs_names];xy=a.obsm['spatial']*a.uns['scalefactors']['tissue_lowres_scalef'];image=a.uns['spatial'][sample]['images']['lowres'];high=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy();e=a.obsp['spatial_connectivities'].tocoo()
    for j,k in enumerate(keys):
     target=o['gate_'+k].to_numpy(bool);use=high[e.row]&target[e.col];ax=axs[row,j];ax.imshow(image);ax.scatter(xy[target,0],xy[target,1],s=12,color=COL[0],alpha=.8)
     if use.any():ax.add_collection(LineCollection(np.stack([xy[e.row[use]],xy[e.col[use]]],axis=1),colors='#222222',linewidths=.7,alpha=.6))
     ax.scatter(xy[high,0],xy[high,1],s=20,facecolors='none',edgecolors='#B33232',linewidths=.8);ax.set_title(sample+' | '+NAMES[k]+f' (n={int(target.sum())})',fontsize=9);ax.axis('off')
   f.suptitle(patient+' | True Visium H&E coordinates; red rings: high Mac RNA; blue: CD8 programs; lines: direct neighbors',fontsize=13);f.tight_layout();save(f,f'{patient}_CD8_fine_spatial_part{part}')
if __name__=='__main__':run()
