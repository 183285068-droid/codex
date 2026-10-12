"""Focused T-program spatial analysis; mixed Visium spots, not resolved cells."""
from pathlib import Path
import os, json, hashlib, shutil
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE281759_Tcell_mpl')
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
R=Path(__file__).resolve().parent
S=R.parent/'GSE281759_PGAM5_immune_subgroup_neighbors'
B=R.parent/'GSE281759_spatial_HCC_PGAM5'
K=['All_T','Cytotoxic_T','Exhausted_T','Activated_T']
LABEL=['All T lineage','Cytotoxic program','Exhaustion-supported','Early activation-supported']
P=['A009','B013','C004','C008']
COL=['#285A8E','#B5822C','#A8536A','#607C44']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
def save(f,n):
 for ext in ['png','pdf']:f.savefig(R/'figures'/f'{n}.{ext}',dpi=180,bbox_inches='tight')
 plt.close(f)
def heat(m,title,name,cmap='Blues',lo=0,hi=100):
 f,a=plt.subplots(figsize=(10,8));im=a.imshow(m,aspect='auto',cmap=cmap,vmin=lo,vmax=hi)
 a.set_xticks(range(len(m.columns)),m.columns,rotation=30,ha='right');a.set_yticks(range(len(m)),m.index)
 for i in range(len(m)):
  for j in range(len(m.columns)):
   v=m.iloc[i,j]
   a.text(j,i,f'{v:.1f}' if np.isfinite(v) else 'NA',ha='center',va='center',fontsize=9)
 a.set_title(title,pad=15);f.colorbar(im,ax=a,shrink=.6);f.tight_layout();save(f,name)
manifest=[]
for name in ['neighborhood_counts.csv','spatial_statistics.csv','candidate_LR_codetection.csv','patient_neighborhood_descriptive.csv','patient_equal_section_effects.csv','analysis_config.json','audit_results.json']:
 p=S/name;manifest.append({'path':str(p.relative_to(R.parent)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
 if name.endswith('.csv'):
  t=pd.read_csv(p);t=t[t.target.isin(K)];t.to_csv(R/name,index=False)
 else:shutil.copyfile(p,R/('source_'+name))
counts=pd.read_csv(R/'neighborhood_counts.csv');stats=pd.read_csv(R/'spatial_statistics.csv');lr=pd.read_csv(R/'candidate_LR_codetection.csv');patient=pd.read_csv(R/'patient_neighborhood_descriptive.csv')
summary=counts[counts.source_group.eq('PGAM5_high')].groupby('target')[['n_sources','with_neighbor','n_target_spots']].sum().reindex(K)
summary['neighbor_pct']=100*summary.with_neighbor/summary.n_sources;summary.to_csv(R/'pooled_descriptive.csv')
f,ax=plt.subplots(figsize=(10,5));x=np.arange(4)
for i,p in enumerate(P):
 z=patient[(patient.patient==p)&patient.source_group.eq('PGAM5_high')].set_index('target').reindex(K)
 ax.barh(x+(i-1.5)*.18,z.neighbor_pct,height=.17,color=COL[i],label=f'{p}: {int(z.n_sources.iloc[0])} high source spots')
ax.set_yticks(x,LABEL);ax.invert_yaxis();ax.set_xlim(0,100);ax.set_xlabel('High source spots with >=1 direct target neighbor (%)');ax.legend(fontsize=8,loc='lower right');ax.set_title('Four-patient T-RNA neighborhood coverage\nPrograms overlap; spatial proximity is candidate evidence');f.tight_layout();save(f,'patient_T_program_coverage')
main=stats[stats.comparison.eq('high_vs_other_detected')]
matrix=main.pivot(index='sample',columns='target',values='delta_pp').reindex(columns=K);matrix.columns=LABEL
heat(matrix,'High minus other PGAM5-detected sources (percentage points)\nNA: not testable; original preplanned 144-test FDR retained','section_T_program_difference','RdBu_r',-50,50)
pool=lr.groupby(['target','pair'])[['eligible_edges','supported_edges']].sum();pool['pct']=100*pool.supported_edges/pool.eligible_edges.replace(0,np.nan);pool.to_csv(R/'pooled_LR_descriptive.csv')
heat(pool.pct.unstack().reindex(K),'Adjacent ligand-receptor RNA co-detection (%)\nDifferent edge denominators; no communication significance','T_candidate_LR_codetection')
checks=[];comb=[];edges=[]
for patient_id in P:
 items=[i for i in json.loads((B/'sample_scope.json').read_text()) if i['patient']==patient_id]
 f,axs=plt.subplots(len(items),3,figsize=(13,3.2*len(items)),squeeze=False)
 for row,item in enumerate(items):
  sample=item['sample'];file=B/f'{sample}_spatial_marker_subset.h5ad';a=ad.read_h5ad(file)
  manifest.append({'path':str(file.relative_to(R.parent)),'sha256':hashlib.sha256(file.read_bytes()).hexdigest()})
  annotation=S/f'{sample}_annotations.csv.gz'
  manifest.append({'path':str(annotation.relative_to(R.parent)),'sha256':hashlib.sha256(annotation.read_bytes()).hexdigest()})
  o=pd.read_csv(annotation,index_col=0).loc[a.obs_names]
  W=a.obsp['spatial_connectivities'];high=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy();mac=o.Mac_RNA_enriched.to_numpy(bool)&(np.asarray(W.sum(1)).ravel()>0)
  gates={k:o['gate_'+k].to_numpy(bool) for k in K};Y={k:np.asarray(W@gates[k].astype(int)).ravel()>0 for k in K}
  for group,mask in [('PGAM5_high',mac&high),('Other_detected',mac&~high&(o.PGAM5_raw.to_numpy()>0)),('Undetected',mac&(o.PGAM5_raw.to_numpy()==0))]:
   for k in K:
    t=counts[(counts['sample']==sample)&(counts.source_group==group)&(counts.target==k)].iloc[0]
    assert int(mask.sum())==t.n_sources and int(Y[k][mask].sum())==t.with_neighbor
    checks.append(f'{sample}/{group}/{k}: independent graph recount')
  code=gates[K[1]].astype(int)+2*gates[K[2]]+4*gates[K[3]]
  for bit in range(8):
   target=gates['All_T']&(code==bit);near=np.asarray(W@target.astype(int)).ravel()>0
   comb.append(dict(**item,program_combination=bit,cytotoxic=bool(bit&1),exhaustion=bool(bit&2),activation=bool(bit&4),n_target_spots=int(target.sum()),n_high_sources=int((mac&high).sum()),high_sources_with_neighbor=int(near[mac&high].sum())))
  assert sum(v['n_target_spots'] for v in comb if v['sample']==sample)==int(gates['All_T'].sum());checks.append(sample+': exclusive target partition')
  e=W.tocoo();xy=a.obsm['spatial']*a.uns['scalefactors']['tissue_lowres_scalef'];image=a.uns['spatial'][sample]['images']['lowres']
  for j,k in enumerate(K[1:]):
   use=high[e.row]&gates[k][e.col];src,tgt=e.row[use],e.col[use]
   edges.extend(dict(sample=sample,patient=patient_id,target=k,source_barcode=a.obs_names[s],target_barcode=a.obs_names[t]) for s,t in zip(src,tgt))
   ax=axs[row,j];ax.imshow(image);ax.scatter(xy[gates[k],0],xy[gates[k],1],s=9,color=COL[j],alpha=.7,label=k)
   if len(src):ax.add_collection(LineCollection(np.stack([xy[src],xy[tgt]],axis=1),colors='#444444',linewidths=.6,alpha=.5))
   ax.scatter(xy[high,0],xy[high,1],s=20,facecolors='none',edgecolors='#B33232',linewidths=.8,label='PGAM5-high Mac RNA')
   ax.set_title(f'{sample} | {k}',fontsize=10);ax.axis('off')
   if row==0:ax.legend(fontsize=7,loc='upper right')
 f.suptitle(patient_id+' | True Visium coordinates on H&E; lines: direct neighbors; mixed RNA spots',fontsize=11);f.tight_layout();save(f,patient_id+'_T_spatial_maps')
pd.DataFrame(comb).to_csv(R/'exclusive_T_program_combinations.csv',index=False)
pd.DataFrame(edges).to_csv(R/'high_source_T_neighbor_edges.csv.gz',index=False)
ct=pd.DataFrame(comb).groupby('program_combination').n_target_spots.sum()
f,ax=plt.subplots(figsize=(11,5));names=['No supported state','Cytotoxic only','Exhaustion only','Cytotoxic + exhaustion','Activation only','Cytotoxic + activation','Exhaustion + activation','All three']
ax.barh(names,ct.reindex(range(8)),color=COL[0]);ax.invert_yaxis();ax.set_xlabel('T-lineage RNA target spots (mutually exclusive combinations)');ax.set_title('State overlap within T-lineage RNA gates\nNo supported state does not establish naive or memory identity');f.tight_layout();save(f,'T_program_overlap')
for name in ['primary_q_144','joint_q_432']:
 assert not (stats[name].dropna()<.05).any();checks.append(name+': no FDR-significant result')
(R/'audit_results.json').write_text(json.dumps({'passed':len(checks),'checks':checks,'scope':'Independent graph recount and target partition; source full audit copied separately'},indent=2))
(R/'input_manifest.json').write_text(json.dumps(manifest,indent=2))
print(summary.to_string());print('PASS',len(checks));print(ct.to_string())
