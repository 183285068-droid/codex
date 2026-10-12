"""Source-backed review figures; reuse audited fixed gates and planned test families."""
from pathlib import Path
import os,json,hashlib
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE245908_review_mpl')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/GSE245908_review_cache')
import numpy as np,pandas as pd,anndata as ad,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent
SAMPLES=['CHC20','CHC23'];COLORS=['#285A8E','#B5822C']
NAMES={'All_T':'All T-lineage (baseline)','Plasma_related':'Plasma-related','Inflammatory_monocyte':'Inflammatory myeloid-related','Naive_memory_T':'Naive/memory T-related','cDC2_related':'cDC2-related','CD8_cytotoxic_T':'CD8/cytotoxic T-related','Neutrophil_related':'Granulocyte-related','NK_related':'NK-related','Mast_related':'Mast-related','B_related':'B-related','T_exhaustion':'Exhaustion-supported T','T_activation':'Activation-supported T','Treg_related':'Treg-related','cDC1_related':'cDC1-related','pDC_related':'pDC-related'}
def save(fig,name):
 fig.savefig(ROOT/f'{name}.png',dpi=180,bbox_inches='tight');fig.savefig(ROOT/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
def run():
 files=['subgroup_neighborhood_counts.csv','strict_marker_descriptive_sensitivity.csv','spatial_conditional_tests.csv','source_group_summary.csv','candidate_LR_by_immune_subgroup.csv','audit_results.json']+[s+'_immune_marker_subset.h5ad' for s in SAMPLES]
 manifest=[dict(path=f,sha256=hashlib.sha256((BASE/f).read_bytes()).hexdigest(),bytes=(BASE/f).stat().st_size) for f in files]
 (ROOT/'reference_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 counts=pd.read_csv(BASE/files[0]);strict=pd.read_csv(BASE/files[1]);stats=pd.read_csv(BASE/files[2]);sources=pd.read_csv(BASE/files[3]).set_index('sample');high=counts[counts.source_group.eq('PGAM5_high')];checks=[]
 for sample in SAMPLES:
  a=ad.read_h5ad(BASE/f'{sample}_immune_marker_subset.h5ad');o=a.obs;h=o.source_group.eq('PGAM5_high').to_numpy();G=a.obsp['spatial_connectivities'];assert int(h.sum())==int(sources.loc[sample,'PGAM5_high'])
  for row in high[high['sample']==sample].itertuples():
   mask=o['gate_'+row.target].to_numpy(bool);y=np.asarray(G@mask.astype(int)).ravel();assert int((y[h]>0).sum())==row.sources_with_neighbor;assert int(y[h].sum())==row.directed_edges;assert int(mask.sum())==row.n_target_spots
   checks.append(sample+' '+row.target+' direct source count/edge/target')
  raw=lambda gene:a[:,[gene]].layers['counts'].toarray().ravel()
  masks={'Treg_add_IKZF2_and_IL2RA_or_CTLA4':o.gate_Treg_related.to_numpy(bool)&(raw('IKZF2')>0)&((raw('IL2RA')>0)|(raw('CTLA4')>0)), 'CD8_cytotoxic_require_CD8B':o.gate_CD8_cytotoxic_T.to_numpy(bool)&(raw('CD8B')>0),'Monocyte_require_FCN1':o.gate_Inflammatory_monocyte.to_numpy(bool)&(raw('FCN1')>0),'Neutrophil_require_FCGR3B':o.gate_Neutrophil_related.to_numpy(bool)&(raw('FCGR3B')>0)}
  for row in strict[strict['sample']==sample].itertuples():
   mask=masks[row.rule];y=np.asarray(G@mask.astype(int)).ravel();assert int((y[h]>0).sum())==row.PGAM5_high_sources_with_neighbor;assert int(mask.sum())==row.target_spots;checks.append(sample+' '+row.rule+' strict count')
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
 order=list(NAMES);fig,ax=plt.subplots(figsize=(11,8));y=np.arange(len(order))
 for j,sample in enumerate(SAMPLES):
  d=high[high['sample']==sample].set_index('target').reindex(order);vals=d.source_adjacency_pct.to_numpy();ax.barh(y+(j-.5)*.30,vals,height=.28,color=COLORS[j],label=sample+f' (n={int(sources.loc[sample,"PGAM5_high"])} high sources)')
  for yy,value,row in zip(y+(j-.5)*.30,vals,d.itertuples()):ax.text(value+1,yy,f'{int(row.sources_with_neighbor)}/{int(row.n_sources)}',va='center',fontsize=8)
 ax.set_yticks(y,[NAMES[k] for k in order]);ax.invert_yaxis();ax.set_xlim(0,105);ax.set_xlabel('High source spots with >=1 direct target neighbor (%)');ax.set_title('GSE245908 | PGAM5-high macrophage-RNA immune neighborhoods\nSpatial RNA candidates; no subgroup passes planned FDR correction');ax.legend(loc='lower right');ax.spines[['top','right']].set_visible(False);fig.tight_layout();save(fig,'candidate_neighborhood_overview')
 mappings=[('Inflammatory_monocyte','Monocyte_require_FCN1','Myeloid: add FCN1'),('CD8_cytotoxic_T','CD8_cytotoxic_require_CD8B','Cytotoxic T: add CD8B'),('Treg_related','Treg_add_IKZF2_and_IL2RA_or_CTLA4','Treg: add IKZF2 + IL2RA/CTLA4'),('Neutrophil_related','Neutrophil_require_FCGR3B','Granulocyte: add FCGR3B')]
 fig,axes=plt.subplots(1,2,figsize=(13,5),sharey=True);review=[]
 for ax,sample in zip(axes,SAMPLES):
  denominator=int(sources.loc[sample,'PGAM5_high']);broad=[];narrow=[]
  for k,r,label in mappings:
   a=high[(high['sample']==sample)&high.target.eq(k)].iloc[0];b=strict[(strict['sample']==sample)&strict.rule.eq(r)].iloc[0];broad.append(int(a.sources_with_neighbor));narrow.append(int(b.PGAM5_high_sources_with_neighbor));review.append(dict(sample=sample,target=k,strict_rule=r,n_high=denominator,broad_sources=int(a.sources_with_neighbor),strict_sources=int(b.PGAM5_high_sources_with_neighbor),post_hoc_descriptive_only=True))
  y=np.arange(len(mappings));ax.barh(y-.15,np.array(broad)/denominator*100,height=.28,color='#285A8E',label='Fixed main RNA gate');ax.barh(y+.15,np.array(narrow)/denominator*100,height=.28,color='#B5822C',label='Stricter marker support')
  for yy,b,n in zip(y,broad,narrow):
   ax.text(100*b/denominator+1,yy-.15,f'{b}/{denominator}',va='center',fontsize=9);ax.text(100*n/denominator+1,yy+.15,f'{n}/{denominator}',va='center',fontsize=9)
  ax.set_xlim(0,105);ax.set_title(sample);ax.set_xlabel('High sources with target neighbor (%)');ax.spines[['top','right']].set_visible(False)
 axes[0].set_yticks(range(len(mappings)),[x[2] for x in mappings]);axes[0].invert_yaxis();axes[1].legend(loc='lower right');fig.suptitle('Marker specificity matters | Post-hoc descriptive checks; no new significance tests');fig.tight_layout();save(fig,'strict_marker_comparison')
 pd.DataFrame(review).to_csv(ROOT/'marker_sensitivity_source_counts.csv',index=False);high.to_csv(ROOT/'high_source_neighborhood_counts.csv',index=False)
 report=dict(status='passed',checks=len(checks),details=checks,main_min_q=float(stats.primary_BH_q_30.min()),joint_min_q=float(stats.BH_q_90_planned_tests.min()),significant_primary=int(stats.primary_BH_q_30.lt(.05).sum()),significant_joint=int(stats.BH_q_90_planned_tests.lt(.05).sum()),original_audit=json.loads((BASE/'audit_results.json').read_text())['checks'])
 (ROOT/'review_checks.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items() if k!='details'},indent=2))
if __name__=='__main__':run()
