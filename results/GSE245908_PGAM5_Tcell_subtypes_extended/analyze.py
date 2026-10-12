"""Fixed T RNA subtype gates and conditional spatial comparisons."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
import frozen_spatial_methods as f
R=Path(__file__).resolve().parent;OLD=R.parent/'GSE245908_PGAM5_Tcell_spatial_interactions';f.ROOT=R;f.WORK=R/'unused_cache'
PANELS={'CD8_identity':['CD8A','CD8B'],'CD4_identity':['CD4'],'Treg':['FOXP3','IKZF2','IL2RA','CTLA4'],'Naive_memory':['CCR7','SELL','TCF7','LEF1','IL7R'],**{k:f.PANELS[k] for k in ['Activation_early','Effector_cytotoxic','Exhaustion_core']}}
TYPES=['All_CD8','CD8_exhaustion','CD8_activation','CD8_effector','CD8_naive_memory','CD4_other','Treg','T_unresolved']
GROUPS=['PGAM5_high','Other_detected','Undetected'];STATES={'CD8_exhaustion':'Exhaustion_core','CD8_activation':'Activation_early','CD8_effector':'Effector_cytotoxic','CD8_naive_memory':'Naive_memory'}
def run():
 stats=[];support=[];summary=[];neighbors=[];lr=[];coverage=[];thresholds=[];sensitivity=[];manifest=[]
 for s in f.SAMPLES:
  print('Processing',s,flush=True);a=f.load_sample(s);o=pd.read_csv(OLD/f'{s}_spot_states.csv',index_col=0).loc[a.obs_names].copy();W=f.graph(o);degree=np.asarray(W.sum(1)).ravel();t=o.T_enriched.to_numpy(bool);high=o.source_group.eq('PGAM5_high').to_numpy();pos=o.source_group.isin(GROUPS[:2]).to_numpy();mac=o.Mac_enriched.to_numpy(bool)
  for k,gs in PANELS.items():
   assert all(g in a.var_names for g in gs),(s,k);z=f.genes(a,gs);o['score_'+k]=z.mean(1);o['detected_'+k]=(z>0).sum(1);coverage.append(dict(sample=s,panel=k,requested=';'.join(gs),n_present=len(gs)))
  cd8=t&(f.genes(a,['CD8A']).ravel()>0)&(o.score_CD8_identity.to_numpy()>o.score_CD4_identity.to_numpy());cd4=t&(o.score_CD4_identity.to_numpy()>0)&(o.score_CD4_identity.to_numpy()>o.score_CD8_identity.to_numpy());unresolved=t&~(cd8|cd4)
  treg=cd4&(f.genes(a,['FOXP3']).ravel()>0)&(o.detected_Treg.to_numpy()>=2)
  masks={'All_CD8':cd8,'CD4_other':cd4&~treg,'Treg':treg,'T_unresolved':unresolved}
  for label,panel in STATES.items():
   cut=float(o.loc[cd8,'score_'+panel].quantile(.75)) if cd8.any() else np.nan
   masks[label]=cd8&(o['score_'+panel].to_numpy()>=cut)&(o['score_'+panel].to_numpy()>0)&(o['detected_'+panel].to_numpy()>=2);thresholds.append(dict(sample=s,subtype=label,Q75_cutoff=cut,reference_n=int(cd8.sum())))
  for label in TYPES:o['gate_'+label]=masks[label]
  o['identity']=np.select([treg,cd4&~treg,cd8,unresolved],['Treg_candidate','CD4_other_candidate','CD8_candidate','T_unresolved'],default='Outside_T_gate')
  rng=np.random.default_rng(24590842+(s=='CHC23'))
  Y=np.column_stack([np.divide(W@masks[k].astype(float),degree,out=np.zeros(len(o)),where=degree>0) for k in TYPES]);rows,info=f.permutation_results(o,pos,high,Y,TYPES,'subtype_coverage',rng,s);stats+=rows;support.append(info)
  n8=W@cd8.astype(float);valid=n8>0;Yscore=np.column_stack([np.divide(W@(o['score_'+p].to_numpy()*cd8),n8,out=np.zeros(len(o)),where=n8>0) for p in STATES.values()]);rows,info=f.permutation_results(o,pos&valid,high,Yscore,list(STATES.values()),'CD8_neighbor_programs',rng,s);stats+=rows;support.append(info)
  for j,label in enumerate(TYPES):
   o['neighbor_fraction_'+label]=Y[:,j]
   for g in GROUPS:
    keep=o.source_group.eq(g).to_numpy();nt=W@masks[label].astype(float);summary.append(dict(sample=s,subtype=label,n_target_spots=int(masks[label].sum()),source_group=g,n_sources=int(keep.sum()),n_sources_with_target=int((keep&(nt>0)).sum()),directed_edges=int(nt[keep].sum()),mean_neighbor_fraction=float(Y[keep,j].mean())))
  for j,panel in enumerate(STATES.values()):o['CD8_neighbor_'+panel]=np.where(valid,Yscore[:,j],np.nan)
  # Gate sensitivity is descriptive, never used to choose the primary gate.
  for name,receiver in [('CD8A_only_broad',t&(f.genes(a,['CD8A']).ravel()>0)),('CD8A_CD8B_both_strict',cd8&(f.genes(a,['CD8B']).ravel()>0)),('Treg_FOXP3_three_markers_strict',treg&(o.detected_Treg.to_numpy()>=3))]:
   nt=W@receiver.astype(float)
   for g in GROUPS:
    keep=o.source_group.eq(g).to_numpy();sensitivity.append(dict(sample=s,gate=name,source_group=g,n_targets=int(receiver.sum()),n_sources=int(keep.sum()),n_sources_with_target=int((keep&(nt>0)).sum()),n_edges=int(nt[keep].sum()),mean_neighbor_fraction=float((nt[keep]/degree[keep]).mean())))
  ri,ci=W.nonzero();keep=mac[ri]&t[ci];ri=ri[keep];ci=ci[keep];edges=pd.DataFrame({'sample':s,'source_barcode':o.index[ri],'target_barcode':o.index[ci],'source_group':o.source_group.iloc[ri].to_numpy(),'target_identity':o.identity.iloc[ci].to_numpy()})
  for k in TYPES:edges['target_'+k]=masks[k][ci]
  for ligand,receptor,context in f.PAIRS:
   l=f.genes(a,[ligand],layer='counts').ravel()/o.total_counts.to_numpy()*1e4;r=np.min(f.genes(a,receptor.split('_'),layer='counts')/o.total_counts.to_numpy()[:,None]*1e4,axis=1);z=np.sqrt(l[ri]*r[ci]);edges[ligand+'--'+receptor]=z
   for k in TYPES:
    for g in GROUPS:
     select=(o.source_group.iloc[ri].to_numpy()==g)&masks[k][ci];v=z[select];lr.append(dict(sample=s,subtype=k,source_group=g,ligand=ligand,receptor=receptor,n_edges=len(v),n_sources=len(np.unique(ri[select])),n_targets=len(np.unique(ci[select])),n_codetected=int((v>0).sum()),codetection_pct=100*(v>0).mean() if len(v) else np.nan,mean_edge_compatibility=float(v.mean()) if len(v) else np.nan))
  edges.to_csv(R/f'{s}_subtype_LR_edges.csv.gz',index=False);o.to_csv(R/f'{s}_spot_subtypes.csv',index_label='barcode',float_format='%.17g')
  for p in [f.BASE/f'{s}_spatial_edges.csv.gz',OLD/f'{s}_spot_states.csv',*sorted((f.BASE/'inputs'/s).glob('*'))]:
   manifest.append(dict(path=str(p.relative_to(R.parent)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
 df=pd.DataFrame(stats);df['BH_q']=np.nan;ok=df.spatial_conditional_p.notna();df.loc[ok,'BH_q']=multipletests(df.loc[ok,'spatial_conditional_p'],method='fdr_bh')[1];df.to_csv(R/'spatial_subtype_tests.csv',index=False)
 pd.DataFrame(summary).to_csv(R/'subtype_neighbor_summary.csv',index=False);pd.DataFrame(lr).to_csv(R/'LR_by_subtype.csv',index=False);pd.DataFrame(coverage).to_csv(R/'marker_coverage.csv',index=False);pd.DataFrame(thresholds).to_csv(R/'state_thresholds.csv',index=False);pd.DataFrame(sensitivity).to_csv(R/'gate_sensitivity.csv',index=False);pd.DataFrame(support).to_csv(R/'permutation_support.csv',index=False)
 (R/'input_manifest.json').write_text(json.dumps(manifest,indent=2));(R/'analysis_config.json').write_text(json.dumps(dict(panels=PANELS,types=TYPES,seed=24590842,permutations=f.B,multiple_testing='Single BH family: 2 sections x (8 coverage + 4 continuous CD8 programs) = 24 fixed tests; no LR subtype p-values',states_nonexclusive=True,identity_mutually_exclusive=True,LR_pairs=f.PAIRS),indent=2));print(df.to_string(index=False),flush=True)
if __name__=='__main__':run()
