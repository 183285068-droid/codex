from pathlib import Path
import os,json,hashlib
import numpy as np,pandas as pd,h5py
from scipy import sparse
from scipy.stats import false_discovery_control
import frozen_spatial_methods as f
R=Path(__file__).resolve().parent;BASE=R.parent/'GSE238264_PGAM5_spatial_response';RAW=Path(os.environ.get('GSE238264_INPUT','/workspace/scratch/GSE238264/inputs'));f.ROOT=R
SAMPLES=['HCC1R','HCC2R','HCC3R','HCC4R','HCC5NR','HCC6NR','HCC7NR'];GROUPS=['PGAM5_high','Other_detected','Undetected']
PANELS={'CD8_identity':['CD8A','CD8B'],'CD4_identity':['CD4'],'Treg_program':['FOXP3','IKZF2','IL2RA','CTLA4'],'Exhaustion':['PDCD1','LAG3','TIGIT','CTLA4','TOX','CXCL13'],'Activation':['CD69','IL2RA','TNFRSF9','ICOS','CD40LG','IL2'],'Effector':['IFNG','PRF1','GZMB','GZMA','CCL5','TBX21'],'Naive_memory':['CCR7','SELL','TCF7','LEF1','IL7R'],'Exhaustion_sensitivity':['PDCD1','LAG3','TOX','CXCL13','HAVCR2','ENTPD1']}
TYPES=['All_T','All_CD8','CD8_exhaustion','CD8_activation','CD8_effector','CD8_naive_memory','CD4_other','Treg','T_unresolved'];STATES={'CD8_exhaustion':'Exhaustion','CD8_activation':'Activation','CD8_effector':'Effector','CD8_naive_memory':'Naive_memory'}
PAIRS=[('CD274','PDCD1'),('PDCD1LG2','PDCD1'),('PVR','TIGIT'),('NECTIN2','TIGIT'),('LGALS9','HAVCR2'),('CD80','CTLA4'),('CD86','CTLA4'),('CD80','CD28'),('CD86','CD28'),('ICOSLG','ICOS'),('TNFSF9','TNFRSF9'),('TNFSF4','TNFRSF4'),('TGFB1','TGFBR1_TGFBR2'),('IL10','IL10RA_IL10RB'),('CXCL9','CXCR3'),('CXCL10','CXCR3'),('CXCL11','CXCR3'),('CCL3','CCR5'),('CCL4','CCR5'),('CCL5','CCR5'),('SPP1','CD44'),('CXCL12','CXCR4')]
def run():
 stats=[];summary=[];lrs=[];coverage=[];support=[];manifest=[];thresholds=[];sens=[]
 resource=pd.read_csv(R.parent/'GSE245908_PGAM5_Tcell_spatial_interactions/source_inputs/omni_resource.csv');selected=[]
 for l,r in PAIRS:
  hit=resource[(resource.source_genesymbol==l)&(resource.target_genesymbol==r)];assert len(hit);selected.append(dict(ligand=l,receptor=r,resources=';'.join(sorted(hit.resource.unique())),snapshot='LIANA 1.10.0 omni_resource.csv'))
 pd.DataFrame(selected).to_csv(R/'selected_LR_pairs.csv',index=False)
 for si,sample in enumerate(SAMPLES):
  bp=BASE/f'{sample}_spot_annotations.csv.gz';o=pd.read_csv(bp,index_col=0);e0=pd.read_csv(BASE/f'{sample}_spatial_edges.csv.gz');ix={b:i for i,b in enumerate(o.index)};G=sparse.csr_matrix((np.ones(len(e0)),([ix[b] for b in e0.source_barcode],[ix[b] for b in e0.target_barcode])),shape=(len(o),len(o)));degree=G.sum(1).A1;n=len(o);hp=RAW/sample/'filtered_feature_bc_matrix.h5'
  with h5py.File(hp) as h:
   m=h['matrix'];raw=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=m['shape'][:]);names=np.array([g.decode() for g in m['features']['name'][:]]);bc=[g.decode() for g in m['barcodes'][:]]
  bix={b:i for i,b in enumerate(bc)};raw=raw[:,[bix[b] for b in o.index]];assert np.allclose(raw.sum(0).A1,o.total_counts)
  genes=list(dict.fromkeys([g for gs in PANELS.values() for g in gs]+[l for l,r in PAIRS]+[g for l,r in PAIRS for g in r.split('_')]));count=np.column_stack([raw[names==g].sum(0).A1 for g in genes]);cp=count/o.total_counts.to_numpy()[:,None]*1e4;x=np.log1p(cp)
  for k,gs in PANELS.items():
   assert all(g in names for g in gs);idx=[genes.index(g) for g in gs];o['score_'+k]=x[:,idx].mean(1);o['detected_'+k]=(count[:,idx]>0).sum(1);coverage.append(dict(sample=sample,panel=k,genes=';'.join(gs),n_present=len(gs)))
  t=o.T_RNA.to_numpy(bool);cd8=t&(count[:,genes.index('CD8A')]>0)&(o.score_CD8_identity.to_numpy()>o.score_CD4_identity.to_numpy());cd4=t&(o.score_CD4_identity.to_numpy()>0)&(o.score_CD4_identity.to_numpy()>o.score_CD8_identity.to_numpy());treg=cd4&(count[:,genes.index('FOXP3')]>0)&(o.detected_Treg_program.to_numpy()>=2);un=t&~(cd8|cd4)
  masks={'All_T':t,'All_CD8':cd8,'CD4_other':cd4&~treg,'Treg':treg,'T_unresolved':un}
  for k,panel in STATES.items():
   cut=float(o.loc[cd8,'score_'+panel].quantile(.75)) if cd8.any() else np.nan;mask=cd8&(o['score_'+panel].to_numpy()>=cut)&(o['score_'+panel].to_numpy()>0)&(o['detected_'+panel].to_numpy()>=2);masks[k]=mask;thresholds.append(dict(sample=sample,subtype=k,n_CD8=int(cd8.sum()),cutoff=cut))
  for k in TYPES:o['gate_'+k]=masks[k]
  o['T_identity']=np.select([treg,cd4&~treg,cd8,un],['Treg_candidate','CD4_other','CD8_candidate','T_unresolved'],default='Outside_T_gate');mac=o.Mac_RNA.to_numpy(bool);pos=mac&(o.PGAM5_raw.to_numpy()>0);high=o.source_group.eq('PGAM5_high').to_numpy();assert np.array_equal(high,pos&(o.PGAM5_log1p_CP10k.to_numpy()>=o.loc[pos,'PGAM5_log1p_CP10k'].quantile(.75)))
  Y=np.column_stack([np.divide(G@masks[k].astype(float),degree,out=np.zeros(n),where=degree>0) for k in TYPES]);rng=np.random.default_rng(238264751+si);rows,info=f.permutation_results(o,pos,high,Y,TYPES,'T_subtype_coverage',rng,sample)
  if info['movable_sources']==0:
   for row in rows:row['spatial_conditional_p']=np.nan;row['status']='not_testable: no exchangeable labels'
  stats+=rows;support.append(info)
  n8=G@cd8.astype(float);valid=n8>0;programs=['Exhaustion','Activation','Effector','Naive_memory'];Ys=np.column_stack([np.divide(G@(o['score_'+k].to_numpy()*cd8),n8,out=np.zeros(n),where=n8>0) for k in programs]);rows,info=f.permutation_results(o,pos&valid,high,Ys,programs,'CD8_neighbor_programs',rng,sample)
  if info['movable_sources']==0:
   for row in rows:row['spatial_conditional_p']=np.nan;row['status']='not_testable: no exchangeable labels'
  stats+=rows;support.append(info)
  for j,k in enumerate(TYPES):
   nt=G@masks[k].astype(float);o['neighbor_fraction_'+k]=Y[:,j]
   for group in GROUPS:
    use=o.source_group.eq(group).to_numpy();summary.append(dict(sample=sample,subtype=k,source_group=group,n_targets=int(masks[k].sum()),n_sources=int(use.sum()),n_sources_with_target=int((use&(nt>0)).sum()),n_edges=int(nt[use].sum()),mean_neighbor_fraction=float(Y[use,j].mean())))
  for j,k in enumerate(programs):o['CD8_neighbor_'+k]=np.where(valid,Ys[:,j],np.nan)
  strict=cd8&(count[:,genes.index('CD8B')]>0);broad=t&(count[:,genes.index('CD8A')]>0);sensreg=treg&(o.detected_Treg_program.to_numpy()>=3)
  for k,mask in [('CD8A_only_broad',broad),('CD8A_CD8B_both_strict',strict),('Treg_3_markers_strict',sensreg)]:
   for group in GROUPS:
    use=o.source_group.eq(group).to_numpy();nt=G@mask.astype(float);sens.append(dict(sample=sample,gate=k,source_group=group,n_targets=int(mask.sum()),n_edges=int(nt[use].sum())))
  noTG=cd8&(o.score_Exhaustion_sensitivity>=o.loc[cd8,'score_Exhaustion_sensitivity'].quantile(.75))&(o.score_Exhaustion_sensitivity>0)&(o.detected_Exhaustion_sensitivity>=2);o['gate_CD8_exhaustion_sensitivity']=noTG
  ri,ci=G.nonzero();sel=mac[ri]&t[ci];ri=ri[sel];ci=ci[sel];e=pd.DataFrame({'source_barcode':o.index[ri],'target_barcode':o.index[ci],'source_group':o.source_group.iloc[ri].to_numpy(),'target_identity':o.T_identity.iloc[ci].to_numpy()})
  for k in TYPES:e['target_'+k]=masks[k][ci]
  for l,r in PAIRS:
   rec=cp[:,[genes.index(g) for g in r.split('_')]].min(1);value=np.sqrt(cp[ri,genes.index(l)]*rec[ci]);e[l+'--'+r]=value
   for k in TYPES:
    target=masks[k][ci];nt=G@masks[k].astype(float);agg=np.bincount(ri[target],weights=value[target],minlength=n);mean=np.divide(agg,nt,out=np.full(n,np.nan),where=nt>0)
    for group in GROUPS:
     sel=e.source_group.eq(group).to_numpy()&target;v=value[sel];use=o.source_group.eq(group).to_numpy()&(nt>0);lrs.append(dict(sample=sample,subtype=k,source_group=group,ligand=l,receptor=r,n_edges=len(v),n_sources=len(np.unique(ri[sel])),n_targets=len(np.unique(ci[sel])),codetected_edges=int((v>0).sum()),codetection_pct=100*(v>0).mean() if len(v) else np.nan,mean_source_compatibility=float(mean[use].mean()) if use.any() else np.nan))
  o.to_csv(R/f'{sample}_T_subtypes.csv.gz',index_label='barcode',float_format='%.17g');e.to_csv(R/f'{sample}_T_LR_edges.csv.gz',index=False)
  for path in [bp,BASE/f'{sample}_spatial_edges.csv.gz',hp]:manifest.append(dict(path=str(path),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
  print(sample,'CD8',cd8.sum(),'Treg',treg.sum(),'High',high.sum(),flush=True)
 stats=pd.DataFrame(stats);assert len(stats)==91;stats['BH_q_91_planned']=false_discovery_control(stats.spatial_conditional_p.fillna(1).to_numpy());stats.loc[stats.spatial_conditional_p.isna(),'BH_q_91_planned']=np.nan;stats.to_csv(R/'spatial_T_subtype_tests.csv',index=False);pd.DataFrame(summary).to_csv(R/'T_subtype_neighbors.csv',index=False);pd.DataFrame(lrs).to_csv(R/'LR_by_T_subtype.csv',index=False);pd.DataFrame(sens).to_csv(R/'gate_sensitivity.csv',index=False);pd.DataFrame(thresholds).to_csv(R/'CD8_state_thresholds.csv',index=False);pd.DataFrame(support).to_csv(R/'permutation_support.csv',index=False);pd.DataFrame(coverage).to_csv(R/'marker_coverage.csv',index=False);(R/'input_manifest.json').write_text(json.dumps(manifest,indent=2));(R/'analysis_config.json').write_text(json.dumps(dict(PGAM5_Q=.75,source='frozen original baseline',panels=PANELS,states=STATES,types=TYPES,states_overlap=True,Treg_FOXP3_required=True,statistics='7*(9 coverage+4 continuous CD8 programs)=91 jointly corrected planned tests',seed=238264751,permutations=4999,LR='22 resource-defined pairs, descriptive only; no LIANA algorithm',minimum_sources_per_group=3),indent=2));print('Significant corrected tests:',int(stats.BH_q_91_planned.lt(.05).sum()))
if __name__=='__main__':run()
