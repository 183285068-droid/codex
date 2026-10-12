"""Spatial RNA compatibility analysis, not identified single-cell communication."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE245908_Tcell_mpl')
import numpy as np,pandas as pd,scanpy as sc,scipy.sparse as sp,json,hashlib
import itertools,math
from statsmodels.stats.multitest import multipletests
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent/'GSE245908_spatial_PGAM5';WORK=Path(os.environ.get('GSE245908_WORK','/workspace/scratch/GSE245908'))
SAMPLES=['CHC20','CHC23'];SEED=2459083;B=4999
PANELS={'T_lineage':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2'],
'Activation_early':['CD69','IL2RA','TNFRSF9','ICOS','CD40LG','IL2'],
'Effector_cytotoxic':['IFNG','PRF1','GZMB','GZMA','CCL5','TBX21'],
'Exhaustion_core':['PDCD1','LAG3','TIGIT','CTLA4','TOX','CXCL13'],
'Exhaustion_expanded':['PDCD1','LAG3','TIGIT','CTLA4','TOX','CXCL13','HAVCR2','ENTPD1'],
'Macrophage_core':['C1QA','C1QB','C1QC','CD68','CSF1R','MERTK','MSR1']}
PAIRS=[('CD274','PDCD1','checkpoint'),('PDCD1LG2','PDCD1','checkpoint'),('PVR','TIGIT','checkpoint'),('NECTIN2','TIGIT','checkpoint'),('LGALS9','HAVCR2','context-dependent'),('CD80','CTLA4','checkpoint'),('CD86','CTLA4','checkpoint'),('CD80','CD28','costimulation'),('CD86','CD28','costimulation'),('ICOSLG','ICOS','costimulation'),('TNFSF9','TNFRSF9','costimulation'),('TNFSF4','TNFRSF4','costimulation'),('TGFB1','TGFBR1_TGFBR2','context-dependent'),('IL10','IL10RA_IL10RB','context-dependent'),('CXCL9','CXCR3','chemotaxis'),('CXCL10','CXCR3','chemotaxis'),('CXCL11','CXCR3','chemotaxis'),('CCL3','CCR5','chemotaxis'),('CCL4','CCR5','chemotaxis'),('CCL5','CCR5','chemotaxis'),('SPP1','CD44','context-dependent')]
def genes(a,names,layer=None):
 x=a[:,names].layers[layer] if layer else a[:,names].X
 return x.toarray().astype(float)
def load_sample(s):
 obs=pd.read_csv(BASE/f'{s}_spot_annotations.csv',index_col=0)
 if (WORK/f'{s}_full_processed.h5ad').exists():a=sc.read_h5ad(WORK/f'{s}_full_processed.h5ad')
 else:
  gsm={'CHC20':'GSM7850822','CHC23':'GSM7850823'}[s];a=sc.read_10x_h5(BASE/'inputs'/s/f'{gsm}_{s}_filtered_feature_bc_matrix.h5');a.var_names_make_unique();a=a[obs.index].copy();a.layers['counts']=a.X.copy();sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
 a.obs=obs.loc[a.obs_names];return a
def graph(obs):
 ed=pd.read_csv(BASE/f'{obs["sample"].iloc[0]}_spatial_edges.csv.gz');ix={v:i for i,v in enumerate(obs.index)}
 return sp.csr_matrix((np.ones(len(ed)),([ix[x] for x in ed.barcode_from],[ix[x] for x in ed.barcode_to])),shape=(len(obs),len(obs)))
def perm_matrix(obs,eligible,high,rng):
 ix=np.flatnonzero(eligible);d=obs.iloc[ix];xy=obs[['array_col','array_row']].to_numpy();lo=xy.min(0);span=np.ptp(xy,axis=0)+1
 tile=np.minimum(((xy[ix]-lo)/span*3).astype(int),2);blocks=tile[:,0]+3*tile[:,1];strata=[]
 for tileid in sorted(set(blocks)):
  idx=np.flatnonzero(blocks==tileid);depth=d.iloc[idx].total_counts.to_numpy();groups=(depth>np.median(depth)).astype(int)
  for level in [0,1]:
   members=idx[groups==level]
   if len(members):strata.append(members)
 h=high[ix].astype(float);space=math.prod(math.comb(len(m),int(h[m].sum())) for m in strata);movable=sum(len(m) for m in strata if 0<h[m].sum()<len(m))
 if space<=B:
  choices=[list(itertools.combinations(m,int(h[m].sum()))) for m in strata];mat=np.zeros((space,len(h)))
  for row,parts in enumerate(itertools.product(*choices)):
   for part in parts:mat[row,list(part)]=1
  mode='exact_conditional_enumeration'
 else:
  mat=np.tile(h,(B,1));mode='Monte_Carlo_conditional_permutation'
 for members in strata:
  if h[members].sum()>0 and h[members].sum()<len(members):
   if mode=='Monte_Carlo_conditional_permutation':
    for k in range(B):mat[k,members]=rng.permutation(h[members])
 return ix,mat,{'n_eligible_sources':len(ix),'n_high':int(h.sum()),'n_other_detected':int(len(h)-h.sum()),'n_strata':len(strata),'movable_sources':movable,'distinct_null_assignments':len(np.unique(mat,axis=0)),'configuration_space':str(space),'null_rows':len(mat),'mode':mode}
def permutation_results(obs,eligible,high,Y,labels,kind,rng,s):
 ix,mat,info=perm_matrix(obs,eligible,high,rng);h=high[ix].astype(float);y=Y[ix];rows=[]
 if h.sum()<3 or (len(h)-h.sum())<3:
  for label in labels:rows.append({'sample':s,'kind':kind,'feature':label,'n_high':int(h.sum()),'n_other_detected':int(len(h)-h.sum()),'mean_high':np.nan,'mean_other_detected':np.nan,'delta_high_minus_other':np.nan,'spatial_conditional_p':np.nan,'status':'insufficient sources (<3 in a group)'})
 else:
  obsdelta=h@y/h.sum()-(1-h)@y/(len(h)-h.sum());nulldelta=mat@y/h.sum()-(1-mat)@y/(len(h)-h.sum())
  # Conditional permutation null may have a nonzero expected global delta.
  # Center on the exact conditional expectation, avoiding a zero-centered false null.
  expected_h=np.zeros(len(h));
  xy=obs[['array_col','array_row']].to_numpy();lo=xy.min(0);span=np.ptp(xy,axis=0)+1;tile=np.minimum(((xy[ix]-lo)/span*3).astype(int),2);blocks=tile[:,0]+3*tile[:,1]
  for block in set(blocks):
   members=np.flatnonzero(blocks==block);level=(obs.iloc[ix[members]].total_counts.to_numpy()>obs.iloc[ix[members]].total_counts.median()).astype(int)
   for lv in [0,1]:
    part=members[level==lv]
    if len(part):expected_h[part]=h[part].mean()
  nullcenter=expected_h@y/h.sum()-(1-expected_h)@y/(len(h)-h.sum())
  for j,label in enumerate(labels):
   exceed=(abs(nulldelta[:,j]-nullcenter[j])>=abs(obsdelta[j]-nullcenter[j])-1e-14).sum()
   p=exceed/len(mat) if info['mode']=='exact_conditional_enumeration' else (1+exceed)/(len(mat)+1)
   rows.append({'sample':s,'kind':kind,'feature':label,'n_high':int(h.sum()),'n_other_detected':int(len(h)-h.sum()),'mean_high':float(y[h>0,j].mean()),'mean_other_detected':float(y[h==0,j].mean()),'delta_high_minus_other':float(obsdelta[j]),'conditional_null_expected_delta':float(nullcenter[j]),'delta_minus_conditional_expectation':float(obsdelta[j]-nullcenter[j]),'spatial_conditional_p':float(p),'status':'all scores zero' if np.all(y[:,j]==0) else 'tested'})
  np.savez_compressed(ROOT/f'{s}_{kind}_permutation_audit.npz',source_indices=ix,high=h,outcomes=y,labels=np.array(labels),null_deltas=nulldelta,conditional_null_center=nullcenter,observed_deltas=obsdelta,mode=np.array(info['mode']))
 return rows,dict(sample=s,kind=kind,**info)
if __name__=='__main__':
 ROOT.mkdir(exist_ok=True);allstats=[];allsummary=[];edgesout=[];groupsum=[];coverage=[];perm_info=[];lr_detail=[];sens=[]
 resource=pd.read_csv(ROOT/'source_inputs/omni_resource.csv');selected=[]
 for l,r,c in PAIRS:
  d=resource[(resource.source_genesymbol==l)&(resource.target_genesymbol==r)];assert len(d)
  selected.append({'ligand':l,'receptor_complex':r,'context':c,'database_resources':';'.join(sorted(set(d.resource))),'database_snapshot':'LIANA wheel 1.10.0 omni_resource.csv'})
 pd.DataFrame(selected).to_csv(ROOT/'selected_LR_pairs.csv',index=False)
 (ROOT/'programs.json').write_text(json.dumps(PANELS,indent=2))
 for s in SAMPLES:
  print('Processing',s,flush=True);a=load_sample(s);o=a.obs.copy();W=graph(o);deg=np.asarray(W.sum(1)).ravel();wr=sp.diags(np.divide(1,deg,out=np.zeros_like(deg),where=deg>0))@W
  for k,gs in PANELS.items():
   available=[g for g in gs if g in a.var_names];coverage.append({'sample':s,'program':k,'n_requested':len(gs),'n_present':len(available),'genes_present':';'.join(available)});z=genes(a,available);o[k]=z.mean(1);o[k+'_detected_genes']=(z>0).sum(1)
  maccut=float(o.Macrophage_core.quantile(.75));mycut=float(o.Pan_myeloid.median());tcut=float(o.T_lineage.quantile(.75))
  mac=(o.Macrophage_core.to_numpy()>=maccut)&(o.Pan_myeloid.to_numpy()>=mycut)&(o.Macrophage_core_detected_genes.to_numpy()>=3)&(deg>0)
  t=(o.T_lineage.to_numpy()>=tcut)&(o.T_lineage_detected_genes.to_numpy()>=2)&(deg>0)
  pos=mac&(o.PGAM5_raw.to_numpy()>0);cut=float(o.loc[pos,'PGAM5_log1p_CP10k'].quantile(.75));high=pos&(o.PGAM5_log1p_CP10k.to_numpy()>=cut)
  group=np.where(mac,np.where(high,'PGAM5_high',np.where(pos,'Other_detected','Undetected')),'Outside_mac_gate');o['source_group']=group;o['T_enriched']=t;o['Mac_enriched']=mac
  states={}
  for k in ['Activation_early','Exhaustion_core','Effector_cytotoxic']:
   cutoff=float(o.loc[t,k].quantile(.75));mask=t&(o[k].to_numpy()>=cutoff)&(o[k].to_numpy()>0)&(o[k+'_detected_genes'].to_numpy()>=2);o[k+'_supported_high']=mask;states[k]=mask
  o['T_state']=np.select([t&states['Activation_early']&states['Exhaustion_core'],t&states['Exhaustion_core'],t&states['Activation_early'],t],['Both_programs','Exhaustion_supported','Activation_supported','T_other'],default='Outside_T_gate')
  nT=W@t.astype(float);o['T_neighbor_count']=nT;o['T_neighbor_fraction']=np.divide(nT,deg,out=np.zeros_like(nT),where=deg>0);valid=nT>0
  outcomes=[]
  for k in ['Activation_early','Effector_cytotoxic','Exhaustion_core','Exhaustion_expanded']:
   z=np.divide(W@(o[k].to_numpy()*t),nT,out=np.zeros_like(nT),where=nT>0);o['T_neighbor_'+k]=np.where(valid,z,np.nan);outcomes.append(z)
  rng=np.random.default_rng(SEED+(s=='CHC23'))
  r,info=permutation_results(o,pos,high,o.T_neighbor_fraction.to_numpy()[:,None],['T_neighbor_fraction'],'T_coverage',rng,s);allstats+=r;perm_info.append(info)
  Y=np.column_stack(outcomes[:3]);r,info=permutation_results(o,pos&valid,high,Y,['Activation_early','Effector_cytotoxic','Exhaustion_core'],'T_programs',rng,s);allstats+=r;perm_info.append(info)
  # Source-specific LR means across adjacent T-enriched targets; zero edges retained.
  ri,ci=W.nonzero();edges=pd.DataFrame({'source_index':ri,'target_index':ci});edges=edges[mac[ri]&t[ci]];ri=edges.source_index.to_numpy();ci=edges.target_index.to_numpy()
  edge_base=pd.DataFrame({'sample':s,'source_barcode':o.index[ri],'target_barcode':o.index[ci],'source_group':group[ri],'target_T_state':o.T_state.iloc[ci].to_numpy()});edgesout.append(edge_base)
  per_source_lr=[]
  for l,r,c in PAIRS:
   lig=genes(a,[l],layer='counts').ravel()/o.total_counts.to_numpy()*1e4;rec=np.min(genes(a,r.split('_'),layer='counts')/o.total_counts.to_numpy()[:,None]*1e4,axis=1);score=np.sqrt(lig[ri]*rec[ci]);aggregate=np.bincount(ri,weights=score,minlength=len(o));mean=np.divide(aggregate,nT,out=np.zeros(len(o)),where=nT>0);per_source_lr.append(mean)
   for g in ['PGAM5_high','Other_detected','Undetected']:
    for state in ['All_T','Exhaustion_supported','Activation_supported','Both_programs','T_other']:
     mask=group[ri]==g
     if state!='All_T':mask&=o.T_state.iloc[ci].to_numpy()==state
     vals=score[mask];lr_detail.append({'sample':s,'source_group':g,'receiver_T_state':state,'ligand':l,'receptor_complex':r,'n_edges':len(vals),'n_source_spots':len(np.unique(ri[mask])),'n_receiver_spots':len(np.unique(ci[mask])),'n_codetected_edges':int((vals>0).sum()),'codetection_pct':float((vals>0).mean()*100) if len(vals) else np.nan,'mean_edge_RNA_compatibility':float(vals.mean()) if len(vals) else np.nan})
   edge_base[l+'--'+r]=score
  LRY=np.column_stack(per_source_lr);labels=[l+'--'+r for l,r,_ in PAIRS]
  r,info=permutation_results(o,pos&valid,high,LRY,labels,'LR',rng,s);allstats+=r;perm_info.append(info)
  edge_base.to_csv(ROOT/f'{s}_directed_LR_edges.csv.gz',index=False)
  for g in ['PGAM5_high','Other_detected','Undetected']:
   mask=group==g
   for k in ['T_neighbor_fraction','T_neighbor_Activation_early','T_neighbor_Effector_cytotoxic','T_neighbor_Exhaustion_core','T_neighbor_Exhaustion_expanded']:
    vals=o.loc[mask,k].dropna();groupsum.append({'sample':s,'source_group':g,'outcome':k,'n_sources_with_outcome':len(vals),'mean':float(vals.mean()) if len(vals) else np.nan,'median':float(vals.median()) if len(vals) else np.nan})
  # Narrow/broad exhaustion and T gate >=3 markers are effect-size sensitivity only.
  strict=t&(o.T_lineage_detected_genes.to_numpy()>=3);ns=W@strict.astype(float)
  for receiver,mask in [('primary_T_gate',t),('strict_3_T_markers',strict)]:
   n=W@mask.astype(float)
   for k in ['Activation_early','Effector_cytotoxic','Exhaustion_core','Exhaustion_expanded']:
    z=np.divide(W@(o[k].to_numpy()*mask),n,out=np.zeros(len(o)),where=n>0)
    for g in ['PGAM5_high','Other_detected','Undetected']:
     keep=(group==g)&(n>0);sens.append({'sample':s,'receiver_gate':receiver,'program':k,'source_group':g,'n_sources':int(keep.sum()),'mean_neighbor_score':float(z[keep].mean()) if keep.any() else np.nan})
  allsummary.append({'sample':s,'QC_spots':len(o),'mac_enriched_spots':int(mac.sum()),'PGAM5_detected_mac_spots':int(pos.sum()),'PGAM5_high_mac_spots':int(high.sum()),'PGAM5_other_detected_mac_spots':int((pos&~high).sum()),'PGAM5_undetected_mac_spots':int((mac&~pos).sum()),'PGAM5_high_cutoff_log1p_CP10k':cut,'mac_program_Q75_cutoff':maccut,'myeloid_Q50_cutoff':mycut,'T_program_Q75_cutoff':tcut,'T_enriched_spots':int(t.sum()),'T_activation_supported':int(states['Activation_early'].sum()),'T_exhaustion_supported':int(states['Exhaustion_core'].sum()),'T_effector_supported':int(states['Effector_cytotoxic'].sum()),'directed_mac_T_edges':len(edge_base),'high_sources_with_T_neighbors':int((high&valid).sum())})
  o.to_csv(ROOT/f'{s}_spot_states.csv',index_label='barcode',float_format='%.17g')
  pd.DataFrame({'barcode':o.index[pos&valid],**{lab:LRY[pos&valid,j] for j,lab in enumerate(labels)}}).to_csv(ROOT/f'{s}_per_source_LR_scores.csv',index=False)
 df=pd.DataFrame(allstats);df['BH_q']=np.nan;valid=df.spatial_conditional_p.notna();df.loc[valid,'BH_q']=multipletests(df.loc[valid,'spatial_conditional_p'],method='fdr_bh')[1];df.to_csv(ROOT/'spatial_conditional_tests.csv',index=False)
 pd.DataFrame(allsummary).to_csv(ROOT/'sample_state_summary.csv',index=False);pd.DataFrame(groupsum).to_csv(ROOT/'source_group_neighbor_program_summary.csv',index=False);pd.DataFrame(lr_detail).to_csv(ROOT/'LR_by_receiver_state.csv',index=False);pd.DataFrame(coverage).to_csv(ROOT/'program_coverage.csv',index=False);pd.DataFrame(perm_info).to_csv(ROOT/'permutation_support.csv',index=False);pd.DataFrame(sens).to_csv(ROOT/'receiver_gate_and_exhaustion_panel_sensitivity.csv',index=False)
 manifest=[]
 for p in [*sorted((BASE/'inputs').rglob('*')),*[BASE/f'{s}_{suffix}' for s in SAMPLES for suffix in ['spot_annotations.csv','spatial_edges.csv.gz']],ROOT/'source_inputs/omni_resource.csv',ROOT/'source_inputs/base_requirements.json',*sorted((ROOT/'source_inputs').glob('*.soft'))]:
  if p.is_file():manifest.append({'path':str(p.relative_to(ROOT.parent)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
 (ROOT/'input_manifest.json').write_text(json.dumps(manifest,indent=2));print(pd.DataFrame(allsummary).to_string(index=False),flush=True);print(df[df.kind!='LR'].to_string(index=False),flush=True)
