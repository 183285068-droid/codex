"""Requested five-group focus: four overlapping CD8 programs and separate Treg."""
from pathlib import Path
import json,hashlib,importlib.util
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.stats import false_discovery_control
R=Path(__file__).resolve().parent;BASE=R.parent/'GSE281759_spatial_HCC_PGAM5';FINE=R.parent/'GSE281759_PGAM5_CD8_fine_programs';EXT=R.parent/'GSE281759_PGAM5_Tcell_extended';OLD=R.parent/'GSE281759_PGAM5_immune_subgroup_neighbors';SAMPLES=json.loads((BASE/'sample_scope.json').read_text());SEED=281759502
MAP={'CD8_exhaustion':['Exhaustion_supported','Progenitor_exhaustion_like','Terminal_exhaustion_like'],'CD8_activation':['Early_activation'],'CD8_memory':['Central_memory_like','GZMK_memory_like','Resident_memory_like'],'CD8_effector':['Cytotoxic_effector','Terminal_effector_like']};K=list(MAP)+['Treg_supported']
spec=importlib.util.spec_from_file_location('prior',OLD/'analyze.py');helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper);helper.ROOT=R;helper.TARGETS=K
def dense(x):return np.asarray(x.toarray() if sparse.issparse(x) else x)
def run():
 des=[];stats=[];lr=[];edges=[];zones=[];sensitivity=[];comb=[];refs=[];checks=[]
 pairs=pd.read_csv(OLD/'candidate_LR_resource.csv');pairs.to_csv(R/'candidate_LR_resource.csv',index=False)
 config=dict(CD8_union_rules=MAP,Treg='fixed T-lineage gate, FOXP3 positive, >=2 of FOXP3 IL2RA CTLA4 IKZF2 TIGIT, score>=T-gate Q75 and >0',activation_does_not_use_cycling_or_IFN_gate=True,memory_union_does_not_use_naive_gate=True,overlapping=True,seed=SEED,permutations=4999,primary_tests=80,joint_tests=240,primary_only_tests=75)
 (R/'analysis_config.json').write_text(json.dumps(config,indent=2))
 for item in SAMPLES:
  name=item['sample'];paths=[FINE/f'{name}_CD8_fine_markers.h5ad',EXT/f'{name}_supplemental_T_markers.h5ad']
  for p in paths:refs.append(dict(path=str(p.relative_to(R.parent)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
  a,e=[ad.read_h5ad(p) for p in paths];e=e[a.obs_names].copy();o=a.obs.copy();check=np.array_equal(o.total_counts,e.obs.total_counts);assert check;checks.append(name+' full-library and barcode alignment')
  shared=a.var_names.intersection(e.var_names);assert np.array_equal(dense(a[:,shared].layers['counts']),dense(e[:,shared].layers['counts']));checks.append(name+' shared raw genes agree')
  missing=e.var_names.difference(a.var_names);genes=a.var_names.to_list()+missing.to_list();raw=np.column_stack([dense(a.layers['counts']),dense(e[:,missing].layers['counts'])]);x=np.column_stack([dense(a.X),dense(e[:,missing].X)]);assert np.allclose(x,np.log1p(raw/o.total_counts.to_numpy()[:,None]*1e4));checks.append(name+' full-library normalization')
  gates={k:o[['gate_'+v for v in values]].any(axis=1).to_numpy() for k,values in MAP.items()};gates['Treg_supported']=e.obs.gate_Treg_supported.to_numpy(bool);cd8=o.gate_All_CD8_reference.to_numpy(bool)
  for k in MAP:assert not (gates[k]&~cd8).any()
  W=a.obsp['spatial_connectivities'];Y=np.column_stack([np.asarray(W@gates[k].astype(int)).ravel()>0 for k in K]);mac=o.Mac_RNA_enriched.to_numpy(bool)&(np.asarray(W.sum(1)).ravel()>0);positive=o.PGAM5_raw.to_numpy()>0;high=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy();rng=np.random.default_rng(SEED+int(item['gsm'][3:]))
  for kind,eligible in [('high_vs_other_detected',mac&positive),('high_vs_all_other_mac',mac),('nonT_high_vs_other_detected',mac&positive&~o.T_RNA_enriched.to_numpy(bool))]:stats.extend(helper.evaluate(o,Y,eligible,item,kind,rng))
  nh=np.asarray(W@(mac&high).astype(int)).ravel()>0;nd=np.asarray(W@(mac&positive&~high).astype(int)).ravel()>0;o['near_high_source']=nh;o['near_other_detected_source']=nd
  for j,k in enumerate(K):
   o['gate_'+k]=gates[k];o['neighbor_'+k]=Y[:,j]
   for group,mask in [('PGAM5_high',mac&high),('Other_detected',mac&positive&~high),('Undetected',mac&~positive)]:
    n=int(mask.sum());count=int(Y[mask,j].sum());manual=sum(any(gates[k][t] for t in W.indices[W.indptr[s]:W.indptr[s+1]]) for s in np.flatnonzero(mask));assert manual==count;checks.append(name+'/'+group+'/'+k+' explicit neighbors');des.append(dict(**item,target=k,source_group=group,n_sources=n,with_neighbor=count,neighbor_pct=100*count/n if n else np.nan,n_targets=int(gates[k].sum())))
   for zone,zmask in [('High_near_including_shared',nh),('Other_only_near',~nh&nd),('Neither_high_nor_other_near',~nh&~nd)]:zones.append(dict(**item,target=k,zone=zone,n_targets=int(gates[k].sum()),n_in_zone=int((gates[k]&zmask).sum())))
  for label,mask in [('Treg_CD4_supported',gates['Treg_supported']&e.obs.gate_CD4_related.to_numpy(bool)),('Treg_CD8_overlap',gates['Treg_supported']&cd8)]:
   near=np.asarray(W@mask.astype(int)).ravel()>0;sensitivity.append(dict(**item,label=label,n_targets=int(mask.sum()),high_near_unique_targets=int((mask&nh).sum()),high_sources_with_neighbor=int(near[high&mac].sum())))
  sig=np.array(['+'.join(k for k in MAP if gates[k][i]) or 'No_four_program_supported' for i in range(len(o))],dtype=object);o['CD8_four_program_signature']=np.where(cd8,sig,'Outside_CD8_gate')
  for label in sorted(set(sig[cd8])):
   mask=cd8&(sig==label);comb.append(dict(**item,signature=label,n_CD8_targets=int(mask.sum()),n_high_near_targets=int((mask&nh).sum())))
  edge=W.tocoo()
  for k in K:
   use=high[edge.row]&gates[k][edge.col];src,tgt=edge.row[use],edge.col[use];edges.extend(dict(**item,target=k,source_barcode=o.index[s],target_barcode=o.index[t]) for s,t in zip(src,tgt))
   for p in pairs.itertuples():
    li=genes.index(p.ligand);ri=[genes.index(g) for g in p.receptor_complex.split('_')];support=(raw[src,li]>0)&(raw[np.ix_(tgt,ri)]>0).all(1);lr.append(dict(**item,target=k,pair=p.ligand+'--'+p.receptor_complex,eligible_edges=len(src),unique_sources=len(np.unique(src)),unique_targets=len(np.unique(tgt)),supported_edges=int(support.sum()),supported_sources=len(np.unique(src[support])),supported_targets=len(np.unique(tgt[support])),codetection_pct=100*support.mean() if len(src) else np.nan))
  o.to_csv(R/f'{name}_annotations.csv.gz',float_format='%.17g');out=ad.AnnData(sparse.csr_matrix(x),obs=o,var=pd.DataFrame(index=genes));out.layers['counts']=sparse.csr_matrix(raw);out.obsp['spatial_connectivities']=W.copy();out.obsm['spatial']=a.obsm['spatial'].copy();out.write_h5ad(R/f'{name}_focused_T_markers.h5ad',compression='gzip')
  print(name,{k:int(gates[k].sum()) for k in K},flush=True)
 dd=pd.DataFrame(des);dd.to_csv(R/'neighborhood_counts.csv',index=False);pd.DataFrame(lr).to_csv(R/'candidate_LR_codetection.csv',index=False);pd.DataFrame(edges).to_csv(R/'high_source_T_directed_edges.csv.gz',index=False);pd.DataFrame(zones).to_csv(R/'unique_target_zones.csv',index=False);pd.DataFrame(sensitivity).to_csv(R/'Treg_identity_sensitivity.csv',index=False);pd.DataFrame(comb).to_csv(R/'CD8_four_program_signatures.csv',index=False)
 st=pd.DataFrame(stats);st['joint_q_240']=false_discovery_control(st.p.fillna(1).to_numpy());st['primary_q_80']=np.nan;sel=st.comparison.eq('high_vs_other_detected');st.loc[sel,'primary_q_80']=false_discovery_control(st.loc[sel,'p'].fillna(1).to_numpy());st['primary_only_q_75']=np.nan;sel2=sel&st.tumor_status.eq('primary');st.loc[sel2,'primary_only_q_75']=false_discovery_control(st.loc[sel2,'p'].fillna(1).to_numpy());st.to_csv(R/'spatial_statistics.csv',index=False)
 p=dd.groupby(['patient','source_group','target'],as_index=False).agg(n_sources=('n_sources','sum'),with_neighbor=('with_neighbor','sum'),n_targets=('n_targets','sum'));p['neighbor_pct']=100*p.with_neighbor/p.n_sources;p.to_csv(R/'patient_descriptive.csv',index=False);st[st.status.eq('tested')].groupby(['patient','comparison','target'],as_index=False).agg(median_delta_pp=('delta_pp','median'),n_testable_sections=('sample','size')).to_csv(R/'patient_equal_section_effects.csv',index=False)
 h=dd[dd.source_group.eq('PGAM5_high')].groupby('target')[['n_sources','with_neighbor','n_targets']].sum();h['coverage_pct']=100*h.with_neighbor/h.n_sources;h.to_csv(R/'pooled_descriptive.csv');(R/'reference_manifest.json').write_text(json.dumps(refs,indent=2));(R/'audit_results.json').write_text(json.dumps(dict(passed=len(checks),checks=checks),indent=2));print(h.to_string());print('min primary q',st.primary_q_80.min(),'CHECKS',len(checks))
if __name__=='__main__':run()
