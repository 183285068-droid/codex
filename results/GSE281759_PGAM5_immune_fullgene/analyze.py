"""Broad overlapping immune-RNA programs from verified full-gene Visium counts."""
from pathlib import Path
import os,json,gzip,hashlib,importlib.util
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np,pandas as pd,anndata as ad
from scipy import io,sparse
from scipy.stats import false_discovery_control
R=Path(__file__).resolve().parent;BASE=R.parent/'GSE281759_spatial_HCC_PGAM5';OLD=R.parent/'GSE281759_PGAM5_immune_subgroup_neighbors';EXT=R.parent/'GSE281759_PGAM5_Tcell_extended';SAMPLES=json.loads((BASE/'sample_scope.json').read_text());B=4999;SEED=281759402
spec=importlib.util.spec_from_file_location('previous',OLD/'analyze.py');helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
PANELS={
 'All_T':helper.PANELS['All_T'],
 'CD4_related':['CD4','IL7R','LTB'], 'CD8_supported':['CD8A','CD8B'],
 'Treg_supported':['FOXP3','IL2RA','CTLA4','IKZF2','TIGIT'],
 'Naive_T_like':['CCR7','SELL','TCF7','LEF1','MAL'],
 'Central_memory_T_like':['CCR7','IL7R','LTB','TCF7'],
 'GZMK_memory_T_like':['GZMK','CCL5','IL7R'],
 'Exhausted_T':helper.PANELS['Exhausted_T'], 'Activated_T':helper.PANELS['Activated_T'], 'Cytotoxic_T':helper.PANELS['Cytotoxic_T'],
 'NK_related':['NCR1','NCR3','NCAM1','KLRD1','NKG7','GNLY'],
 'Cytotoxic_NK_like':['FGFBP2','FCGR3A','PRF1','GNLY','NKG7','KLRD1'],
 'XCL_NK_like':['XCL1','XCL2','KLRC1','NCAM1','KLRD1'],
 'B_related':helper.PANELS['B_related'],
 'Naive_B_like':['TCL1A','IGHD','IL4R','FCER2','BANK1'],
 'Memory_B_like':['CD27','TNFRSF13B','CD37','AIM2','MS4A1','CD79A'],
 'Plasma_related':helper.PANELS['Plasma_related'],
 'IgA_plasma_related':['IGHA1','IGHA2','JCHAIN','MZB1'],
 'IgG_plasma_related':['IGHG1','IGHG2','IGHG3','IGHG4','MZB1'],
 'IgM_plasma_related':['IGHM','JCHAIN','MZB1'],
 'cDC1_related':helper.PANELS['cDC1_related'], 'cDC2_related':helper.PANELS['cDC2_related'],
 'pDC_related':['TCF4','IL3RA','CLEC4C','IRF7','LILRA4','GZMB'],
 'Classical_monocyte_related':['FCN1','CD14','VCAN','S100A8','S100A9','LST1'],
 'FCGR3A_myeloid_related':['FCGR3A','MS4A7','LILRB1','LST1','CTSS'],
 'Neutrophil_related':['FCGR3B','CSF3R','CXCR2','FPR1','S100A8','S100A9','S100A12'],
 'Mast_related':['KIT','TPSAB1','TPSB2','CPA3','MS4A2','HDC'],
 'Mac_RNA_reference':['C1QA','C1QB','C1QC','CD68','CSF1R','MERTK','MSR1'],
 'SPP1_mac_related':['SPP1','TREM2','APOE','LPL','GPNMB','LGALS3'],
 'TREM2_mac_related':['TREM2','APOE','LPL','GPNMB','CST7']}
K=list(PANELS);helper.ROOT=R;helper.TARGETS=K
NEW_RULES={
 'NK_related':dict(parent='all_QC',minimum=3,anchor='NCR1 OR NCR3 OR NCAM1'),
 'Cytotoxic_NK_like':dict(parent='NK_related',minimum=3,anchor='PRF1 OR GNLY'),
 'XCL_NK_like':dict(parent='NK_related',minimum=2,anchor='XCL1 OR XCL2'),
 'Naive_B_like':dict(parent='B_related',minimum=2,anchor='TCL1A OR IGHD'),
 'Memory_B_like':dict(parent='B_related',minimum=3,anchor='CD27 AND TNFRSF13B'),
 'IgA_plasma_related':dict(parent='Plasma_related',minimum=2,anchor='IGHA1 OR IGHA2'),
 'IgG_plasma_related':dict(parent='Plasma_related',minimum=2,anchor='IGHG1 OR IGHG2 OR IGHG3 OR IGHG4'),
 'IgM_plasma_related':dict(parent='Plasma_related',minimum=2,anchor='IGHM'),
 'pDC_related':dict(parent='all_QC',minimum=3,anchor='(CLEC4C OR LILRA4) AND (TCF4 OR IL3RA)'),
 'Classical_monocyte_related':dict(parent='all_QC',minimum=3,anchor='(FCN1 OR CD14) AND Pan-myeloid >= QC median'),
 'FCGR3A_myeloid_related':dict(parent='all_QC',minimum=3,anchor='FCGR3A AND (MS4A7 OR LILRB1) AND Pan-myeloid >= QC median'),
 'Neutrophil_related':dict(parent='all_QC',minimum=3,anchor='FCGR3B OR CSF3R OR CXCR2'),
 'Mast_related':dict(parent='all_QC',minimum=3,anchor='(TPSAB1 OR TPSB2) AND (CPA3 OR KIT)'),
 'SPP1_mac_related':dict(parent='Mac_RNA_reference',minimum=3,anchor='SPP1'),
 'TREM2_mac_related':dict(parent='Mac_RNA_reference',minimum=3,anchor='TREM2')}
EXT_MAP={'CD4_related':'CD4_related','CD8_supported':'CD8_supported','Treg_supported':'Treg_supported','Naive_T_like':'Naive_like','Central_memory_T_like':'Central_memory_like','GZMK_memory_T_like':'Effector_memory_like'}
OLD_KEYS=['All_T','Exhausted_T','Activated_T','Cytotoxic_T','B_related','Plasma_related','cDC1_related','cDC2_related']
def dense(x):return np.asarray(x.toarray() if sparse.issparse(x) else x)
def locate(name):
 p=list(Path(os.environ.get('GSE281759_INPUT_DIR','/workspace/scratch/GSE281759_uploaded/inputs')).glob(name))+list(Path('/workspace/attachments').rglob(name))
 if not p:raise FileNotFoundError(name)
 return p[0]
def anchor(k,d,o):
 if k=='NK_related':return d('NCR1')|d('NCR3')|d('NCAM1')
 if k=='Cytotoxic_NK_like':return d('PRF1')|d('GNLY')
 if k=='XCL_NK_like':return d('XCL1')|d('XCL2')
 if k=='Naive_B_like':return d('TCL1A')|d('IGHD')
 if k=='Memory_B_like':return d('CD27')&d('TNFRSF13B')
 if k=='IgA_plasma_related':return d('IGHA1')|d('IGHA2')
 if k=='IgG_plasma_related':return d('IGHG1')|d('IGHG2')|d('IGHG3')|d('IGHG4')
 if k=='IgM_plasma_related':return d('IGHM')
 if k=='pDC_related':return (d('CLEC4C')|d('LILRA4'))&(d('TCF4')|d('IL3RA'))
 if k=='Classical_monocyte_related':return (d('FCN1')|d('CD14'))&(o.Pan_myeloid.to_numpy()>=o.Pan_myeloid.median())
 if k=='FCGR3A_myeloid_related':return d('FCGR3A')&(d('MS4A7')|d('LILRB1'))&(o.Pan_myeloid.to_numpy()>=o.Pan_myeloid.median())
 if k=='Neutrophil_related':return d('FCGR3B')|d('CSF3R')|d('CXCR2')
 if k=='Mast_related':return (d('TPSAB1')|d('TPSB2'))&(d('CPA3')|d('KIT'))
 if k=='SPP1_mac_related':return d('SPP1')
 if k=='TREM2_mac_related':return d('TREM2')
def run():
 config=dict(targets=PANELS,new_rules=NEW_RULES,reused_T_mapping=EXT_MAP,reused_original_gates=OLD_KEYS,quantile=.75,quantile_population='new target parent gate or all QC per new_rules',score='mean log1p(CP10k), full-library denominator, zeros included',seed=SEED,permutations=B,primary_tests=16*len(K),joint_tests=16*len(K)*3,primary_only_tests=15*len(K),labels_overlap=True,single_cell_identity=False)
 (R/'analysis_config.json').write_text(json.dumps(config,indent=2));manifest={v['name']:v for v in json.loads((BASE/'source_manifest.json').read_text())};inputs=[];refs=[];checks=[];des=[];stats=[];thresholds=[];markers=[];lr=[];edges=[];target_zones=[];macro_nonhigh=[]
 pairs=pd.read_csv(OLD/'candidate_LR_resource.csv');pairs.to_csv(R/'candidate_LR_resource.csv',index=False);genes=list(dict.fromkeys(['PGAM5']+[g for gs in PANELS.values() for g in gs]+pairs.ligand.to_list()+[g for c in pairs.receptor_complex for g in c.split('_')]))
 for item in SAMPLES:
  s=item['sample'];paths={}
  for suffix in ['matrix.mtx.gz','features.tsv.gz','barcodes.tsv.gz']:
   name=f'{item["gsm"]}_{s}_{suffix}';p=locate(name);v=manifest[name];assert p.stat().st_size==v['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==v['sha256'];inputs.append(dict(v,local_path=str(p)));checks.append(s+'/'+suffix+' SHA');paths[suffix]=p
  for p in [BASE/f'{s}_spatial_marker_subset.h5ad',OLD/f'{s}_annotations.csv.gz',EXT/f'{s}_supplemental_T_markers.h5ad']:refs.append(dict(path=str(p.relative_to(R.parent)),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size))
  f=pd.read_csv(paths['features.tsv.gz'],sep='\t',header=None);bc=pd.Index(pd.read_csv(paths['barcodes.tsv.gz'],header=None)[0].astype(str));assert bc.is_unique and f[0].is_unique
  with gzip.open(paths['matrix.mtx.gz'],'rb') as stream:full=io.mmread(stream,spmatrix=True).T.tocsr()
  assert full.shape==(len(bc),len(f)) and np.all(full.data>=0) and np.all(full.data==np.floor(full.data));checks.append(s+': matrix integer counts/dimensions')
  a=ad.read_h5ad(BASE/f'{s}_spatial_marker_subset.h5ad');ix=bc.get_indexer(a.obs_names);assert (ix>=0).all();full=full[ix];total=np.asarray(full.sum(1)).ravel();assert np.array_equal(total,a.obs.total_counts.to_numpy());oldix=pd.Index(f[0]).get_indexer(a.var.gene_id);assert np.array_equal(dense(full[:,oldix]),dense(a.layers['counts']));checks.append(s+': barcode/full UMI/old markers match')
  for g in genes:assert int((f[1]==g).sum())==1,g
  symbols=ad.utils.make_index_unique(pd.Index(f[1]));gi=symbols.get_indexer(genes);raw=dense(full[:,gi]);x=np.log1p(raw/total[:,None]*10000);o=a.obs.copy();previous=pd.read_csv(OLD/f'{s}_annotations.csv.gz',index_col=0).loc[o.index];ext=ad.read_h5ad(EXT/f'{s}_supplemental_T_markers.h5ad').obs.loc[o.index]
  gates={k:previous['gate_'+k].to_numpy(bool) for k in OLD_KEYS};gates.update({k:ext['gate_'+v].to_numpy(bool) for k,v in EXT_MAP.items()});gates['Mac_RNA_reference']=o.Mac_RNA_enriched.to_numpy(bool)
  det=lambda g:raw[:,genes.index(g)]>0
  for k in K:
   gs=PANELS[k];j=[genes.index(g) for g in gs];score=x[:,j].mean(1);n=(raw[:,j]>0).sum(1)
   cut=np.nan;parent_n=np.nan
   if k in NEW_RULES:
    rule=NEW_RULES[k];parent=np.ones(len(o),bool) if rule['parent']=='all_QC' else gates[rule['parent']];parent_n=int(parent.sum());cut=float(np.quantile(score[parent],.75)) if parent.any() else np.nan;gates[k]=parent&(score>=cut)&(score>0)&(n>=rule['minimum'])&anchor(k,det,o)
   gate=gates[k];o['gate_'+k]=gate;o['score_'+k]=score;o['detected_'+k]=n;thresholds.append(dict(**item,target=k,threshold=cut,parent_n=parent_n,n_target=int(gate.sum()),rule_type='new' if k in NEW_RULES else 'reused'))
   for g in gs:markers.append(dict(**item,target=k,gene=g,qc_detected=int(det(g).sum()),target_detected=int(det(g)[gate].sum()),n_qc=len(o),n_target=int(gate.sum())))
  o=o.copy()
  W=a.obsp['spatial_connectivities'];Y=np.column_stack([np.asarray(W@gates[k].astype(int)).ravel()>0 for k in K]);mac=o.Mac_RNA_enriched.to_numpy(bool)&(np.asarray(W.sum(1)).ravel()>0);positive=o.PGAM5_raw.to_numpy()>0;high=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy();rng=np.random.default_rng(SEED+int(item['gsm'][3:]))
  for kind,eligible in [('high_vs_other_detected',mac&positive),('high_vs_all_other_mac',mac),('nonT_high_vs_other_detected',mac&positive&~o.T_RNA_enriched.to_numpy(bool))]:stats.extend(helper.evaluate(o,Y,eligible,item,kind,rng))
  near_high=np.asarray(W@(mac&high).astype(int)).ravel()>0;near_other=np.asarray(W@(mac&positive&~high).astype(int)).ravel()>0
  for group,mask in [('PGAM5_high',mac&high),('Other_detected',mac&positive&~high),('Undetected',mac&~positive)]:
   for j,k in enumerate(K):
    n=int(mask.sum());v=int(Y[mask,j].sum());des.append(dict(**item,source_group=group,target=k,n_sources=n,with_neighbor=v,neighbor_pct=100*v/n if n else np.nan,n_target_spots=int(gates[k].sum())))
    count=sum(any(gates[k][z] for z in W.indices[W.indptr[src]:W.indptr[src+1]]) for src in np.flatnonzero(mask));assert count==v;checks.append(s+'/'+group+'/'+k+' explicit neighbors')
  edge=W.tocoo()
  for k in K:
   use=high[edge.row]&gates[k][edge.col];src,tgt=edge.row[use],edge.col[use];edges.extend(dict(**item,target=k,source_barcode=o.index[z],target_barcode=o.index[t]) for z,t in zip(src,tgt))
   for zone,zmask in [('High_near_including_shared',near_high),('Other_only_near',~near_high&near_other),('Neither_high_nor_other_near',~near_high&~near_other)]:target_zones.append(dict(**item,target=k,zone=zone,n_target=int(gates[k].sum()),n_targets_in_zone=int((gates[k]&zmask).sum())))
   for p in pairs.itertuples():
    li=genes.index(p.ligand);ri=[genes.index(g) for g in p.receptor_complex.split('_')];support=(raw[src,li]>0)&(raw[np.ix_(tgt,ri)]>0).all(1);lr.append(dict(**item,target=k,pair=p.ligand+'--'+p.receptor_complex,eligible_edges=len(src),unique_sources=len(np.unique(src)),unique_targets=len(np.unique(tgt)),supported_edges=int(support.sum()),supported_sources=len(np.unique(src[support])),supported_targets=len(np.unique(tgt[support])),codetection_pct=100*support.mean() if len(src) else np.nan))
  for k in ['Mac_RNA_reference','SPP1_mac_related','TREM2_mac_related']:
   mask=gates[k]&~high;near=np.asarray(W@mask.astype(int)).ravel()>0;macro_nonhigh.append(dict(**item,target=k,n_targets_without_high=int(mask.sum()),n_high_sources=int((mac&high).sum()),high_sources_with_nonhigh_neighbor=int(near[mac&high].sum())))
  for j,k in enumerate(K):o['neighbor_'+k]=Y[:,j]
  o['near_high_source']=near_high;o['near_other_detected_source']=near_other;o.to_csv(R/f'{s}_annotations.csv.gz',float_format='%.17g');extra=ad.AnnData(sparse.csr_matrix(x),obs=o,var=pd.DataFrame(index=genes));extra.layers['counts']=sparse.csr_matrix(raw);extra.obsp['spatial_connectivities']=W.copy();extra.obsm['spatial']=a.obsm['spatial'].copy();extra.write_h5ad(R/f'{s}_immune_markers.h5ad',compression='gzip');print(s,{k:int(gates[k].sum()) for k in ['NK_related','Naive_B_like','Memory_B_like','pDC_related','Neutrophil_related','Mast_related']},flush=True)
 dd=pd.DataFrame(des);dd.to_csv(R/'neighborhood_counts.csv',index=False);pd.DataFrame(thresholds).to_csv(R/'thresholds.csv',index=False);pd.DataFrame(markers).to_csv(R/'marker_detection.csv',index=False);pd.DataFrame(lr).to_csv(R/'candidate_LR_codetection.csv',index=False);pd.DataFrame(edges).to_csv(R/'high_source_immune_directed_edges.csv.gz',index=False);pd.DataFrame(target_zones).to_csv(R/'unique_target_neighborhood_zones.csv',index=False);pd.DataFrame(macro_nonhigh).to_csv(R/'macrophage_targets_excluding_high_descriptive.csv',index=False)
 st=pd.DataFrame(stats);st['joint_q_1440']=false_discovery_control(st.p.fillna(1).to_numpy());st['primary_q_480']=np.nan;sel=st.comparison.eq('high_vs_other_detected');st.loc[sel,'primary_q_480']=false_discovery_control(st.loc[sel,'p'].fillna(1).to_numpy());st['primary_only_q_450']=np.nan;sel2=sel&st.tumor_status.eq('primary');st.loc[sel2,'primary_only_q_450']=false_discovery_control(st.loc[sel2,'p'].fillna(1).to_numpy());st.to_csv(R/'spatial_statistics.csv',index=False)
 patient=dd.groupby(['patient','source_group','target'],as_index=False).agg(n_sources=('n_sources','sum'),with_neighbor=('with_neighbor','sum'),n_target_spots=('n_target_spots','sum'));patient['neighbor_pct']=100*patient.with_neighbor/patient.n_sources;patient.to_csv(R/'patient_descriptive.csv',index=False);st[st.status.eq('tested')].groupby(['patient','comparison','target'],as_index=False).agg(median_delta_pp=('delta_pp','median'),n_testable_sections=('sample','size')).to_csv(R/'patient_equal_section_effects.csv',index=False)
 h=dd[dd.source_group.eq('PGAM5_high')].groupby('target')[['n_sources','with_neighbor','n_target_spots']].sum();h['coverage_pct']=100*h.with_neighbor/h.n_sources;h.to_csv(R/'pooled_descriptive.csv');(R/'input_manifest.json').write_text(json.dumps(inputs,indent=2));(R/'reference_manifest.json').write_text(json.dumps(refs,indent=2));(R/'audit_results.json').write_text(json.dumps(dict(passed=len(checks),checks=checks),indent=2));print(h.to_string());print('TARGETS',len(K),'min primary q',st.primary_q_480.min(),'CHECKS',len(checks))
if __name__=='__main__':run()
