"""Fine, overlapping CD8-supported spatial RNA programs; not single-cell types."""
from pathlib import Path
import os,json,hashlib,gzip,importlib.util
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse,io
from scipy.stats import false_discovery_control
R=Path(__file__).resolve().parent;BASE=R.parent/'GSE281759_spatial_HCC_PGAM5';EXT=R.parent/'GSE281759_PGAM5_Tcell_extended';OLD=R.parent/'GSE281759_PGAM5_immune_subgroup_neighbors'
SAMPLES=json.loads((BASE/'sample_scope.json').read_text());SEED=281759302;B=4999
PANELS={
 'Naive_like':['CCR7','SELL','TCF7','LEF1','MAL'],
 'Central_memory_like':['CCR7','IL7R','LTB','TCF7'],
 'GZMK_memory_like':['GZMK','IL7R','CCL5'],
 'Cytotoxic_effector':['NKG7','PRF1','GZMB','GNLY','CTSW','CCL5'],
 'Terminal_effector_like':['FGFBP2','CX3CR1','GZMH','PRF1','NKG7','KLRG1'],
 'Exhaustion_supported':['PDCD1','LAG3','TIGIT','HAVCR2','TOX','ENTPD1','CXCL13','CTLA4'],
 'Progenitor_exhaustion_like':['TCF7','SLAMF6','PDCD1','TOX'],
 'Terminal_exhaustion_like':['HAVCR2','ENTPD1','CXCL13','TOX','LAG3','PDCD1'],
 'Resident_memory_like':['CD69','ITGAE','CXCR6','ZNF683','ITGA1','RGS1'],
 'Early_activation':['CD69','IL2RA','TNFRSF9','ICOS','CD40LG','IL2'],
 'Cycling':['MKI67','TOP2A','STMN1','TYMS','PCNA'],
 'IFN_response':['STAT1','ISG15','IFIT1','IFIT3','MX1','OAS1','IFI6']}
K=list(PANELS);ALL=['All_CD8_reference']+K
MIN={k:(2 if k in ['GZMK_memory_like','Early_activation'] else 3) for k in K}
ANCHOR={'Naive_like':'(CCR7 OR SELL) AND (TCF7 OR LEF1)','Central_memory_like':'CCR7 AND IL7R','GZMK_memory_like':'GZMK AND (IL7R OR CCL5)','Cytotoxic_effector':'PRF1 OR GZMB','Terminal_effector_like':'FGFBP2 OR CX3CR1','Exhaustion_supported':'PDCD1 OR TOX','Progenitor_exhaustion_like':'TCF7 AND PDCD1','Terminal_exhaustion_like':'(HAVCR2 OR ENTPD1) AND (TOX OR PDCD1)','Resident_memory_like':'ITGAE OR ZNF683 OR ITGA1','Early_activation':'CD69 OR IL2RA OR TNFRSF9 OR ICOS','Cycling':'MKI67 OR TOP2A','IFN_response':'IFIT1 OR IFIT3 OR MX1'}
spec=importlib.util.spec_from_file_location('prior',OLD/'analyze.py');helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper);helper.ROOT=R;helper.TARGETS=K
def dense(x):return np.asarray(x.toarray() if sparse.issparse(x) else x)
def anchor_mask(k,d):
 if k=='Naive_like':return (d('CCR7')|d('SELL'))&(d('TCF7')|d('LEF1'))
 if k=='Central_memory_like':return d('CCR7')&d('IL7R')
 if k=='GZMK_memory_like':return d('GZMK')&(d('IL7R')|d('CCL5'))
 if k=='Cytotoxic_effector':return d('PRF1')|d('GZMB')
 if k=='Terminal_effector_like':return d('FGFBP2')|d('CX3CR1')
 if k=='Exhaustion_supported':return d('PDCD1')|d('TOX')
 if k=='Progenitor_exhaustion_like':return d('TCF7')&d('PDCD1')
 if k=='Terminal_exhaustion_like':return (d('HAVCR2')|d('ENTPD1'))&(d('TOX')|d('PDCD1'))
 if k=='Resident_memory_like':return d('ITGAE')|d('ZNF683')|d('ITGA1')
 if k=='Early_activation':return d('CD69')|d('IL2RA')|d('TNFRSF9')|d('ICOS')
 if k=='Cycling':return d('MKI67')|d('TOP2A')
 if k=='IFN_response':return d('IFIT1')|d('IFIT3')|d('MX1')
def locate(name):
 p=list(Path(os.environ.get('GSE281759_INPUT_DIR','/workspace/scratch/GSE281759_uploaded/inputs')).glob(name))+list(Path('/workspace/attachments').rglob(name))
 if not p:raise FileNotFoundError(name)
 return p[0]
def reference_manifest():
 paths=[OLD/'analyze.py',OLD/'candidate_LR_resource.csv',BASE/'sample_scope.json']
 for item in SAMPLES:
  paths.extend([BASE/f'{item["sample"]}_spatial_marker_subset.h5ad',EXT/f'{item["sample"]}_supplemental_T_markers.h5ad'])
 (R/'reference_manifest.json').write_text(json.dumps([dict(path=str(p.relative_to(R.parent)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths],indent=2))
def run():
 reference_manifest()
 config=dict(targets=PANELS,minimum_detected=MIN,anchors=ANCHOR,quantile=.75,quantile_population='fixed CD8-supported gate in each section',score='mean log1p(CP10k), full-library denominator, zeros included',cd8='previous fixed T-lineage gate AND CD8A>0 AND CD8B>0, CD8A/B mean>=T-gate Q75 and >0',seed=SEED,permutations=B,primary_tests=192,joint_tests=576,primary_only_tests=180,labels_overlap=True,signature_inference=False)
 (R/'analysis_config.json').write_text(json.dumps(config,indent=2))
 manifest={v['name']:v for v in json.loads((BASE/'source_manifest.json').read_text())};inputs=[];checks=[];des=[];stats=[];thresholds=[];markers=[];lr=[];comb=[];inter=[];target_groups=[];edge_rows=[]
 pairs=pd.read_csv(OLD/'candidate_LR_resource.csv');pairs.to_csv(R/'candidate_LR_resource.csv',index=False)
 genes=list(dict.fromkeys(['CD8A','CD8B','PGAM5']+[g for gs in PANELS.values() for g in gs]+pairs.ligand.to_list()+[g for c in pairs.receptor_complex for g in c.split('_')]))
 for item in SAMPLES:
  sample=item['sample'];paths={}
  for suffix in ['matrix.mtx.gz','features.tsv.gz','barcodes.tsv.gz']:
   name=f'{item["gsm"]}_{sample}_{suffix}';p=locate(name);v=manifest[name];assert p.stat().st_size==v['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==v['sha256'];inputs.append(dict(v,local_path=str(p)));checks.append(sample+'/'+suffix+' source SHA');paths[suffix]=p
  f=pd.read_csv(paths['features.tsv.gz'],sep='\t',header=None);bc=pd.Index(pd.read_csv(paths['barcodes.tsv.gz'],header=None)[0].astype(str));assert bc.is_unique and f[0].is_unique
  with gzip.open(paths['matrix.mtx.gz'],'rb') as h:full=io.mmread(h,spmatrix=True).T.tocsr()
  assert full.shape==(len(bc),len(f)) and np.all(full.data>=0) and np.all(full.data==np.floor(full.data));checks.append(sample+': matrix integer counts/dimensions')
  a=ad.read_h5ad(BASE/f'{sample}_spatial_marker_subset.h5ad');ix=bc.get_indexer(a.obs_names);assert (ix>=0).all();full=full[ix];total=np.asarray(full.sum(1)).ravel();assert np.array_equal(total,a.obs.total_counts.to_numpy());checks.append(sample+': barcode/full-library UMI')
  ids=pd.Index(f[0]);oldix=ids.get_indexer(a.var.gene_id);assert np.array_equal(dense(full[:,oldix]),dense(a.layers['counts']));checks.append(sample+': all original marker counts')
  for g in genes:assert int((f[1]==g).sum())==1,g
  symbols=ad.utils.make_index_unique(pd.Index(f[1]));gi=symbols.get_indexer(genes);assert (gi>=0).all();raw=dense(full[:,gi]);x=np.log1p(raw/total[:,None]*10000);o=a.obs.copy()
  previous=ad.read_h5ad(EXT/f'{sample}_supplemental_T_markers.h5ad');cd8=previous.obs.loc[o.index,'gate_CD8_supported'].to_numpy(bool)
  cd8_ix=[genes.index(g) for g in ['CD8A','CD8B']];t=o.T_RNA_enriched.to_numpy(bool);cd8score=x[:,cd8_ix].mean(1);expected=t&(cd8score>=np.quantile(cd8score[t],.75))&(raw[:,cd8_ix]>0).all(1)&(cd8score>0);assert np.array_equal(cd8,expected);checks.append(sample+': fixed CD8 gate reconstructed')
  gates={'All_CD8_reference':cd8};o['gate_All_CD8_reference']=cd8
  d=lambda g:raw[:,genes.index(g)]>0
  for k,gs in PANELS.items():
   j=[genes.index(g) for g in gs];score=x[:,j].mean(1);n=(raw[:,j]>0).sum(1);cut=float(np.quantile(score[cd8],.75)) if cd8.any() else np.nan
   gate=cd8&(score>=cut)&(score>0)&(n>=MIN[k])&anchor_mask(k,d);gates[k]=gate;o['gate_'+k]=gate;o['score_'+k]=score;o['detected_'+k]=n
   thresholds.append(dict(**item,target=k,n_CD8=int(cd8.sum()),threshold=cut,n_target=int(gate.sum()),status='defined' if cd8.any() else 'no_CD8_reference_spots'))
   for g in gs:markers.append(dict(**item,target=k,gene=g,qc_detected=int(d(g).sum()),CD8_detected=int(d(g)[cd8].sum()),target_detected=int(d(g)[gate].sum()),n_qc=len(o),n_CD8=int(cd8.sum()),n_target=int(gate.sum()),CD8_mean_log1p=float(x[cd8,genes.index(g)].mean()) if cd8.any() else np.nan))
  W=a.obsp['spatial_connectivities'];Y=np.column_stack([np.asarray(W@gates[k].astype(int)).ravel()>0 for k in K]);mac=o.Mac_RNA_enriched.to_numpy(bool)&(np.asarray(W.sum(1)).ravel()>0);positive=o.PGAM5_raw.to_numpy()>0;high=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy();rng=np.random.default_rng(SEED+int(item['gsm'][3:]))
  for kind,eligible in [('high_vs_other_detected',mac&positive),('high_vs_all_other_mac',mac),('nonT_high_vs_other_detected',mac&positive&~t)]:stats.extend(helper.evaluate(o,Y,eligible,item,kind,rng))
  near_high=np.asarray(W@(mac&high).astype(int)).ravel()>0;near_other=np.asarray(W@(mac&positive&~high).astype(int)).ravel()>0
  o['CD8_near_high']=cd8&near_high;o['CD8_near_other_detected']=cd8&near_other
  for zone,zone_mask in [('High_near_including_shared',cd8&near_high),('Other_only_near',cd8&~near_high&near_other),('Neither_high_nor_other_near',cd8&~near_high&~near_other)]:
   for k in ALL:
    ng=int(zone_mask.sum());nt=int((zone_mask&gates[k]).sum());target_groups.append(dict(**item,zone=zone,target=k,n_CD8=ng,n_program=nt,program_pct=100*nt/ng if ng else np.nan,mean_program_score=float(o.loc[zone_mask,'score_'+k].mean()) if ng and k in K else np.nan))
  for group,mask in [('PGAM5_high',mac&high),('Other_detected',mac&positive&~high),('Undetected',mac&~positive)]:
   for k in ALL:
    near=np.asarray(W@gates[k].astype(int)).ravel()>0;n=int(mask.sum());v=int(near[mask].sum());des.append(dict(**item,source_group=group,target=k,n_sources=n,with_neighbor=v,neighbor_pct=100*v/n if n else np.nan,n_target_spots=int(gates[k].sum())))
    count=sum(any(gates[k][j] for j in W.indices[W.indptr[src]:W.indptr[src+1]]) for src in np.flatnonzero(mask));assert count==v;checks.append(sample+'/'+group+'/'+k+' explicit neighbors')
  # A target can have >1 supporting program; exact signatures are the finest descriptive labels.
  signatures=np.array(['+'.join(k for k in K if gates[k][i]) or 'No_supported_program' for i in range(len(o))],dtype=object);o['CD8_program_signature']=np.where(cd8,signatures,'Outside_CD8_gate')
  for label in sorted(set(signatures[cd8])):
   target=cd8&(signatures==label);near=np.asarray(W@target.astype(int)).ravel()>0;comb.append(dict(**item,signature=label,n_target=int(target.sum()),high_sources_with_neighbor=int(near[high&mac].sum()),n_high_sources=int((high&mac).sum())))
  assert sum(v['n_target'] for v in comb if v['sample']==sample)==int(cd8.sum());checks.append(sample+': exact signature partition')
  for k in K:
   for k2 in K:inter.append(dict(**item,first=k,second=k2,both=int((gates[k]&gates[k2]).sum())))
  edge=W.tocoo()
  cd8_edge=high[edge.row]&cd8[edge.col]
  edge_rows.extend(dict(**item,source_barcode=o.index[src],target_barcode=o.index[tgt],signature=signatures[tgt]) for src,tgt in zip(edge.row[cd8_edge],edge.col[cd8_edge]))
  for k in ALL:
   use=high[edge.row]&gates[k][edge.col];src,tgt=edge.row[use],edge.col[use]
   for p in pairs.itertuples():
    li=genes.index(p.ligand);ri=[genes.index(g) for g in p.receptor_complex.split('_')];support=(raw[src,li]>0)&(raw[np.ix_(tgt,ri)]>0).all(1)
    lr.append(dict(**item,target=k,pair=p.ligand+'--'+p.receptor_complex,eligible_edges=len(src),unique_sources=len(np.unique(src)),unique_targets=len(np.unique(tgt)),supported_edges=int(support.sum()),supported_sources=len(np.unique(src[support])),supported_targets=len(np.unique(tgt[support])),codetection_pct=100*support.mean() if len(src) else np.nan))
  for j,k in enumerate(K):o['neighbor_'+k]=Y[:,j]
  o.to_csv(R/f'{sample}_annotations.csv.gz',float_format='%.17g');extra=ad.AnnData(sparse.csr_matrix(x),obs=o,var=pd.DataFrame(index=genes));extra.layers['counts']=sparse.csr_matrix(raw);extra.obsp['spatial_connectivities']=W.copy();extra.obsm['spatial']=a.obsm['spatial'].copy();extra.write_h5ad(R/f'{sample}_CD8_fine_markers.h5ad',compression='gzip')
  print(sample,'CD8',int(cd8.sum()),{k:int(gates[k].sum()) for k in K},flush=True)
 dd=pd.DataFrame(des);dd.to_csv(R/'neighborhood_counts.csv',index=False);pd.DataFrame(thresholds).to_csv(R/'thresholds.csv',index=False);pd.DataFrame(markers).to_csv(R/'marker_detection.csv',index=False);pd.DataFrame(lr).to_csv(R/'candidate_LR_codetection.csv',index=False);pd.DataFrame(comb).to_csv(R/'exact_program_signatures.csv',index=False);pd.DataFrame(inter).to_csv(R/'program_overlap.csv',index=False);pd.DataFrame(target_groups).to_csv(R/'CD8_target_neighborhood_programs.csv',index=False);pd.DataFrame(edge_rows).to_csv(R/'high_source_CD8_directed_edges.csv.gz',index=False)
 st=pd.DataFrame(stats);st['joint_q_576']=false_discovery_control(st.p.fillna(1).to_numpy());st['primary_q_192']=np.nan;sel=st.comparison.eq('high_vs_other_detected');st.loc[sel,'primary_q_192']=false_discovery_control(st.loc[sel,'p'].fillna(1).to_numpy());st['primary_only_q_180']=np.nan;sel2=sel&st.tumor_status.eq('primary');st.loc[sel2,'primary_only_q_180']=false_discovery_control(st.loc[sel2,'p'].fillna(1).to_numpy());st.to_csv(R/'spatial_statistics.csv',index=False)
 patient=dd.groupby(['patient','source_group','target'],as_index=False).agg(n_sources=('n_sources','sum'),with_neighbor=('with_neighbor','sum'),n_target_spots=('n_target_spots','sum'));patient['neighbor_pct']=100*patient.with_neighbor/patient.n_sources;patient.to_csv(R/'patient_descriptive.csv',index=False)
 st[st.status.eq('tested')].groupby(['patient','comparison','target'],as_index=False).agg(median_delta_pp=('delta_pp','median'),n_testable_sections=('sample','size')).to_csv(R/'patient_equal_section_effects.csv',index=False)
 h=dd[dd.source_group.eq('PGAM5_high')].groupby('target')[['n_sources','with_neighbor','n_target_spots']].sum();h['coverage_pct']=100*h.with_neighbor/h.n_sources;h.to_csv(R/'pooled_descriptive.csv')
 (R/'input_manifest.json').write_text(json.dumps(inputs,indent=2));(R/'audit_results.json').write_text(json.dumps(dict(passed=len(checks),checks=checks),indent=2));print(h.to_string());print('min primary q',st.primary_q_192.min(),'CHECKS',len(checks))
if __name__=='__main__':run()
