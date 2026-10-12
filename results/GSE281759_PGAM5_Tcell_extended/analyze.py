"""Extended mixed-spot T-RNA programs from verified full-gene GEO matrices."""
from pathlib import Path
import json,hashlib,gzip,importlib.util,os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse,io
from scipy.stats import false_discovery_control
R=Path(__file__).resolve().parent;BASE=R.parent/'GSE281759_spatial_HCC_PGAM5';OLD=R.parent/'GSE281759_PGAM5_immune_subgroup_neighbors'
PANELS={'CD4_related':['CD4','IL7R','LTB'],'CD8_supported':['CD8A','CD8B'],'Treg_supported':['FOXP3','IL2RA','CTLA4','IKZF2','TIGIT'],'Naive_like':['CCR7','SELL','TCF7','LEF1','MAL'],'Central_memory_like':['CCR7','IL7R','LTB','TCF7'],'Effector_memory_like':['GZMK','CCL5','IL7R']}
K=list(PANELS);SAMPLES=json.loads((BASE/'sample_scope.json').read_text());B=4999;SEED=281759202
spec=importlib.util.spec_from_file_location('prior_analysis',OLD/'analyze.py');helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper);helper.ROOT=R;helper.TARGETS=K
def dense(x):return np.asarray(x.toarray() if sparse.issparse(x) else x)
def locate(name):
 input_dir=Path(os.environ.get('GSE281759_INPUT_DIR','/workspace/scratch/GSE281759_uploaded/inputs'))
 paths=list(input_dir.glob(name))+list(Path('/workspace/attachments').rglob(name))
 if not paths:raise FileNotFoundError(name)
 return paths[0]
def run():
 manifest={v['name']:v for v in json.loads((BASE/'source_manifest.json').read_text())};inputs=[];checks=[];des=[];stats=[];markers=[];thresholds=[];lr=[]
 genes=list(dict.fromkeys(g for gs in PANELS.values() for g in gs));pairs=pd.read_csv(OLD/'candidate_LR_resource.csv');pairs.to_csv(R/'candidate_LR_resource.csv',index=False)
 for item in SAMPLES:
  s=item['sample'];prefix=f'{item["gsm"]}_{s}';paths={}
  for suffix in ['matrix.mtx.gz','features.tsv.gz','barcodes.tsv.gz']:
   name=prefix+'_'+suffix;p=locate(name);expected=manifest[name]
   assert p.stat().st_size==expected['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==expected['sha256'];inputs.append(dict(expected,local_path=str(p)));checks.append(s+'/'+suffix+': source SHA256')
   paths[suffix]=p
  f=pd.read_csv(paths['features.tsv.gz'],sep='\t',header=None);bc=pd.Index(pd.read_csv(paths['barcodes.tsv.gz'],header=None)[0].astype(str));assert bc.is_unique
  with gzip.open(paths['matrix.mtx.gz'],'rb') as stream:full=io.mmread(stream,spmatrix=True).T.tocsr()
  assert full.shape==(len(bc),len(f)) and np.all(full.data>=0) and np.all(full.data==np.floor(full.data));checks.append(s+': matrix dimensions/integer counts')
  a=ad.read_h5ad(BASE/f'{s}_spatial_marker_subset.h5ad');ix=bc.get_indexer(a.obs_names);assert (ix>=0).all()
  full=full[ix];total=np.asarray(full.sum(1)).ravel();assert np.array_equal(total,a.obs.total_counts.to_numpy());checks.append(s+': full-library UMI and barcode alignment')
  symbols=ad.utils.make_index_unique(pd.Index(f[1]));oldix=pd.Index(f[0]).get_indexer(a.var.gene_id) if 'gene_id' in a.var else symbols.get_indexer(a.var_names)
  assert (oldix>=0).all() and np.array_equal(dense(full[:,oldix]),dense(a.layers['counts']));checks.append(s+': all existing marker counts match full matrix')
  gi=symbols.get_indexer(genes);assert (gi>=0).all();raw=dense(full[:,gi]);norm=np.log1p(raw/total[:,None]*10000);o=a.obs.copy();t=o.T_RNA_enriched.to_numpy(bool);gates={}
  def det(g):return raw[:,genes.index(g)]>0
  for k,gs in PANELS.items():
   j=[genes.index(g) for g in gs];score=norm[:,j].mean(1);n=(raw[:,j]>0).sum(1);cut=float(np.quantile(score[t],.75));minimum=3 if k in ['Naive_like','Central_memory_like'] else 2
   gate=t&(score>=cut)&(score>0)&(n>=minimum)
   if k=='CD4_related':gate&=det('CD4')
   if k=='CD8_supported':gate&=det('CD8A')&det('CD8B')
   if k=='Treg_supported':gate&=det('FOXP3')
   if k=='Naive_like':gate&=(det('CCR7')|det('SELL'))&(det('TCF7')|det('LEF1'))
   if k=='Central_memory_like':gate&=det('CCR7')&det('IL7R')
   if k=='Effector_memory_like':gate&=det('GZMK')&(det('CCL5')|det('IL7R'))
   gates[k]=gate;o['gate_'+k]=gate;o['score_'+k]=score;o['detected_'+k]=n;thresholds.append(dict(**item,target=k,score_q75_in_T=cut,n_target=int(gate.sum())))
   for g in gs:markers.append(dict(**item,target=k,gene=g,in_feature_list=True,qc_spots_detected=int(det(g).sum()),T_spots_detected=int(det(g)[t].sum()),target_spots_detected=int(det(g)[gate].sum()),n_qc=len(o),n_T=int(t.sum()),n_target=int(gate.sum())))
  W=a.obsp['spatial_connectivities'];Y=np.column_stack([np.asarray(W@gates[k].astype(int)).ravel()>0 for k in K]);degree=np.asarray(W.sum(1)).ravel();mac=o.Mac_RNA_enriched.to_numpy(bool)&(degree>0);positive=o.PGAM5_raw.to_numpy()>0;high=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy();rng=np.random.default_rng(SEED+int(item['gsm'][3:]))
  for kind,eligible in [('high_vs_other_detected',mac&positive),('high_vs_all_other_mac',mac),('nonT_high_vs_other_detected',mac&positive&~t)]:stats.extend(helper.evaluate(o,Y,eligible,item,kind,rng))
  for group,mask in [('PGAM5_high',mac&high),('Other_detected',mac&positive&~high),('Undetected',mac&~positive)]:
   for j,k in enumerate(K):
    n=int(mask.sum());v=int(Y[mask,j].sum());des.append(dict(**item,source_group=group,target=k,n_sources=n,with_neighbor=v,neighbor_pct=100*v/n if n else np.nan,n_target_spots=int(gates[k].sum())))
    manual=sum(any(gates[k][z] for z in W.indices[W.indptr[src]:W.indptr[src+1]]) for src in np.flatnonzero(mask));assert manual==v;checks.append(s+'/'+group+'/'+k+': explicit neighbor recount')
  edge=W.tocoo()
  for k in K:
   use=high[edge.row]&gates[k][edge.col];src,tgt=edge.row[use],edge.col[use]
   for p in pairs.itertuples():
    ligand=symbols.get_loc(p.ligand);ri=symbols.get_indexer(p.receptor_complex.split('_'));assert (ri>=0).all()
    support=(dense(full[src,ligand]).ravel()>0)&(dense(full[tgt][:,ri])>0).all(1)
    lr.append(dict(**item,target=k,pair=p.ligand+'--'+p.receptor_complex,eligible_edges=len(src),unique_sources=len(np.unique(src)),unique_targets=len(np.unique(tgt)),supported_edges=int(support.sum()),supported_sources=len(np.unique(src[support])),supported_targets=len(np.unique(tgt[support])),codetection_pct=100*support.mean() if len(src) else np.nan))
  for j,k in enumerate(K):o['neighbor_'+k]=Y[:,j]
  o.to_csv(R/f'{s}_annotations.csv.gz',float_format='%.17g')
  # Raw/norm supplemental genes suffice to reconstruct new target gates; base supplies source/geometry.
  extra=ad.AnnData(sparse.csr_matrix(norm),obs=o,var=pd.DataFrame(index=genes));extra.layers['counts']=sparse.csr_matrix(raw);extra.obsp['spatial_connectivities']=W.copy();extra.obsm['spatial']=a.obsm['spatial'].copy();extra.write_h5ad(R/f'{s}_supplemental_T_markers.h5ad',compression='gzip')
  print(s,{k:int(gates[k].sum()) for k in K},flush=True)
 dd=pd.DataFrame(des);dd.to_csv(R/'neighborhood_counts.csv',index=False);pd.DataFrame(markers).to_csv(R/'marker_detection.csv',index=False);pd.DataFrame(thresholds).to_csv(R/'thresholds.csv',index=False);pd.DataFrame(lr).to_csv(R/'candidate_LR_codetection.csv',index=False)
 st=pd.DataFrame(stats);st['joint_q_288']=false_discovery_control(st.p.fillna(1).to_numpy());st['primary_q_96']=np.nan;sel=st.comparison.eq('high_vs_other_detected');st.loc[sel,'primary_q_96']=false_discovery_control(st.loc[sel,'p'].fillna(1).to_numpy());st['primary_only_q_90']=np.nan;sel2=sel&st.tumor_status.eq('primary');st.loc[sel2,'primary_only_q_90']=false_discovery_control(st.loc[sel2,'p'].fillna(1).to_numpy());st.to_csv(R/'spatial_statistics.csv',index=False)
 patient=dd.groupby(['patient','source_group','target'],as_index=False).agg(n_sources=('n_sources','sum'),with_neighbor=('with_neighbor','sum'),n_target_spots=('n_target_spots','sum'));patient['neighbor_pct']=100*patient.with_neighbor/patient.n_sources;patient.to_csv(R/'patient_descriptive.csv',index=False)
 st[st.status.eq('tested')].groupby(['patient','comparison','target'],as_index=False).agg(median_delta_pp=('delta_pp','median'),n_testable_sections=('sample','size')).to_csv(R/'patient_equal_section_effects.csv',index=False)
 (R/'input_manifest.json').write_text(json.dumps(inputs,indent=2));(R/'audit_results.json').write_text(json.dumps(dict(passed=len(checks),checks=checks),indent=2))
 (R/'analysis_config.json').write_text(json.dumps(dict(targets=PANELS,permutations=B,seed=SEED,primary_planned_tests=96,joint_planned_tests=288,labels_overlap=True,primary_only_tests=90,source_gate_unchanged=True,minimum_detected={'Naive_like':3,'Central_memory_like':3,'others':2},anchors={'CD4_related':'CD4','CD8_supported':'CD8A AND CD8B','Treg_supported':'FOXP3','Naive_like':'(CCR7 OR SELL) AND (TCF7 OR LEF1)','Central_memory_like':'CCR7 AND IL7R','Effector_memory_like':'GZMK AND (CCL5 OR IL7R)'},T_state_cutoff='within existing T gate Q75, score>0',not_resolved_cell_identity=True),indent=2))
 print('CHECKS',len(checks));print('min primary q',st.primary_q_96.min());print('significant',st[st.primary_q_96<.05].to_string())
if __name__=='__main__':run()
