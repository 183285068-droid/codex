"""Tumor-section expression-only T candidates. No coordinates or cell pairing."""
from pathlib import Path
import json,hashlib,shutil
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
R=Path(__file__).resolve().parent;SOURCE=R.parent/'GSE283406_PGAM5_immune_expression_candidates';K=['All_T','Exhaustion_T','Activation_T','Cytotoxic_T']
PANELS={'All_T':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2'],'Exhaustion_T':['PDCD1','LAG3','TIGIT','CTLA4','TOX','CXCL13'],'Activation_T':['CD69','IL2RA','TNFRSF9','ICOS','CD40LG','IL2'],'Cytotoxic_T':['IFNG','PRF1','GZMB','GZMA','CCL5','TBX21']}
FIELDS=['T_RNA_enriched','Exhaustion_core_supported','Activation_early_supported','Effector_cytotoxic_supported'];PAIRS=[('SPP1','CD44'),('PVR','TIGIT'),('NECTIN2','TIGIT'),('LGALS9','HAVCR2'),('CD274','PDCD1')]
def dense(x):return np.asarray(x.toarray() if sparse.issparse(x) else x)
def input_path():return R/'inputs/J259X_AVITI_expression_marker_subset.h5ad'
def run():
 path=input_path();path.parent.mkdir(exist_ok=True)
 if not path.exists():shutil.copyfile(SOURCE/'inputs/J259X_AVITI_expression_marker_subset.h5ad',path)
 if not (R/'candidate_LR_resource_snapshot.csv').exists():shutil.copyfile(SOURCE/'candidate_LR_resource_snapshot.csv',R/'candidate_LR_resource_snapshot.csv')
 a=ad.read_h5ad(path);o=a.obs.copy();assert 'spatial' not in a.obsm and 'spatial_connectivities' not in a.obsp
 gates={k:o[field].to_numpy(bool) for k,field in zip(K,FIELDS)};counts=dense(a.layers['counts']);x=dense(a.X);assert np.allclose(x,np.log1p(counts/o.total_counts.to_numpy()[:,None]*10000),rtol=1e-6,atol=1e-6)
 high=o.PGAM5_mac_RNA_group.eq('PGAM5_high').to_numpy();other=o.PGAM5_mac_RNA_group.eq('Other_detected').to_numpy();und=o.PGAM5_mac_RNA_group.eq('Undetected').to_numpy();groups={'PGAM5_high':high,'Other_detected':other,'Undetected':und};over=[];conditional=[];lr=[];marker=[];signatures=[]
 for k in K:
  ix=a.var_names.get_indexer(PANELS[k]);assert (ix>=0).all();score=x[:,ix].mean(1);o['gate_'+k]=gates[k];o['score_'+k]=score
  for group,mask in groups.items():
   n=int(mask.sum());v=int((mask&gates[k]).sum());over.append(dict(target=k,source_group=group,n_sources=n,n_target=int(gates[k].sum()),same_barcode_overlap=v,same_barcode_pct=100*v/n,mean_program_score=float(score[mask].mean())))
   tg=mask&gates['All_T'];conditional.append(dict(target=k,source_group=group,n_T_coexpressing_sources=int(tg.sum()),state_supported=int((tg&gates[k]).sum()),state_pct=100*(tg&gates[k]).sum()/tg.sum() if tg.any() else np.nan))
  for g in PANELS[k]:
   j=a.var_names.get_loc(g);mask=gates[k];marker.append(dict(target=k,gene=g,n_targets=int(mask.sum()),detection_pct=100*(counts[mask,j]>0).mean(),mean_log1p_CP10k=float(x[mask,j].mean())))
  for ligand,receptor in PAIRS:
   li=a.var_names.get_loc(ligand);ri=a.var_names.get_loc(receptor);mask=gates[k]
   for group,source in groups.items():lr.append(dict(target=k,source_group=group,pair=ligand+'--'+receptor,n_sources=int(source.sum()),ligand_detected=int((counts[source,li]>0).sum()),ligand_pct=100*(counts[source,li]>0).mean(),n_target=int(mask.sum()),receptor_detected=int((counts[mask,ri]>0).sum()),receptor_pct=100*(counts[mask,ri]>0).mean(),source_target_pairing=False))
 sig=np.array(['+'.join(k for k in K[1:] if gates[k][i]) or 'No_supported_state' for i in range(len(o))]);o['T_program_signature']=np.where(gates['All_T'],sig,'Outside_T_gate')
 for label in sorted(set(sig[gates['All_T']])):
  t=gates['All_T']&(sig==label);signatures.append(dict(signature=label,n_T_targets=int(t.sum()),high_same_barcode=int((t&high).sum()),other_same_barcode=int((t&other).sum()),undetected_same_barcode=int((t&und).sum())))
 coverage=[]
 for panel,gs in {'CD8':['CD8A','CD8B'],'CD4':['CD4','IL7R','LTB'],'Treg':['FOXP3','IKZF2','IL2RA','CTLA4','TIGIT'],'Naive_memory':['CCR7','SELL','TCF7','LEF1','MAL','IL7R','GZMK','ITGAE','ZNF683'],**PANELS}.items():
  for g in gs:coverage.append(dict(panel=panel,gene=g,in_current_marker_input=g in a.var_names))
 pd.DataFrame(over).to_csv(R/'same_barcode_T_program_overlap.csv',index=False);pd.DataFrame(conditional).to_csv(R/'conditional_T_program_descriptive.csv',index=False);pd.DataFrame(lr).to_csv(R/'ligand_receptor_independent_expression.csv',index=False);pd.DataFrame(marker).to_csv(R/'T_marker_evidence.csv',index=False);pd.DataFrame(signatures).to_csv(R/'T_program_signatures.csv',index=False);pd.DataFrame(coverage).to_csv(R/'marker_coverage.csv',index=False);o.to_csv(R/'J259X_T_annotations.csv.gz',float_format='%.17g')
 (R/'reference_manifest.json').write_text(json.dumps([dict(path=str(path.relative_to(R)),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest())],indent=2));(R/'analysis_config.json').write_text(json.dumps(dict(primary='J259X tumor section AVITI GSM8661792',patients=1,QC='UMI>=500 and detected genes>=200, FFPE probe panel; no polyA MT25 gate',panels=PANELS,source_gate_fixed=True,high_Q75_among_detected_Mac_RNA=True,state_cutoff='T-gate Q75, score>0, >=2 state genes',coordinates_available=False,cell_pairing=False,LR='separate source ligand and target receptor expression only',no_p_values=True,context_J259B_excluded_from_tumor_primary=True),indent=2));print(pd.DataFrame(over).to_string(index=False));print(pd.DataFrame(signatures).to_string(index=False))
if __name__=='__main__':run()
