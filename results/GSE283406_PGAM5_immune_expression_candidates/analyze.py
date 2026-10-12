"""Expression-only immune RNA candidates. No coordinates, edges or confirmed interactions."""
from pathlib import Path
import json,hashlib,shutil
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
ROOT=Path(__file__).resolve().parent;REF=ROOT.parent/'GSE283406_HCC_expression_preanalysis'
PANELS={'All_T':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2'],
 'Exhaustion_T':['PDCD1','LAG3','TIGIT','CTLA4','TOX','CXCL13'],
 'Activation_T':['CD69','IL2RA','TNFRSF9','ICOS','CD40LG','IL2'],
 'Cytotoxic_T':['IFNG','PRF1','GZMB','GZMA','CCL5','TBX21'],
 'B_related':['MS4A1','CD79A','CD79B','CD19','CD22'],
 'Plasma_related':['MZB1','JCHAIN','SDC1','IGKC'],
 'cDC2_related':['CD1C','FCER1A','CLEC10A'],
 'cDC1_related':['CLEC9A','XCR1'],
 'Inflammatory_myeloid':['S100A8','S100A9','VCAN','IL1B']}
PAIRS=[('SPP1','CD44'),('PVR','TIGIT'),('NECTIN2','TIGIT'),('LGALS9','HAVCR2'),('CD274','PDCD1')]
def dense(x):return np.asarray(x.toarray() if sparse.issparse(x) else x)
def annotate(a):
 o=a.obs
 for k,gs in PANELS.items():
  assert all(g in a.var_names for g in gs);norm=dense(a[:,gs].X);counts=dense(a[:,gs].layers['counts']);o['score_'+k]=norm.mean(1);o['detected_'+k]=(counts>0).sum(1)
  prev={'All_T':'T_RNA_enriched','Exhaustion_T':'Exhaustion_core_supported','Activation_T':'Activation_early_supported','Cytotoxic_T':'Effector_cytotoxic_supported'}
  if k in prev:gate=o[prev[k]].to_numpy(bool)
  else:
   score=o['score_'+k].to_numpy();gate=(score>=np.quantile(score,.75))&(score>0)&((counts>0).sum(1)>=(3 if k=='Plasma_related' else 2))
   def has(gs):return (dense(a[:,gs].layers['counts'])>0).any(1)
   if k=='B_related':gate&=has(['MS4A1','CD19'])
   if k=='Plasma_related':gate&=has(['MZB1','JCHAIN'])
   if k=='cDC2_related':gate&=has(['CD1C','CLEC10A'])
   if k=='cDC1_related':gate&=has(['CLEC9A'])&has(['XCR1'])
   if k=='Inflammatory_myeloid':gate&=has(['VCAN'])&o.Pan_myeloid.ge(o.Pan_myeloid.median()).to_numpy()
  o['gate_'+k]=gate
 return a
def run():
 cfg=dict(kind='expression-only proxy analysis, not spatial interaction',primary='J259X tumor section AVITI',context='J259B mixed tumor/adjacent section, not tumor ROI confirmed',biological_patients=1,coordinates_available=False,subgroup_panels=PANELS,missing_subgroup_markers=['NK','Treg','naive/memory T','CD4/CD8 identity','pDC','granulocyte','mast'],LR_pairs=PAIRS,LR_method='source ligand and target receptor detected/mean RNA independently summarized, no source-target pairing, no interaction score/p',cooccurrence='immune gate overlaps PGAM5-high Mac RNA SAME barcode; mixed RNA not cell neighbors',no_UMAP_adjacency=True,no_patient_p_values=True)
 (ROOT/'analysis_config.json').write_text(json.dumps(cfg,indent=2)+'\n');manifest=[];summary=[];overlap=[];lr=[];coverage=[];marker=[]
 for section in ['J259X','J259B']:
  dst=ROOT/'inputs'/f'{section}_AVITI_expression_marker_subset.h5ad'
  if not dst.exists():shutil.copyfile(REF/dst.name,dst)
  a=annotate(ad.read_h5ad(dst));o=a.obs;check=not ('spatial' in a.obsm or 'spatial_connectivities' in a.obsp);assert check
  manifest.append(dict(path=str(dst.relative_to(ROOT)),bytes=dst.stat().st_size,sha256=hashlib.sha256(dst.read_bytes()).hexdigest()));high=o.PGAM5_mac_RNA_group.eq('PGAM5_high').to_numpy();other=o.PGAM5_mac_RNA_group.eq('Other_detected').to_numpy();und=o.PGAM5_mac_RNA_group.eq('Undetected').to_numpy();scope='tumor_primary_analysis' if section=='J259X' else 'mixed_tissue_context_only'
  summary.append(dict(section=section,scope=scope,qc_spots=len(a),mac_RNA=int(o.Mac_RNA_enriched.sum()),PGAM5_high=int(high.sum()),other_detected=int(other.sum()),undetected=int(und.sum()),high_cutoff=float(a.uns['PGAM5_high_threshold_log1p_CP10k'])))
  for k,gs in PANELS.items():
   target=o['gate_'+k].to_numpy(bool)
   for group,mask in [('PGAM5_high',high),('Other_detected',other),('Undetected',und)]:overlap.append(dict(section=section,scope=scope,target=k,source_group=group,n_sources=int(mask.sum()),total_target_spots=int(target.sum()),same_spot_overlap=int((mask&target).sum()),same_spot_overlap_pct=100*(mask&target).sum()/mask.sum() if mask.any() else np.nan,mean_target_score_in_sources=float(o.loc[mask,'score_'+k].mean()) if mask.any() else np.nan))
   for gene in gs:
    v=dense(a[:,[gene]].layers['counts']).ravel();norm=dense(a[:,[gene]].X).ravel();marker.append(dict(section=section,target=k,gene=gene,n_target=int(target.sum()),detected_pct=100*(v[target]>0).mean() if target.any() else np.nan,mean_expression=norm[target].mean() if target.any() else np.nan))
   for ligand,receptor in PAIRS:
    lig=dense(a[:,[ligand]].layers['counts']).ravel()>0;rec=dense(a[:,[receptor]].layers['counts']).ravel()>0
    lr.append(dict(section=section,scope=scope,target=k,pair=ligand+'--'+receptor,ligand=ligand,receptor=receptor,n_high_sources=int(high.sum()),high_ligand_detected=int((high&lig).sum()),high_ligand_detect_pct=100*(high&lig).sum()/high.sum() if high.any() else np.nan,n_target_spots=int(target.sum()),target_receptor_detected=int((target&rec).sum()),target_receptor_detect_pct=100*(target&rec).sum()/target.sum() if target.any() else np.nan,source_target_gate_overlap=int((high&target).sum()),spatial_pairing_available=False))
  for panel,gs in {'NK':['NKG7','GNLY','KLRD1'],'Treg':['FOXP3','IKZF2','IL2RA'],'Naive_memory_T':['CCR7','TCF7','LEF1','IL7R'],'CD8':['CD8A','CD8B'],'MHCII':['HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','CD74'],'Key_LR':['CXCL12','CXCR4','TGFBR1','TGFBR2','IL10','IL10RA','IL10RB'],**PANELS}.items():
   for g in gs:coverage.append(dict(section=section,panel=panel,gene=g,in_marker_input=g in a.var_names))
  o.to_csv(ROOT/f'{section}_immune_candidate_annotations.csv.gz',float_format='%.17g');print(section,summary[-1],flush=True)
 pd.DataFrame(summary).to_csv(ROOT/'section_summary.csv',index=False);pd.DataFrame(overlap).to_csv(ROOT/'same_spot_RNA_overlap.csv',index=False);pd.DataFrame(lr).to_csv(ROOT/'candidate_LR_independent_expression.csv',index=False);pd.DataFrame(marker).to_csv(ROOT/'target_marker_evidence.csv',index=False);pd.DataFrame(coverage).to_csv(ROOT/'marker_coverage.csv',index=False);(ROOT/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
if __name__=='__main__':run()
