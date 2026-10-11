"""Provisional spot-region labels reviewed against descriptive marker tables."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE245908_mpl')
import pandas as pd,scanpy as sc
from preprocess import ROOT,SAMPLES,WORK
LABELS={
 'CHC20':{
 'S0':('CYP2A6/ALDOB','Metabolic hepatocytic RNA','代谢相关肝细胞RNA区域'),
 'S1':('RBP4/FGA','Secretory hepatocytic RNA','分泌相关肝细胞RNA区域'),
 'S2':('COL1A1/COMP','Matrix / myeloid mixed RNA','基质及髓系混合RNA区域'),
 'S3':('MALAT1/HMGB1','Nuclear / chromatin RNA','核及染色质相关RNA区域'),
 'S4':('LYZ/AGR2','Myeloid / epithelial mixed RNA','髓系及上皮混合RNA区域'),
 'S5':('CD74/C1QA','Antigen-presentation / myeloid RNA','抗原呈递及髓系RNA区域')},
 'CHC23':{
 'S0':('CYP2A6/BHMT','Metabolic hepatocytic RNA','代谢相关肝细胞RNA区域'),
 'S1':('TTR/GSTA1','Secretory / detoxification RNA','分泌及解毒相关RNA区域'),
 'S2':('IGFBP7/CD74','Stromal / immune mixed RNA','基质及免疫混合RNA区域'),
 'S3':('CYP2E1/IGFBP2','Metabolic / MT-high RNA','代谢及线粒体高表达RNA区域'),
 'S4':('MTRNR2L8/SAA1','MT / acute-phase mixed RNA','线粒体及急性期混合RNA区域'),
 'S5':('CRP/ORM1','Acute-phase hepatocytic RNA','急性期肝细胞RNA区域')}
}
if __name__=='__main__':
 rows=[]
 for s in SAMPLES:
  p=ROOT/f'{s}_spot_annotations.csv';d=pd.read_csv(p,index_col=0)
  d['region_name']=d.cluster.map({k:f'{k}: {v[0]} | {v[1]}' for k,v in LABELS[s].items()});d['region_name_zh']=d.cluster.map({k:v[2] for k,v in LABELS[s].items()});assert d.region_name.notna().all();d.to_csv(p,index_label='barcode',float_format='%.17g')
  a=sc.read_h5ad(ROOT/f'{s}_spatial_marker_subset.h5ad');a.obs=d.loc[a.obs_names];a.write_h5ad(ROOT/f'{s}_spatial_marker_subset.h5ad',compression='gzip')
  profiles=pd.read_csv(ROOT/f'{s}_cluster_profiles.csv');profiles['region_name']=profiles.cluster.map({k:f'{k}: {v[0]} | {v[1]}' for k,v in LABELS[s].items()});profiles['region_name_zh']=profiles.cluster.map({k:v[2] for k,v in LABELS[s].items()});profiles.to_csv(ROOT/f'{s}_cluster_profiles.csv',index=False)
  markers=pd.read_csv(ROOT/f'{s}_cluster_top100_markers.csv');full=sc.read_h5ad(WORK/f'{s}_full_processed.h5ad')
  for c,v in LABELS[s].items():
   m=markers[markers.group==c];e=m[m.names.isin(v[0].split('/'))]
   for gene in v[0].split('/'):
    # Marker may be below top100; use full normalized reference to verify expression.
    import numpy as np
    idx=np.flatnonzero(full.var_names==gene)[0];x=np.asarray(full.X[:,idx].toarray()).ravel();ingroup=full.obs.cluster==c
    rows.append({'sample':s,'cluster':c,'marker':gene,'mean_log1p_CP10k_in_region':float(x[ingroup].mean()),'mean_log1p_CP10k_outside':float(x[~ingroup].mean()),'detection_pct_in_region':float((x[ingroup]>0).mean()*100),'name_zh':v[2]})
 pd.DataFrame(rows).to_csv(ROOT/'region_naming_marker_evidence.csv',index=False)
 pd.DataFrame([{'sample':s,'cluster':c,'marker_pair':v[0],'RNA_program':v[1],'name_zh':v[2],'annotation_level':'provisional mixed-spot region; not a single-cell identity'} for s in SAMPLES for c,v in LABELS[s].items()]).to_csv(ROOT/'region_names.csv',index=False)
