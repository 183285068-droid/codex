"""Independent raw normalization, new and reused gates, lattice, LR, p/FDR audit."""
import json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy.spatial import cKDTree
from scipy.stats import false_discovery_control
from analyze import R,BASE,OLD,EXT,K,PANELS,NEW_RULES,EXT_MAP,OLD_KEYS,SAMPLES,dense
checks=[]
def check(name,v):
 if not v:raise AssertionError(name)
 checks.append(name)
def run():
 s=pd.read_csv(R/'spatial_statistics.csv');d=pd.read_csv(R/'neighborhood_counts.csv');lr=pd.read_csv(R/'candidate_LR_codetection.csv');pairs=pd.read_csv(R/'candidate_LR_resource.csv');edges=pd.read_csv(R/'high_source_immune_directed_edges.csv.gz');zones=pd.read_csv(R/'unique_target_neighborhood_zones.csv');mn=pd.read_csv(R/'macrophage_targets_excluding_high_descriptive.csv')
 for v in json.loads((R/'input_manifest.json').read_text()):check(v['name']+' SHA',hashlib.sha256(open(v['local_path'],'rb').read()).hexdigest()==v['sha256'])
 for v in json.loads((R/'reference_manifest.json').read_text()):check(v['path']+' reference SHA',hashlib.sha256((R.parent/v['path']).read_bytes()).hexdigest()==v['sha256'])
 for item in SAMPLES:
  name=item['sample'];a=ad.read_h5ad(R/f'{name}_immune_markers.h5ad');o=a.obs;raw=dense(a.layers['counts']);x=dense(a.X);W=a.obsp['spatial_connectivities'];check(name+' full library normalization',np.allclose(x,np.log1p(raw/o.total_counts.to_numpy()[:,None]*10000),rtol=1e-6,atol=1e-6))
  xy=o[['array_row','array_col']].to_numpy(float);xy[:,0]*=np.sqrt(3);ij=cKDTree(xy).query_pairs(2.000001,output_type='ndarray');expected=set(map(tuple,ij))|set(map(tuple,ij[:,::-1]));e=W.tocoo();check(name+' independent lattice',expected==set(zip(e.row,e.col)));check(name+' no self loop',np.all(e.row!=e.col))
  old=pd.read_csv(OLD/f'{name}_annotations.csv.gz',index_col=0).loc[o.index];ext=ad.read_h5ad(EXT/f'{name}_supplemental_T_markers.h5ad').obs.loc[o.index];gates={}
  for k in OLD_KEYS:gates[k]=old['gate_'+k].to_numpy(bool);check(name+'/'+k+' fixed prior gate',np.array_equal(gates[k],o['gate_'+k]))
  for k,v in EXT_MAP.items():gates[k]=ext['gate_'+v].to_numpy(bool);check(name+'/'+k+' fixed T gate',np.array_equal(gates[k],o['gate_'+k]))
  gates['Mac_RNA_reference']=o.Mac_RNA_enriched.to_numpy(bool);det=lambda g:raw[:,a.var_names.get_loc(g)]>0
  for k,rule in NEW_RULES.items():
   parent=np.ones(len(o),bool) if rule['parent']=='all_QC' else gates[rule['parent']];ix=a.var_names.get_indexer(PANELS[k]);score=x[:,ix].mean(1);n=(raw[:,ix]>0).sum(1);g=parent&(score>=np.quantile(score[parent],.75))&(score>0)&(n>=rule['minimum']) if parent.any() else np.zeros(len(o),bool)
   if k=='NK_related':g&=det('NCR1')|det('NCR3')|det('NCAM1')
   elif k=='Cytotoxic_NK_like':g&=det('PRF1')|det('GNLY')
   elif k=='XCL_NK_like':g&=det('XCL1')|det('XCL2')
   elif k=='Naive_B_like':g&=det('TCL1A')|det('IGHD')
   elif k=='Memory_B_like':g&=det('CD27')&det('TNFRSF13B')
   elif k=='IgA_plasma_related':g&=det('IGHA1')|det('IGHA2')
   elif k=='IgG_plasma_related':g&=det('IGHG1')|det('IGHG2')|det('IGHG3')|det('IGHG4')
   elif k=='IgM_plasma_related':g&=det('IGHM')
   elif k=='pDC_related':g&=(det('CLEC4C')|det('LILRA4'))&(det('TCF4')|det('IL3RA'))
   elif k=='Classical_monocyte_related':g&=(det('FCN1')|det('CD14'))&(o.Pan_myeloid.to_numpy()>=o.Pan_myeloid.median())
   elif k=='FCGR3A_myeloid_related':g&=det('FCGR3A')&(det('MS4A7')|det('LILRB1'))&(o.Pan_myeloid.to_numpy()>=o.Pan_myeloid.median())
   elif k=='Neutrophil_related':g&=det('FCGR3B')|det('CSF3R')|det('CXCR2')
   elif k=='Mast_related':g&=(det('TPSAB1')|det('TPSB2'))&(det('CPA3')|det('KIT'))
   elif k=='SPP1_mac_related':g&=det('SPP1')
   elif k=='TREM2_mac_related':g&=det('TREM2')
   gates[k]=g;check(name+'/'+k+' independently reconstructed gate',np.array_equal(g,o['gate_'+k]))
  mac=(o.Macrophage_core>=o.Macrophage_core.quantile(.75))&(o.Macrophage_core_n_detected>=3)&(o.Pan_myeloid>=o.Pan_myeloid.median());positive=o.PGAM5_raw>0;high=(mac&positive&(o.PGAM5_log1p_CP10k>=o.loc[mac&positive,'PGAM5_log1p_CP10k'].quantile(.75))).to_numpy();check(name+' fixed source gate',np.array_equal(high,o.PGAM5_mac_group.eq('PGAM5_high')));mac=mac.to_numpy()&(np.asarray(W.sum(1)).ravel()>0);positive=positive.to_numpy()
  near_high=np.asarray(W@(mac&high).astype(int)).ravel()>0;near_other=np.asarray(W@(mac&positive&~high).astype(int)).ravel()>0
  for k in ['Mac_RNA_reference','SPP1_mac_related','TREM2_mac_related']:
   target=gates[k]&~high;near=np.asarray(W@target.astype(int)).ravel()>0;row=mn[(mn['sample']==name)&(mn.target==k)].iloc[0];check(name+'/'+k+' targets excluding high',target.sum()==row.n_targets_without_high and near[mac&high].sum()==row.high_sources_with_nonhigh_neighbor)
  for k in K:
   use=high[e.row]&gates[k][e.col];src,tgt=e.row[use],e.col[use];ed=edges[(edges['sample']==name)&(edges.target==k)];check(name+'/'+k+' directed edges',set(zip(o.index[src],o.index[tgt]))==set(zip(ed.source_barcode,ed.target_barcode)))
   for group,mask in [('PGAM5_high',mac&high),('Other_detected',mac&positive&~high),('Undetected',mac&~positive)]:
    row=d[(d['sample']==name)&(d.target==k)&(d.source_group==group)].iloc[0];near=np.asarray(W@gates[k].astype(int)).ravel()>0;check(name+'/'+k+'/'+group+' counts',mask.sum()==row.n_sources and near[mask].sum()==row.with_neighbor)
   for zone,zmask in [('High_near_including_shared',near_high),('Other_only_near',~near_high&near_other),('Neither_high_nor_other_near',~near_high&~near_other)]:
    row=zones[(zones['sample']==name)&(zones.target==k)&(zones.zone==zone)].iloc[0];check(name+'/'+k+'/'+zone+' unique targets',gates[k].sum()==row.n_target and (gates[k]&zmask).sum()==row.n_targets_in_zone)
   for p in pairs.itertuples():
    li=a.var_names.get_loc(p.ligand);ri=a.var_names.get_indexer(p.receptor_complex.split('_'));support=(raw[src,li]>0)&(raw[np.ix_(tgt,ri)]>0).all(1);row=lr[(lr['sample']==name)&(lr.target==k)&(lr.pair==p.ligand+'--'+p.receptor_complex)].iloc[0];check(name+'/'+k+'/'+p.ligand+'--'+p.receptor_complex+' LR',len(src)==row.eligible_edges and support.sum()==row.supported_edges and len(np.unique(src[support]))==row.supported_sources and len(np.unique(tgt[support]))==row.supported_targets)
  for kind in s.comparison.unique():
   path=R/f'{name}_{kind}_audit.npz'
   if not path.exists():continue
   z=np.load(path);idx=o.index.get_indexer(z['barcodes']);h=z['high'];y=z['outcomes'];check(name+'/'+kind+' high labels',np.array_equal(h,high[idx]));check(name+'/'+kind+' outcomes',np.array_equal(y,np.column_stack([np.asarray(W@gates[k].astype(int)).ravel()>0 for k in K])[idx]));obs=100*(h@y/h.sum()-(1-h)@y/(len(h)-h.sum()));check(name+'/'+kind+' observed statistic',np.allclose(obs,z['observed']))
   rows=s[(s['sample']==name)&(s.comparison==kind)].set_index('target')
   for j,k in enumerate(K):
    row=rows.loc[k];null=z['null'][:,j]
    if row.status=='tested':pval=min(1,2*min((1+(null>=obs[j]-1e-10).sum())/5000,(1+(null<=obs[j]+1e-10).sum())/5000));check(name+'/'+kind+'/'+k+' equal-tail p',abs(pval-row.p)<1e-12)
 for col,mask in [('joint_q_1440',np.ones(len(s),bool)),('primary_q_480',s.comparison.eq('high_vs_other_detected')),('primary_only_q_450',s.comparison.eq('high_vs_other_detected')&s.tumor_status.eq('primary'))]:check(col+' BH',np.allclose(false_discovery_control(s.loc[mask,'p'].fillna(1).to_numpy()),s.loc[mask,col]))
 (R/'independent_audit_results.json').write_text(json.dumps(dict(passed=len(checks),checks=checks),indent=2));print('PASS',len(checks))
if __name__=='__main__':run()
