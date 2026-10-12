"""Independent reconstructed gates, lattice, LR edges, randomization p and FDR."""
import json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy.spatial import cKDTree
from scipy.stats import false_discovery_control
from analyze import R,BASE,EXT,K,ALL,PANELS,SAMPLES,dense
checks=[]
def check(name,v):
 if not v:raise AssertionError(name)
 checks.append(name)
def run():
 d=pd.read_csv(R/'neighborhood_counts.csv');s=pd.read_csv(R/'spatial_statistics.csv');lr=pd.read_csv(R/'candidate_LR_codetection.csv');pairs=pd.read_csv(R/'candidate_LR_resource.csv');zones=pd.read_csv(R/'CD8_target_neighborhood_programs.csv');ed=pd.read_csv(R/'high_source_CD8_directed_edges.csv.gz')
 for v in json.loads((R/'input_manifest.json').read_text()):check(v['name']+' SHA',hashlib.sha256(open(v['local_path'],'rb').read()).hexdigest()==v['sha256'])
 for item in SAMPLES:
  name=item['sample'];a=ad.read_h5ad(R/f'{name}_CD8_fine_markers.h5ad');o=a.obs;raw=dense(a.layers['counts']);x=dense(a.X);W=a.obsp['spatial_connectivities']
  check(name+' full-library normalization',np.allclose(x,np.log1p(raw/o.total_counts.to_numpy()[:,None]*10000),rtol=1e-6,atol=1e-6))
  cd8=o.gate_All_CD8_reference.to_numpy(bool);check(name+' unchanged CD8 reference',np.array_equal(cd8,ad.read_h5ad(EXT/f'{name}_supplemental_T_markers.h5ad').obs.loc[o.index,'gate_CD8_supported']))
  xy=o[['array_row','array_col']].to_numpy(float);xy[:,0]*=np.sqrt(3);ij=cKDTree(xy).query_pairs(2.000001,output_type='ndarray');expected=set(map(tuple,ij))|set(map(tuple,ij[:,::-1]));coo=W.tocoo();check(name+' true lattice reconstruction',expected==set(zip(coo.row,coo.col)))
  core=o.Macrophage_core.to_numpy();mac=(core>=np.quantile(core,.75))&(o.Macrophage_core_n_detected.to_numpy()>=3)&(o.Pan_myeloid.to_numpy()>=o.Pan_myeloid.median());positive=o.PGAM5_raw.to_numpy()>0;pg=o.PGAM5_log1p_CP10k.to_numpy();high=mac&positive&(pg>=np.quantile(pg[mac&positive],.75));check(name+' PGAM5 source gate',np.array_equal(high,o.PGAM5_mac_group.eq('PGAM5_high')));mac&=np.asarray(W.sum(1)).ravel()>0
  det=lambda g:raw[:,a.var_names.get_loc(g)]>0
  gates={'All_CD8_reference':cd8}
  for k,gs in PANELS.items():
   ix=a.var_names.get_indexer(gs);score=x[:,ix].mean(1);n=(raw[:,ix]>0).sum(1);minimum=2 if k in ['GZMK_memory_like','Early_activation'] else 3
   g=cd8&(score>=np.quantile(score[cd8],.75))&(score>0)&(n>=minimum) if cd8.any() else np.zeros(len(o),bool)
   if k=='Naive_like':g&=(det('CCR7')|det('SELL'))&(det('TCF7')|det('LEF1'))
   elif k=='Central_memory_like':g&=det('CCR7')&det('IL7R')
   elif k=='GZMK_memory_like':g&=det('GZMK')&(det('IL7R')|det('CCL5'))
   elif k=='Cytotoxic_effector':g&=det('PRF1')|det('GZMB')
   elif k=='Terminal_effector_like':g&=det('FGFBP2')|det('CX3CR1')
   elif k=='Exhaustion_supported':g&=det('PDCD1')|det('TOX')
   elif k=='Progenitor_exhaustion_like':g&=det('TCF7')&det('PDCD1')
   elif k=='Terminal_exhaustion_like':g&=(det('HAVCR2')|det('ENTPD1'))&(det('TOX')|det('PDCD1'))
   elif k=='Resident_memory_like':g&=det('ITGAE')|det('ZNF683')|det('ITGA1')
   elif k=='Early_activation':g&=det('CD69')|det('IL2RA')|det('TNFRSF9')|det('ICOS')
   elif k=='Cycling':g&=det('MKI67')|det('TOP2A')
   elif k=='IFN_response':g&=det('IFIT1')|det('IFIT3')|det('MX1')
   check(name+'/'+k+' independently reconstructed gate',np.array_equal(g,o['gate_'+k]));gates[k]=g
  near_high=np.asarray(W@(mac&high).astype(int)).ravel()>0;near_other=np.asarray(W@(mac&positive&~high).astype(int)).ravel()>0
  for zone,mask in [('High_near_including_shared',cd8&near_high),('Other_only_near',cd8&~near_high&near_other),('Neither_high_nor_other_near',cd8&~near_high&~near_other)]:
   for k in ALL:
    row=zones[(zones['sample']==name)&(zones.zone==zone)&(zones.target==k)].iloc[0];check(name+'/'+zone+'/'+k+' target-zone counts',int(mask.sum())==row.n_CD8 and int((mask&gates[k]).sum())==row.n_program)
  target_edge=high[coo.row]&cd8[coo.col];actual=set(zip(o.index[coo.row[target_edge]],o.index[coo.col[target_edge]]));ee=ed[ed['sample']==name];check(name+' actual CD8 source-target edges',actual==set(zip(ee.source_barcode,ee.target_barcode)))
  sig=pd.read_csv(R/'exact_program_signatures.csv');ss=sig[sig['sample']==name];check(name+' exact signatures partition',ss.n_target.sum()==cd8.sum())
  for k in ALL:
   use=high[coo.row]&gates[k][coo.col];src,tgt=coo.row[use],coo.col[use]
   for p in pairs.itertuples():
    li=a.var_names.get_loc(p.ligand);ri=a.var_names.get_indexer(p.receptor_complex.split('_'));support=(raw[src,li]>0)&(raw[np.ix_(tgt,ri)]>0).all(1);row=lr[(lr['sample']==name)&(lr.target==k)&(lr.pair==p.ligand+'--'+p.receptor_complex)].iloc[0]
    check(name+'/'+k+'/'+p.ligand+'--'+p.receptor_complex+' LR edges',len(src)==row.eligible_edges and support.sum()==row.supported_edges and len(np.unique(src[support]))==row.supported_sources and len(np.unique(tgt[support]))==row.supported_targets)
  for kind in s.comparison.unique():
   path=R/f'{name}_{kind}_audit.npz'
   if not path.exists():continue
   z=np.load(path);idx=o.index.get_indexer(z['barcodes']);h=z['high'];y=z['outcomes'];check(name+'/'+kind+' high labels',np.array_equal(h,high[idx]));check(name+'/'+kind+' outcomes',np.array_equal(y,np.column_stack([np.asarray(W@gates[k].astype(int)).ravel()>0 for k in K])[idx]));obs=100*(h@y/h.sum()-(1-h)@y/(len(h)-h.sum()));check(name+'/'+kind+' observed statistic',np.allclose(obs,z['observed']))
   rows=s[(s['sample']==name)&(s.comparison==kind)].set_index('target')
   for j,k in enumerate(K):
    row=rows.loc[k];null=z['null'][:,j]
    if row.status=='tested':
     pval=min(1,2*min((1+(null>=obs[j]-1e-10).sum())/5000,(1+(null<=obs[j]+1e-10).sum())/5000));check(name+'/'+kind+'/'+k+' equal-tail p',abs(pval-row.p)<1e-12)
 for col,mask in [('joint_q_576',np.ones(len(s),bool)),('primary_q_192',s.comparison.eq('high_vs_other_detected')),('primary_only_q_180',s.comparison.eq('high_vs_other_detected')&s.tumor_status.eq('primary'))]:check(col+' BH',np.allclose(false_discovery_control(s.loc[mask,'p'].fillna(1).to_numpy()),s.loc[mask,col]))
 (R/'independent_audit_results.json').write_text(json.dumps(dict(passed=len(checks),checks=checks),indent=2));print('PASS',len(checks))
if __name__=='__main__':run()
