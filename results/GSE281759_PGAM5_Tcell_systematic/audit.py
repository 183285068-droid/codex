import json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy.spatial import cKDTree
from scipy.stats import false_discovery_control
from analyze import R,FINE,EXT,K,MAP,SAMPLES,dense
checks=[]
def check(name,v):
 assert v,name
 checks.append(name)
def run():
 st=pd.read_csv(R/'spatial_statistics.csv');lr=pd.read_csv(R/'candidate_LR_codetection.csv');pairs=pd.read_csv(R/'candidate_LR_resource.csv');sig=pd.read_csv(R/'CD8_four_program_signatures.csv');sens=pd.read_csv(R/'Treg_identity_sensitivity.csv')
 for v in json.loads((R/'reference_manifest.json').read_text()):check(v['path']+' SHA',hashlib.sha256((R.parent/v['path']).read_bytes()).hexdigest()==v['sha256'])
 for item in SAMPLES:
  name=item['sample'];a=ad.read_h5ad(R/f'{name}_focused_T_markers.h5ad');f=ad.read_h5ad(FINE/f'{name}_CD8_fine_markers.h5ad');e=ad.read_h5ad(EXT/f'{name}_supplemental_T_markers.h5ad');o=a.obs;raw=dense(a.layers['counts']);x=dense(a.X);W=a.obsp['spatial_connectivities']
  check(name+' full-library normalization',np.allclose(x,np.log1p(raw/o.total_counts.to_numpy()[:,None]*10000),rtol=1e-6,atol=1e-6));xy=o[['array_row','array_col']].to_numpy(float);xy[:,0]*=np.sqrt(3);ij=cKDTree(xy).query_pairs(2.000001,output_type='ndarray');expected=set(map(tuple,ij))|set(map(tuple,ij[:,::-1]));ed=W.tocoo();check(name+' lattice',expected==set(zip(ed.row,ed.col)))
  gates={k:f.obs[['gate_'+v for v in values]].any(axis=1).to_numpy() for k,values in MAP.items()};t=o.T_RNA_enriched.to_numpy(bool);ix=a.var_names.get_indexer(['FOXP3','IL2RA','CTLA4','IKZF2','TIGIT']);score=x[:,ix].mean(1);gates['Treg_supported']=t&(score>=np.quantile(score[t],.75))&(score>0)&((raw[:,ix]>0).sum(1)>=2)&(raw[:,a.var_names.get_loc('FOXP3')]>0)
  for k in K:check(name+'/'+k+' union or Treg reconstruction',np.array_equal(gates[k],o['gate_'+k]))
  check(name+' Treg unchanged',np.array_equal(gates['Treg_supported'],e.obs.gate_Treg_supported))
  mac=(o.Macrophage_core>=o.Macrophage_core.quantile(.75))&(o.Macrophage_core_n_detected>=3)&(o.Pan_myeloid>=o.Pan_myeloid.median());positive=o.PGAM5_raw.gt(0);high=(mac&positive&o.PGAM5_log1p_CP10k.ge(o.loc[mac&positive,'PGAM5_log1p_CP10k'].quantile(.75))).to_numpy();check(name+' source gate',np.array_equal(high,o.PGAM5_mac_group.eq('PGAM5_high')));nh=np.asarray(W@high.astype(int)).ravel()>0
  check(name+' signature partition',sig[sig['sample']==name].n_CD8_targets.sum()==o.gate_All_CD8_reference.sum())
  for label,mask in [('Treg_CD4_supported',gates['Treg_supported']&e.obs.gate_CD4_related.to_numpy(bool)),('Treg_CD8_overlap',gates['Treg_supported']&o.gate_All_CD8_reference.to_numpy(bool))]:
   row=sens[(sens['sample']==name)&(sens.label==label)].iloc[0];near=np.asarray(W@mask.astype(int)).ravel()>0;check(name+'/'+label+' sensitivity',mask.sum()==row.n_targets and (mask&nh).sum()==row.high_near_unique_targets and near[high].sum()==row.high_sources_with_neighbor)
  for k in K:
   use=high[ed.row]&gates[k][ed.col];src,tgt=ed.row[use],ed.col[use]
   for p in pairs.itertuples():
    li=a.var_names.get_loc(p.ligand);ri=a.var_names.get_indexer(p.receptor_complex.split('_'));support=(raw[src,li]>0)&(raw[np.ix_(tgt,ri)]>0).all(1);row=lr[(lr['sample']==name)&(lr.target==k)&(lr.pair==p.ligand+'--'+p.receptor_complex)].iloc[0];check(name+'/'+k+'/'+p.ligand+'--'+p.receptor_complex+' LR',len(src)==row.eligible_edges and support.sum()==row.supported_edges and len(np.unique(src[support]))==row.supported_sources and len(np.unique(tgt[support]))==row.supported_targets)
  for kind in st.comparison.unique():
   path=R/f'{name}_{kind}_audit.npz'
   if not path.exists():continue
   z=np.load(path);idx=o.index.get_indexer(z['barcodes']);h=z['high'];y=z['outcomes'];check(name+'/'+kind+' labels',np.array_equal(h,high[idx]));check(name+'/'+kind+' outcomes',np.array_equal(y,np.column_stack([np.asarray(W@gates[k].astype(int)).ravel()>0 for k in K])[idx]));obs=100*(h@y/h.sum()-(1-h)@y/(len(h)-h.sum()));check(name+'/'+kind+' observed',np.allclose(obs,z['observed']));rows=st[(st['sample']==name)&(st.comparison==kind)].set_index('target')
   for j,k in enumerate(K):
    row=rows.loc[k];null=z['null'][:,j]
    if row.status=='tested':pval=min(1,2*min((1+(null>=obs[j]-1e-10).sum())/5000,(1+(null<=obs[j]+1e-10).sum())/5000));check(name+'/'+kind+'/'+k+' p',abs(pval-row.p)<1e-12)
 for col,mask in [('joint_q_240',np.ones(len(st),bool)),('primary_q_80',st.comparison.eq('high_vs_other_detected')),('primary_only_q_75',st.comparison.eq('high_vs_other_detected')&st.tumor_status.eq('primary'))]:check(col+' BH',np.allclose(false_discovery_control(st.loc[mask,'p'].fillna(1).to_numpy()),st.loc[mask,col]))
 (R/'independent_audit_results.json').write_text(json.dumps(dict(passed=len(checks),checks=checks),indent=2));print('PASS',len(checks))
if __name__=='__main__':run()
