import json,hashlib
import numpy as np,pandas as pd,anndata as ad
from analyze import R,SOURCE,K,PANELS,PAIRS,dense
checks=[]
def check(n,v):
 assert v,n
 checks.append(n)
def run():
 ref=json.loads((R/'reference_manifest.json').read_text())[0];p=R/ref['path'];check('Input SHA',hashlib.sha256(p.read_bytes()).hexdigest()==ref['sha256']);a=ad.read_h5ad(p);o=pd.read_csv(R/'J259X_T_annotations.csv.gz',index_col=0).loc[a.obs_names];raw=dense(a.layers['counts']);x=dense(a.X);check('Full library normalization',np.allclose(x,np.log1p(raw/o.total_counts.to_numpy()[:,None]*10000),rtol=1e-6,atol=1e-6));check('No spatial coordinates or graph invented','spatial' not in a.obsm and 'spatial_connectivities' not in a.obsp)
 pos=o.Mac_RNA_enriched.to_numpy(bool)&(o.PGAM5_raw.to_numpy()>0);high=pos&(o.PGAM5_log1p_CP10k.to_numpy()>=np.quantile(o.PGAM5_log1p_CP10k.to_numpy()[pos],.75));check('High Q75 fixed',np.array_equal(high,o.PGAM5_mac_RNA_group.eq('PGAM5_high')))
 gates={};tix=a.var_names.get_indexer(PANELS['All_T']);ts=x[:,tix].mean(1);gates['All_T']=(ts>=np.quantile(ts,.75))&((raw[:,tix]>0).sum(1)>=2)
 for k in K[1:]:
  ix=a.var_names.get_indexer(PANELS[k]);score=x[:,ix].mean(1);gates[k]=gates['All_T']&(score>=np.quantile(score[gates['All_T']],.75))&(score>0)&((raw[:,ix]>0).sum(1)>=2)
 d=pd.read_csv(R/'same_barcode_T_program_overlap.csv');c=pd.read_csv(R/'conditional_T_program_descriptive.csv');lr=pd.read_csv(R/'ligand_receptor_independent_expression.csv');resource=pd.read_csv(R/'candidate_LR_resource_snapshot.csv')
 for k in K:
  check(k+' independently reconstructed gate',np.array_equal(gates[k],o['gate_'+k]))
  for group in ['PGAM5_high','Other_detected','Undetected']:
   mask=o.PGAM5_mac_RNA_group.eq(group).to_numpy();row=d[d.target.eq(k)&d.source_group.eq(group)].iloc[0];check(k+group+' same barcode counts',row.n_sources==mask.sum() and row.n_target==gates[k].sum() and row.same_barcode_overlap==(mask&gates[k]).sum());row=c[c.target.eq(k)&c.source_group.eq(group)].iloc[0];check(k+group+' conditional denominator',row.n_T_coexpressing_sources==(mask&gates['All_T']).sum() and row.state_supported==(mask&gates[k]).sum())
   for l,r in PAIRS:
    row=lr[lr.target.eq(k)&lr.source_group.eq(group)&lr.pair.eq(l+'--'+r)].iloc[0];check(k+group+l+r+' independent expression',row.ligand_detected==(raw[mask,a.var_names.get_loc(l)]>0).sum() and row.receptor_detected==(raw[gates[k],a.var_names.get_loc(r)]>0).sum())
 for l,r in PAIRS:check(l+r+' resource definition',((resource.source_genesymbol==l)&(resource.target_genesymbol==r)).any())
 signatures=pd.read_csv(R/'T_program_signatures.csv');check('Signature partition',signatures.n_T_targets.sum()==575 and signatures.high_same_barcode.sum()==7);check('No LR pairing',not lr.source_target_pairing.any());check('One biological patient',json.loads((R/'analysis_config.json').read_text())['patients']==1)
 (R/'audit_results.json').write_text(json.dumps(dict(passed=len(checks),checks=checks),indent=2));print('PASS',len(checks))
if __name__=='__main__':run()
