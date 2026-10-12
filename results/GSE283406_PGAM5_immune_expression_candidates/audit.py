import json,hashlib
import numpy as np,pandas as pd,anndata as ad
from analyze import ROOT,PANELS,PAIRS,dense
checks=[]
def check(n,v):
 checks.append(dict(check=n,passed=bool(v)))
 if not v:raise AssertionError(n)
def run():
 d=pd.read_csv(ROOT/'same_spot_RNA_overlap.csv');lr=pd.read_csv(ROOT/'candidate_LR_independent_expression.csv');resource=pd.read_csv(ROOT/'candidate_LR_resource_snapshot.csv')
 for r in json.loads((ROOT/'input_manifest.json').read_text()):check(r['path']+' SHA',hashlib.sha256((ROOT/r['path']).read_bytes()).hexdigest()==r['sha256'])
 for l,r in PAIRS:check(l+r+' resource definition',((resource.source_genesymbol==l)&(resource.target_genesymbol==r)).any())
 for section in ['J259X','J259B']:
  a=ad.read_h5ad(ROOT/'inputs'/f'{section}_AVITI_expression_marker_subset.h5ad');o=pd.read_csv(ROOT/f'{section}_immune_candidate_annotations.csv.gz',index_col=0).loc[a.obs_names];raw=dense(a.layers['counts']);norm=np.log1p(raw/a.obs.total_counts.to_numpy()[:,None]*1e4);check(section+' full panel normalization',np.allclose(norm,dense(a.X),rtol=3e-6,atol=3e-6));check(section+' no invented spatial geometry',not('spatial' in a.obsm or 'spatial_connectivities' in a.obsp));check(section+' original PGAM5 raw',np.array_equal(raw[:,a.var_names.get_loc('PGAM5')],o.PGAM5_raw))
  m=o.Mac_RNA_enriched.to_numpy(bool);pos=m&o.PGAM5_raw.gt(0).to_numpy();cut=o.loc[pos,'PGAM5_log1p_CP10k'].quantile(.75);h=pos&o.PGAM5_log1p_CP10k.ge(cut).to_numpy();check(section+' fixed high Q75',np.array_equal(h,o.PGAM5_mac_RNA_group.eq('PGAM5_high')))
  for k,gs in PANELS.items():
   vals=dense(a[:,gs].X);counts=dense(a[:,gs].layers['counts']);check(section+k+' score and detection',np.allclose(vals.mean(1),o['score_'+k]) and np.array_equal((counts>0).sum(1),o['detected_'+k]));gate=o['gate_'+k].to_numpy(bool)
   prev={'All_T':'T_RNA_enriched','Exhaustion_T':'Exhaustion_core_supported','Activation_T':'Activation_early_supported','Cytotoxic_T':'Effector_cytotoxic_supported'}
   if k in prev:expected=a.obs[prev[k]].to_numpy(bool)
   else:
    score=o['score_'+k].to_numpy();expected=(score>=np.quantile(score,.75))&(score>0)&((counts>0).sum(1)>=(3 if k=='Plasma_related' else 2))
    def has(g):return (raw[:,a.var_names.get_indexer(g)]>0).any(1)
    if k=='B_related':expected&=has(['MS4A1','CD19'])
    if k=='Plasma_related':expected&=has(['MZB1','JCHAIN'])
    if k=='cDC2_related':expected&=has(['CD1C','CLEC10A'])
    if k=='cDC1_related':expected&=has(['CLEC9A'])&has(['XCR1'])
    if k=='Inflammatory_myeloid':expected&=has(['VCAN'])&o.Pan_myeloid.ge(o.Pan_myeloid.median()).to_numpy()
   check(section+k+' independently reconstructed target gate',np.array_equal(expected,gate))
   for group in ['PGAM5_high','Other_detected','Undetected']:
    mask=o.PGAM5_mac_RNA_group.eq(group).to_numpy();row=d[d.section.eq(section)&d.target.eq(k)&d.source_group.eq(group)].iloc[0];check(section+k+group+' same-spot counts',row.n_sources==mask.sum() and row.total_target_spots==gate.sum() and row.same_spot_overlap==(mask&gate).sum())
   for l,r in PAIRS:
    row=lr[lr.section.eq(section)&lr.target.eq(k)&lr.pair.eq(l+'--'+r)].iloc[0];lig=raw[:,a.var_names.get_loc(l)]>0;rec=raw[:,a.var_names.get_loc(r)]>0;check(section+k+l+r+' independent ligand/receptor counts',row.high_ligand_detected==(h&lig).sum() and row.target_receptor_detected==(gate&rec).sum() and row.source_target_gate_overlap==(h&gate).sum())
 check('no spatial pairs claimed',not lr.spatial_pairing_available.any());check('9 candidate classes per section',len(d)==54);check('1 biological patient documented',json.loads((ROOT/'analysis_config.json').read_text())['biological_patients']==1)
 (ROOT/'audit_results.json').write_text(json.dumps(dict(checks=len(checks),passed=sum(c['passed'] for c in checks),details=checks),indent=2)+'\n');print('PASS',len(checks))
if __name__=='__main__':run()
