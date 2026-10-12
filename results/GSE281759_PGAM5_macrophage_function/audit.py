"""Independently reconstruct endpoints and full OLS fits; verify permutation/FDR outputs."""
import json,hashlib,itertools
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import false_discovery_control
from analyze import ROOT,REF,SAMPLES,PANELS,PROGRAMS,COMPS,dense,design
checks=[]
def check(name,value):
 checks.append(dict(check=name,passed=bool(value)))
 if not value:raise AssertionError(name)
def run():
 stat=pd.read_csv(ROOT/'program_statistics.csv');desc=pd.read_csv(ROOT/'program_group_means.csv');meta=pd.read_csv(ROOT/'patient_direction_summary.csv');patients=pd.read_csv(ROOT/'patient_equal_section_effects.csv');genes=pd.read_csv(ROOT/'marker_gene_descriptive_contrasts.csv')
 for row in json.loads((ROOT/'reference_manifest.json').read_text()):check(row['path']+' SHA',hashlib.sha256((ROOT.parent/row['path']).read_bytes()).hexdigest()==row['sha256'])
 for item in SAMPLES:
  sample=item['sample'];a=ad.read_h5ad(REF/f'{sample}_spatial_marker_subset.h5ad');o=a.obs;raw=dense(a.layers['counts']);norm=np.log1p(raw/o.total_counts.to_numpy()[:,None]*1e4);saved=pd.read_csv(ROOT/f'{sample}_spot_program_scores.csv.gz',index_col=0)
  check(sample+' normalization',np.allclose(norm,dense(a.X),rtol=3e-6,atol=3e-6));check(sample+' barcode alignment',np.array_equal(saved.index,a.obs_names));check(sample+' PGAM5 excluded from endpoints',all('PGAM5' not in gs for gs in PANELS.values()))
  mac=o.Mac_RNA_enriched.to_numpy(bool);high=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy();pos=mac&o.PGAM5_raw.gt(0).to_numpy();expected=pos&o.PGAM5_log1p_CP10k.ge(o.loc[pos,'PGAM5_log1p_CP10k'].quantile(.75)).to_numpy();check(sample+' original high definition',np.array_equal(high,expected))
  for k,gs in PANELS.items():
   score=norm[:,a.var_names.get_indexer(gs)].mean(1);check(sample+k+' score',np.allclose(score,saved[k],rtol=3e-6,atol=3e-6))
   for group in ['PGAM5_high','Other_detected','Undetected']:
    mask=o.PGAM5_mac_group.eq(group).to_numpy();row=desc[(desc['sample']==sample)&desc.program.eq(k)&desc.source_group.eq(group)].iloc[0];check(sample+k+group+' count and mean',row.n_spots==mask.sum() and np.isclose(row.mean_score,score[mask].mean(),rtol=3e-6,atol=3e-6))
  for comp in COMPS:
   eligible=mac.copy()
   if comp==COMPS[1]:eligible&=o.PGAM5_raw.gt(0).to_numpy()
   if comp==COMPS[2]:eligible&=~o.T_RNA_enriched.to_numpy(bool)
   if comp==COMPS[3]:eligible&=o.pass_qc_mt5_author.to_numpy(bool)
   ix=np.flatnonzero(eligible);h=high[ix].astype(float);path=ROOT/f'{sample}_{comp}_audit.npz'
   if comp in COMPS[:2]:
    for row in genes[(genes['sample']==sample)&genes.comparison.eq(comp)].itertuples():
     li=a.var_names.get_loc(row.gene);rr=raw[ix,li];xx=dense(a.X)[ix,li];u=o.total_counts.to_numpy()[ix];pb1=1e4*rr[h>0].sum()/u[h>0].sum();pb0=1e4*rr[h==0].sum()/u[h==0].sum()
     check(sample+comp+row.gene+' gene means/detection',np.isclose(row.high_mean,xx[h>0].mean()) and np.isclose(row.comparator_mean,xx[h==0].mean()) and np.isclose(row.high_detect_pct,100*(rr[h>0]>0).mean()) and np.isclose(row.comparator_detect_pct,100*(rr[h==0]>0).mean()))
     check(sample+comp+row.gene+' full-library pseudobulk',np.isclose(row.pseudobulk_CP10k_high,pb1) and np.isclose(row.pseudobulk_CP10k_other,pb0) and np.isclose(row.descriptive_log2_ratio_offset_0_1,np.log2((pb1+.1)/(pb0+.1))))
   for k in PROGRAMS:
    row=stat[(stat['sample']==sample)&stat.comparison.eq(comp)&stat.program.eq(k)].iloc[0];check(sample+comp+k+' counts',row.n_high==h.sum() and row.n_comparator==len(h)-h.sum())
   if not path.exists():continue
   v=np.load(path);check(sample+comp+' audit cohort',np.array_equal(v['barcodes'],a.obs_names[ix]));check(sample+comp+' high labels',np.array_equal(v['high'],h))
   H=v['assignments_first20'];code=v['groups']
   for c in np.unique(code):check(sample+comp+str(c)+' permutation group counts',np.all(H[:,code==c].sum(1)==h[code==c].sum()))
   for j,k in enumerate(PROGRAMS):
    y=dense(a[:,PANELS[k]].X).mean(1)[ix];C,Q,rank,names=design(a,ix,set(PANELS[k]));row=stat[(stat['sample']==sample)&stat.comparison.eq(comp)&stat.program.eq(k)].iloc[0]
    beta=np.linalg.lstsq(np.column_stack([C,h]),y,rcond=1e-10)[0][-1];check(sample+comp+k+' independent full OLS',np.isclose(beta,row.adjusted_beta,rtol=1e-6,atol=1e-7))
    for b in range(3):
     pb=np.linalg.lstsq(np.column_stack([C,H[b]]),y,rcond=1e-10)[0][-1];check(sample+comp+k+str(b)+' full permutation OLS',np.isclose(pb,v['null'][b,j],rtol=1e-6,atol=1e-7))
    if row.status=='tested':
     null=v['null'][:,j];obs=v['observed'][j];p=min(1,2*min((1+(null>=obs-1e-12).sum())/5000,(1+(null<=obs+1e-12).sum())/5000));check(sample+comp+k+' equal-tail p',np.isclose(p,row.p))
 check('832 planned tests',len(stat)==832);check('joint BH',np.allclose(false_discovery_control(stat.p.fillna(1).to_numpy()),stat.joint_q_832));ix=stat.comparison.eq(COMPS[0]);check('208 primary BH',np.allclose(false_discovery_control(stat.loc[ix,'p'].fillna(1).to_numpy()),stat.loc[ix,'primary_q_208']))
 for row in meta.itertuples():
  v=patients[patients.comparison.eq(COMPS[0])&patients.program.eq(row.program)].median_standardized_beta.to_numpy();null=np.array(list(itertools.product([-1,1],repeat=len(v))))@v/len(v);p=(np.abs(null)>=abs(v.mean())-1e-12).mean();check(row.program+' patient signflip',np.isclose(p,row.patient_signflip_p))
 (ROOT/'audit_results.json').write_text(json.dumps(dict(checks=len(checks),passed=sum(c['passed'] for c in checks),details=checks),indent=2)+'\n');print('PASS',len(checks))
if __name__=='__main__':run()
