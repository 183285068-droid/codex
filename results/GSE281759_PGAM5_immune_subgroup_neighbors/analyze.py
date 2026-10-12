"""Exploratory mixed-spot spatial neighborhoods; no resolved cell interaction claims."""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.stats import false_discovery_control
ROOT=Path(__file__).resolve().parent
REF=ROOT.parent/'GSE281759_spatial_HCC_PGAM5'
SAMPLES=json.loads((REF/'sample_scope.json').read_text())
B=4999; SEED=281759102
PANELS={'All_T':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2'],
 'Exhausted_T':['PDCD1','LAG3','TIGIT','CTLA4','TOX','CXCL13'],
 'Activated_T':['CD69','IL2RA','TNFRSF9','ICOS','CD40LG','IL2'],
 'Cytotoxic_T':['IFNG','PRF1','GZMB','GZMA','CCL5','TBX21'],
 'B_related':['MS4A1','CD79A','CD79B','CD19','CD22'],
 'Plasma_related':['MZB1','JCHAIN','SDC1','IGKC'],
 'cDC2_related':['CD1C','FCER1A','CLEC10A'],
 'cDC1_related':['CLEC9A','XCR1'],
 'Inflammatory_myeloid':['FCN1','VCAN','S100A8','S100A9','IL1B']}
TARGETS=list(PANELS)
def dense(x):return np.asarray(x.toarray() if sparse.issparse(x) else x)
def annotate(a):
 o=a.obs; raw=dense(a.layers['counts']); x=dense(a.X); gates={}; cuts={}
 for k,gs in PANELS.items():
  ix=a.var_names.get_indexer(gs);assert (ix>=0).all()
  score=x[:,ix].mean(1);det=(raw[:,ix]>0).sum(1)
  o['score_'+k]=score;o['detected_'+k]=det
  if k=='All_T':gates[k]=o.T_RNA_enriched.to_numpy(bool);continue
  prev={'Exhausted_T':'Exhaustion_core_supported','Activated_T':'Activation_early_supported','Cytotoxic_T':'Effector_cytotoxic_supported'}
  if k in prev:gates[k]=o[prev[k]].to_numpy(bool);continue
  cut=float(np.quantile(score,.75));cuts[k]=cut
  mask=(score>=cut)&(score>0)&(det>=(3 if k=='Plasma_related' else 2))
  def has(gs):return (raw[:,a.var_names.get_indexer(gs)]>0).any(1)
  if k=='B_related':mask&=has(['MS4A1','CD19'])
  if k=='Plasma_related':mask&=has(['MZB1','JCHAIN'])
  if k=='cDC2_related':mask&=has(['CD1C','CLEC10A'])
  if k=='cDC1_related':mask&=has(['CLEC9A'])&has(['XCR1'])
  if k=='Inflammatory_myeloid':mask&=has(['FCN1','VCAN'])&(o.Pan_myeloid>=o.Pan_myeloid.median()).to_numpy()
  gates[k]=mask
 for k,v in gates.items():o['gate_'+k]=v
 return gates,cuts
def strata(o,ix):
 xy=o[['array_col','array_row']].to_numpy();tile=np.minimum(((xy[ix]-xy.min(0))/(np.ptp(xy,axis=0)+1)*3).astype(int),2)
 code=tile[:,0]+3*tile[:,1];groups=[]
 for c in np.unique(code):
  z=np.flatnonzero(code==c);depth=(o.total_counts.to_numpy()[ix[z]]>np.median(o.total_counts.to_numpy()[ix[z]])).astype(int)
  for d in [0,1]:
   zz=z[depth==d]
   if not len(zz):continue
   core=(o.Macrophage_core.to_numpy()[ix[zz]]>np.median(o.Macrophage_core.to_numpy()[ix[zz]])).astype(int)
   for m in [0,1]:
    g=zz[core==m]
    if len(g):groups.append(g)
 return groups
def evaluate(o,Y,eligible,item,kind,rng):
 ix=np.flatnonzero(eligible);h=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy()[ix].astype(float);y=Y[ix].astype(float);n1=int(h.sum());n0=len(h)-n1
 base=dict(**item,comparison=kind,n_high=n1,n_comparator=n0)
 if n1<5 or n0<10:return [dict(**base,target=k,status='not_testable',reason='high<5 or comparator<10',p=np.nan,delta_pp=np.nan) for k in TARGETS]
 groups=strata(o,ix);mat=np.tile(h,(B,1));code=np.zeros(len(ix),int)
 for c,g in enumerate(groups):mat[:,g]=rng.permuted(mat[:,g],axis=1);code[g]=c
 null=100*(mat@y/n1-(1-mat)@y/n0);obs=100*(h@y/n1-(1-h)@y/n0)
 np.savez_compressed(ROOT/f'{item["sample"]}_{kind}_audit.npz',barcodes=o.index[ix].to_numpy(str),high=h,outcomes=y,groups=code,observed=obs,null=null)
 rows=[]
 for j,k in enumerate(TARGETS):
  ok=np.ptp(null[:,j])>1e-10
  p=min(1.,2*min((1+(null[:,j]>=obs[j]-1e-10).sum())/(B+1),(1+(null[:,j]<=obs[j]+1e-10).sum())/(B+1))) if ok else np.nan
  rows.append(dict(**base,target=k,status='tested' if ok else 'not_testable',reason='' if ok else 'constant conditional null',delta_pp=obs[j],high_pct=100*y[h>0,j].mean(),comparator_pct=100*y[h==0,j].mean(),null_mean_pp=null[:,j].mean(),p=p))
 return rows
def run():
 config=dict(seed=SEED,permutations=B,targets=PANELS,comparisons=['high_vs_other_detected','high_vs_all_other_mac','nonT_high_vs_other_detected'],planned_tests=16*9*3,primary_tests=16*9,source='previous fixed Mac RNA gate and section-specific PGAM5-positive Q75, ties retained',null='within 3x3 spatial tiles x source UMI halves x macrophage-core halves; two-sided equal-tail +1 p',target_overlap=True,missing_markers=['NK','Treg','naive/memory T','pDC','neutrophil','mast','CD4/CD8 identity'],causal_claim=False)
 (ROOT/'analysis_config.json').write_text(json.dumps(config,indent=2)+'\n')
 rows=[];des=[];summ=[];lr=[];manifest=[];cuts=[];markers=[]
 pairs=pd.read_csv(REF/'selected_LR_pairs.csv');pairs.to_csv(ROOT/'candidate_LR_resource.csv',index=False)
 for item in SAMPLES:
  s=item['sample'];path=REF/f'{s}_spatial_marker_subset.h5ad';manifest.append(dict(path=str(path.relative_to(ROOT.parent)),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
  a=ad.read_h5ad(path);o=a.obs;gates,threshold=annotate(a);cuts.append(dict(**item,**threshold,PGAM5_high=float(a.uns['PGAM5_high_threshold_log1p_CP10k'])))
  W=a.obsp['spatial_connectivities'];degree=np.asarray(W.sum(1)).ravel();Y=np.column_stack([W@gates[k].astype(int)>0 for k in TARGETS]);mac=o.Mac_RNA_enriched.to_numpy(bool)&(degree>0);positive=o.PGAM5_raw.to_numpy()>0
  high=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy();rng=np.random.default_rng(SEED+int(item['gsm'][3:]))
  for kind,eligible in [('high_vs_other_detected',mac&positive),('high_vs_all_other_mac',mac),('nonT_high_vs_other_detected',mac&positive&~o.T_RNA_enriched.to_numpy(bool))]:rows.extend(evaluate(o,Y,eligible,item,kind,rng))
  summ.append(dict(**item,qc_spots=len(a),mac_spots=int(o.Mac_RNA_enriched.sum()),n_high_all=int(high.sum()),n_high_connected=int((high&mac).sum()),n_other_detected=int((mac&positive&~high).sum()),n_undetected=int((mac&~positive).sum())))
  for group,mask in [('PGAM5_high',mac&high),('Other_detected',mac&positive&~high),('Undetected',mac&~positive)]:
   for j,k in enumerate(TARGETS):des.append(dict(**item,source_group=group,target=k,n_sources=int(mask.sum()),with_neighbor=int(Y[mask,j].sum()),neighbor_pct=100*Y[mask,j].mean() if mask.any() else np.nan,n_target_spots=int(gates[k].sum())))
  raw=dense(a.layers['counts']);norm=np.expm1(dense(a.X));edge=W.tocoo();ss,tt=edge.row,edge.col
  for k,gs in PANELS.items():
   for g in gs:
    v=raw[:,a.var_names.get_loc(g)];markers.append(dict(**item,target=k,gene=g,all_detection_pct=100*(v>0).mean(),target_detection_pct=100*(v[gates[k]]>0).mean() if gates[k].any() else np.nan))
   use=high[ss]&gates[k][tt];src,tgt=ss[use],tt[use]
   for p in pairs.itertuples():
    li=a.var_names.get_loc(p.ligand);ri=a.var_names.get_indexer(p.receptor_complex.split('_'));assert (ri>=0).all()
    support=(raw[src,li]>0)&(raw[np.ix_(tgt,ri)]>0).all(1)
    score=np.sqrt(norm[src,li]*norm[np.ix_(tgt,ri)].min(1))
    lr.append(dict(**item,target=k,pair=p.ligand+'--'+p.receptor_complex,eligible_edges=len(src),unique_sources=len(np.unique(src)),unique_targets=len(np.unique(tgt)),supported_edges=int(support.sum()),supported_sources=len(np.unique(src[support])),supported_targets=len(np.unique(tgt[support])),codetection_pct=100*support.mean() if len(src) else np.nan,mean_geometric_CP10k=score.mean() if len(src) else np.nan))
  for j,k in enumerate(TARGETS):o['neighbor_'+k]=Y[:,j]
  o.to_csv(ROOT/f'{s}_annotations.csv.gz',float_format='%.17g')
  print(s,'high',high.sum(),flush=True)
 stat=pd.DataFrame(rows);stat['joint_q_432']=false_discovery_control(stat.p.fillna(1).to_numpy());stat['primary_q_144']=np.nan
 ix=stat.comparison.eq('high_vs_other_detected');stat.loc[ix,'primary_q_144']=false_discovery_control(stat.loc[ix,'p'].fillna(1).to_numpy())
 primary=ix&stat.tumor_status.eq('primary');stat['primary_only_q_135']=np.nan;stat.loc[primary,'primary_only_q_135']=false_discovery_control(stat.loc[primary,'p'].fillna(1).to_numpy())
 stat.to_csv(ROOT/'spatial_statistics.csv',index=False)
 dd=pd.DataFrame(des);dd.to_csv(ROOT/'neighborhood_counts.csv',index=False)
 pd.DataFrame(summ).to_csv(ROOT/'section_summary.csv',index=False);pd.DataFrame(lr).to_csv(ROOT/'candidate_LR_codetection.csv',index=False);pd.DataFrame(cuts).to_csv(ROOT/'thresholds.csv',index=False);pd.DataFrame(markers).to_csv(ROOT/'marker_detection.csv',index=False)
 patient=dd.groupby(['patient','source_group','target'],as_index=False).agg(n_sources=('n_sources','sum'),with_neighbor=('with_neighbor','sum'),n_target_spots=('n_target_spots','sum'));patient['neighbor_pct']=100*patient.with_neighbor/patient.n_sources;patient.to_csv(ROOT/'patient_neighborhood_descriptive.csv',index=False)
 effects=stat[stat.status.eq('tested')].groupby(['patient','comparison','target'],as_index=False).agg(median_delta_pp=('delta_pp','median'),n_testable_sections=('sample','size'));effects.to_csv(ROOT/'patient_equal_section_effects.csv',index=False)
 (ROOT/'reference_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
if __name__=='__main__':run()
