"""Patient is the independent unit. Mixed-spot plasma-related RNA functional exploration."""
from pathlib import Path
import json,itertools,hashlib,shutil
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import false_discovery_control,ttest_1samp,t as student_t
from scipy import sparse
ROOT=Path(__file__).resolve().parent;REF=ROOT.parent/'GSE281759_spatial_HCC_PGAM5';ANN=ROOT.parent/'GSE281759_PGAM5_immune_subgroup_neighbors'
SAMPLES=json.loads((ROOT/'sample_scope.json').read_text());PATIENTS=['A009','B013','C004','C008']
PANELS={'MHCII':['HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','CD74'],
 'Regulatory_mediators':['IL10','TGFB1','LGALS9'],
 'Checkpoint_ligands':['CD274','PDCD1LG2','PVR','NECTIN2'],
 'Costimulation':['CD80','CD86','ICOSLG','TNFSF9','TNFSF4'],
 'Inflammatory_chemokines':['CCL2','CCL3','CCL4','CCL8','CXCL8'],
 'CXCR3_ligands':['CXCL9','CXCL10','CXCL11'],
 'Cycling':['MKI67','TOP2A','UBE2C','CENPF','BIRC5','PCNA'],
 'Hypoxia':['CA9','VEGFA','SLC2A1','LDHA','BNIP3','ADM']}
MIX={'Plasma':['MZB1','JCHAIN','SDC1','IGKC'],
 'Mac':['C1QA','C1QB','C1QC','CD68','CSF1R','MERTK','MSR1'],
 'T':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2'],
 'B':['MS4A1','CD79A','CD79B','CD19','CD22'],
 'Epithelial':['EPCAM','KRT8','KRT18','KRT19'],
 'Hepatocyte':['ALB','APOA1','APOA2','TTR','ASGR1','CYP2E1']}
COMPS=['high_only_vs_other_detected_only','high_any_vs_no_high','strict_high_only_vs_other_detected_only']
def dense(x):return np.asarray(x.toarray() if sparse.issparse(x) else x)
def score(a,gs):return dense(a[:,gs].X).mean(1)
def nuisance(a,ix):
 o=a.obs;cols=[np.ones(len(ix)),np.log1p(o.total_counts.to_numpy()[ix]),o.pct_mt.to_numpy()[ix]]
 for gs in MIX.values():cols.append(score(a,gs)[ix])
 xy=o[['array_col','array_row']].to_numpy();tiles=np.minimum(((xy-xy.min(0))/(np.ptp(xy,axis=0)+1)*3).astype(int),2);tile=(tiles[:,0]+3*tiles[:,1])[ix]
 for t in np.unique(tile)[1:]:cols.append((tile==t).astype(float))
 C=np.column_stack(cols)
 for j in range(1,C.shape[1]):
  sd=C[:,j].std();C[:,j]=(C[:,j]-C[:,j].mean())/sd if sd>1e-10 else 0
 u,s,v=np.linalg.svd(C,full_matrices=False);r=int((s>s[0]*1e-10).sum());return C,u[:,:r],r
def run():
 cfg=dict(unit='one independent paired effect per patient; four patients',panels=PANELS,primary_comparison=COMPS[0],primary_estimator='raw mean log1p(CP10k) difference, equal section means within each patient using the same eligible sections on both sides',primary_tests=8,minimum_targets_per_group_per_section=3,adjusted='mean nuisance-residual contrasts, fit to all fixed plasma targets per section',covariates=MIX,other_covariates=['log UMI','MT percent','3x3 spatial tile'],sensitivity=['any-high vs no-high','MZB1 and JCHAIN co-detected targets','primary tumors only (exclude C004 recurrence)'],patient_test='exact two-sided sign-flip of patient paired differences, 2^n assignments, no pooled-spot p',minimum_patient_n_for_test=4,minimum_two_sided_p_n4=.125,causal_claim=False,cell_subtypes='JCHAIN/MZB1 co-detection RNA states are not IgA/IgG or resolved plasma-cell subtypes')
 (ROOT/'analysis_config.json').write_text(json.dumps(cfg,indent=2)+'\n');sections=[];counts=[];states=[];coverage=[];manifest=[]
 for item in SAMPLES:
  s=item['sample'];dst=ROOT/'inputs'/f'{s}_plasma_marker_inputs.h5ad'
  if not dst.exists():
   a=ad.read_h5ad(REF/f'{s}_spatial_marker_subset.h5ad');o=pd.read_csv(ANN/f'{s}_annotations.csv.gz',index_col=0).loc[a.obs_names];a.obs['gate_Plasma_related']=o.gate_Plasma_related.to_numpy(bool);a.write_h5ad(dst,compression='gzip')
  a=ad.read_h5ad(dst);o=a.obs;G=a.obsp['spatial_connectivities'];pl=o.gate_Plasma_related.to_numpy(bool);pi=np.flatnonzero(pl);h=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy();other=o.PGAM5_mac_group.eq('Other_detected').to_numpy();nh=np.asarray(G@h.astype(int)).ravel()[pi];no=np.asarray(G@other.astype(int)).ravel()[pi];strict=(dense(a[pi,['MZB1','JCHAIN']].layers['counts'])>0).all(1)
  group=np.where(nh>0,np.where(no>0,'Shared','High_only'),np.where(no>0,'Other_only','No_detected_mac_neighbor'))
  raw=dense(a[pi,['JCHAIN','MZB1']].layers['counts']);state=np.where((raw>0).all(1),'JCHAIN_MZB1_both',np.where(raw[:,0]>0,'JCHAIN_only',np.where(raw[:,1]>0,'MZB1_only','Neither')))
  Y=np.column_stack([score(a,gs)[pi] for gs in PANELS.values()]);C,Q,rank=nuisance(a,pi) if len(pi) else (None,None,0);res=Y-Q@(Q.T@Y) if len(pi)-rank>=15 else np.full_like(Y,np.nan)
  q=pd.DataFrame(Y,index=a.obs_names[pi],columns=PANELS);q['group']=group;q['state']=state;q['strict_plasma']=strict;q['n_high_neighbors']=nh;q['n_other_neighbors']=no
  for k,v in MIX.items():q['mix_'+k]=score(a,v)[pi]
  q.to_csv(ROOT/f'{s}_target_profiles.csv.gz',float_format='%.17g')
  counts.append(dict(**item,plasma_targets=len(pi),high_any=int((nh>0).sum()),high_only=int((group=='High_only').sum()),other_only=int((group=='Other_only').sum()),shared=int((group=='Shared').sum()),remote=int((group=='No_detected_mac_neighbor').sum()),strict_targets=int(strict.sum())))
  for gp in ['High_only','Shared','Other_only','No_detected_mac_neighbor']:
   for st in ['JCHAIN_MZB1_both','JCHAIN_only','MZB1_only','Neither']:
    use=(group==gp)&(state==st);states.append(dict(**item,group=gp,state=st,n_targets=int(use.sum()),group_total=int((group==gp).sum())))
  for comp in COMPS:
   if comp=='high_any_vs_no_high':left=nh>0;right=nh==0
   else:left=group=='High_only';right=group=='Other_only'
   if comp.startswith('strict'):left&=strict;right&=strict
   n1=int(left.sum());n0=int(right.sum());eligible=n1>=3 and n0>=3
   for j,k in enumerate(PANELS):
    for method,z in [('raw',Y[:,j]),('adjusted',res[:,j])]:
     valid=eligible and np.isfinite(z).all();sections.append(dict(**item,comparison=comp,method=method,program=k,n_high=n1,n_other=n0,status='eligible' if valid else 'not_eligible',high_mean=float(z[left].mean()) if n1 and np.isfinite(z[left]).all() else np.nan,other_mean=float(z[right].mean()) if n0 and np.isfinite(z[right]).all() else np.nan,delta=float(z[left].mean()-z[right].mean()) if valid else np.nan,covariate_rank=rank))
  for panel,gs in {'Plasma_differentiation':['PRDM1','XBP1','IRF4','TNFRSF17'],'Secretory_ER':['HSPA5','HSP90B1','SEC11C','FKBP11'],'Isotypes':['IGHG1','IGHG2','IGHG3','IGHG4','IGHA1','IGHA2','IGHM'],'Trafficking':['CXCR4','ITGA4','ITGB1'],**PANELS}.items():
   for g in gs:coverage.append(dict(**item,panel=panel,gene=g,present=g in a.var_names))
  manifest.append(dict(path=str(dst.relative_to(ROOT)),bytes=dst.stat().st_size,sha256=hashlib.sha256(dst.read_bytes()).hexdigest()));np.savez_compressed(ROOT/f'{s}_residual_audit.npz',target_barcodes=a.obs_names[pi].to_numpy(str),scores=Y,residuals=res,covariate_rank=rank);print(s,'groups',counts[-1]['high_only'],counts[-1]['other_only'],flush=True)
 sec=pd.DataFrame(sections);sec.to_csv(ROOT/'section_effects_descriptive.csv',index=False);pd.DataFrame(counts).to_csv(ROOT/'section_target_counts.csv',index=False);pd.DataFrame(states).to_csv(ROOT/'RNA_state_counts.csv',index=False);pd.DataFrame(coverage).to_csv(ROOT/'functional_marker_coverage.csv',index=False);(ROOT/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 patient_rows=[]
 for sensitivity,subset in [('all_tumor_sections',sec),('primary_only',sec[sec.tumor_status.eq('primary')])]:
  for comp in COMPS:
   for method in ['raw','adjusted']:
    for program in PANELS:
     for patient in PATIENTS:
      sub=subset[subset.comparison.eq(comp)&subset.method.eq(method)&subset.program.eq(program)&subset.patient.eq(patient)&subset.status.eq('eligible')]
      patient_rows.append(dict(patient=patient,comparison=comp,method=method,program=program,sensitivity=sensitivity,n_eligible_sections=len(sub),n_high_targets=int(sub.n_high.sum()),n_other_targets=int(sub.n_other.sum()),high_section_equal_mean=sub.high_mean.mean(),other_section_equal_mean=sub.other_mean.mean(),patient_delta=sub.delta.mean(),sections=';'.join(sub['sample'])))
 pat=pd.DataFrame(patient_rows);pat.to_csv(ROOT/'patient_paired_effects.csv',index=False);tests=[]
 for (sens,comp,method,k),sub in pat.groupby(['sensitivity','comparison','method','program'],sort=False):
  v=sub.patient_delta.dropna().to_numpy();n=len(v);obs=v.mean() if n else np.nan;null=np.array(list(itertools.product([-1,1],repeat=n)))@v/n if n else np.array([]);p=float((np.abs(null)>=abs(obs)-1e-12).mean()) if n==4 else np.nan
  se=v.std(ddof=1)/np.sqrt(n) if n==4 else np.nan;tp=float(ttest_1samp(v,0).pvalue) if n==4 and se>1e-12 else np.nan;half=student_t.ppf(.975,3)*se if n==4 else np.nan
  tests.append(dict(sensitivity=sens,comparison=comp,method=method,program=k,n_patients=n,mean_patient_delta=obs,median_patient_delta=np.median(v) if n else np.nan,positive_patients=int((v>0).sum()),negative_patients=int((v<0).sum()),exact_signflip_p=p,parametric_paired_t_p=tp,t_model_mean_CI95_low=obs-half,t_model_mean_CI95_high=obs+half,status='tested_4_patients' if n==4 else 'insufficient_paired_patients'))
 tests=pd.DataFrame(tests);tests['family_BH_q_8']=np.nan
 tests['parametric_family_BH_q_8']=np.nan
 for key,ix in tests.groupby(['sensitivity','comparison','method']).groups.items():
  tests.loc[ix,'family_BH_q_8']=false_discovery_control(tests.loc[ix,'exact_signflip_p'].fillna(1).to_numpy());tests.loc[ix,'parametric_family_BH_q_8']=false_discovery_control(tests.loc[ix,'parametric_paired_t_p'].fillna(1).to_numpy())
 tests['joint_BH_q_96']=false_discovery_control(tests.exact_signflip_p.fillna(1).to_numpy());tests.to_csv(ROOT/'patient_level_statistics.csv',index=False)
if __name__=='__main__':run()
