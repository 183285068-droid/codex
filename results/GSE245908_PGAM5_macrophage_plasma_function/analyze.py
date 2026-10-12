"""Plasma-related mixed-spot functional profiles; source-label graph permutations."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,anndata as ad,h5py
from scipy import sparse
from scipy.stats import false_discovery_control
ROOT=Path(__file__).resolve().parent;REF=ROOT.parent/'GSE245908_PGAM5_immune_subgroup_neighbors';RAW=ROOT.parent/'GSE245908_spatial_PGAM5'
SAMPLES={'CHC20':'GSM7850822','CHC23':'GSM7850823'};B=4999;SEED=245908105
PANELS={
 'Plasma_differentiation':['PRDM1','XBP1','IRF4','TNFRSF17'],
 'Secretory_ER':['HSPA5','HSP90B1','SEC11C','SSR1','SSR2','SEC61A1','FKBP11','P4HB'],
 'UPR_stress':['ERN1','ATF4','ATF6','DDIT3','DNAJB9','HERPUD1'],
 'Antibody_heavy_chain':['IGHG1','IGHG2','IGHG3','IGHG4','IGHA1','IGHA2','IGHM','IGHD','IGHE'],
 'MHCII_presentation':['HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','CD74'],
 'Immunoregulation':['IL10','TGFB1','CD274','PDCD1LG2'],
 'Survival':['MCL1','BCL2','BCL2L1','TNFRSF17'],
 'Trafficking_adhesion':['CXCR4','ITGA4','ITGB1','CD44'],
 'Cycling':['MKI67','TOP2A','UBE2C','CENPF','BIRC5','PCNA'],
 'Hypoxia':['CA9','VEGFA','SLC2A1','LDHA','BNIP3','ADM'],
 'IFN_response':['STAT1','ISG15','IFIT1','IFIT2','IFIT3','MX1','IFI6','OAS1']}
PROGRAMS=list(PANELS);COMPS=['high_only_vs_other_detected_only','high_any_vs_no_high','strict_high_only_vs_other_detected_only']
MIX={'Plasma_identity':['MZB1','JCHAIN','SDC1','IGKC','DERL3'],
 'Macrophage':['C1QA','C1QB','C1QC','CD68','CSF1R','MERTK','MSR1'],
 'T':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2'],
 'B':['MS4A1','CD79A','CD79B','CD19','CD22'],
 'Hepatocyte':['ALB','APOA1','APOA2','TTR','ASGR1','CYP2E1'],
 'Epithelial':['EPCAM','KRT8','KRT18','KRT19'],
 'Endothelial':['PECAM1','VWF','CDH5','KDR','EMCN'],
 'Fibroblast':['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA']}
ISOTYPES={'IgG':['IGHG1','IGHG2','IGHG3','IGHG4'],'IgA':['IGHA1','IGHA2'],'IgM':['IGHM'],'IgD':['IGHD'],'IgE':['IGHE']}
def dense(x):return np.asarray(x.toarray() if sparse.issparse(x) else x)
def score(a,genes):return dense(a[:,genes].X).mean(1)
def tiles(o):
 xy=o[['array_col','array_row']].to_numpy();v=np.minimum(((xy-xy.min(0))/(np.ptp(xy,axis=0)+1)*3).astype(int),2);return v[:,0]+3*v[:,1]
def load(sample,gsm):
 dst=ROOT/'inputs'/f'{sample}_functional_marker_inputs.h5ad'
 if dst.exists():return ad.read_h5ad(dst),None
 a=ad.read_h5ad(REF/f'{sample}_immune_marker_subset.h5ad');path=RAW/'inputs'/sample/f'{gsm}_{sample}_filtered_feature_bc_matrix.h5'
 with h5py.File(path) as f:
  m=f['matrix'];raw=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:]));names=[v.decode() for v in m['features/name'][:]];ids=[v.decode() for v in m['features/id'][:]];barcodes=[v.decode() for v in m['barcodes'][:]]
 genes=list(dict.fromkeys(list(a.var_names)+sum(PANELS.values(),[])+sum(MIX.values(),[])));assert all(names.count(g)==1 for g in genes)
 ci=np.array([barcodes.index(b) for b in a.obs_names]);gi=np.array([names.index(g) for g in genes]);counts=raw[gi][:,ci].T.tocsr();total=np.asarray(raw[:,ci].sum(0)).ravel();assert np.array_equal(total,a.obs.total_counts)
 norm=counts.astype(float).multiply(1e4/total[:,None]).tocsr();norm.data=np.log1p(norm.data)
 out=ad.AnnData(norm,obs=a.obs.copy(),var=pd.DataFrame({'gene_ids':[ids[i] for i in gi]},index=genes));out.layers['counts']=counts;out.obsp['spatial_connectivities']=a.obsp['spatial_connectivities'].copy();out.obsm['spatial']=a.obsm['spatial'].copy();out.uns['spatial']=a.uns['spatial'];out.write_h5ad(dst,compression='gzip')
 return out,dict(path=str(path.relative_to(ROOT.parent)),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),full_features=raw.shape[0])
def nuisance(a,ix,exclude):
 o=a.obs;cols=[np.ones(len(ix)),np.log1p(o.total_counts.to_numpy()[ix]),o.pct_mt.to_numpy()[ix]];names=['intercept','log_UMI','pct_MT']
 for k,gs in MIX.items():
  use=[g for g in gs if g not in exclude]
  if use:cols.append(score(a,use)[ix]);names.append(k)
 tile=tiles(o)[ix]
 for t in np.unique(tile)[1:]:cols.append((tile==t).astype(float));names.append('tile_'+str(t))
 C=np.column_stack(cols)
 for j in range(1,C.shape[1]):
  sd=C[:,j].std();C[:,j]=(C[:,j]-C[:,j].mean())/sd if sd>1e-12 else 0
 u,s,v=np.linalg.svd(C,full_matrices=False);rank=int((s>s[0]*1e-10).sum());return C,u[:,:rank],rank,names
def strata(o,ix):
 tile=tiles(o)[ix];out=[]
 for t in np.unique(tile):
  z=np.flatnonzero(tile==t);d=o.total_counts.to_numpy()[ix[z]]>np.median(o.total_counts.to_numpy()[ix[z]])
  for lv in [False,True]:
   zz=z[d==lv]
   if not len(zz):continue
   m=o.score_Macrophage_core.to_numpy()[ix[zz]]>np.median(o.score_Macrophage_core.to_numpy()[ix[zz]])
   for mv in [False,True]:
    g=zz[m==mv]
    if len(g):out.append(g)
 return out
def masks(nh,no,strict,comparison):
 if comparison=='high_any_vs_no_high':return nh>0,nh==0
 left=(nh>0)&(no==0);right=(nh==0)&(no>0)
 if comparison.startswith('strict'):left&=strict;right&=strict
 return left,right
def run():
 (ROOT/'inputs').mkdir(exist_ok=True)
 config=dict(seed=SEED,permutations=B,programs=PANELS,primary=COMPS[0],comparisons=COMPS,planned_primary_tests=22,planned_joint_tests=66,statistical_unit='plasma-related Visium target spot, not plasma cell',covariates=MIX,other_covariates=['log UMI','MT percent','3x3 spatial tile'],adjusted_effect='difference in mean nuisance-regression residual scores; regression fitted to all fixed plasma-related targets without group label',null='permute source high labels among PGAM5-detected Mac-RNA sources inside 3x3 tile x UMI halves x Mac-core halves, then recompute target adjacency groups on the SAME graph',minimum_groups=[5,10],minimum_valid_permutations=1000,invalid_permutations='group sizes below minimum => NaN; conditional inference given analyzable groups, valid denominator reported',strict_gate='original plasma gate plus MZB1,JCHAIN,DERL3 all raw count>0',causal_claim=False)
 (ROOT/'analysis_config.json').write_text(json.dumps(config,indent=2)+'\n');stats=[];summary=[];desc=[];genes=[];iso=[];sources=[];expanded=[]
 for sample,gsm in SAMPLES.items():
  a,source=load(sample,gsm)
  if source:sources.append(source)
  dst=ROOT/'inputs'/f'{sample}_functional_marker_inputs.h5ad';expanded.append(dict(path=str(dst.relative_to(ROOT)),bytes=dst.stat().st_size,sha256=hashlib.sha256(dst.read_bytes()).hexdigest()))
  o=a.obs;G=a.obsp['spatial_connectivities'];p=o.gate_Plasma_related.to_numpy(bool);pi=np.flatnonzero(p);source_mask=o.source_group.isin(['PGAM5_high','Other_detected']).to_numpy();si=np.flatnonzero(source_mask);h=o.source_group.eq('PGAM5_high').to_numpy()[si].astype(float);E=G[pi][:,si];nh=np.asarray(E@h).ravel();nt=np.asarray(E.sum(1)).ravel();no=nt-nh;strict=(dense(a[pi,['MZB1','JCHAIN','DERL3']].layers['counts'])>0).all(1)
  group=np.where(nh>0,np.where(no>0,'Shared_high_other','High_only'),np.where(no>0,'Other_detected_only','No_detected_mac_neighbor'));q=o.iloc[pi].copy();q['n_high_neighbors']=nh.astype(int);q['n_other_detected_neighbors']=no.astype(int);q['neighbor_group']=group;q['strict_plasma']=strict
  for k in PROGRAMS:q[k]=score(a,PANELS[k])[pi]
  q.to_csv(ROOT/f'{sample}_plasma_target_scores.csv.gz',float_format='%.17g')
  summary.append(dict(sample=sample,qc_spots=len(a),plasma_spots=len(pi),high_sources=int(h.sum()),other_detected_sources=int(len(h)-h.sum()),high_any=int((nh>0).sum()),high_only=int((group=='High_only').sum()),other_only=int((group=='Other_detected_only').sum()),shared=int((group=='Shared_high_other').sum()),remote=int((group=='No_detected_mac_neighbor').sum()),strict_plasma=int(strict.sum())))
  for gp in ['High_only','Shared_high_other','Other_detected_only','No_detected_mac_neighbor']:
   use=group==gp
   for k in PROGRAMS:
    z=q[k].to_numpy();desc.append(dict(sample=sample,group=gp,program=k,n_targets=int(use.sum()),mean=z[use].mean() if use.any() else np.nan,median=np.median(z[use]) if use.any() else np.nan))
   igs={k:dense(a[pi,gs].layers['counts']).sum(1) for k,gs in ISOTYPES.items()};den=sum(v[use].sum() for v in igs.values())
   for k,v in igs.items():iso.append(dict(sample=sample,group=gp,isotype=k,n_targets=int(use.sum()),raw_constant_counts=int(v[use].sum()),total_heavy_constant_counts=int(den),count_fraction=v[use].sum()/den if den else np.nan))
  rng=np.random.default_rng(SEED+int(gsm[3:]));H=np.tile(h,(B,1));code=np.zeros(len(si),int)
  for c,g in enumerate(strata(o,si)):H[:,g]=rng.permuted(H[:,g],axis=1);code[g]=c
  NH=np.asarray(E@H.T).T;NO=nt[None,:]-NH
  for comp in COMPS:
   left,right=masks(nh,no,strict,comp);L,R=masks(NH,NO,strict[None,:],comp);nl=L.sum(1);nr=R.sum(1);valid=(nl>=5)&(nr>=10);n1=int(left.sum());n0=int(right.sum());resids=[];nulls=[]
   for k,gs in PANELS.items():
    y=q[k].to_numpy();C,Q,rank,names=nuisance(a,pi,set(gs));res=y-Q@(Q.T@y);sd=y.std(ddof=1);obs=res[left].mean()-res[right].mean() if n1 and n0 else np.nan;null=np.divide(L@res,nl,out=np.full(B,np.nan),where=nl>0)-np.divide(R@res,nr,out=np.full(B,np.nan),where=nr>0);null[~valid]=np.nan;vv=null[np.isfinite(null)]
    ok=n1>=5 and n0>=10 and len(vv)>=1000 and len(pi)-rank>=15 and sd>1e-10 and np.ptp(vv)>1e-10
    pv=min(1.,2*min((1+(vv>=obs-1e-12).sum())/(len(vv)+1),(1+(vv<=obs+1e-12).sum())/(len(vv)+1))) if ok else np.nan
    stats.append(dict(sample=sample,comparison=comp,program=k,n_high_targets=n1,n_comparator_targets=n0,high_mean=y[left].mean() if n1 else np.nan,comparator_mean=y[right].mean() if n0 else np.nan,raw_delta=y[left].mean()-y[right].mean() if n1 and n0 else np.nan,adjusted_residual_delta=obs,standardized_delta=obs/sd if sd>1e-10 else np.nan,covariate_rank=rank,valid_permutations=len(vv),invalid_permutations=B-len(vv),null_mean=vv.mean() if len(vv) else np.nan,p=pv,status='tested' if ok else 'not_testable',reason='' if ok else 'insufficient groups/valid permutations or degenerate outcome'))
    resids.append(res);nulls.append(null)
   np.savez_compressed(ROOT/f'{sample}_{comp}_audit.npz',target_barcodes=a.obs_names[pi].to_numpy(str),source_barcodes=a.obs_names[si].to_numpy(str),source_high=h,source_strata=code,assignments_first20=H[:20].astype(np.uint8),target_high_neighbor_counts_first20=NH[:20],null_high_n=nl,null_comparator_n=nr,strict=strict,residuals=np.column_stack(resids),null=np.column_stack(nulls))
   if comp==COMPS[0]:
    for gene in sorted(set(sum(PANELS.values(),[])+['MZB1','JCHAIN','DERL3','SDC1','IGKC','PAX5','MS4A1'])):
     if gene not in a.var_names:continue
     vals=dense(a[pi,[gene]].X).ravel();rr=dense(a[pi,[gene]].layers['counts']).ravel();genes.append(dict(sample=sample,gene=gene,n_high=n1,n_other=n0,high_mean=vals[left].mean() if n1 else np.nan,other_mean=vals[right].mean() if n0 else np.nan,high_detect_pct=100*(rr[left]>0).mean() if n1 else np.nan,other_detect_pct=100*(rr[right]>0).mean() if n0 else np.nan))
  print(sample,summary[-1],flush=True)
 stat=pd.DataFrame(stats);stat['joint_q_66']=false_discovery_control(stat.p.fillna(1).to_numpy());ix=stat.comparison.eq(COMPS[0]);stat['primary_q_22']=np.nan;stat.loc[ix,'primary_q_22']=false_discovery_control(stat.loc[ix,'p'].fillna(1).to_numpy());stat.to_csv(ROOT/'program_statistics.csv',index=False);pd.DataFrame(summary).to_csv(ROOT/'target_group_summary.csv',index=False);pd.DataFrame(desc).to_csv(ROOT/'program_group_profiles.csv',index=False);pd.DataFrame(genes).to_csv(ROOT/'marker_gene_descriptive_profiles.csv',index=False);pd.DataFrame(iso).to_csv(ROOT/'immunoglobulin_isotype_counts.csv',index=False)
 (ROOT/'input_manifest.json').write_text(json.dumps(expanded,indent=2)+'\n')
 if sources:(ROOT/'raw_source_manifest.json').write_text(json.dumps(sources,indent=2)+'\n')
if __name__=='__main__':run()
