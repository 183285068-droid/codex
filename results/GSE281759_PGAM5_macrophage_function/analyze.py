"""Targeted RNA-program contrasts on mixed Visium macrophage-associated spots."""
from pathlib import Path
import json,hashlib,itertools
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.stats import false_discovery_control
ROOT=Path(__file__).resolve().parent;REF=ROOT.parent/'GSE281759_spatial_HCC_PGAM5'
SAMPLES=json.loads((REF/'sample_scope.json').read_text());B=4999;SEED=281759103
PANELS={
 'MHCII_presentation':['HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','CD74'],
 'Inflammatory_myeloid':['S100A8','S100A9','FCN1','VCAN','IL1B'],
 'CCL_CXCL8_chemokines':['CCL2','CCL3','CCL4','CCL8','CXCL8'],
 'CXCR3_ligands':['CXCL9','CXCL10','CXCL11'],
 'Checkpoint_ligands':['CD274','PDCD1LG2','PVR','NECTIN2'],
 'Regulatory_mediators':['IL10','TGFB1','LGALS9'],
 'Costimulation_ligands':['CD80','CD86','ICOSLG','TNFSF9','TNFSF4'],
 'Scavenger_LPL_state':['LPL','CD163','STAB1','MARCO'],
 'C1Q_complement':['C1QA','C1QB','C1QC'],
 'Uptake_clearance_proxy':['MERTK','MSR1','CD68','CTSS'],
 'SPP1_MMP12_remodeling':['SPP1','MMP12'],
 'Hypoxia_response':['CA9','VEGFA','SLC2A1','LDHA','BNIP3','ADM'],
 'Cycling':['MKI67','TOP2A','UBE2C','CENPF','BIRC5','PCNA']}
PROGRAMS=list(PANELS);GENES=sorted(set(sum(PANELS.values(),[])))
IDENTITY={'Mac_core':['C1QA','C1QB','C1QC','CD68','CSF1R','MERTK','MSR1'],
 'Pan_myeloid':['LST1','TYROBP','FCER1G','AIF1','LYZ','CTSS']}
MIX={'Hepatocyte':['ALB','APOA1','APOA2','TTR','ASGR1','CYP2E1'],
 'Epithelial':['EPCAM','KRT8','KRT18','KRT19'],
 'T':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2'],
 'B':['MS4A1','CD79A','CD79B','CD19','CD22'],
 'Plasma':['JCHAIN','MZB1','IGKC','SDC1'],
 'Endothelial':['PECAM1','VWF','CDH5','KDR','EMCN'],
 'Fibroblast':['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA']}
COMPS=['high_vs_all_other_mac','high_vs_other_detected','nonT_high_vs_all_other_mac','mt5_high_vs_all_other_mac']
def dense(x):return np.asarray(x.toarray() if sparse.issparse(x) else x)
def scores(a,genes):return dense(a[:,genes].X).mean(1)
def tiles(o):
 xy=o[['array_col','array_row']].to_numpy();z=np.minimum(((xy-xy.min(0))/(np.ptp(xy,axis=0)+1)*3).astype(int),2)
 return z[:,0]+3*z[:,1]
def groups(o,ix):
 tile=tiles(o)[ix];parts=[]
 for t in np.unique(tile):
  z=np.flatnonzero(tile==t);v=o.total_counts.to_numpy()[ix[z]];d=v>np.median(v)
  for level in [False,True]:
   zz=z[d==level]
   if not len(zz):continue
   v=o.Macrophage_core.to_numpy()[ix[zz]];m=v>np.median(v)
   for state in [False,True]:
    part=zz[m==state]
    if len(part):parts.append(part)
 return parts
def design(a,ix,excluded):
 o=a.obs;cols={'intercept':np.ones(len(ix)),'log_UMI':np.log1p(o.total_counts.to_numpy()[ix]),'pct_MT':o.pct_mt.to_numpy()[ix]}
 for k,gs in {**IDENTITY,**MIX}.items():
  use=[g for g in gs if g not in excluded]
  if use:cols[k]=scores(a,use)[ix]
 tile=tiles(o)[ix]
 for t in np.unique(tile)[1:]:cols['tile_'+str(t)]=(tile==t).astype(float)
 names=list(cols);C=np.column_stack(list(cols.values()))
 for j in range(1,C.shape[1]):
  sd=C[:,j].std();C[:,j]=(C[:,j]-C[:,j].mean())/sd if sd>1e-12 else 0
 u,s,v=np.linalg.svd(C,full_matrices=False);rank=int((s>s[0]*1e-10).sum());Q=u[:,:rank]
 return C,Q,rank,names
def fit(a,ix,h,y,excluded):
 C,Q,rank,names=design(a,ix,excluded);hr=h-Q@(Q.T@h);yr=y-Q@(Q.T@y);den=float(hr@hr)
 return C,Q,rank,names,hr,yr,den,float(hr@yr/den) if den>1e-8 else np.nan
def run():
 config=dict(seed=SEED,permutations=B,programs=PANELS,comparisons=COMPS,primary='high_vs_all_other_mac',primary_tests=208,joint_tests=832,primary_only_tests=195,program_score='arithmetic mean log1p(full-library CP10k), zeros retained',covariates={**IDENTITY,**MIX},other_covariates=['log1p(total UMI)','MT percent','3x3 spatial tile indicators'],overlap='remove outcome genes from identity/mixing covariate programs, retain original source gates',null='permute high labels within 3x3 tile x UMI halves x macrophage-core halves; refit FWL adjusted coefficient',minimums={'high':5,'comparator':10,'residual_df':15},sensitivity='fixed original high/source gates; restrict T-negative or author MT<=5% QC without rethresholding',patient='equal-section median standardized coefficient; 4-patient exact sign-flip exploratory check',limitations=['mixed Visium spots, no resolved macrophage identity','179 selected markers, no genome-wide DEG/GO/KEGG/GSEA','no mitochondrial/mitophagy/OXPHOS complete gene sets','proxy scores are not functional assays'])
 (ROOT/'analysis_config.json').write_text(json.dumps(config,indent=2)+'\n')
 stats=[];des=[];gene_rows=[];summ=[];manifest=[];QC=[]
 for item in SAMPLES:
  sample=item['sample'];path=REF/f'{sample}_spatial_marker_subset.h5ad';a=ad.read_h5ad(path);o=a.obs;raw=dense(a.layers['counts']);norm=dense(a.X);assert all(g in a.var_names for g in GENES)
  manifest.append(dict(path=str(path.relative_to(ROOT.parent)),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),bytes=path.stat().st_size))
  mac=o.Mac_RNA_enriched.to_numpy(bool);high=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy();det=o.PGAM5_raw.gt(0).to_numpy()
  yall=np.column_stack([scores(a,PANELS[k]) for k in PROGRAMS]);out=pd.DataFrame(yall,index=a.obs_names,columns=PROGRAMS)
  for col in ['PGAM5_mac_group','Mac_RNA_enriched','total_counts','pct_mt','T_RNA_enriched','pass_qc_mt5_author','PGAM5_raw','PGAM5_log1p_CP10k']:out[col]=o[col]
  out.to_csv(ROOT/f'{sample}_spot_program_scores.csv.gz',float_format='%.17g')
  summ.append(dict(**item,qc_spots=len(a),mac_spots=int(mac.sum()),high=int(high.sum()),other_detected=int((mac&det&~high).sum()),undetected=int((mac&~det).sum()),mt5_high=int((high&o.pass_qc_mt5_author.to_numpy(bool)).sum())))
  for group,mask in [('PGAM5_high',high),('Other_detected',mac&det&~high),('Undetected',mac&~det)]:
   for j,k in enumerate(PROGRAMS):des.append(dict(**item,source_group=group,program=k,n_spots=int(mask.sum()),mean_score=yall[mask,j].mean() if mask.any() else np.nan,median_score=np.median(yall[mask,j]) if mask.any() else np.nan))
   for name,genes in MIX.items():QC.append(dict(**item,source_group=group,variable=name,mean=float(scores(a,genes)[mask].mean()) if mask.any() else np.nan))
   for name in ['total_counts','pct_mt','Macrophage_core','Pan_myeloid']:QC.append(dict(**item,source_group=group,variable=name,mean=float(o[name].to_numpy()[mask].mean()) if mask.any() else np.nan))
  rng=np.random.default_rng(SEED+int(item['gsm'][3:]))
  for comp in COMPS:
   eligible=mac.copy()
   if comp=='high_vs_other_detected':eligible&=det
   if comp.startswith('nonT'):eligible&=~o.T_RNA_enriched.to_numpy(bool)
   if comp.startswith('mt5'):eligible&=o.pass_qc_mt5_author.to_numpy(bool)
   ix=np.flatnonzero(eligible);h=high[ix].astype(float);n1=int(h.sum());n0=len(h)-n1;valid=n1>=5 and n0>=10
   base=dict(**item,comparison=comp,n_high=n1,n_comparator=n0)
   if not valid:
    stats.extend(dict(**base,program=k,status='not_testable',reason='high<5 or comparator<10',p=np.nan) for k in PROGRAMS);continue
   strata=groups(o,ix);H=np.tile(h,(B,1));code=np.zeros(len(ix),int)
   for c,g in enumerate(strata):H[:,g]=rng.permuted(H[:,g],axis=1);code[g]=c
   nulllist=[];betas=[];ranks=[];ys=[]
   for j,k in enumerate(PROGRAMS):
    y=yall[ix,j];C,Q,rank,names,hr,yr,den,beta=fit(a,ix,h,y,set(PANELS[k]));sd=y.std(ddof=1);df=len(h)-rank-1
    resH=H-(H@Q)@Q.T;denoms=(resH*resH).sum(1);null=np.divide(resH@yr,denoms,out=np.full(B,np.nan),where=denoms>1e-8)
    ok=df>=15 and den>1e-8 and sd>1e-10 and np.isfinite(null).all() and np.ptp(null)>1e-10
    p=min(1.,2*min((1+(null>=beta-1e-12).sum())/(B+1),(1+(null<=beta+1e-12).sum())/(B+1))) if ok else np.nan
    stats.append(dict(**base,program=k,status='tested' if ok else 'not_testable',reason='' if ok else 'insufficient residual df / degenerate adjusted null',high_mean=float(y[h>0].mean()),comparator_mean=float(y[h==0].mean()),unadjusted_delta=float(y[h>0].mean()-y[h==0].mean()),adjusted_beta=beta,standardized_beta=beta/sd if sd>1e-10 else np.nan,df=df,covariate_rank=rank,high_residual_ss=den,p=p,null_mean=float(np.nanmean(null)),null_q025=float(np.nanquantile(null,.025)),null_q975=float(np.nanquantile(null,.975))))
    nulllist.append(null);betas.append(beta);ranks.append(rank);ys.append(y)
   np.savez_compressed(ROOT/f'{sample}_{comp}_audit.npz',barcodes=a.obs_names[ix].to_numpy(str),high=h,outcomes=np.column_stack(ys),assignments_first20=H[:20].astype(np.uint8),groups=code,observed=betas,covariate_rank=ranks,null=np.column_stack(nulllist))
   if comp in ['high_vs_all_other_mac','high_vs_other_detected']:
    for gene in GENES:
     z=norm[ix,a.var_names.get_loc(gene)];C,Q,rank,names,hr,yr,den,beta=fit(a,ix,h,z,{gene});li=a.var_names.get_loc(gene);counts=raw[ix,li]
     pb1=1e4*counts[h>0].sum()/o.total_counts.to_numpy()[ix[h>0]].sum();pb0=1e4*counts[h==0].sum()/o.total_counts.to_numpy()[ix[h==0]].sum()
     gene_rows.append(dict(**base,gene=gene,high_mean=float(z[h>0].mean()),comparator_mean=float(z[h==0].mean()),high_detect_pct=float(100*(counts[h>0]>0).mean()),comparator_detect_pct=float(100*(counts[h==0]>0).mean()),adjusted_beta=beta,standardized_beta=beta/z.std(ddof=1) if z.std(ddof=1)>1e-10 else np.nan,pseudobulk_CP10k_high=pb1,pseudobulk_CP10k_other=pb0,descriptive_log2_ratio_offset_0_1=np.log2((pb1+.1)/(pb0+.1))))
  print(sample,'done',flush=True)
 stat=pd.DataFrame(stats);stat['joint_q_832']=false_discovery_control(stat.p.fillna(1).to_numpy());stat['primary_q_208']=np.nan;ix=stat.comparison.eq(COMPS[0]);stat.loc[ix,'primary_q_208']=false_discovery_control(stat.loc[ix,'p'].fillna(1).to_numpy());ix2=ix&stat.tumor_status.eq('primary');stat['primary_only_q_195']=np.nan;stat.loc[ix2,'primary_only_q_195']=false_discovery_control(stat.loc[ix2,'p'].fillna(1).to_numpy());stat.to_csv(ROOT/'program_statistics.csv',index=False)
 pd.DataFrame(des).to_csv(ROOT/'program_group_means.csv',index=False);pd.DataFrame(gene_rows).to_csv(ROOT/'marker_gene_descriptive_contrasts.csv',index=False);pd.DataFrame(summ).to_csv(ROOT/'section_summary.csv',index=False);pd.DataFrame(QC).to_csv(ROOT/'group_quality_mixing.csv',index=False)
 patient=stat[stat.status.eq('tested')].groupby(['patient','comparison','program'],as_index=False).agg(median_standardized_beta=('standardized_beta','median'),median_adjusted_beta=('adjusted_beta','median'),n_sections=('sample','size'));patient.to_csv(ROOT/'patient_equal_section_effects.csv',index=False)
 meta=[]
 for k in PROGRAMS:
  pp=patient[patient.comparison.eq(COMPS[0])&patient.program.eq(k)];v=pp.median_standardized_beta.to_numpy();n=len(v);obs=v.mean() if n else np.nan;null=np.array(list(itertools.product([-1,1],repeat=n)))@v/n if n else np.array([]);p=(np.abs(null)>=abs(obs)-1e-12).mean() if n==4 else np.nan
  meta.append(dict(program=k,n_patients=n,mean_patient_effect=obs,median_patient_effect=np.median(v) if n else np.nan,positive_patients=int((v>0).sum()),negative_patients=int((v<0).sum()),patient_signflip_p=p))
 meta=pd.DataFrame(meta);meta['patient_BH_q_13']=false_discovery_control(meta.patient_signflip_p.fillna(1).to_numpy());meta.to_csv(ROOT/'patient_direction_summary.csv',index=False)
 (ROOT/'reference_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
if __name__=='__main__':run()
