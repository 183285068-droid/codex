"""Section-specific conditional spatial association; no causal or cell attribution claim."""
from pathlib import Path
import json,math
import anndata as ad
import numpy as np,pandas as pd
from scipy.stats import rankdata,false_discovery_control
ROOT=Path(__file__).resolve().parent
REF=ROOT.parent/'GSE281759_spatial_HCC_PGAM5'
SAMPLES=json.loads((REF/'sample_scope.json').read_text())
SEED=28175942;B=4999
PROGRAMS=['Exhaustion_core','Exhaustion_expanded','Activation_early','Effector_cytotoxic']
COV=['source_log_UMI','source_pct_mt','source_Macrophage_core','source_Hepatocyte','source_Epithelial','source_T_lineage','neighbor_T_lineage','neighbor_log_UMI','neighbor_pct_mt']

def residualizer(d):
    C=np.column_stack([np.ones(len(d)),rankdata(d[COV].to_numpy(float),axis=0)])
    u,s,_=np.linalg.svd(C,full_matrices=False);rank=int(np.sum(s>max(C.shape)*np.finfo(float).eps*s[0]))
    return u[:,:rank],rank

def partial(x,y,q):
    x=rankdata(np.asarray(x),axis=-1);y=rankdata(np.asarray(y))
    xr=x-(x@q)@q.T;yr=y-q@(q.T@y)
    return (xr@yr)/np.sqrt(np.sum(xr*xr,axis=-1)*np.dot(yr,yr))

def strata(d):
    # Spatial tiles fixed by full-section coordinates, plus within-cohort UMI halves.
    return (d.spatial_tile.to_numpy(int)*2+pd.qcut(d.source_log_UMI,2,labels=False,duplicates='drop').to_numpy(int))

def two_sided(obs,null):
    return min(1.,2*min((1+np.sum(null>=obs-1e-12))/(len(null)+1),(1+np.sum(null<=obs+1e-12))/(len(null)+1)))

def evaluate(d,item,program,cohort,kind,rng):
    d=d.reset_index(drop=True);y=d['neighbor_'+program].to_numpy(float);x=d.PGAM5_log1p_CP10k.to_numpy(float)
    groups=strata(d) if len(d)>1 else np.zeros(len(d),int)
    ids=[np.flatnonzero(groups==g) for g in np.unique(groups)]
    high=d.PGAM5_mac_group.eq('PGAM5_high').to_numpy();q,rank=residualizer(d) if len(d)>1 else (np.ones((len(d),1)),1)
    result=dict(**item,program=program,cohort=cohort,contrast=kind,n_sources=len(d),n_high=int(high.sum()),n_other_detected=int(d.PGAM5_mac_group.eq('Other_detected').sum()),covariate_rank=rank,
                strata=len(ids),permutable_sources=sum(len(i) for i in ids if len(i)>1),statistic=np.nan,raw_spearman=np.nan,
                raw_mean_high=np.nan,raw_mean_other=np.nan,raw_mean_difference=np.nan,permutation_p=np.nan,permutations=B,status='not_testable')
    if len(d)<30 or len(d)-rank<15 or np.ptp(y)==0 or np.ptp(x)==0:return result
    if kind=='adjusted_high_minus_other' and (high.sum()<5 or (~high).sum()<10):return result
    yr=rankdata(y);yr=yr-q@(q.T@yr)
    if np.dot(yr,yr)<1e-10:return result
    result['raw_spearman']=float(np.corrcoef(rankdata(x),rankdata(y))[0,1])
    if kind=='adjusted_high_minus_other':
        # Rank-outcome residual difference in units of residual SD; NOT RNA fold change.
        sd=float(np.std(yr,ddof=1));obs=float((yr[high].mean()-yr[~high].mean())/sd)
        result.update(raw_mean_high=float(y[high].mean()),raw_mean_other=float(y[~high].mean()),raw_mean_difference=float(y[high].mean()-y[~high].mean()))
    else:obs=float(partial(x,y,q))
    if not np.isfinite(obs):return result
    null=[];example=[]
    for start in range(0,B,128):
        count=min(128,B-start);perms=np.tile(np.arange(len(d)),(count,1))
        for index in ids:perms[:,index]=rng.permuted(perms[:,index],axis=1)
        if start==0:example=perms[:3].copy()
        if kind=='adjusted_high_minus_other':
            h=high[perms];n=h.sum(1);v=((h@yr)/n-((~h)@yr)/(len(d)-n))/sd
        else:v=partial(x[perms],y,q)
        null.extend(v)
    null=np.asarray(null);assert np.isfinite(null).all()
    result.update(statistic=obs,permutation_p=two_sided(obs,null),status='tested',null_mean=float(null.mean()),null_q025=float(np.quantile(null,.025)),null_q975=float(np.quantile(null,.975)))
    name=item['sample']+'_'+cohort+'_'+program+'_'+kind
    np.savez_compressed(ROOT/(name+'_audit.npz'),source_barcodes=d.barcode.to_numpy(str),x=x,y=y,covariates=d[COV].to_numpy(float),groups=groups,high=high,
                        observed=obs,null=null,example_permutations=example)
    return result

def prepare(a,item):
    G=a.obsp['spatial_connectivities'];o=a.obs;t=o.T_RNA_enriched.to_numpy(bool);mac=o.Mac_RNA_enriched.to_numpy(bool)
    n=np.asarray(G@t.astype(int)).ravel()
    denom=np.where(n>0,n,np.nan)
    d=o.loc[mac,['PGAM5_raw','PGAM5_log1p_CP10k','PGAM5_mac_group','array_row','array_col','T_neighbor_count','T_exhaustion_neighbor_count','T_activation_neighbor_count']].copy()
    for k,v in item.items():d[k]=v
    for k,v in [('source_log_UMI',np.log1p(o.total_counts)),('source_pct_mt',o.pct_mt),('source_Macrophage_core',o.Macrophage_core),('source_Hepatocyte',o.Hepatocyte),('source_Epithelial',o.Epithelial),('source_T_lineage',o.T_lineage)]:d[k]=np.asarray(v)[mac]
    for name,values in [('T_lineage',o.T_lineage),('log_UMI',np.log1p(o.total_counts)),('pct_mt',o.pct_mt)]+[(p,o[p]) for p in PROGRAMS]:
        d['neighbor_'+name]=(np.asarray(G@(np.asarray(values,float)*t)).ravel()/denom)[mac]
    d['neighbor_Tex_fraction']=(o.T_exhaustion_neighbor_count.to_numpy()/denom)[mac]
    ar=o[['array_row','array_col']].to_numpy();tile=np.minimum(2,np.floor((ar-ar.min(0))/np.maximum(ar.max(0)-ar.min(0),1)*3).astype(int));d['spatial_tile']=(tile[:,0]*3+tile[:,1])[mac]
    d.index.name='barcode';d=d.reset_index();d.to_csv(ROOT/f'{item["sample"]}_source_neighbor_RNA.csv',index=False,float_format='%.17g')
    return d
if __name__=='__main__':
    rows=[];descriptive=[]
    for item in SAMPLES:
        a=ad.read_h5ad(REF/f'{item["sample"]}_spatial_marker_subset.h5ad');d=prepare(a,item)
        valid=d[d.T_neighbor_count>0].copy();positive=valid[valid.PGAM5_raw>0].copy()
        rng=np.random.default_rng(SEED+int(item['gsm'][3:]))
        for p in PROGRAMS:rows.append(evaluate(valid,item,p,'all_mac_including_zero','partial_Spearman',rng))
        rows.append(evaluate(positive,item,'Exhaustion_core','PGAM5_detected_only','partial_Spearman',rng))
        rows.append(evaluate(positive,item,'Exhaustion_core','PGAM5_detected_only','adjusted_high_minus_other',rng))
        for group in ['PGAM5_high','Other_detected','Undetected']:
            dd=valid[valid.PGAM5_mac_group.eq(group)]
            for p in PROGRAMS:
                descriptive.append(dict(**item,source_group=group,program=p,n_sources=len(dd),mean_neighbor_score=dd['neighbor_'+p].mean(),median_neighbor_score=dd['neighbor_'+p].median(),mean_neighbor_Tex_fraction=dd.neighbor_Tex_fraction.mean()))
        print('Analyzed',item['sample'],'Mac-with-T',len(valid),'PGAM5 detected',len(positive),flush=True)
    table=pd.DataFrame(rows);assert len(table)==96
    # Maintain all predeclared tests; untestable tests assigned p=1 for family correction only.
    table['BH_q_96_planned_tests']=false_discovery_control(table.permutation_p.fillna(1).to_numpy())
    table['primary_BH_q_16_sections']=np.nan
    primary=(table.cohort=='all_mac_including_zero')&(table.program=='Exhaustion_core')&(table.contrast=='partial_Spearman')
    table.loc[primary,'primary_BH_q_16_sections']=false_discovery_control(table.loc[primary,'permutation_p'].fillna(1).to_numpy())
    table.to_csv(ROOT/'association_statistics.csv',index=False)
    pd.DataFrame(descriptive).to_csv(ROOT/'PGAM5_group_neighbor_program_summary.csv',index=False)
    pd.DataFrame([dict(patient=pat,program=p,sections_tested=int(len(d)),median_section_partial_rho=d.statistic.median(),positive_sections=int(d.statistic.gt(0).sum()),negative_sections=int(d.statistic.lt(0).sum()),q_significant_sections=int(d.BH_q_96_planned_tests.lt(.05).sum())) for (pat,p),d in table[(table.contrast=='partial_Spearman')&(table.cohort=='all_mac_including_zero')&(table.status=='tested')].groupby(['patient','program'])]).to_csv(ROOT/'patient_descriptive_direction_summary.csv',index=False)
    config=dict(seed=SEED,permutations=B,covariates=COV,programs=PROGRAMS,source='prior QC and fixed RNA gates, unchanged',null='3x3 full-section array-coordinate tiles x within-analysis-source UMI median; shuffle PGAM5 RNA / high labels only within strata',p='two-sided equal-tail randomization p with +1 correction; no zero-centered null assumption',primary_FDR_family='16 main exhaustion-core section tests; untestable entry p=1 for BH only',FDR_family='96 planned tests (16 sections x 6); untestable entries p=1 for BH only',minimum_n=30,minimum_residual_df=15,minimum_high_other=[5,10],causal_claim=False)
    (ROOT/'analysis_config.json').write_text(json.dumps(config,indent=2)+'\n')
