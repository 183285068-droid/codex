"""SPP1 macrophage-associated mixed-spot neighborhoods; never cell or causal attribution."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import rankdata, false_discovery_control

ROOT=Path(__file__).resolve().parent
REF=ROOT.parent/'GSE281759_spatial_HCC_PGAM5'
SAMPLES=json.loads((ROOT/'sample_scope.json').read_text())
SEED=28175981; B=4999
PROGRAMS=['Exhaustion_core','Exhaustion_expanded','Activation_early','Effector_cytotoxic','Exhaustion_no_CTLA4_TIGIT']
NO_CHECKPOINT=['PDCD1','LAG3','TOX','CXCL13','HAVCR2','ENTPD1']
COV=['source_log_UMI','source_pct_mt','source_Macrophage_core','source_Hepatocyte','source_Epithelial','source_T_lineage','neighbor_T_lineage','neighbor_log_UMI','neighbor_pct_mt']
STATES={'All_T':'T_RNA_enriched','Exhaustion':'Exhaustion_core_supported','Activation':'Activation_early_supported'}
LR_NAMES=['SPP1--CD44','CD274--PDCD1','PVR--TIGIT','LGALS9--HAVCR2','TGFB1--TGFBR1_TGFBR2','IL10--IL10RA_IL10RB','CXCL9--CXCR3','CD80--CD28']

def dense(x): return np.asarray(x.toarray() if hasattr(x,'toarray') else x)

def augment(a):
    o=a.obs
    o['SPP1_raw']=dense(a[:,['SPP1']].layers['counts']).ravel().astype(int)
    o['SPP1_log1p_CP10k']=dense(a[:,['SPP1']].X).ravel()
    m=o.Mac_RNA_enriched.to_numpy(bool); pos=m&o.SPP1_raw.gt(0).to_numpy()
    cutoff=float(o.loc[pos,'SPP1_log1p_CP10k'].quantile(.75)) if pos.any() else np.nan
    o['SPP1_mac_group']='Outside_mac_RNA_gate'
    o.loc[m,'SPP1_mac_group']='Undetected';o.loc[pos,'SPP1_mac_group']='Other_detected'
    o.loc[pos&o.SPP1_log1p_CP10k.ge(cutoff).to_numpy(),'SPP1_mac_group']='SPP1_high'
    o['Mac_RNA_nonT']=m&~o.T_RNA_enriched.to_numpy(bool)
    o['Exhaustion_no_CTLA4_TIGIT']=dense(a[:,NO_CHECKPOINT].X).mean(1)
    ar=o[['array_row','array_col']].to_numpy()
    tile=np.minimum(2,np.floor((ar-ar.min(0))/np.maximum(ar.max(0)-ar.min(0),1)*3).astype(int))
    o['spatial_tile']=tile[:,0]*3+tile[:,1]
    a.uns['SPP1_high_threshold_log1p_CP10k']=cutoff
    return a

def source_table(a,item):
    o=a.obs;G=a.obsp['spatial_connectivities'];m=o.Mac_RNA_enriched.to_numpy(bool);t=o.T_RNA_enriched.to_numpy(bool)
    n=np.asarray(G@t.astype(int)).ravel();den=np.where(n>0,n,np.nan)
    keys=['SPP1_raw','SPP1_log1p_CP10k','SPP1_mac_group','Mac_RNA_nonT','array_row','array_col','spatial_tile','spatial_degree','T_neighbor_count','T_exhaustion_neighbor_count','T_activation_neighbor_count']
    d=o.loc[m,keys].copy()
    for k,v in item.items():d[k]=v
    for k,v in [('source_log_UMI',np.log1p(o.total_counts)),('source_pct_mt',o.pct_mt),('source_Macrophage_core',o.Macrophage_core),('source_Hepatocyte',o.Hepatocyte),('source_Epithelial',o.Epithelial),('source_T_lineage',o.T_lineage)]: d[k]=np.asarray(v)[m]
    for name,v in [('T_lineage',o.T_lineage),('log_UMI',np.log1p(o.total_counts)),('pct_mt',o.pct_mt)]+[(p,o[p]) for p in PROGRAMS]:
        d['neighbor_'+name]=(np.asarray(G@(np.asarray(v,float)*t)).ravel()/den)[m]
    for state,col in STATES.items():
        count=np.asarray(G@o[col].to_numpy(int)).ravel()[m]
        d[state+'_neighbor_count']=count;d[state+'_neighbor_present']=count>0
    d.index.name='barcode';return d.reset_index()

def strata(d):
    def split(v):
        bins=pd.qcut(v,2,labels=False,duplicates='drop').to_numpy(float)
        return np.nan_to_num(bins,nan=0).astype(int)
    return d.spatial_tile.to_numpy(int)*4+split(d.source_log_UMI)*2+split(d.source_Macrophage_core)

def equal_tail(obs,null):
    return min(1.,2*min((1+(null>=obs-1e-12).sum())/(len(null)+1),(1+(null<=obs+1e-12).sum())/(len(null)+1)))

def projector(d):
    C=np.column_stack([np.ones(len(d)),rankdata(d[COV].to_numpy(float),axis=0)])
    u,s,_=np.linalg.svd(C,full_matrices=False)
    r=int(np.sum(s>max(C.shape)*np.finfo(float).eps*s[0]));return u[:,:r],r

def partial(x,y,q):
    x=rankdata(x,axis=-1);y=rankdata(y)
    xr=x-(x@q)@q.T;yr=y-q@(q.T@y)
    return (xr@yr)/np.sqrt(np.sum(xr*xr,axis=-1)*np.dot(yr,yr))

def evaluate(d,item,family,endpoint,cohort,rng):
    d=d.reset_index(drop=True);g=strata(d) if len(d) else np.array([],int);ids=[np.flatnonzero(g==v) for v in np.unique(g)]
    result=dict(**item,family=family,endpoint=endpoint,cohort=cohort,n_sources=len(d),statistic=np.nan,permutation_p=np.nan,status='not_testable',reason='',permutations=B,strata=len(ids),permutable_sources=sum(len(ix) for ix in ids if len(ix)>1))
    if family=='adjacency':
        h=(d.SPP1_mac_group.eq('SPP1_high') if cohort=='high_vs_other_detected' else d.SPP1_raw.gt(0)).to_numpy()
        y=d[endpoint+'_neighbor_present'].to_numpy(float)
        n1=int(h.sum());n0=len(d)-n1
        result.update(n_positive=n1,n_comparator=n0,positive_neighbor_count=int(y[h].sum()),comparator_neighbor_count=int(y[~h].sum()),positive_neighbor_pct=100*y[h].mean() if n1 else np.nan,comparator_neighbor_pct=100*y[~h].mean() if n0 else np.nan)
        if n1 and n0:result['statistic']=100*(y[h].mean()-y[~h].mean())
        minimum_positive=5 if cohort=='high_vs_other_detected' else 10
        if len(d)<30 or n1<minimum_positive or n0<10:result['reason']='Insufficient sources or group sizes';return result
        if np.ptp(y)==0:result['reason']='Constant neighbor presence';return result
        # A conditional randomization needs a stratum where both labels occur.
        if not any(np.ptp(h[ix].astype(int))>0 for ix in ids):result['reason']='No exchangeable labels within strata';return result
        x=h.astype(float);obs=result['statistic']
    else:
        x=d.SPP1_log1p_CP10k.to_numpy(float);y=d['neighbor_'+endpoint].to_numpy(float)
        result['n_positive']=int(d.SPP1_raw.gt(0).sum());result['n_comparator']=int(d.SPP1_raw.eq(0).sum())
        if len(d)<30 or np.ptp(x)==0 or np.ptp(y)==0:result['reason']='Insufficient sources or constant expression';return result
        q,rank=projector(d);result['covariate_rank']=rank
        if len(d)-rank<15:result['reason']='Insufficient residual degrees of freedom';return result
        obs=float(partial(x,y,q))
        if not np.isfinite(obs):result['reason']='Constant covariate residuals';return result
        if not any(np.ptp(x[ix])>0 for ix in ids):result['reason']='No exchangeable expression within strata';return result
        result['raw_spearman']=float(np.corrcoef(rankdata(x),rankdata(y))[0,1])
    null=[];examples=None
    for start in range(0,B,128):
        count=min(128,B-start);perms=np.tile(np.arange(len(d)),(count,1))
        for ix in ids:perms[:,ix]=rng.permuted(perms[:,ix],axis=1)
        if start==0:examples=perms[:3].copy()
        if family=='adjacency':
            h=x[perms];v=100*((h@y)/n1-((1-h)@y)/n0)
        else:v=partial(x[perms],y,q)
        null.extend(v)
    null=np.asarray(null);assert np.isfinite(null).all()
    # Save degenerate randomization as untestable, not false certainty.
    if np.ptp(null)<1e-10:result['reason']='Degenerate conditional null';return result
    result.update(statistic=float(obs),permutation_p=equal_tail(obs,null),status='tested',null_mean=float(null.mean()),null_q025=float(np.quantile(null,.025)),null_q975=float(np.quantile(null,.975)))
    np.savez_compressed(ROOT/f'{item["sample"]}_{family}_{endpoint}_{cohort}_audit.npz',source_barcodes=d.barcode.to_numpy(str),x=x,y=y,groups=g,covariates=d[COV].to_numpy(float) if family=='correlation' else np.empty((len(d),0)),observed=obs,null=null,example_permutations=examples)
    return result

def lr_edges(a,item,pairs):
    o=a.obs;coo=a.obsp['spatial_connectivities'].tocoo();keep=o.Mac_RNA_enriched.to_numpy(bool)[coo.row]&o.T_RNA_enriched.to_numpy(bool)[coo.col]
    s,t=coo.row[keep],coo.col[keep]
    d=pd.DataFrame({'source_barcode':a.obs_names[s],'target_barcode':a.obs_names[t],'SPP1_source_group':o.SPP1_mac_group.to_numpy()[s],'source_SPP1_raw':o.SPP1_raw.to_numpy()[s],'source_nonT':o.Mac_RNA_nonT.to_numpy()[s],'target_exhaustion_supported':o.Exhaustion_core_supported.to_numpy()[t],'target_activation_supported':o.Activation_early_supported.to_numpy()[t]})
    expression={g:np.expm1(dense(a[:,[g]].X).ravel().astype(float)) for g in a.var_names}
    for r in pairs.itertuples():d[r.pair]=np.sqrt(expression[r.ligand][s]*np.min(np.vstack([expression[g][t] for g in r.receptor_complex.split('_')]),axis=0))
    d.to_csv(ROOT/f'{item["sample"]}_mac_T_candidate_edges.csv.gz',index=False)
    rows=[]
    for group in ['All_mac','SPP1_detected','SPP1_high','Other_detected','Undetected']:
        sub=d if group=='All_mac' else d[d.source_SPP1_raw.gt(0)] if group=='SPP1_detected' else d[d.SPP1_source_group.eq(group)]
        for state,col in [('All_T',None),('Exhaustion','target_exhaustion_supported'),('Activation','target_activation_supported')]:
            e=sub if col is None else sub[sub[col]]
            for r in pairs.itertuples():
                pos=e[r.pair].gt(0)
                rows.append(dict(**item,source_group=group,target_state=state,pair=r.pair,eligible_edges=len(e),unique_sources=e.source_barcode.nunique(),unique_targets=e.target_barcode.nunique(),supported_edges=int(pos.sum()),supported_sources=e.loc[pos,'source_barcode'].nunique(),supported_targets=e.loc[pos,'target_barcode'].nunique(),edge_codetection_pct=100*pos.mean() if len(e) else np.nan,mean_geometric_CP10k=e[r.pair].mean() if len(e) else np.nan))
    return rows

def run():
    config=dict(seed=SEED,permutations=B,SPP1_positive='raw count >0 in prior Mac-RNA gate; NOT cell-specific expression',SPP1_high='within-section Q75 of positive Mac-RNA normalized SPP1, ties included',strata='3x3 full-section Visium coordinate tile x source UMI halves x Macrophage-core halves; bins recomputed within each analysis source set',adjacency='Difference in % source spots with >=1 directly adjacent target; no self-edges; sources with degree=0 excluded',correlation='Partial Spearman on all Mac-RNA sources with >=1 T-RNA neighbor, including SPP1 zeros; neighbor program mean over T-RNA neighbors only',covariates=COV,planned_tests=176,minimum_sources=30,minimum_group_sizes=[10,10],high_group_sizes=[5,10],minimum_residual_df=15,programs=PROGRAMS,no_checkpoint_program=NO_CHECKPOINT,primary_family='16 main All_T detection-vs-undetected adjacency tests',joint_FDR='176 planned tests; untestable p stays missing and contributes 1 only for BH',causal_claim=False)
    (ROOT/'analysis_config.json').write_text(json.dumps(config,indent=2)+'\n')
    pairs=pd.read_csv(ROOT/'selected_LR_pairs.csv');pairs['pair']=pairs.ligand+'--'+pairs.receptor_complex;pairs=pairs[pairs.pair.isin(LR_NAMES)];assert len(pairs)==8
    pairs.to_csv(ROOT/'tested_candidate_LR_pairs.csv',index=False)
    rows=[];summaries=[];lr=[];descriptive=[]
    for item in SAMPLES:
        a=augment(ad.read_h5ad(REF/f'{item["sample"]}_spatial_marker_subset.h5ad'));o=a.obs;d=source_table(a,item)
        d.to_csv(ROOT/f'{item["sample"]}_source_neighbor_RNA.csv',index=False,float_format='%.17g')
        o.to_csv(ROOT/f'{item["sample"]}_spot_annotations.csv.gz',float_format='%.17g')
        base=d[d.spatial_degree>0];strict=base[base.Mac_RNA_nonT];positive=base[base.SPP1_raw>0];valid=d[d.T_neighbor_count>0]
        rng=np.random.default_rng(SEED+int(item['gsm'][3:]))
        for state in STATES:rows.append(evaluate(base,item,'adjacency',state,'detected_vs_undetected',rng))
        rows.append(evaluate(strict,item,'adjacency','Exhaustion','nonT_detected_vs_undetected',rng))
        rows.append(evaluate(positive,item,'adjacency','Exhaustion','high_vs_other_detected',rng))
        for p in PROGRAMS:rows.append(evaluate(valid,item,'correlation',p,'all_mac_including_zero',rng))
        rows.append(evaluate(valid[valid.Mac_RNA_nonT],item,'correlation','Exhaustion_core','nonT_mac_including_zero',rng))
        for group,mask in [('SPP1_detected',d.SPP1_raw.gt(0)),('SPP1_high',d.SPP1_mac_group.eq('SPP1_high')),('Undetected',d.SPP1_raw.eq(0))]:
            dd=d[mask&d.spatial_degree.gt(0)]
            for state in STATES:descriptive.append(dict(**item,source_group=group,target_state=state,n_sources=len(dd),with_neighbor=int(dd[state+'_neighbor_present'].sum()),neighbor_pct=100*dd[state+'_neighbor_present'].mean() if len(dd) else np.nan))
        summaries.append(dict(**item,qc_spots=len(a),SPP1_all_spots=int(o.SPP1_raw.gt(0).sum()),mac_RNA_spots=len(d),SPP1_mac_detected=int(d.SPP1_raw.gt(0).sum()),SPP1_mac_undetected=int(d.SPP1_raw.eq(0).sum()),SPP1_mac_high=int(d.SPP1_mac_group.eq('SPP1_high').sum()),SPP1_high_threshold=a.uns['SPP1_high_threshold_log1p_CP10k'],nonT_mac_RNA_spots=int(d.Mac_RNA_nonT.sum()),T_RNA_spots=int(o.T_RNA_enriched.sum()),Tex_RNA_spots=int(o.Exhaustion_core_supported.sum()),Tact_RNA_spots=int(o.Activation_early_supported.sum()),SPP1_detected_with_T_neighbor=int((d.SPP1_raw.gt(0)&d.All_T_neighbor_present).sum()),SPP1_detected_with_Tex_neighbor=int((d.SPP1_raw.gt(0)&d.Exhaustion_neighbor_present).sum()),SPP1_detected_with_Tact_neighbor=int((d.SPP1_raw.gt(0)&d.Activation_neighbor_present).sum())))
        lr.extend(lr_edges(a,item,pairs))
        print('Analyzed',item['sample'],'Mac',len(d),'SPP1+',d.SPP1_raw.gt(0).sum(),flush=True)
    table=pd.DataFrame(rows);assert len(table)==176
    table['BH_q_176_planned_tests']=false_discovery_control(table.permutation_p.fillna(1).to_numpy())
    primary=(table.family=='adjacency')&(table.endpoint=='All_T')&(table.cohort=='detected_vs_undetected')
    table['primary_BH_q_16_sections']=np.nan;table.loc[primary,'primary_BH_q_16_sections']=false_discovery_control(table.loc[primary,'permutation_p'].fillna(1).to_numpy())
    primonly=primary&table.tumor_status.eq('primary');table['primary_only_BH_q_15_sections']=np.nan;table.loc[primonly,'primary_only_BH_q_15_sections']=false_discovery_control(table.loc[primonly,'permutation_p'].fillna(1).to_numpy())
    table.to_csv(ROOT/'interaction_statistics.csv',index=False)
    summary=pd.DataFrame(summaries);summary.to_csv(ROOT/'section_summary.csv',index=False)
    counts=[k for k in summary.columns if k not in SAMPLES[0] and k!='SPP1_high_threshold']
    pt=summary.groupby('patient')[counts].sum().reset_index();pt['sections']=summary.groupby('patient').size().reindex(pt.patient).to_numpy();pt.to_csv(ROOT/'patient_counts.csv',index=False)
    pd.DataFrame(lr).to_csv(ROOT/'candidate_LR_codetection_summary.csv',index=False)
    pd.DataFrame(descriptive).to_csv(ROOT/'source_group_neighbor_summary.csv',index=False)
    patient=[]
    for scope,ss in [('all_sections',table),('primary_only',table[table.tumor_status.eq('primary')])]:
        for (pat,fam,ep,coh),dd in ss.groupby(['patient','family','endpoint','cohort']):
            test=dd[dd.status.eq('tested')]
            patient.append(dict(scope=scope,patient=pat,family=fam,endpoint=ep,cohort=coh,sections_planned=len(dd),sections_tested=len(test),median_tested_section_statistic=test.statistic.median(),positive_tested_sections=int(test.statistic.gt(0).sum()),negative_tested_sections=int(test.statistic.lt(0).sum()),joint_q_significant_sections=int(test.BH_q_176_planned_tests.lt(.05).sum())))
    pd.DataFrame(patient).to_csv(ROOT/'patient_descriptive_effect_summary.csv',index=False)

if __name__=='__main__':run()
