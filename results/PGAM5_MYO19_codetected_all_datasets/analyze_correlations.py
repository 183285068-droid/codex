"""Per-cohort correlations restricted to the two-gene co-detected cells."""
from pathlib import Path
import itertools
import json
import platform
import importlib.metadata as md
import numpy as np
import pandas as pd
from scipy.stats import spearmanr,pearsonr,rankdata
from statsmodels.stats.multitest import multipletests

R=Path(__file__).resolve().parent
B=200000
specs=json.loads((R/'cohorts.json').read_text())
d=pd.read_csv(R/'codetected_macrophage_expression.csv.gz')
summary=pd.read_csv(R/'cohort_detection_summary.csv')
assert d[['PGAM5_raw_count','MYO19_raw_count']].gt(0).all().all()
protocol={'scope':'All8completed single-cell datasets; GSE140228 split by platform gives9strata. Preserve prior tissue/histology scope and latest frozen macrophage identities. GSE154906 cancelled','selection':'Strictly PGAM5 raw symbol-summed count>0 AND MYO19 raw symbol-summed count>0; no zero/undetected cells in correlations or scatter plots','primary':'Spearman on natural log1p(full-library CP10k), separate cohort/platform. Minimum3co-detected cells and both vectors nonconstant','primary_p':'Two-sided absolute-statistic permutation p. n<=8: exact enumeration of all n! labeled MYO19 permutations, including identity, extreme/total. n>8: fixed200000MonteCarlo shuffles, (extreme+1)/(B+1), fixed cohort/method/scale seeds','multiple_testing':'Four separate families: Spearman/log1p_CP10k(primary),Spearman/raw_counts,Pearson/log1p_CP10k,Pearson/raw_counts(sensitivity). Each BH-adjusts permutation p over estimable strata only; family size reported. Unestimable n<3 or constant vector =>NA, never p=0','analytical_p':'Also export standard scipy two-sided Spearman/Pearson asymptotic p for inspection; BH families use permutation p, not analytical p','sensitivity':'Raw counts preserve original UMI/read distinctions. Normalized correlation can arise through common library denominator even if raw counts are constant; raw non-estimability must not be hidden','dependence':'Cell exchangeability assumed by pooled-cell tests; source patients/libraries are not independent cell replicates. No cross-cohort pooled p, no treating overlapping125449/151530 or140228platforms as independent validation','adjustments':'No depth matching/regression, partial correlation, new annotation, donor pairing or gene-specific expression threshold beyond raw>0'}
(R/'protocol.json').write_text(json.dumps(protocol,indent=2))

def permutation_test(xx,yy,seed):
    x=xx-xx.mean();y=yy-yy.mean()
    den=np.linalg.norm(x)*np.linalg.norm(y)
    observed=float(np.dot(x,y)/den)
    extreme=0;n=len(x)
    if n<=8:
        total=0
        for perm in itertools.permutations(y):
            value=np.dot(x,perm)/den
            extreme+=abs(value)>=abs(observed)-1e-12
            total+=1
        return observed,float(extreme/total),int(extreme),int(total),'exact_labeled_permutation'
    rng=np.random.default_rng(seed)
    for k in range(B):
        value=np.dot(x,rng.permutation(y))/den
        extreme+=abs(value)>=abs(observed)-1e-12
    return observed,float((extreme+1)/(B+1)),int(extreme),B,'MonteCarlo_permutation'

rows=[]
for i,spec in enumerate(specs):
    t=d[d.cohort.eq(spec['cohort'])]
    assert len(t)==spec['expected_codetected']
    for j,scale in enumerate(['log1p_CP10k','raw_counts']):
        suffix='_log1p_CP10k' if scale=='log1p_CP10k' else '_raw_count'
        x=t['PGAM5'+suffix].to_numpy(dtype=float);y=t['MYO19'+suffix].to_numpy(dtype=float)
        for k,method in enumerate(['Spearman','Pearson']):
            row={'cohort':spec['cohort'],'dataset':spec['dataset'],'platform':spec['platform'],'count_type':spec['count_type'],'scope':spec['scope'],'n_codetected':len(t),'method':method,'expression_scale':scale,'family':method+'/'+scale,'role':'primary' if method=='Spearman' and scale=='log1p_CP10k' else 'sensitivity','PGAM5_unique_values':len(np.unique(x)),'MYO19_unique_values':len(np.unique(y)),'correlation':None,'analytical_p':None,'permutation_p':None,'permutation_type':None,'permutations':None,'extreme_permutations':None,'permutation_seed':None,'status':'estimable'}
            if len(t)<3:row['status']='not_estimable_n<3'
            elif np.ptp(x)==0 or np.ptp(y)==0:row['status']='not_estimable_constant_vector'
            else:
                stat=spearmanr(x,y) if method=='Spearman' else pearsonr(x,y)
                xx=rankdata(x) if method=='Spearman' else x
                yy=rankdata(y) if method=='Spearman' else y
                seed=20261010+i*100+j*10+k
                value,p,extreme,total,kind=permutation_test(xx,yy,seed)
                assert np.isclose(value,stat.statistic,atol=1e-12)
                row.update(correlation=float(stat.statistic),analytical_p=float(stat.pvalue),permutation_p=p,permutation_type=kind,permutations=total,extreme_permutations=extreme,permutation_seed=seed if kind.startswith('MonteCarlo') else None)
            rows.append(row)
    print(spec['cohort'],'selected',len(t),flush=True)
statistics=pd.DataFrame(rows)
statistics['BH_q_permutation']=np.nan
statistics['BH_family_estimable_tests']=0
for family,t in statistics.groupby('family',sort=False):
    valid=t.permutation_p.notna()
    ix=t.index[valid]
    statistics.loc[t.index,'BH_family_estimable_tests']=len(ix)
    if len(ix):statistics.loc[ix,'BH_q_permutation']=multipletests(statistics.loc[ix,'permutation_p'],method='fdr_bh')[1]
statistics.to_csv(R/'correlation_statistics_full.csv',index=False)
d.groupby(['cohort','PGAM5_raw_count','MYO19_raw_count']).size().rename('cells').reset_index().to_csv(R/'raw_count_patterns.csv',index=False)
wide=summary.copy()
for method in ['Spearman','Pearson']:
    for scale in ['log1p_CP10k','raw_counts']:
        t=statistics[(statistics.method==method)&(statistics.expression_scale==scale)].set_index('cohort')
        prefix=method+'_'+('log1p' if scale=='log1p_CP10k' else 'raw')
        for field,label in [('correlation','r'),('permutation_p','p'),('BH_q_permutation','q'),('status','status')]:wide[prefix+'_'+label]=wide.cohort.map(t[field])
wide.to_csv(R/'correlation_summary.csv',index=False)
primary=statistics[(statistics.method=='Spearman')&statistics.expression_scale.eq('log1p_CP10k')]
result={'datasets':8,'platform_strata':len(specs),'selected_co_detected_rows':len(d),'primary_estimable_strata':int(primary.permutation_p.notna().sum()),'primary_BH_significant_strata':primary.loc[primary.BH_q_permutation.lt(.05),'cohort'].tolist(),'unestimable_primary_strata':primary.loc[primary.status.ne('estimable'),['cohort','n_codetected','status']].to_dict('records'),'selection_verified_all_both_raw_positive':True,'no_cross_cohort_pooled_test':True,'primary_unit':'Co-detected cell; labels/overlap do not establish independent patient validation','interpretation':'Selected-subset expression correlations, not all-macrophage or causal coexpression. Shared CP10k denominator and low raw-count variance can dominate.'}
(R/'result.json').write_text(json.dumps(result,indent=2))
versions={'python':platform.python_version(),**{name:md.version(name) for name in ['numpy','pandas','scipy','statsmodels','matplotlib','anndata','h5py']}}
(R/'environment_versions.json').write_text(json.dumps(versions,indent=2));(R/'requirements.txt').write_text('\n'.join(k+'=='+v for k,v in versions.items() if k!='python')+'\n')
from plot_correlations import render
render(R)
print(wide[['cohort','both_detected','Spearman_log1p_r','Spearman_log1p_p','Spearman_log1p_q','Spearman_raw_r','Spearman_raw_p','Spearman_raw_q']].to_string(index=False),flush=True)
print(json.dumps(result,indent=2),flush=True)
