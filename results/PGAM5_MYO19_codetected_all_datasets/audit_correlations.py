"""Independent source CSR selection, rank/linear coefficients, p formulas and BH."""
from pathlib import Path
import argparse
import itertools
import json
import hashlib
import math
import h5py
import numpy as np
import pandas as pd
from scipy.stats import rankdata,t,beta

R=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--saved-results-only',action='store_true',help='Skip large source h5ad reads; check portable selected-cell results only')
args=parser.parse_args()
specs=json.loads((R/'cohorts.json').read_text())
d=pd.read_csv(R/'codetected_macrophage_expression.csv.gz')
summary=pd.read_csv(R/'cohort_detection_summary.csv').set_index('cohort')
stats=pd.read_csv(R/'correlation_statistics_full.csv')
checks=[]

def check(name,value):
    assert value,name
    checks.append(name)

def close(a,b):return np.allclose(a,b,atol=1e-12,rtol=1e-9,equal_nan=True)

def column(node):
    if isinstance(node,h5py.Group):
        if 'values' in node:
            values=column(node['values']);assert not node['mask'][:].any();return values
        values=column(node['categories']);codes=node['codes'][:];assert (codes>=0).all();return values[codes]
    if h5py.check_string_dtype(node.dtype) is not None:return node.asstr()[:]
    return node[:]

check('8datasets,9platform strata,224selected rows, unique within cohort',len(specs)==9 and len(set(z['dataset'] for z in specs))==8 and len(d)==224 and not d[['cohort','cell_id']].duplicated().any())
check('Every selected cell strictly positive for both raw genes',d[['PGAM5_raw_count','MYO19_raw_count']].gt(0).all().all())
for g in ['PGAM5','MYO19']:
    cp=d[g+'_raw_count']/d.total_counts*10000
    check(g+' CP10k/log1p reconstructed from full-library denominator',close(cp,d[g+'_CP10k']) and close(np.log1p(cp),d[g+'_log1p_CP10k']))
source_pairs=0
for spec in specs:
    cohort=spec['cohort'];z=d[d.cohort.eq(cohort)]
    check(cohort+' selected count and cohort summary agree',len(z)==spec['expected_codetected'] and int(summary.loc[cohort,'both_detected'])==len(z))
    if args.saved_results_only:continue
    with h5py.File(spec['path'],'r') as f:
        n=len(f['X/indptr'])-1
        mask=column(f['obs/harmonized_group'])=='Macrophage' if 'harmonized_group' in f['obs'] else np.ones(n,bool)
        check(cohort+' source macrophage denominator fixed',int(mask.sum())==spec['expected_macrophages'])
        genes=column(f['var/gene']) if 'gene' in f['var'] else column(f['var'][f['var'].attrs['_index']])
        assert f['X'].attrs['encoding-type']=='csr_matrix'
        pointer=f['X/indptr'][:];indices=f['X/indices'][:];values=f['X/data'][:]
        recovered=[]
        for g in ['PGAM5','MYO19']:
            ix=np.flatnonzero(genes==g);check(cohort+' '+g+' source gene found',len(ix)>0)
            where=np.flatnonzero(np.isin(indices,ix))
            obs_rows=np.searchsorted(pointer,where,side='right')-1
            recovered.append(np.bincount(obs_rows,weights=values[where],minlength=n))
        counts=np.column_stack(recovered)
        expected=np.flatnonzero(mask & (counts[:,0]>0) & (counts[:,1]>0))
        check(cohort+' exact double-detected membership independently reconstructed',set(expected)==set(z.source_obs_row) and len(expected)==len(z))
        order=z.source_obs_row.to_numpy(dtype=int)
        check(cohort+' all selected raw gene counts match independent CSR parse',np.array_equal(counts[order],z[['PGAM5_raw_count','MYO19_raw_count']].to_numpy()))
        check(cohort+' full-library denominator and source row identity match',close(column(f['obs/total_counts'])[order],z.total_counts) and np.array_equal(column(f['obs'][f['obs'].attrs['_index']])[order],z.source_obs_name))
        check(cohort+' gene detection summary includes complete macrophage pool',int((counts[mask,0]>0).sum())==summary.loc[cohort,'PGAM5_detected'] and int((counts[mask,1]>0).sum())==summary.loc[cohort,'MYO19_detected'])
        source_pairs+=int(mask.sum())*2
        del indices,values,counts
    print('Audited source',cohort,'macrophages',spec['expected_macrophages'],'both',len(z),flush=True)
for row in stats.itertuples():
    z=d[d.cohort.eq(row.cohort)]
    suffix='_log1p_CP10k' if row.expression_scale=='log1p_CP10k' else '_raw_count'
    x=z['PGAM5'+suffix].to_numpy(dtype=float);y=z['MYO19'+suffix].to_numpy(dtype=float)
    name=row.cohort+' '+row.method+'/'+row.expression_scale
    if len(z)<3:
        check(name+' n<3 reportsNA not p=0',row.status=='not_estimable_n<3' and pd.isna(row.correlation) and pd.isna(row.permutation_p) and pd.isna(row.BH_q_permutation));continue
    if np.ptp(x)==0 or np.ptp(y)==0:
        check(name+' constant vector reportsNA',row.status=='not_estimable_constant_vector' and pd.isna(row.correlation) and pd.isna(row.permutation_p) and pd.isna(row.BH_q_permutation));continue
    if row.method=='Spearman':x=rankdata(x);y=rankdata(y)
    x=x-x.mean();y=y-y.mean();den=np.sqrt(np.sum(x*x)*np.sum(y*y))
    r=float(np.sum(x*y)/den)
    if row.method=='Spearman':
        p=0. if abs(r)>=1-1e-15 else 2*t.sf(abs(r)*np.sqrt((len(x)-2)/((1+r)*(1-r))),len(x)-2)
    else:p=2*beta(len(x)/2-1,len(x)/2-1,loc=-1,scale=2).cdf(-abs(r))
    check(name+' independent coefficient and analytical p',close(r,row.correlation) and close(p,row.analytical_p))
    if row.permutation_type=='exact_labeled_permutation':
        unique=list(set(itertools.permutations(y)))
        vals=np.array([np.sum(x*np.asarray(q))/den for q in unique])
        p_exact=float((np.abs(vals)>=abs(r)-1e-12).mean())
        check(name+' exact p independently using distinct rank/value permutations',close(p_exact,row.permutation_p) and row.permutations==math.factorial(len(x)) and close(row.permutation_p,row.extreme_permutations/row.permutations))
    else:check(name+' fixed200000MonteCarlo formula',row.permutations==200000 and close(row.permutation_p,(row.extreme_permutations+1)/(row.permutations+1)))
for family,z in stats.groupby('family'):
    valid=z[z.permutation_p.notna()];p=valid.permutation_p.to_numpy();order=np.argsort(p)
    q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    check(family+' BH independently reconstructed with recorded family size',close(q,valid.BH_q_permutation) and z.BH_family_estimable_tests.eq(len(p)).all())
primary=stats[(stats.method=='Spearman')&stats.expression_scale.eq('log1p_CP10k')]
check('Primary family7estimable,2insufficient,5BH-significant strata',primary.status.eq('estimable').sum()==7 and primary.status.eq('not_estimable_n<3').sum()==2 and primary.BH_q_permutation.lt(.05).sum()==5)
check('GSE140228Droplet includes exactly5HCC+1CC co-detected cells',d[d.cohort.eq('GSE140228_Droplet')].histology.value_counts().to_dict()=={'HCC':5,'CC':1})
prev=pd.read_csv('/workspace/codex/results/GSE151530_PGAM5_MYO19_correlation/macrophage_expression.csv.gz') if Path('/workspace/codex/results/GSE151530_PGAM5_MYO19_correlation/macrophage_expression.csv.gz').exists() else None
if prev is not None:
    prev=prev[(prev.PGAM5_raw_count>0)&(prev.MYO19_raw_count>0)].set_index('cell_id').sort_index()
    now=d[d.cohort.eq('GSE151530')].set_index('cell_id').sort_index()
    fields=['PGAM5_raw_count','MYO19_raw_count','PGAM5_CP10k','MYO19_CP10k']
    check('GSE151530matches22co-detected cells from immediately prior analysis',prev.index.equals(now.index) and close(prev[fields],now[fields]))
provenance=json.loads((R/'source_provenance.json').read_text())
check('Selected expression SHA256 matches extraction record',hashlib.sha256((R/'codetected_macrophage_expression.csv.gz').read_bytes()).hexdigest()==provenance['source_expression_sha256'])
result={'status':'PASS','checks_passed':len(checks),'source_HDF5_audited':not args.saved_results_only,'source_macrophage_gene_cell_counts_checked':source_pairs,'selected_gene_cell_counts_checked':len(d)*2,'checks':checks,'limits':'Implementation audit only; does not establish independent cell/patient replicates or biological coexpression'}
name='audit_saved_results.json' if args.saved_results_only else 'audit.json'
(R/name).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
