"""Independently check saved expression, source CSR counts and correlations."""
from pathlib import Path
import argparse
import hashlib
import json
import h5py
import numpy as np
import pandas as pd
from scipy.stats import rankdata, t, beta

R = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('--h5ad',type=Path,default=Path('/workspace/scratch/PGAM5_myeloid_reference_v2/GSE151530_harmonized_counts_QC.h5ad'))
parser.add_argument('--saved-results-only',action='store_true',help='Skip source h5ad counts; reproduce saved-statistic checks from the portable result bundle')
args = parser.parse_args()
d = pd.read_csv(R/'macrophage_expression.csv.gz')
s = pd.read_csv(R/'correlation_statistics.csv').set_index('analysis')
checks = []

def check(name,condition):
    assert condition,name
    checks.append(name)

def close(x,y):
    return np.allclose(x,y,rtol=1e-10,atol=1e-12,equal_nan=True)

check('Exactly3493existing HCC macrophages with unique cell IDs',len(d)==3493 and d.cell_id.is_unique and d.diagnosis.eq('Hepatocellular carcinoma').all())
for gene in ['PGAM5','MYO19']:
    raw = d[gene+'_raw_count'].to_numpy()
    cp = raw/d.total_counts.to_numpy()*10000
    check(gene+' raw counts integer and nonnegative',np.isfinite(raw).all() and (raw>=0).all() and (raw==np.floor(raw)).all())
    check(gene+' CP10k/log1p independently reconstructed',close(cp,d[gene+'_CP10k']) and close(np.log1p(cp),d[gene+'_log1p_CP10k']))
    check(gene+' detection label is raw_count>0',np.array_equal(raw>0,d[gene+'_detected']))
check('PGAM5=161,MYO19=189,both22,neither3165',d.PGAM5_detected.sum()==161 and d.MYO19_detected.sum()==189 and (d.PGAM5_detected&d.MYO19_detected).sum()==22 and (~d.PGAM5_detected&~d.MYO19_detected).sum()==3165)
check('27libraries mapped to20existing donor proxies',d.Sample.nunique()==27 and d.donor_id.nunique()==20 and d.groupby('Sample').donor_id.nunique().eq(1).all())
if not args.saved_results_only:
    with h5py.File(args.h5ad,'r') as f:
        def column(node):
            if isinstance(node,h5py.Group):
                if 'values' in node:
                    assert not node['mask'][:].any()
                    return column(node['values'])
                cats=column(node['categories']);codes=node['codes'][:]
                assert (codes>=0).all()
                return cats[codes]
            if h5py.check_string_dtype(node.dtype) is not None:
                return node.asstr()[:]
            return node[:]
        mask=column(f['obs/harmonized_group'])=='Macrophage'
        source_rows=np.flatnonzero(mask)
        source_columns=column(f['obs/matrix_column'])[mask]
        genes=column(f['var'][f['var'].attrs['_index']])
        gene_index=[int(np.flatnonzero(genes==g)[0]) for g in ['PGAM5','MYO19']]
        by_column=d.set_index('matrix_column').loc[source_columns]
        check('All3493source macrophage memberships/IDs/library totals match',len(source_rows)==3493 and set(source_columns)==set(d.matrix_column) and np.array_equal(column(f['obs/cell_id'])[mask],by_column.cell_id) and close(column(f['obs/total_counts'])[mask],by_column.total_counts))
        pointer=f['X/indptr'][:]
        recovered=np.zeros((len(source_rows),2),dtype=float)
        for k,row in enumerate(source_rows):
            start,end=pointer[row:row+2]
            indices=f['X/indices'][start:end]
            values=f['X/data'][start:end]
            for j,g in enumerate(gene_index):
                recovered[k,j]=values[indices==g].sum()
        check('All6986gene-cell raw counts independently parsed from HDF5 CSR',np.array_equal(recovered,by_column[['PGAM5_raw_count','MYO19_raw_count']].to_numpy()))

donors=pd.read_csv(R/'donor_proxy_mean_expression.csv').set_index('donor_id')
means={}
for donor in sorted(d.donor_id.unique()):
    z=d[d.donor_id.eq(donor)]
    means[donor]=[z.PGAM5_CP10k.sum()/len(z),z.MYO19_CP10k.sum()/len(z)]
    check('Donor '+donor+' count and mean expression independently reconstructed',int(donors.loc[donor,'macrophages'])==len(z) and close(means[donor],donors.loc[donor,['PGAM5_mean_CP10k','MYO19_mean_CP10k']].to_numpy(dtype=float)))
positive=d.PGAM5_detected.to_numpy()
both=positive & d.MYO19_detected.to_numpy()
x=d.PGAM5_log1p_CP10k.to_numpy();y=d.MYO19_log1p_CP10k.to_numpy()
spec=[('Spearman_all_cells_log1p_CP10k',x,y,True),('Pearson_all_cells_log1p_CP10k',x,y,False),('Spearman_PGAM5_detected_log1p_CP10k',x[positive],y[positive],True),('Spearman_both_detected_log1p_CP10k',x[both],y[both],True),('Spearman_both_detected_raw_counts',d.PGAM5_raw_count.to_numpy()[both],d.MYO19_raw_count.to_numpy()[both],True),('Spearman_donor_proxy_mean_CP10k',np.array(list(means.values()))[:,0],np.array(list(means.values()))[:,1],True)]
for name,xx,yy,ranks in spec:
    n=len(xx)
    if ranks:xx=rankdata(xx);yy=rankdata(yy)
    xx=xx-xx.mean();yy=yy-yy.mean()
    r=np.sum(xx*yy)/np.sqrt(np.sum(xx*xx)*np.sum(yy*yy))
    p=2*t.sf(abs(r)*np.sqrt((n-2)/((1+r)*(1-r))),n-2) if ranks else 2*beta(n/2-1,n/2-1,loc=-1,scale=2).cdf(-abs(r))
    row=s.loc[name]
    check(name+' manual coefficient and analytical p',int(row.n)==n and close([r,p],[row.correlation,row.asymptotic_p]))
z=s[s.family.eq('secondary')]
pvalues=z.asymptotic_p.to_numpy();order=np.argsort(pvalues)
adjusted=np.empty(len(order));adjusted[order]=np.minimum(1,np.minimum.accumulate((pvalues[order]*len(order)/np.arange(1,len(order)+1))[::-1])[::-1])
check('One primary test and five-secondary BH family independently checked',s.family.eq('primary').sum()==1 and len(z)==5 and close(adjusted,z.secondary_BH_q))
primary=s.loc['Spearman_all_cells_log1p_CP10k']
check('Fixed200000permutation p formula and extreme count',int(primary.permutations)==200000 and close(primary.permutation_p,(primary.extreme_permutations+1)/(primary.permutations+1)))
provenance=json.loads((R/'source_provenance.json').read_text())
check('Extracted expression file matches recorded SHA256',hashlib.sha256((R/'macrophage_expression.csv.gz').read_bytes()).hexdigest()==provenance['extracted_expression_sha256'])
result={'status':'PASS','checks_passed':len(checks),'source_HDF5_counts_audited':not args.saved_results_only,'checks':checks,'scope':'Implementation checks; not an independent biological validation or proof of cell/patient exchangeability'}
output='audit_saved_results.json' if args.saved_results_only else 'audit.json'
(R/output).write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
