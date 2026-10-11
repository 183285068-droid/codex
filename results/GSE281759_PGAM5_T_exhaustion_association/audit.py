"""Independent matrix-edge/outcome/rank-regression/permutation/FDR checks."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import rankdata,false_discovery_control
from PIL import Image
ROOT=Path(__file__).resolve().parent;REF=ROOT.parent/'GSE281759_spatial_HCC_PGAM5'
config=json.loads((ROOT/'analysis_config.json').read_text());samples=json.loads((REF/'sample_scope.json').read_text());tab=pd.read_csv(ROOT/'association_statistics.csv',float_precision='round_trip');COV=config['covariates'];checks=[]
def check(name,ok):
    assert bool(ok),name;checks.append(name)
def adjusted_stat(x,y,cov,high,kind):
    C=np.column_stack([np.ones(len(x)),rankdata(cov,axis=0)]);yr=rankdata(y);yr-=C@np.linalg.lstsq(C,yr,rcond=None)[0]
    if kind=='adjusted_high_minus_other':return (yr[high].mean()-yr[~high].mean())/np.std(yr,ddof=1)
    xr=rankdata(x);xr-=C@np.linalg.lstsq(C,xr,rcond=None)[0];return np.corrcoef(xr,yr)[0,1]
for it in samples:
    s=it['sample'];a=ad.read_h5ad(REF/f'{s}_spatial_marker_subset.h5ad');o=a.obs;d=pd.read_csv(ROOT/f'{s}_source_neighbor_RNA.csv',float_precision='round_trip');mac=o.Mac_RNA_enriched.to_numpy(bool);t=o.T_RNA_enriched.to_numpy(bool);G=a.obsp['spatial_connectivities'];adj=G.tocoo();neighbors={i:[] for i in range(len(a))}
    for i,j in zip(adj.row,adj.col):
        if t[j]:neighbors[i].append(j)
    check(s+' unique source barcodes',d.barcode.is_unique);check(s+' all Mac-RNA sources exactly',list(d.barcode)==list(a.obs_names[mac]))
    ix=pd.Index(a.obs_names).get_indexer(d.barcode);check(s+' source exposure with zeros',np.array_equal(d.PGAM5_raw,o.PGAM5_raw.iloc[ix]) and np.allclose(d.PGAM5_log1p_CP10k,o.PGAM5_log1p_CP10k.iloc[ix],atol=1e-14))
    for i,row in enumerate(d.itertuples()):
        ns=neighbors[ix[i]];check(s+' source T-neighbor count '+row.barcode,len(ns)==row.T_neighbor_count)
        for p in config['programs']+['T_lineage','log_UMI','pct_mt']:
            values=np.log1p(o.total_counts.iloc[ns].to_numpy()) if p=='log_UMI' else o[p].iloc[ns].to_numpy(dtype=float)
            expected=np.mean(values) if ns else np.nan;actual=getattr(row,'neighbor_'+p)
            check(s+' neighbor '+p+' '+row.barcode,(np.isnan(expected) and np.isnan(actual)) or np.isclose(expected,actual,atol=1e-13))
    valid=d[d.T_neighbor_count>0]
    for row in tab[tab['sample']==s].itertuples():
        dd=valid if row.cohort=='all_mac_including_zero' else valid[valid.PGAM5_raw>0]
        check(s+' source selection '+row.program+row.contrast,row.n_sources==len(dd))
        if row.status!='tested':
            check(s+' no p for untestable '+row.program+row.contrast,np.isnan(row.permutation_p));continue
        name=s+'_'+row.cohort+'_'+row.program+'_'+row.contrast
        v=np.load(ROOT/(name+'_audit.npz'))
        check(name+' saved source IDs',np.array_equal(v['source_barcodes'],dd.barcode.to_numpy(str)))
        check(name+' full source exposure',np.array_equal(v['x'],dd.PGAM5_log1p_CP10k.to_numpy()))
        check(name+' zero outcome retained',np.array_equal(v['y'],dd['neighbor_'+row.program].to_numpy()))
        check(name+' nine covariates',np.array_equal(v['covariates'],dd[COV].to_numpy(float)))
        bins=dd.spatial_tile.to_numpy(int)*2+pd.qcut(dd.source_log_UMI,2,labels=False,duplicates='drop').to_numpy(int)
        check(name+' spatial/UMI strata',np.array_equal(bins,v['groups']))
        stat=adjusted_stat(v['x'],v['y'],v['covariates'],v['high'],row.contrast);check(name+' independent OLS statistic',np.isclose(stat,row.statistic,atol=1e-10))
        for k,perm in enumerate(v['example_permutations']):
            check(name+' permutation '+str(k),np.array_equal(np.sort(perm),np.arange(len(dd))) and np.array_equal(v['groups'][perm],v['groups']))
            calc=adjusted_stat(v['x'][perm],v['y'],v['covariates'],v['high'][perm],row.contrast)
            check(name+' independently recomputed null '+str(k),np.isclose(calc,v['null'][k],atol=1e-10))
        null=v['null'];check(name+' 4999 finite null statistics',len(null)==4999 and np.isfinite(null).all())
        p=min(1,2*min((1+(null>=stat-1e-12).sum())/5000,(1+(null<=stat+1e-12).sum())/5000))
        check(name+' equal-tail p',np.isclose(p,row.permutation_p))
        if row.contrast=='adjusted_high_minus_other':
            high=dd.PGAM5_mac_group.eq('PGAM5_high').to_numpy();y=dd.neighbor_Exhaustion_core.to_numpy()
            check(name+' observed group sizes',high.sum()==row.n_high and (~high).sum()==row.n_other_detected)
            check(name+' raw mean difference',np.isclose(y[high].mean()-y[~high].mean(),row.raw_mean_difference,atol=1e-12))
    print('Audited',s,flush=True)
check('planned 96 tests',len(tab)==96)
check('joint BH96',np.allclose(false_discovery_control(tab.permutation_p.fillna(1).to_numpy()),tab.BH_q_96_planned_tests))
main=tab.primary_BH_q_16_sections.notna();check('main planned 16 tests',main.sum()==16);check('main BH16',np.allclose(false_discovery_control(tab.loc[main,'permutation_p'].fillna(1).to_numpy()),tab.loc[main,'primary_BH_q_16_sections']))
check('four patients',tab.patient.nunique()==4)
figs=list((ROOT/'figures').glob('*.png'));check('8 standalone PNG figures',len(figs)==8)
for p in figs:
    with Image.open(p) as im:im.verify()
    check(p.name+' PDF companion',(p.with_suffix('.pdf')).read_bytes().startswith(b'%PDF-'))
(ROOT/'audit_results.json').write_text(json.dumps(dict(passed=True,checks_passed=len(checks),checks=checks),ensure_ascii=False,indent=2)+'\n');print('Passed',len(checks),'checks')
