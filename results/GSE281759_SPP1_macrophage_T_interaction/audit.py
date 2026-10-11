"""Independent QC, adjacency, OLS/permutation, FDR and export verification."""
from pathlib import Path
import hashlib,json
import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import rankdata
from PIL import Image

ROOT=Path(__file__).resolve().parent;REF=ROOT.parent/'GSE281759_spatial_HCC_PGAM5'
SAMPLES=json.loads((ROOT/'sample_scope.json').read_text());PROGRAMS=json.loads((ROOT/'programs.json').read_text());PROGRAMS['Exhaustion_no_CTLA4_TIGIT']=['PDCD1','LAG3','TOX','CXCL13','HAVCR2','ENTPD1']
stats=pd.read_csv(ROOT/'interaction_statistics.csv');summary=pd.read_csv(ROOT/'section_summary.csv');lrsummary=pd.read_csv(ROOT/'candidate_LR_codetection_summary.csv');pairs=pd.read_csv(ROOT/'tested_candidate_LR_pairs.csv');cfg=json.loads((ROOT/'analysis_config.json').read_text())
checks=[]
def check(ok,label):
    if not ok:raise AssertionError(label)
    checks.append(label)
def eq(a,b,label,tol=1e-6):check(np.allclose(np.asarray(a),np.asarray(b),atol=tol,rtol=tol,equal_nan=True),label)
def dense(x):return x.toarray() if hasattr(x,'toarray') else np.asarray(x)
def bh(p):
    p=np.asarray(p);ix=np.argsort(p);s=p[ix]*len(p)/np.arange(1,len(p)+1);out=np.empty(len(p));out[ix]=np.minimum(1,np.minimum.accumulate(s[::-1])[::-1]);return out
def residual(v,C):
    r=rankdata(v);return r-C@np.linalg.lstsq(C,r,rcond=None)[0]

for entry in json.loads((ROOT/'reference_manifest.json').read_text()):
    p=REF/entry['name'];check(p.stat().st_size==entry['bytes'],entry['name']+' reference size');check(hashlib.sha256(p.read_bytes()).hexdigest()==entry['sha256'],entry['name']+' reference SHA256')

for item in SAMPLES:
    s=item['sample'];a=ad.read_h5ad(REF/f'{s}_spatial_marker_subset.h5ad');o=a.obs
    ann=pd.read_csv(ROOT/f'{s}_spot_annotations.csv.gz',index_col=0);d=pd.read_csv(ROOT/f'{s}_source_neighbor_RNA.csv').set_index('barcode');r=summary.set_index('sample').loc[s]
    check(list(ann.index)==list(a.obs_names),s+' all annotation barcodes')
    cnt=dense(a[:,['SPP1']].layers['counts']).ravel();log=np.log1p(cnt*1e4/o.total_counts.to_numpy())
    eq(ann.SPP1_raw,cnt,s+' raw SPP1');eq(ann.SPP1_log1p_CP10k,log,s+' full-library normalized SPP1')
    for p,genes in PROGRAMS.items():
        x=dense(a[:,genes].X);raw=dense(a[:,genes].layers['counts']);score=np.log1p(raw*1e4/o.total_counts.to_numpy()[:,None]).mean(1)
        eq(x,np.log1p(raw*1e4/o.total_counts.to_numpy()[:,None]),s+' full-library normalization '+p)
        eq(ann[p],score,s+' zero-inclusive program '+p)
    m=(o.Macrophage_core>=o.Macrophage_core.quantile(.75))&(o.Macrophage_core_n_detected>=3)&(o.Pan_myeloid>=o.Pan_myeloid.median())
    t=(o.T_lineage>=o.T_lineage.quantile(.75))&(o.T_lineage_n_detected>=2)
    check(np.array_equal(ann.Mac_RNA_enriched,m),s+' original macrophage gate');check(np.array_equal(ann.T_RNA_enriched,t),s+' original T gate')
    check(np.array_equal(ann.Mac_RNA_nonT,m&~t),s+' nonT source sensitivity')
    check(list(d.index)==list(o.index[m]),s+' source barcodes')
    for p in ['Exhaustion_core','Activation_early','Effector_cytotoxic']:
        gate=t&(o[p]>=o.loc[t,p].quantile(.75))&(o[p]>0)&(o[p+'_n_detected']>=2)
        check(np.array_equal(ann[p+'_supported'],gate),s+' target state '+p)
    pos=m.to_numpy()&(cnt>0);threshold=np.quantile(log[pos],.75)
    eq(threshold,r.SPP1_high_threshold,s+' positive Q75');check(np.array_equal(ann.SPP1_mac_group.eq('SPP1_high'),pos&(log>=threshold-1e-7)),s+' high groups / ties')
    check(int(pos.sum())==r.SPP1_mac_detected,s+' detection count');check(int(m.sum())==r.mac_RNA_spots,s+' macrophage count');check(len(a)==r.qc_spots,s+' QC spots')
    # Independent KD-tree reconstruction of the hexagonal lattice.
    xy=np.column_stack([o.array_col.to_numpy(),np.sqrt(3)*o.array_row.to_numpy()]);links=cKDTree(xy).query_pairs(2.000001)
    expect={(i,j) for i,j in links}|{(j,i) for i,j in links};g=a.obsp['spatial_connectivities'].tocoo();actual=set(zip(g.row,g.col))
    check(expect==actual,s+' independently reconstructed genuine hex graph');check(all(i!=j for i,j in actual),s+' no self edges');check(np.array_equal(np.bincount(g.row,minlength=len(a)),o.spatial_degree),s+' spatial degrees');check(o.spatial_degree.max()<=6,s+' six-neighbor bound')
    xyz=a.obsm['spatial']*a.uns['scalefactors']['tissue_lowres_scalef'];im=a.uns['spatial'][s]['images']['lowres'];check(np.all((xyz[:,0]>=0)&(xyz[:,0]<im.shape[1])&(xyz[:,1]>=0)&(xyz[:,1]<im.shape[0])),s+' real pixel/image bounds')
    bysource={i:[] for i in range(len(a))}
    for i,j in expect:bysource[i].append(j)
    allnums=[];exnums=[];actnums=[];neighbor={k:[] for k in ['T_lineage','log_UMI','pct_mt']+cfg['programs']}
    values={'T_lineage':ann.T_lineage.to_numpy(),'log_UMI':np.log1p(ann.total_counts.to_numpy()),'pct_mt':ann.pct_mt.to_numpy(),**{p:ann[p].to_numpy() for p in cfg['programs']}}
    for i in np.flatnonzero(m):
        ix=np.array(bysource[i],int);T=ix[t.to_numpy()[ix]]
        allnums.append(len(T));exnums.append(int(ann.Exhaustion_core_supported.iloc[ix].sum()));actnums.append(int(ann.Activation_early_supported.iloc[ix].sum()))
        for key in neighbor:neighbor[key].append(values[key][T].mean() if len(T) else np.nan)
    for key,v in [('All_T',allnums),('Exhaustion',exnums),('Activation',actnums)]:
        eq(d[key+'_neighbor_count'],v,s+' explicit neighbor counts '+key);check(np.array_equal(d[key+'_neighbor_present'],np.asarray(v)>0),s+' neighbor presence '+key)
    for key,v in neighbor.items():eq(d['neighbor_'+key],v,s+' explicit T-neighbor mean '+key)
    for key,val in [('source_log_UMI',np.log1p(o.total_counts)),('source_pct_mt',o.pct_mt),('source_Macrophage_core',o.Macrophage_core),('source_Hepatocyte',o.Hepatocyte),('source_Epithelial',o.Epithelial),('source_T_lineage',o.T_lineage)]:eq(d[key],np.asarray(val)[m],s+' source covariate '+key)
    for test in stats[stats['sample']==s].itertuples():
        base=d[d.spatial_degree>0] if test.family=='adjacency' else d[d.T_neighbor_count>0]
        if test.cohort.startswith('nonT'):base=base[base.Mac_RNA_nonT]
        if test.cohort=='high_vs_other_detected':base=base[base.SPP1_raw>0]
        check(len(base)==test.n_sources,s+' analysis denominator '+test.family+test.endpoint+test.cohort)
        if test.status!='tested':check(pd.isna(test.permutation_p),s+' untestable p missing '+test.family+test.endpoint+test.cohort);continue
        path=ROOT/f'{s}_{test.family}_{test.endpoint}_{test.cohort}_audit.npz'
        with np.load(path,allow_pickle=False) as z:
            label=s+' '+test.family+' '+test.endpoint+' '+test.cohort
            check(list(z['source_barcodes'])==list(base.index),label+' ordered barcodes')
            x=z['x'];y=z['y'];groups=z['groups'];perms=z['example_permutations'];null=z['null'];obs=float(z['observed'])
            bins=[]
            for k in ['source_log_UMI','source_Macrophage_core']:bins.append(np.nan_to_num(pd.qcut(base[k],2,labels=False,duplicates='drop').to_numpy(float),nan=0).astype(int))
            eq(groups,base.spatial_tile.to_numpy()*4+bins[0]*2+bins[1],label+' stratification')
            check(len(null)==4999 and np.isfinite(null).all(),label+' complete finite null')
            for perm in perms:
                check(np.array_equal(np.sort(perm),np.arange(len(base))),label+' permutation bijection');check(np.array_equal(groups[perm],groups),label+' within-stratum permutation')
            if test.family=='adjacency':
                h=(base.SPP1_mac_group.eq('SPP1_high') if test.cohort=='high_vs_other_detected' else base.SPP1_raw.gt(0)).to_numpy();yy=base[test.endpoint+'_neighbor_present'].to_numpy(float)
                eq(x,h.astype(float),label+' source labels');eq(y,yy,label+' fixed target-neighbor outcome')
                value=100*(yy[h].mean()-yy[~h].mean());v=np.array([100*(yy[h[p]].mean()-yy[~h[p]].mean()) for p in perms])
                eq(test.positive_neighbor_pct,100*yy[h].mean(),label+' positive proportion');eq(test.comparator_neighbor_pct,100*yy[~h].mean(),label+' comparator proportion')
            else:
                eq(x,base.SPP1_log1p_CP10k,label+' continuous source including zeros');eq(y,base['neighbor_'+test.endpoint],label+' fixed neighboring program')
                eq(z['covariates'],base[cfg['covariates']].to_numpy(float),label+' nine covariates')
                C=np.column_stack([np.ones(len(base)),rankdata(z['covariates'],axis=0)]);rx=residual(x,C);ry=residual(y,C);value=np.corrcoef(rx,ry)[0,1]
                v=np.array([np.corrcoef(residual(x[p],C),ry)[0,1] for p in perms])
            eq(value,obs,label+' independent observed statistic',1e-9);eq(value,test.statistic,label+' exported statistic',1e-9);eq(v,null[:3],label+' independent first null statistics',1e-9)
            p=min(1.,2*min((1+sum(null>=obs-1e-12))/5000,(1+sum(null<=obs+1e-12))/5000))
            eq(p,test.permutation_p,label+' equal-tail p',1e-12)
    e=pd.read_csv(ROOT/f'{s}_mac_T_candidate_edges.csv.gz');lookup={b:i for i,b in enumerate(o.index)};si=np.array([lookup[b] for b in e.source_barcode]);ti=np.array([lookup[b] for b in e.target_barcode])
    edges={(i,j) for i,j in expect if m.iloc[i] and t.iloc[j]};check(set(zip(si,ti))==edges,s+' exact Mac to T edge set')
    eq(e.source_SPP1_raw,cnt[si],s+' edge source SPP1');check(np.array_equal(e.source_nonT,(m&~t).to_numpy()[si]),s+' edge strict source status')
    for p in pairs.itertuples():
        lig=np.expm1(dense(a[:,[p.ligand]].X).ravel().astype(float));rec=np.min(np.expm1(dense(a[:,p.receptor_complex.split('_')].X).astype(float)),axis=1)
        eq(e[p.pair],np.sqrt(lig[si]*rec[ti]),s+' independent LR expression '+p.pair)
        for row in lrsummary[(lrsummary['sample']==s)&(lrsummary.pair==p.pair)].itertuples():
            mask=np.ones(len(e),bool) if row.source_group=='All_mac' else e.source_SPP1_raw.gt(0).to_numpy() if row.source_group=='SPP1_detected' else e.SPP1_source_group.eq(row.source_group).to_numpy()
            if row.target_state!='All_T':mask&=e['target_'+('exhaustion' if row.target_state=='Exhaustion' else 'activation')+'_supported'].to_numpy()
            dd=e[mask];sup=dd[p.pair]>0
            check(len(dd)==row.eligible_edges,s+' LR denominator '+p.pair+row.source_group+row.target_state)
            check(int(sup.sum())==row.supported_edges,s+' LR support '+p.pair+row.source_group+row.target_state)
            eq(100*sup.mean() if len(dd) else np.nan,row.edge_codetection_pct,s+' LR rate '+p.pair+row.source_group+row.target_state)
    desc=pd.read_csv(ROOT/'source_group_neighbor_summary.csv');base=d[d.spatial_degree>0]
    for row in desc[desc['sample']==s].itertuples():
        mask=base.SPP1_raw.gt(0) if row.source_group=='SPP1_detected' else base.SPP1_raw.eq(0) if row.source_group=='Undetected' else base.SPP1_mac_group.eq('SPP1_high')
        dd=base[mask];check(len(dd)==row.n_sources,s+' descriptive denominator '+row.source_group+row.target_state)
        eq(100*dd[row.target_state+'_neighbor_present'].mean() if len(dd) else np.nan,row.neighbor_pct,s+' descriptive rate '+row.source_group+row.target_state)
    print('Audited',s,flush=True)

check(len(stats)==176,'All 176 planned tests present');eq(stats.BH_q_176_planned_tests,bh(stats.permutation_p.fillna(1).to_numpy()),'Independent BH for joint family')
primary=(stats.family=='adjacency')&(stats.endpoint=='All_T')&(stats.cohort=='detected_vs_undetected');check(primary.sum()==16,'16 primary adjacency tests');eq(stats.loc[primary,'primary_BH_q_16_sections'],bh(stats.loc[primary,'permutation_p'].fillna(1).to_numpy()),'Independent primary BH16')
only=primary&stats.tumor_status.eq('primary');check(only.sum()==15,'15 primary-tumor sections');eq(stats.loc[only,'primary_only_BH_q_15_sections'],bh(stats.loc[only,'permutation_p'].fillna(1).to_numpy()),'Independent primary-only BH15')
pc=pd.read_csv(ROOT/'patient_counts.csv').set_index('patient');check(len(pc)==4,'Four independent patients')
for k in ['qc_spots','SPP1_mac_detected','mac_RNA_spots','SPP1_detected_with_T_neighbor','SPP1_detected_with_Tex_neighbor','SPP1_detected_with_Tact_neighbor']:eq(pc[k],summary.groupby('patient')[k].sum(),k+' patient totals')
pe=pd.read_csv(ROOT/'patient_descriptive_effect_summary.csv')
for row in pe.itertuples():
    dd=stats[(stats.patient==row.patient)&(stats.family==row.family)&(stats.endpoint==row.endpoint)&(stats.cohort==row.cohort)]
    if row.scope=='primary_only':dd=dd[dd.tumor_status=='primary']
    dd=dd[dd.status=='tested'];eq(dd.statistic.median(),row.median_tested_section_statistic,'Patient equal-section median '+row.scope+row.patient+row.family+row.endpoint+row.cohort)
pngs=list((ROOT/'figures').glob('*.png'));check(len(pngs)==9,'Nine figure pairs')
for p in pngs:
    with Image.open(p) as im:check(im.width>1000 and im.height>500,p.name+' rendered size');im.verify()
    check(p.with_suffix('.pdf').read_bytes().startswith(b'%PDF'),p.name+' PDF export')
out=dict(status='passed',checks=len(checks),tested=int((stats.status=='tested').sum()),not_testable=int((stats.status!='tested').sum()),checks_passed=checks)
(ROOT/'audit_results.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED',len(checks),'independent checks',flush=True)
