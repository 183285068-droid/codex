"""Independent original-count, gate, neighbor, conditional-null and export checks."""
from pathlib import Path
import json,hashlib,math,gzip
import numpy as np,pandas as pd,anndata as ad,h5py,scipy.sparse as sp
from scipy.spatial import cKDTree
from PIL import Image
from analyze import ROOT,BASE,SAMPLES,PANELS,TARGETS
checks=[]
def check(ok,label):
    if not ok:raise AssertionError(label)
    checks.append(label)
def eq(a,b,label,tol=1e-6):check(np.allclose(a,b,atol=tol,rtol=tol,equal_nan=True),label)
def bh(p):
    p=np.asarray(p);ix=np.argsort(p);v=p[ix]*len(p)/np.arange(1,len(p)+1);result=np.empty(len(p));result[ix]=np.minimum(1,np.minimum.accumulate(v[::-1])[::-1]);return result
def textdecode(v):return [x.decode() if isinstance(x,bytes) else str(x) for x in v]

for r in json.loads((ROOT/'source_inputs/base_requirements.json').read_text()):
    p=BASE/r['path'];check(p.stat().st_size==r['bytes'],r['path']+' frozen bytes');check(hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],r['path']+' frozen SHA256')
original_lr_manifest=json.loads((ROOT.parent/'GSE245908_PGAM5_Tcell_spatial_interactions/input_manifest.json').read_text()) if (ROOT.parent/'GSE245908_PGAM5_Tcell_spatial_interactions/input_manifest.json').exists() else []
for r in original_lr_manifest:
    if r['path'].endswith('/source_inputs/omni_resource.csv'):check(hashlib.sha256((ROOT/'source_inputs/omni_resource.csv').read_bytes()).hexdigest()==r['sha256'],'Same pinned LIANA resource as earlier analysis')
stats=pd.read_csv(ROOT/'spatial_conditional_tests.csv');counts=pd.read_csv(ROOT/'subgroup_neighborhood_counts.csv');lr=pd.read_csv(ROOT/'candidate_LR_by_immune_subgroup.csv');summary=pd.read_csv(ROOT/'source_group_summary.csv');pairs=pd.read_csv(ROOT/'selected_candidate_LR_pairs.csv');resource=pd.read_csv(ROOT/'source_inputs/omni_resource.csv')
for r in pairs.itertuples():check(((resource.source_genesymbol==r.ligand)&(resource.target_genesymbol==r.receptor_complex)).any(),'LR resource definition '+r.pair)

for s,gsm in SAMPLES.items():
    a=ad.read_h5ad(ROOT/f'{s}_immune_marker_subset.h5ad');o=a.obs;ann=pd.read_csv(ROOT/f'{s}_spot_annotations.csv.gz',index_col=0);src=pd.read_csv(ROOT/f'{s}_source_neighbor_subgroups.csv',index_col=0)
    check(list(ann.index)==list(a.obs_names),s+' annotation barcodes')
    with h5py.File(BASE/'inputs'/s/f'{gsm}_{s}_filtered_feature_bc_matrix.h5') as f:
        m=f['matrix'];raw=sp.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:]));barcodes=textdecode(m['barcodes'][:]);ids=textdecode(m['features/id'][:]);names=textdecode(m['features/name'][:])
    bi={b:i for i,b in enumerate(barcodes)};gi={g:i for i,g in enumerate(ids)};ci=np.array([bi[b] for b in a.obs_names]);fi=np.array([gi[g] for g in a.var.gene_ids]);selected=raw[fi][:,ci].T.toarray();total=np.asarray(raw[:,ci].sum(0)).ravel()
    eq(selected,a.layers['counts'].toarray(),s+' full original HDF5 panel counts');eq(total,o.total_counts.to_numpy(),s+' full-library sizes');eq(np.log1p(selected*1e4/total[:,None]),a.X.toarray(),s+' full-library normalized marker expressions')
    check(raw.shape[0]==36601,s+' all original 36601 features');check(len(ids)==len(set(ids)),s+' unique original feature IDs')
    mt=np.array([n.startswith('MT-') for n in names]);pct=np.asarray(raw[mt][:,ci].sum(0)).ravel()/total*100;eq(pct,o.pct_mt.to_numpy(),s+' raw MT fraction')
    check(np.all((total>=500)&(raw[:,ci].getnnz(axis=0)>=200)&(pct<=25)),s+' all main QC conditions')
    pg=np.array([i for i,g in enumerate(names) if g=='PGAM5']);check(len(pg)==1,s+' unique PGAM5 original feature');eq(raw[pg][:,ci].toarray().ravel(),o.PGAM5_raw.to_numpy(),s+' original PGAM5 counts')
    for k,gs in PANELS.items():
        p=a[:,gs];vals=p.X.toarray();cnt=p.layers['counts'].toarray()
        eq(vals.mean(1),o['score_'+k],s+' zero-inclusive program '+k);eq((cnt>0).sum(1),o['detected_'+k],s+' detected marker count '+k)
    # Reconstruct graph by distance on true transformed hex-array geometry, independent of offset lookup.
    xy=np.column_stack([o.array_col.to_numpy(),np.sqrt(3)*o.array_row.to_numpy()]);undirected=cKDTree(xy).query_pairs(2.000001);edges=undirected|{(j,i) for i,j in undirected};G=a.obsp['spatial_connectivities'].tocoo()
    check(edges==set(zip(G.row,G.col)),s+' independently reconstructed hex graph');check(all(i!=j for i,j in edges),s+' no self edges');check(o.spatial_degree.max()<=6,s+' max six neighbors');eq(np.bincount(G.row,minlength=len(a)),o.spatial_degree,s+' genuine degree')
    prior=pd.read_csv(BASE/f'{s}_spatial_edges.csv.gz');lookup={b:i for i,b in enumerate(o.index)};check(set((lookup[r.barcode_from],lookup[r.barcode_to]) for r in prior.itertuples())==edges,s+' preserved prior spatial graph')
    m=(o.score_Macrophage_core>=o.score_Macrophage_core.quantile(.75))&(o.detected_Macrophage_core>=3)&(o.score_Pan_myeloid>=o.score_Pan_myeloid.median())&o.spatial_degree.gt(0)
    t=(o.score_All_T>=o.score_All_T.quantile(.75))&(o.detected_All_T>=2)&o.spatial_degree.gt(0)
    check(np.array_equal(m,o.Mac_RNA_enriched),s+' original Mac-RNA gate');check(np.array_equal(t,o.gate_All_T),s+' original T-RNA gate')
    positive=m&o.PGAM5_raw.gt(0);q=o.loc[positive,'PGAM5_log1p_CP10k'].quantile(.75);high=positive&o.PGAM5_log1p_CP10k.ge(q)
    eq(q,a.uns['PGAM5_high_threshold'],s+' preserved high cutoff');check(np.array_equal(high,o.source_group.eq('PGAM5_high')),s+' preserved Q75 source group');check(list(src.index)==list(o.index[m]),s+' source table barcodes')
    old=ROOT.parent/'GSE245908_PGAM5_Tcell_spatial_interactions'/f'{s}_spot_states.csv'
    if old.exists():
        oldo=pd.read_csv(old,index_col=0);check(np.array_equal(o.source_group,oldo.loc[o.index,'source_group']),s+' exactly matching prior published source labels')
    def rawgene(g):return a[:,[g]].layers['counts'].toarray().ravel()
    def anchor(gs):return np.any(np.vstack([rawgene(g)>0 for g in gs]),axis=0)
    for k in TARGETS[1:]:
        isT=k in ['T_exhaustion','T_activation','CD8_cytotoxic_T','Naive_memory_T','Treg_related'];pop=t if isT else o.spatial_degree.gt(0);cut=float(o.loc[t,'score_'+k].quantile(.75)) if isT else float(o['score_'+k].quantile(.75));minimum=3 if k in ['Plasma_related','pDC_related'] else 2
        target=pop&(o['score_'+k]>=cut)&o['score_'+k].gt(0)&o['detected_'+k].ge(minimum)
        special={'CD8_cytotoxic_T':['CD8A','CD8B'],'NK_related':['KLRD1','NCR1','KLRF1'],'B_related':['MS4A1','CD19'],'Plasma_related':['MZB1','JCHAIN'],'cDC2_related':['CD1C','CLEC10A'],'cDC1_related':['CLEC9A','XCR1'],'Inflammatory_monocyte':['FCN1','VCAN'],'Neutrophil_related':['FCGR3B','CSF3R'],'Mast_related':['TPSB2','CPA3']}
        if k in special:target&=anchor(special[k])
        if k=='Treg_related':target&=rawgene('FOXP3')>0
        if k=='pDC_related':target&=rawgene('CLEC4C')>0
        if k=='NK_related':target&=~t
        if k=='Inflammatory_monocyte':target&=o.score_Pan_myeloid.ge(o.score_Pan_myeloid.median())
        check(np.array_equal(target,o['gate_'+k]),s+' fixed anchor/lineage gate '+k)
    neighbor={i:[] for i in range(len(o))}
    for i,j in edges:neighbor[i].append(j)
    for k in TARGETS:
        gate=o['gate_'+k].to_numpy(bool);n=np.array([int(gate[neighbor[i]].sum()) for i in range(len(o))]);frac=n/np.where(o.spatial_degree.to_numpy()>0,o.spatial_degree.to_numpy(),np.nan)
        eq(src[k+'_neighbor_count'],n[m],s+' explicit neighbor count '+k);eq(src[k+'_neighbor_fraction'],frac[m],s+' explicit neighbor fraction '+k)
        for row in counts[(counts['sample']==s)&(counts.target==k)].itertuples():
            keep=o.source_group.eq(row.source_group).to_numpy();srcidx=np.flatnonzero(keep);tar={j for i in srcidx for j in neighbor[i] if gate[j]}
            check(int(keep.sum())==row.n_sources,s+' descriptive source denominator '+k+row.source_group);check(int(gate.sum())==row.n_target_spots,s+' target denominator '+k+row.source_group);check(int((n[keep]>0).sum())==row.sources_with_neighbor,s+' unique source neighbor support '+k+row.source_group);check(len(tar)==row.unique_neighbor_targets,s+' unique target neighbor support '+k+row.source_group);check(int(n[keep].sum())==row.directed_edges,s+' descriptive directed edges '+k+row.source_group)
            eq(frac[keep].mean(),row.mean_neighbor_fraction,s+' descriptive mean fraction '+k+row.source_group)
    for kind in ['neighbor_fraction','neighbor_present','source_not_target_fraction']:
        labels=TARGETS if kind!='source_not_target_fraction' else TARGETS
        auditsets=[(TARGETS,ROOT/f'{s}_{kind}_permutation_audit.npz')] if kind!='source_not_target_fraction' else [([k],ROOT/f'{s}_{kind}_{k}_permutation_audit.npz') for k in TARGETS]
        for keys,path in auditsets:
            if not path.exists():
                rr=stats[(stats['sample']==s)&(stats.kind==kind)&stats.target.isin(keys)];check(rr.status.eq('not_testable').all(),s+' no-audit insufficient groups '+kind+str(keys));continue
            with np.load(path,allow_pickle=False) as z:
                ix=np.array([lookup[b] for b in z['source_barcodes']]);h=z['high'];groups=z['groups'];mat=z['assignments'].astype(float);Y=z['outcomes'];center=z['conditional_expected'];null=z['null_deltas'];obs=z['observed'];mode=str(z['mode'])
                check(list(z['labels'])==keys,s+' audit target labels '+kind+str(keys));check(np.array_equal(h,high.iloc[ix].astype(float)),s+' high labels '+kind+str(keys));check(np.all(positive.iloc[ix]),s+' positive-only main comparator '+kind+str(keys))
                check(np.allclose(mat.sum(1),h.sum()),s+' global high counts in all assignments '+kind+str(keys));check(np.all((mat==0)|(mat==1)),s+' binary assignments '+kind+str(keys));eh=np.zeros(len(ix));config=1
                # Independently recover exact conditional expectation from per-stratum label counts.
                xy=o[['array_col','array_row']].to_numpy();tile=np.minimum(((xy[ix]-xy.min(0))/(np.ptp(xy,axis=0)+1)*3).astype(int),2);block=tile[:,0]+tile[:,1]*3;expectedcodes=np.zeros(len(ix),int)
                for b in set(block):
                    part=np.flatnonzero(block==b);depth=o.total_counts.to_numpy()[ix[part]];lv=(depth>np.median(depth)).astype(int)
                    for v in [0,1]:expectedcodes[part[lv==v]]=b*2+v
                check(np.array_equal(groups,expectedcodes),s+' spatial/depth strata '+kind+str(keys))
                for g in np.unique(groups):
                    members=np.flatnonzero(groups==g);nh=int(h[members].sum());eh[members]=nh/len(members);config*=math.comb(len(members),nh);check(np.all(mat[:,members].sum(1)==nh),s+' preserved within-stratum high count '+kind+str(keys)+str(g))
                if mode=='exact_enumeration':check(len(mat)==config and len(np.unique(mat,axis=0))==config,s+' complete exact enumeration '+kind+str(keys))
                else:check(len(mat)==4999 and config>4999,s+' full Monte Carlo support '+kind+str(keys))
                Yexpect=[]
                for k in keys:
                    target=o['gate_'+k].to_numpy(bool);n=np.array([int(target[neighbor[i]].sum()) for i in ix]);vals=(n>0).astype(float) if kind=='neighbor_present' else n/o.spatial_degree.to_numpy()[ix];Yexpect.append(vals)
                    if kind=='source_not_target_fraction':check(not target[ix].any(),s+' exclude matching target RNA in source '+k)
                eq(Y,np.array(Yexpect).T,s+' explicit permutation outcomes '+kind+str(keys));observed=Y[h>0].mean(0)-Y[h==0].mean(0);expected=np.average(Y,axis=0,weights=eh)-np.average(Y,axis=0,weights=1-eh)
                eq(obs,observed,s+' independent observed effect '+kind+str(keys));eq(center,expected,s+' independent exact conditional expectation '+kind+str(keys))
                first=np.array([Y[m.astype(bool)].mean(0)-Y[~m.astype(bool)].mean(0) for m in mat[:3]]);eq(null[:3],first,s+' independent first null effects '+kind+str(keys))
                for j,k in enumerate(keys):
                    r=stats[(stats['sample']==s)&(stats.kind==kind)&(stats.target==k)].iloc[0]
                    check(len(ix)==r.n_eligible,s+' statistical denominator '+kind+k);eq(obs[j],r.statistic,s+' exported effect '+kind+k)
                    if r.status=='tested':
                        e=(np.abs(null[:,j]-expected[j])>=abs(obs[j]-expected[j])-1e-12).sum();p=e/len(mat) if mode=='exact_enumeration' else (1+e)/(len(mat)+1);eq(p,r.permutation_p,s+' centered conditional p '+kind+k,1e-12)
                    else:check(pd.isna(r.permutation_p),s+' untestable p missing '+kind+k)
    e=pd.read_csv(ROOT/f'{s}_directed_immune_LR_edges.csv.gz');ri=np.array([lookup[b] for b in e.source_barcode]);ci=np.array([lookup[b] for b in e.target_barcode]);union=o[['gate_'+k for k in TARGETS]].any(axis=1).to_numpy();expected={(i,j) for i,j in edges if m.iloc[i] and union[j]};check(set(zip(ri,ci))==expected,s+' exact Mac to immune-union edge set');check(not e.duplicated(['source_barcode','target_barcode']).any(),s+' no duplicate directed edges')
    for k in TARGETS:check(np.array_equal(e['target_'+k],o['gate_'+k].to_numpy()[ci]),s+' edge target gate '+k)
    for r in pairs.itertuples():
        ligand=rawgene(r.ligand)/total*1e4;rec=(a[:,r.receptor_complex.split('_')].layers['counts'].toarray()/total[:,None]*1e4).min(1);score=np.sqrt(ligand[ri]*rec[ci]);eq(score,e[r.pair].to_numpy(),s+' original counts / all LR subunits '+r.pair)
        for row in lr[(lr['sample']==s)&(lr.pair==r.pair)].itertuples():
            keep=o.source_group.to_numpy()[ri]==row.source_group;keep&=o['gate_'+row.target].to_numpy()[ci];supported=keep&(score>0)
            check(int(keep.sum())==row.eligible_edges,s+' LR eligible denominator '+row.target+row.source_group+r.pair);check(int(supported.sum())==row.codetected_edges,s+' LR co-detection count '+row.target+row.source_group+r.pair);check(len(set(ri[supported]))==row.supported_sources,s+' LR unique supported sources '+row.target+row.source_group+r.pair);check(len(set(ci[supported]))==row.supported_targets,s+' LR unique supported targets '+row.target+row.source_group+r.pair)
    specificity=pd.read_csv(ROOT/'marker_detection_specificity.csv');sensitivity=pd.read_csv(ROOT/'strict_marker_descriptive_sensitivity.csv')
    for row in specificity[specificity['sample']==s].itertuples():
        z=rawgene(row.gene);check(int((z>0).sum())==row.detected_spots,s+' original marker detection '+row.gene);eq(100*(z>0).mean(),row.detected_pct,s+' marker detection percent '+row.gene)
    strict={'Treg_add_IKZF2_and_IL2RA_or_CTLA4':o.gate_Treg_related.to_numpy(bool)&(rawgene('IKZF2')>0)&((rawgene('IL2RA')>0)|(rawgene('CTLA4')>0)),'CD8_cytotoxic_require_CD8B':o.gate_CD8_cytotoxic_T.to_numpy(bool)&(rawgene('CD8B')>0),'Monocyte_require_FCN1':o.gate_Inflammatory_monocyte.to_numpy(bool)&(rawgene('FCN1')>0),'Neutrophil_require_FCGR3B':o.gate_Neutrophil_related.to_numpy(bool)&(rawgene('FCGR3B')>0)}
    for row in sensitivity[sensitivity['sample']==s].itertuples():
        gate=strict[row.rule];nneigh=np.array([int(gate[neighbor[i]].sum()) for i in range(len(o))]);nedge=sum(1 for i,j in edges if high.iloc[i] and gate[j]);check(int(gate.sum())==row.target_spots,s+' strict target count '+row.rule);check(int((nneigh[high]>0).sum())==row.PGAM5_high_sources_with_neighbor,s+' strict source support '+row.rule);check(nedge==row.directed_edges,s+' strict edge support '+row.rule)
    sf=a.uns['spatial'][s]['scalefactors'];xy=a.obsm['spatial']*sf['tissue_lowres_scalef'];im=a.uns['spatial'][s]['images']['lowres'];check(np.all((xy[:,0]>=0)&(xy[:,0]<im.shape[1])&(xy[:,1]>=0)&(xy[:,1]<im.shape[0])),s+' true image coordinate bounds')
    print('Audited',s,flush=True)
check(len(stats)==90,'90 planned tests');eq(stats.BH_q_90_planned_tests,bh(stats.permutation_p.fillna(1).to_numpy()),'Independent joint BH90');primary=stats.kind.eq('neighbor_fraction');check(primary.sum()==30,'30 primary tests');eq(stats.loc[primary,'primary_BH_q_30'],bh(stats.loc[primary,'permutation_p'].fillna(1).to_numpy()),'Independent primary BH30')
for p in (ROOT/'figures').glob('*.png'):
    with Image.open(p) as im:check(im.width>1000 and im.height>500,p.name+' rendered dimensions');im.verify()
    check(p.with_suffix('.pdf').read_bytes().startswith(b'%PDF'),p.name+' PDF exported')
check(len(list((ROOT/'figures').glob('*.png')))==7,'Seven PNG/PDF figure pairs')
(ROOT/'audit_results.json').write_text(json.dumps(dict(status='passed',checks=len(checks),checks_passed=checks),indent=2)+'\n');print('PASSED',len(checks),'independent checks',flush=True)
