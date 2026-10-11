"""Independent raw-matrix, coordinate, graph, normalized RNA, gate and permutation audits."""
from pathlib import Path
import gzip,io,json,hashlib
import numpy as np,pandas as pd,anndata as ad
import scipy.io
from scipy.spatial import cKDTree
from scipy.stats import false_discovery_control
from PIL import Image
from preprocess import ROOT,WORK,SAMPLES,PROGRAMS
checks=[]
def check(name,condition):
    assert bool(condition),name
    checks.append(name)
summary=pd.read_csv(ROOT/'sample_QC_summary.csv').set_index('sample')
stats=pd.read_csv(ROOT/'PGAM5_spatial_autocorrelation.csv')
lr=pd.read_csv(ROOT/'candidate_LR_codetection_summary.csv')
pairs=pd.read_csv(ROOT/'selected_LR_pairs.csv')
for item in json.loads((ROOT/'source_manifest.json').read_text()):
    p=WORK/'inputs'/item['name'];check('source bytes '+p.name,p.stat().st_size==item['bytes']);check('source SHA256 '+p.name,hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256'])
for item in SAMPLES:
    s=item['sample'];pref=WORK/'inputs'/f'{item["gsm"]}_{s}';a=ad.read_h5ad(ROOT/f'{s}_spatial_marker_subset.h5ad')
    f=pd.read_csv(str(pref)+'_features.tsv.gz',sep='\t',header=None);b=pd.read_csv(str(pref)+'_barcodes.tsv.gz',header=None)[0]
    with gzip.open(str(pref)+'_matrix.mtx.gz','rb') as h:raw=scipy.io.mmread(h,spmatrix=True).T.tocsr()
    raw.eliminate_zeros();total=np.asarray(raw.sum(1)).ravel();ng=raw.getnnz(axis=1)
    mt=np.asarray(raw[:,f[1].str.startswith('MT-').to_numpy()].sum(1)).ravel()/total*100
    pos=pd.read_csv(str(pref)+'_tissue_positions_list.csv.gz',header=None,index_col=0)
    pos.columns=['in_tissue','array_row','array_col','pixel_row','pixel_col'];pos=pos.loc[b]
    keep=(pos.in_tissue.to_numpy()==1)&(ng>=200)&(total>=500)&(mt<=25)
    check(s+' exact QC barcodes',list(b[keep])==list(a.obs_names))
    check(s+' raw QC count',keep.sum()==summary.loc[s,'qc_spots'])
    check(s+' total-count normalization library',np.array_equal(a.obs.total_counts.to_numpy(),total[keep]))
    check(s+' raw coordinates',np.array_equal(a.obsm['spatial'],pos.loc[a.obs_names,['pixel_col','pixel_row']].to_numpy()))
    gene_ids=pd.Index(f[0]);idx=gene_ids.get_indexer(a.var.gene_id);check(s+' gene IDs found',np.all(idx>=0))
    counts=raw[keep][:,idx].toarray();check(s+' portable raw counts',np.array_equal(counts,a.layers['counts'].toarray()))
    norm=np.log1p(counts/total[keep,None]*1e4);check(s+' full-library CP10k',np.allclose(norm,a.X.toarray(),rtol=2e-6,atol=2e-6))
    names=pd.Index(f[1]);gx=np.flatnonzero(f[1].eq('PGAM5'));check(s+' unique PGAM5 feature',len(gx)==1)
    pgamraw=raw[keep][:,gx].toarray().ravel();pgam=np.log1p(pgamraw/total[keep]*1e4)
    check(s+' PGAM5 raw',np.array_equal(pgamraw,a.obs.PGAM5_raw));check(s+' PGAM5 expression',np.allclose(pgam,a.obs.PGAM5_log1p_CP10k,rtol=2e-6,atol=2e-6))
    for prog,genes in PROGRAMS.items():
        ix=[np.flatnonzero(f[1].eq(g))[0] for g in genes if g in names]
        x=raw[keep][:,ix].toarray();score=np.log1p(x/total[keep,None]*1e4).mean(1)
        check(s+' program '+prog,np.allclose(score,a.obs[prog],rtol=2e-6,atol=2e-6));check(s+' detections '+prog,np.array_equal((x>0).sum(1),a.obs[prog+'_n_detected']))
    ar=a.obs[['array_row','array_col']].to_numpy().astype(float);ar[:,0]*=np.sqrt(3)
    ij=cKDTree(ar).query_pairs(2.000001,output_type='ndarray');G=a.obsp['spatial_connectivities'];coo=G.tocoo()
    actual=set(zip(coo.row,coo.col));expected=set(map(tuple,ij))|set(map(tuple,ij[:,::-1]))
    check(s+' independent lattice neighbors',actual==expected);check(s+' max degree',np.max(np.asarray(G.sum(1)))<=6);check(s+' no self edges',G.diagonal().sum()==0)
    check(s+' symmetric graph',(G-G.T).nnz==0)
    mac=a.obs.Mac_RNA_enriched.to_numpy();t=a.obs.T_RNA_enriched.to_numpy()
    m=(a.obs.Macrophage_core>=a.obs.Macrophage_core.quantile(.75))&(a.obs.Macrophage_core_n_detected>=3)&(a.obs.Pan_myeloid>=a.obs.Pan_myeloid.median())
    check(s+' Mac gate',np.array_equal(m,mac));check(s+' T gate',np.array_equal(t,(a.obs.T_lineage>=a.obs.T_lineage.quantile(.75))&(a.obs.T_lineage_n_detected>=2)))
    for p in ['Activation_early','Exhaustion_core','Effector_cytotoxic']:
        supported=t&(a.obs[p]>=a.obs.loc[t,p].quantile(.75))&(a.obs[p]>0)&(a.obs[p+'_n_detected']>=2)
        check(s+' T state '+p,np.array_equal(supported,a.obs[p+'_supported']))
    positive=mac&(pgamraw>0);cut=np.quantile(a.obs.PGAM5_log1p_CP10k.to_numpy()[positive],.75)
    check(s+' high PGAM5 gate',np.array_equal(positive&(a.obs.PGAM5_log1p_CP10k>=cut),a.obs.PGAM5_mac_group.eq('PGAM5_high')))
    for p,obs in [('T_RNA_enriched','T_neighbor_count'),('Exhaustion_core_supported','T_exhaustion_neighbor_count'),('Activation_early_supported','T_activation_neighbor_count')]:check(s+' neighbor count '+p,np.array_equal(G@a.obs[p].to_numpy(int),a.obs[obs]))
    im=np.asarray(Image.open(io.BytesIO(gzip.decompress(Path(str(pref)+'_tissue_lowres_image.png.gz').read_bytes()))).convert('RGB'))
    check(s+' original H&E retained',np.array_equal(im,a.uns['spatial'][s]['images']['lowres']))
    sf=json.loads(gzip.decompress(Path(str(pref)+'_scalefactors_json.json.gz').read_bytes()));check(s+' original scale factors',all(a.uns['scalefactors'][k]==v for k,v in sf.items()))
    xy=a.obsm['spatial']*sf['tissue_lowres_scalef'];check(s+' coordinates within H&E',np.all((xy[:,0]>=0)&(xy[:,0]<im.shape[1])&(xy[:,1]>=0)&(xy[:,1]<im.shape[0])))
    for method in ['unrestricted','UMI_MT_stratified']:
        v=np.load(ROOT/f'{s}_{method}_Moran_audit.npz');z=v['z'];groups=v['groups'];scale=len(a)/(G.sum()*(z*z).sum())
        I=scale*sum(z[i]*z[j] for i,j in actual)
        means={g:z[groups==g].mean() for g in np.unique(groups)};var={g:z[groups==g].var() for g in means};n={g:(groups==g).sum() for g in means}
        expected_I=scale*sum(means[groups[i]]*means[groups[j]]-(var[groups[i]]/(n[groups[i]]-1) if groups[i]==groups[j] and n[groups[i]]>1 else 0) for i,j in actual)
        check(s+' Moran I '+method,np.isclose(I,v['observed'],atol=1e-12));check(s+' exact conditional center '+method,np.isclose(expected_I,v['expected'],atol=1e-12))
        for perm in v['example_permutations']:
            check(s+' complete permutation '+method,np.array_equal(np.sort(perm),np.arange(len(a))));check(s+' within-stratum permutation '+method,np.array_equal(groups[perm],groups))
        p=(1+np.sum(np.abs(v['null_I']-expected_I)>=np.abs(I-expected_I)-1e-12))/5000
        row=stats[(stats['sample']==s)&(stats.method==method)].iloc[0];check(s+' permutation p '+method,np.isclose(p,row.permutation_p))
    e=pd.read_csv(ROOT/f'{s}_mac_T_candidate_edges.csv.gz');index=pd.Index(a.obs_names);source=index.get_indexer(e.source_barcode);target=index.get_indexer(e.target_barcode)
    check(s+' edge uniqueness',not e[['source_barcode','target_barcode']].duplicated().any());check(s+' source and target gates',np.all(mac[source]) and np.all(t[target]));check(s+' spatial edges',all((i,j) in actual for i,j in zip(source,target)))
    lookup={g:np.expm1(norm[:,list(a.var_names).index(g)]) for g in a.var_names}
    for pair in pairs.itertuples():
        name=pair.ligand+'--'+pair.receptor_complex;value=np.sqrt(lookup[pair.ligand][source]*np.min(np.vstack([lookup[g][target] for g in pair.receptor_complex.split('_')]),axis=0))
        check(s+' LR all-subunit score '+name,np.allclose(value,e[name],rtol=1e-5,atol=1e-6))
    for row in lr[lr['sample']==s].itertuples():
        d=e if row.source_group=='All_mac_sources' else e[e.source_PGAM5_group.eq(row.source_group)]
        if row.target_state!='All_T_targets':d=d[d['target_activation_supported' if row.target_state=='Activation_supported' else 'target_exhaustion_supported']]
        good=d[row.pair]>0
        check(s+' LR summary '+row.source_group+row.target_state+row.pair,len(d)==row.eligible_edges and good.sum()==row.supported_edges and d.loc[good,'source_barcode'].nunique()==row.supported_unique_sources)
    print('Audited',s,flush=True)
check('32-test BH family',len(stats)==32 and np.allclose(false_discovery_control(stats.permutation_p.to_numpy()),stats.BH_q_32_tests))
check('16 sections / 4 patients',len(SAMPLES)==16 and len({i['patient'] for i in SAMPLES})==4)
check('normal section excluded',all(i['sample']!='A009N' for i in SAMPLES))
check('one recurrent section',sum(i['tumor_status']=='recurrent' for i in SAMPLES)==1)
(ROOT/'audit_results.json').write_text(json.dumps({'passed':True,'checks_passed':len(checks),'checks':checks},ensure_ascii=False,indent=2)+'\n')
print('Passed',len(checks),'checks',flush=True)
