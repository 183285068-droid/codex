"""Independent input, raw count, gene/cell header, normalization, QC and spatial audit."""
from pathlib import Path
import json,gzip,io,hashlib,subprocess,sys
import numpy as np,pandas as pd,anndata as ad,scipy.io,h5py
from scipy.sparse import csc_matrix
from scipy.spatial import cKDTree
from scipy.stats import false_discovery_control
from PIL import Image
from spatial_core import ROOT,WORK,PROGRAMS
checks=[]
def check(name,ok):
    assert bool(ok),name;checks.append(name)
for item in json.loads((ROOT/'source_manifest.json').read_text()):
    f=WORK/'inputs'/item['name'];check('source length '+f.name,f.stat().st_size==item['bytes']);check('source SHA '+f.name,hashlib.sha256(f.read_bytes()).hexdigest()==item['sha256'])
a=ad.read_h5ad(ROOT/'singlecell_marker_subset.h5ad');scobs=[];raws=[];cellids=[];features=None
for lib in range(1,4):
    prefix=WORK/'inputs'/f'GSM6177611_NYU_LIHC1_lib{lib}_lib{lib}'
    f=pd.read_csv(str(prefix)+'_genes.tsv.gz',sep='\t');comments=[]
    with gzip.open(str(prefix)+'_gene_expression.mtx.gz','rt') as h:
        h.readline()
        for line in h:
            if not line.startswith('%'):break
            if '\t' in line:comments.append(line.rstrip()[1:].split('\t'))
    check(f'lib{lib} feature header order',comments[0]==f.name.tolist());check(f'lib{lib} feature IDs identical',features is None or f.equals(features));features=f
    with gzip.open(str(prefix)+'_gene_expression.mtx.gz','rb') as h:x=scipy.io.mmread(h,spmatrix=True).T.tocsr()
    x.eliminate_zeros();total=np.asarray(x.sum(1)).ravel();ng=x.getnnz(axis=1)
    mt=np.asarray(x[:,f.name.str.startswith('MT-').to_numpy()].sum(1)).ravel()/total*100;ribo=np.asarray(x[:,f.name.str.match(r'^RP[SL]\d').to_numpy()].sum(1)).ravel()/total*100
    keep=(total>=500)&(ng>=200)&(mt<=25)&(ribo<=30)
    obs=pd.read_csv(ROOT/f'lib{lib}_all_cell_QC.csv',index_col=0,float_precision='round_trip');check(f'lib{lib} original cell IDs',obs.original_cell_id.tolist()==comments[1]);check(f'lib{lib} QC',np.array_equal(keep,obs.pass_qc));check(f'lib{lib} raw dimensions',x.shape==(2000,19936));check(f'lib{lib} total/genes',np.array_equal(total,obs.total_counts) and np.array_equal(ng,obs.n_genes))
    raws.append(x[keep]);cellids.extend([f'lib{lib}|{c}' for c,m in zip(comments[1],keep) if m])
from scipy.sparse import vstack
raw=vstack(raws).tocsr();check('single-cell stable original IDs',cellids==list(a.obs_names));check('single-cell unique IDs',a.obs_names.is_unique);check('single-cell retained 1783',len(a)==1783)
def expression_checks(obj,raw,fids,fnames,label):
    total=np.asarray(raw.sum(1)).ravel();lookup=pd.Index(fids);ix=lookup.get_indexer(obj.var.gene_id);check(label+' unique matched gene IDs',np.all(ix>=0) and lookup.is_unique)
    count=raw[:,ix].toarray();check(label+' portable raw count',np.array_equal(count,obj.layers['counts'].toarray()))
    norm=np.log1p(count/total[:,None]*1e4);check(label+' full-library normalization',np.allclose(norm,obj.X.toarray(),rtol=2e-6,atol=2e-6));check(label+' total library counts',np.array_equal(total,obj.obs.total_counts))
    for p,genes in PROGRAMS.items():
        indices=[np.flatnonzero(np.asarray(fnames)==g)[0] for g in genes if g in fnames];x=raw[:,indices].toarray();score=np.log1p(x/total[:,None]*1e4).mean(1)
        check(label+' program '+p,np.allclose(score,obj.obs[p],rtol=2e-6,atol=2e-6));check(label+' program detection '+p,np.array_equal((x>0).sum(1),obj.obs[p+'_n_detected']))
expression_checks(a,raw,features.ensembl_id.to_numpy(),features.name.to_numpy(),'single-cell')
flag=a.obs.Mac_RNA_supported&a.obs.expression_cluster.isin(['C3','C6']);check('macrophage cluster/lineage intersection',np.array_equal(flag,a.obs.Mac_RNA_cluster_confirmed));check('246 lineage-supported candidates',flag.sum()==246);check('6 PGAM5 detections all 1 UMI',(a.obs.loc[flag,'PGAM5_raw']>0).sum()==6 and np.all(a.obs.loc[flag&a.obs.PGAM5_raw.gt(0),'PGAM5_raw']==1));check('macrophage related C3/C6 counts',a.obs.expression_cluster.isin(['C3','C6']).sum()==325)
for assay,ob in [('singlecell',a),('spatial',ad.read_h5ad(ROOT/'NYU_LIHC1_Vis_spatial_marker_subset.h5ad'))]:
    rules=json.loads((ROOT/'annotation_rules.json').read_text())[assay]
    for c,r in rules.items():check(assay+' annotation '+c,ob.obs.loc[ob.obs.expression_cluster.eq(c),'RNA_annotation'].eq(c+': '+r[0]+' | '+r[1]).all())
b=ob;pref=WORK/'inputs'/'GSM6177612_NYU_LIHC1_Vis_processed'
with h5py.File(str(pref)+'_filtered_feature_bc_matrix.h5') as h:
    m=h['matrix'];x=csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:])).T.tocsr();ids=np.array([s.decode() for s in m['features']['id'][:]]);names=np.array([s.decode() for s in m['features']['name'][:]]);barcodes=[s.decode() for s in m['barcodes'][:]]
pos=pd.read_csv(str(pref)+'_spatial_tissue_positions_list.csv.gz',header=None,index_col=0);pos.columns=['in_tissue','array_row','array_col','pixel_row','pixel_col'];pos=pos.loc[barcodes]
total=np.asarray(x.sum(1)).ravel();ng=x.getnnz(axis=1);mt=np.asarray(x[:,np.char.startswith(names,'MT-')].sum(1)).ravel()/total*100;ribo=np.asarray(x[:,pd.Series(names).str.match(r'^RP[SL]\d').to_numpy()].sum(1)).ravel()/total*100
keep=(pos.in_tissue.to_numpy()==1)&(total>=500)&(ng>=200)&(mt<=25)&(ribo<=30)
check('spatial QC IDs',list(np.asarray(barcodes)[keep])==list(b.obs_names));check('spatial original coordinates',np.array_equal(b.obsm['spatial'],pos.loc[b.obs_names,['pixel_col','pixel_row']].to_numpy()));expression_checks(b,x[keep],ids,names,'spatial')
array=b.obs[['array_row','array_col']].to_numpy(float);array[:,0]*=np.sqrt(3);ij=cKDTree(array).query_pairs(2.000001,output_type='ndarray');expected=set(map(tuple,ij))|set(map(tuple,ij[:,::-1]));G=b.obsp['spatial_connectivities'];coo=G.tocoo();actual=set(zip(coo.row,coo.col));check('independent six-neighbor lattice',actual==expected);check('no self and max 6',G.diagonal().sum()==0 and np.asarray(G.sum(1)).max()<=6)
sf=json.loads(gzip.decompress(Path(str(pref)+'_spatial_scalefactors_json.json.gz').read_bytes()));im=np.asarray(Image.open(io.BytesIO(gzip.decompress(Path(str(pref)+'_spatial_tissue_lowres_image.png.gz').read_bytes()))).convert('RGB'));check('original H&E',np.array_equal(im,b.uns['spatial']['NYU_LIHC1_Vis']['images']['lowres']));check('original scales',all(sf[k]==v for k,v in b.uns['scalefactors'].items()));xy=b.obsm['spatial']*sf['tissue_lowres_scalef'];check('positions inside image',np.all((xy[:,0]>=0)&(xy[:,0]<im.shape[1])&(xy[:,1]>=0)&(xy[:,1]<im.shape[0])))
stats=pd.read_csv(ROOT/'PGAM5_spatial_autocorrelation.csv')
for method in ['unrestricted','UMI_MT_stratified']:
    v=np.load(ROOT/f'NYU_LIHC1_Vis_{method}_Moran_audit.npz');z=v['z'];g=v['groups'];scale=len(b)/(G.sum()*(z*z).sum());I=scale*sum(z[i]*z[j] for i,j in actual);check(method+' Moran',np.isclose(I,v['observed'],atol=1e-12))
    means={k:z[g==k].mean() for k in np.unique(g)};var={k:z[g==k].var() for k in means};n={k:(g==k).sum() for k in means};center=scale*sum(means[g[i]]*means[g[j]]-(var[g[i]]/(n[g[i]]-1) if g[i]==g[j] and n[g[i]]>1 else 0) for i,j in actual)
    check(method+' exact center',np.isclose(center,v['expected'],atol=1e-12));p=(1+np.sum(np.abs(v['null_I']-center)>=abs(I-center)-1e-12))/5000;check(method+' centered two-sided p',np.isclose(p,stats.loc[stats.method.eq(method),'permutation_p'].iloc[0]))
    for perm in v['example_permutations']:check(method+' permutation strata',np.array_equal(np.sort(perm),np.arange(len(b))) and np.array_equal(g[perm],g))
check('Moran BH2',np.allclose(false_discovery_control(stats.permutation_p.to_numpy()),stats.BH_q_2_tests))
edges=pd.read_csv(ROOT/'NYU_LIHC1_Vis_mac_T_candidate_edges.csv.gz');index=pd.Index(b.obs_names);s=index.get_indexer(edges.source_barcode);t=index.get_indexer(edges.target_barcode)
for pair in pd.read_csv(ROOT/'selected_LR_pairs.csv').itertuples():
    ligand=np.expm1(b[:,[pair.ligand]].X.toarray().ravel())[s];receptor=np.min(np.vstack([np.expm1(b[:,[g]].X.toarray().ravel())[t] for g in pair.receptor_complex.split('_')]),axis=0);check('LR all subunits '+pair.ligand+'--'+pair.receptor_complex,np.allclose(np.sqrt(ligand*receptor),edges[pair.ligand+'--'+pair.receptor_complex],rtol=1e-5,atol=1e-6))
subprocess.run([sys.executable,str(ROOT/'association_audit.py')],check=True)
assoc=json.loads((ROOT/'association_audit_results.json').read_text());checks.extend(assoc['checks']);(ROOT/'audit_results.json').write_text(json.dumps(dict(passed=True,checks_passed=len(checks),checks=checks),indent=2)+'\n');print('Passed',len(checks),'checks')
