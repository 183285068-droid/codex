"""Within-section expression and genuine-coordinate spatial analysis of HCC Visium spots."""
from pathlib import Path
import gzip
import json
import os
os.environ.setdefault('NUMBA_CACHE_DIR','/tmp/GSE203612_numba')
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp

ROOT=Path(__file__).resolve().parent
WORK=Path(os.environ.get('GSE203612_WORK',ROOT/'work'))
WORK.mkdir(parents=True,exist_ok=True)
SAMPLES=json.loads((ROOT/'sample_scope.json').read_text())
PROGRAMS=json.loads((ROOT/'programs.json').read_text())
EXTRA=['PGAM5','MYO19','SPP1','MMP12','LPL','STAB1','MARCO','CD163','AFP','GPC3','GLUL','ALDOB','RBP4','FGA','FGB','FGG','MALAT1','HMGB1','CXCL8','CCL2','CCL8','CD274','PDCD1LG2','PVR','NECTIN2','CD44','LGALS9','HAVCR2','TGFB1','TGFBR1','TGFBR2','CD80','CD86','CD28','ICOSLG','TNFSF9','TNFSF4','TNFRSF4','IL10','IL10RA','IL10RB','CXCL9','CXCL10','CXCL11','CXCR3','CCL3','CCL4','CCL5','CCR5','INTS6','ANGPTL4','SDC1']
SEED=203612
def dense(x):return np.asarray(x.toarray() if sp.issparse(x) else x)

def load(item):
    sample=item['sample'];pref=WORK/'inputs'/f'{item["gsm"]}_{sample}_processed'
    a=sc.read_10x_h5(str(pref)+'_filtered_feature_bc_matrix.h5')
    a.var['gene_id']=a.var['gene_ids'];a.var['original_gene_symbol']=a.var_names.to_numpy();a.var_names_make_unique()
    assert a.obs_names.is_unique and a.var.gene_id.is_unique and np.all(a.X.data>=0) and np.all(a.X.data==np.floor(a.X.data))
    assert list(a.var.original_gene_symbol).count('PGAM5')==1
    pos=pd.read_csv(str(pref)+'_spatial_tissue_positions_list.csv.gz',header=None,names=['barcode','in_tissue','array_row','array_col','pixel_row','pixel_col']).set_index('barcode')
    assert pos.index.is_unique and set(a.obs_names)<=set(pos.index)
    a.obs=pos.loc[a.obs_names].copy()
    for key,value in item.items():a.obs[key]=value
    a.obsm['spatial']=a.obs[['pixel_col','pixel_row']].to_numpy()
    a.uns['scalefactors']=json.loads(gzip.decompress(Path(str(pref)+'_spatial_scalefactors_json.json.gz').read_bytes()))
    a.obs['total_counts']=np.asarray(a.X.astype(np.int64).sum(1)).ravel();a.obs['n_genes']=a.X.getnnz(axis=1)
    a.obs['pct_mt']=dense(a.X[:,a.var_names.str.startswith('MT-')].sum(1)).ravel()/a.obs.total_counts*100
    a.obs['pct_ribo']=dense(a.X[:,a.var_names.str.match(r'^RP[SL]\d')].sum(1)).ravel()/a.obs.total_counts*100
    a.obs['PGAM5_raw']=dense(a[:,['PGAM5']].X).ravel().astype(int)
    base=a.obs.in_tissue.eq(1)&a.obs.n_genes.ge(200)&a.obs.total_counts.ge(500)&a.obs.pct_ribo.le(30)
    a.obs['pass_qc']=base&a.obs.pct_mt.le(25)
    a.obs['pass_qc_author']=a.obs.in_tissue.eq(1)&a.obs.total_counts.ge(500)&a.obs.pct_mt.le(30)&a.obs.pct_ribo.le(30)
    a.obs['pass_qc_mt30']=base&a.obs.pct_mt.le(30)
    a.obs.to_csv(ROOT/f'{sample}_all_spot_QC.csv',float_format='%.17g')
    return a

def graph(a):
    array=a.obs[['array_row','array_col']].to_numpy()
    lookup={tuple(v):i for i,v in enumerate(array)}
    edges=[]
    for i,(row,col) in enumerate(array):
        for dr,dc in [(0,2),(0,-2),(1,1),(1,-1),(-1,1),(-1,-1)]:
            j=lookup.get((row+dr,col+dc))
            if j is not None:edges.append((i,j))
    ij=np.asarray(edges,dtype=int)
    g=sp.csr_matrix((np.ones(len(edges)),(ij[:,0],ij[:,1])),shape=(len(a),len(a))) if len(edges) else sp.csr_matrix((len(a),len(a)))
    a.obsp['spatial_connectivities']=g
    a.obs['spatial_degree']=np.asarray(g.sum(1)).ravel().astype(int)
    pd.DataFrame({'source_barcode':a.obs_names[ij[:,0]],'target_barcode':a.obs_names[ij[:,1]]}).to_csv(ROOT/f'{a.obs["sample"].iloc[0]}_directed_spatial_neighbors.csv.gz',index=False)
    return g

def run(item):
    sample=item['sample'];original=load(item);a=original[original.obs.pass_qc].copy()
    assert len(a)>100,(sample,len(a))
    a.layers['counts']=a.X.copy()
    sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
    a.obs['PGAM5_log1p_CP10k']=dense(a[:,['PGAM5']].X).ravel()
    coverage=[]
    for program,genes in PROGRAMS.items():
        present=[g for g in genes if g in a.var_names]
        assert present,program
        coverage.append(dict(sample=sample,program=program,requested=len(genes),present=len(present),
                             present_genes=';'.join(present),missing_genes=';'.join(g for g in genes if g not in a.var_names)))
        a.obs[program]=dense(a[:,present].X).mean(1)
        a.obs[program+'_n_detected']=(dense(a[:,present].layers['counts'])>0).sum(1)
    g=graph(a)
    a.obs['Mac_RNA_enriched']=(a.obs.Macrophage_core.ge(a.obs.Macrophage_core.quantile(.75))&a.obs.Macrophage_core_n_detected.ge(3)&a.obs.Pan_myeloid.ge(a.obs.Pan_myeloid.median()))
    a.obs['T_RNA_enriched']=a.obs.T_lineage.ge(a.obs.T_lineage.quantile(.75))&a.obs.T_lineage_n_detected.ge(2)
    for program in ['Activation_early','Exhaustion_core','Effector_cytotoxic']:
        t=a.obs.T_RNA_enriched
        cutoff=a.obs.loc[t,program].quantile(.75)
        a.obs[program+'_supported']=t&a.obs[program].ge(cutoff)&a.obs[program].gt(0)&a.obs[program+'_n_detected'].ge(2)
    mac=a.obs.Mac_RNA_enriched;positive=mac&a.obs.PGAM5_raw.gt(0)
    threshold=float(a.obs.loc[positive,'PGAM5_log1p_CP10k'].quantile(.75))
    a.obs['PGAM5_mac_group']='Outside_mac_RNA_gate'
    a.obs.loc[mac,'PGAM5_mac_group']='Undetected'
    a.obs.loc[positive,'PGAM5_mac_group']='Other_detected'
    a.obs.loc[positive&a.obs.PGAM5_log1p_CP10k.ge(threshold),'PGAM5_mac_group']='PGAM5_high'
    a.uns['PGAM5_high_threshold_log1p_CP10k']=threshold
    a.obs['T_neighbor_count']=np.asarray(g@a.obs.T_RNA_enriched.to_numpy(dtype=int)).ravel().astype(int)
    a.obs['T_exhaustion_neighbor_count']=np.asarray(g@a.obs.Exhaustion_core_supported.to_numpy(dtype=int)).ravel().astype(int)
    a.obs['T_activation_neighbor_count']=np.asarray(g@a.obs.Activation_early_supported.to_numpy(dtype=int)).ravel().astype(int)
    a.obs['T_neighbor_fraction']=a.obs.T_neighbor_count/a.obs.spatial_degree.replace(0,np.nan)
    for program in ['Activation_early','Exhaustion_core','Effector_cytotoxic']:
        v=a.obs[program].to_numpy()*a.obs.T_RNA_enriched.to_numpy()
        a.obs['neighbor_T_'+program]=np.asarray(g@v).ravel()/a.obs.T_neighbor_count.replace(0,np.nan)
    h=a.copy();sc.pp.filter_genes(h,min_cells=10)
    h=h[:,~(h.var_names.str.startswith('MT-')|h.var_names.str.match(r'^RP[SL]\d'))].copy()
    sc.pp.highly_variable_genes(h,n_top_genes=2000,flavor='seurat');h=h[:,h.var.highly_variable].copy()
    h.var.to_csv(ROOT/f'{sample}_HVGs.csv',index_label='gene')
    sc.pp.scale(h,max_value=10);sc.tl.pca(h,n_comps=40,random_state=SEED)
    sc.pp.neighbors(h,n_neighbors=15,n_pcs=30,random_state=SEED)
    sc.tl.leiden(h,resolution=.6,key_added='cluster',random_state=SEED,flavor='igraph',n_iterations=2,directed=False)
    sc.tl.umap(h,random_state=SEED)
    a.obs['expression_cluster']=pd.Categorical(['S'+str(x) for x in h.obs.cluster])
    a.obsm['X_umap']=h.obsm['X_umap'];a.obsm['X_pca']=h.obsm['X_pca']
    sc.tl.rank_genes_groups(a,groupby='expression_cluster',method='wilcoxon',tie_correct=True,use_raw=False,n_genes=100,pts=True)
    markers=sc.get.rank_genes_groups_df(a,group=None);markers.to_csv(ROOT/f'{sample}_top100_markers.csv',index=False)
    rows=[]
    for cluster in a.obs.expression_cluster.cat.categories:
        d=a.obs[a.obs.expression_cluster.eq(cluster)]
        rows.append(dict(sample=sample,patient=item['patient'],tumor_status=item['tumor_status'],cluster=cluster,n_spots=len(d),
            PGAM5_detected=int(d.PGAM5_raw.gt(0).sum()),PGAM5_detection_pct=100*d.PGAM5_raw.gt(0).mean(),
            PGAM5_mean_log1p_CP10k=d.PGAM5_log1p_CP10k.mean(),**{k:d[k].mean() for k in PROGRAMS}))
    pd.DataFrame(rows).to_csv(ROOT/f'{sample}_cluster_profiles.csv',index=False)
    a.obs.to_csv(ROOT/f'{sample}_spot_annotations.csv',float_format='%.17g')
    a.write_h5ad(WORK/f'{sample}_full_processed.h5ad',compression='gzip')
    panel=list(dict.fromkeys(EXTRA+[g for values in PROGRAMS.values() for g in values]+markers.groupby('group',observed=True).head(8).names.tolist()))
    portable=a[:,[g for g in panel if g in a.var_names]].copy();portable.uns.pop('rank_genes_groups',None)
    portable.write_h5ad(ROOT/f'{sample}_spatial_marker_subset.h5ad',compression='gzip')
    summary=dict(**item,input_spots=len(original),qc_spots=len(a),qc_author_rule_spots=int(original.obs.pass_qc_author.sum()),qc_mt30_spots=int(original.obs.pass_qc_mt30.sum()),
        genes_in_matrix=a.n_vars,expression_clusters=a.obs.expression_cluster.nunique(),
        median_UMI=a.obs.total_counts.median(),median_genes=a.obs.n_genes.median(),median_pct_mt=a.obs.pct_mt.median(),
        PGAM5_detected=int(a.obs.PGAM5_raw.gt(0).sum()),PGAM5_detection_pct=100*a.obs.PGAM5_raw.gt(0).mean(),
        PGAM5_mean_log1p_CP10k=a.obs.PGAM5_log1p_CP10k.mean(),Mac_RNA_enriched=int(mac.sum()),
        PGAM5_positive_mac_RNA_spots=int(positive.sum()),PGAM5_high_mac_RNA_spots=int(a.obs.PGAM5_mac_group.eq('PGAM5_high').sum()),
        T_RNA_enriched=int(a.obs.T_RNA_enriched.sum()),T_activation_supported=int(a.obs.Activation_early_supported.sum()),
        T_exhaustion_supported=int(a.obs.Exhaustion_core_supported.sum()),mac_sources_with_T_neighbors=int((mac&a.obs.T_neighbor_count.gt(0)).sum()),
        mac_sources_with_Tex_neighbors=int((mac&a.obs.T_exhaustion_neighbor_count.gt(0)).sum()),
        high_mac_sources_with_Tex_neighbors=int((a.obs.PGAM5_mac_group.eq('PGAM5_high')&a.obs.T_exhaustion_neighbor_count.gt(0)).sum()))
    print('Completed',summary,flush=True)
    return summary,coverage

if __name__=='__main__':
    from concurrent.futures import ProcessPoolExecutor
    rows=[];coverage=[]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for summary,c in pool.map(run,SAMPLES):rows.append(summary);coverage.extend(c)
    pd.DataFrame(rows).to_csv(ROOT/'sample_QC_summary.csv',index=False)
    pd.DataFrame(coverage).to_csv(ROOT/'program_gene_coverage.csv',index=False)
