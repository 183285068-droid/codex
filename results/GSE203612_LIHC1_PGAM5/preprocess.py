"""Matched single-cell and Visium preprocessing. Original matrix-header cell identifiers retained."""
from pathlib import Path
import gzip,io,json
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,scipy.io,scipy.sparse as sp
from PIL import Image
from spatial_core import ROOT,WORK,SAMPLES,PROGRAMS,EXTRA,SEED,dense
import spatial_core

def annotate_programs(a):
    coverage=[]
    for p,genes in PROGRAMS.items():
        present=[g for g in genes if g in a.var_names];coverage.append(dict(assay='inDrop',program=p,requested=len(genes),present=len(present),present_genes=';'.join(present),missing_genes=';'.join(g for g in genes if g not in a.var_names)))
        assert present
        a.obs[p]=dense(a[:,present].X).mean(1);a.obs[p+'_n_detected']=(dense(a[:,present].layers['counts'])>0).sum(1)
    return coverage

def singlecell():
    arrays=[];qc=[]
    for lib in range(1,4):
        pref=WORK/'inputs'/f'GSM6177611_NYU_LIHC1_lib{lib}_lib{lib}'
        features=pd.read_csv(str(pref)+'_genes.tsv.gz',sep='\t');comments=[]
        with gzip.open(str(pref)+'_gene_expression.mtx.gz','rt') as h:
            assert h.readline().startswith('%%MatrixMarket')
            for line in h:
                if line.startswith('%'):
                    values=line.rstrip()[1:].split('\t')
                    if len(values)>1:comments.append(values)
                else:break
        with gzip.open(str(pref)+'_gene_expression.mtx.gz','rb') as h:x=scipy.io.mmread(h,spmatrix=True).T.tocsr()
        x.eliminate_zeros();assert len(comments)==2 and comments[0]==features.name.tolist() and x.shape==(len(comments[1]),len(features))
        assert np.all(x.data>=0) and np.all(x.data==np.floor(x.data))
        cells=comments[1];assert len(set(cells))==len(cells) and features.ensembl_id.is_unique
        a=ad.AnnData(x.astype(np.float32),obs=pd.DataFrame({'library':f'lib{lib}','original_cell_id':cells,'patient':'NYU_LIHC1'},index=pd.Index([f'lib{lib}|{c}' for c in cells],name='cell_id')),
            var=pd.DataFrame({'gene_id':features.ensembl_id.to_numpy(),'original_gene_symbol':features.name.to_numpy()},index=pd.Index(features.name.to_numpy(),name='gene')))
        a.var_names_make_unique();a.obs['total_counts']=np.asarray(x.sum(1)).ravel();a.obs['n_genes']=x.getnnz(axis=1)
        a.obs['pct_mt']=dense(a[:,a.var_names.str.startswith('MT-')].X.sum(1)).ravel()/a.obs.total_counts*100
        a.obs['pct_ribo']=dense(a[:,a.var_names.str.match(r'^RP[SL]\d')].X.sum(1)).ravel()/a.obs.total_counts*100
        a.obs['PGAM5_raw']=dense(a[:,['PGAM5']].X).ravel().astype(int)
        a.obs['pass_qc']=a.obs.total_counts.ge(500)&a.obs.n_genes.ge(200)&a.obs.pct_mt.le(25)&a.obs.pct_ribo.le(30)
        a.obs['pass_qc_author']=a.obs.total_counts.ge(500)&a.obs.pct_mt.le(30)&a.obs.pct_ribo.le(30)
        qc.append(dict(library=f'lib{lib}',input_cells=len(a),qc_cells=int(a.obs.pass_qc.sum()),qc_author_rule_cells=int(a.obs.pass_qc_author.sum()),median_input_UMI=a.obs.total_counts.median()))
        a.obs.to_csv(ROOT/f'lib{lib}_all_cell_QC.csv',float_format='%.17g');arrays.append(a[a.obs.pass_qc].copy())
    a=ad.concat(arrays,join='inner',merge='same');assert a.obs_names.is_unique
    a.layers['counts']=a.X.copy();sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a);a.obs['PGAM5_log1p_CP10k']=dense(a[:,['PGAM5']].X).ravel()
    coverage=annotate_programs(a)
    # Conservative RNA-supported candidate labels, not validated cell identity.
    det=lambda gs:(dense(a[:,[g for g in gs if g in a.var_names]].layers['counts'])>0).sum(1)
    a.obs['T_RNA_supported']=a.obs.T_lineage_n_detected.ge(2)
    a.obs['DC_RNA_supported']=det(['CD1C','FCER1A','CLEC10A','CLEC9A','XCR1'])>=2
    a.obs['Mac_RNA_supported']=a.obs.Macrophage_core_n_detected.ge(2)&a.obs.Pan_myeloid_n_detected.ge(2)&~a.obs.T_RNA_supported&a.obs.B_cell_n_detected.lt(2)&~a.obs.DC_RNA_supported&a.obs.Hepatocyte_n_detected.lt(3)
    positive=a.obs.Mac_RNA_supported&a.obs.PGAM5_raw.gt(0)
    a.obs['PGAM5_mac_group']='Outside_mac_RNA_gate';a.obs.loc[a.obs.Mac_RNA_supported,'PGAM5_mac_group']='Undetected';a.obs.loc[positive,'PGAM5_mac_group']='Other_detected'
    if positive.any():a.obs.loc[positive&a.obs.PGAM5_log1p_CP10k.ge(a.obs.loc[positive,'PGAM5_log1p_CP10k'].quantile(.75)),'PGAM5_mac_group']='PGAM5_high'
    h=a.copy();sc.pp.filter_genes(h,min_cells=5);h=h[:,~(h.var_names.str.startswith('MT-')|h.var_names.str.match(r'^RP[SL]\d'))].copy()
    sc.pp.highly_variable_genes(h,n_top_genes=2000,flavor='seurat');h=h[:,h.var.highly_variable].copy();h.var.to_csv(ROOT/'singlecell_HVGs.csv',index_label='gene')
    sc.pp.scale(h,max_value=10);sc.tl.pca(h,n_comps=40,random_state=SEED);sc.pp.neighbors(h,n_neighbors=15,n_pcs=30,random_state=SEED)
    sc.tl.leiden(h,resolution=.6,key_added='cluster',random_state=SEED,flavor='igraph',n_iterations=2,directed=False);sc.tl.umap(h,random_state=SEED)
    a.obs['expression_cluster']=pd.Categorical(['C'+str(x) for x in h.obs.cluster]);a.obsm['X_umap']=h.obsm['X_umap'];a.obsm['X_pca']=h.obsm['X_pca']
    sc.tl.rank_genes_groups(a,groupby='expression_cluster',method='wilcoxon',tie_correct=True,use_raw=False,n_genes=100,pts=True)
    markers=sc.get.rank_genes_groups_df(a,group=None);markers.to_csv(ROOT/'singlecell_top100_markers.csv',index=False)
    profiles=[]
    for c in a.obs.expression_cluster.cat.categories:
        d=a.obs[a.obs.expression_cluster.eq(c)];profiles.append(dict(cluster=c,n_cells=len(d),PGAM5_detected=int(d.PGAM5_raw.gt(0).sum()),PGAM5_detection_pct=100*d.PGAM5_raw.gt(0).mean(),PGAM5_mean_log1p_CP10k=d.PGAM5_log1p_CP10k.mean(),median_UMI=d.total_counts.median(),Mac_RNA_supported=int(d.Mac_RNA_supported.sum()),T_RNA_supported=int(d.T_RNA_supported.sum()),**{p:d[p].mean() for p in PROGRAMS}))
    pd.DataFrame(profiles).to_csv(ROOT/'singlecell_cluster_profiles.csv',index=False)
    a.obs.to_csv(ROOT/'singlecell_annotations.csv',float_format='%.17g');a.write_h5ad(WORK/'singlecell_full_processed.h5ad',compression='gzip')
    panel=list(dict.fromkeys(EXTRA+[g for gs in PROGRAMS.values() for g in gs]+markers.groupby('group',observed=True).head(12).names.tolist()))
    portable=a[:,[g for g in panel if g in a.var_names]].copy();portable.uns.pop('rank_genes_groups',None);portable.write_h5ad(ROOT/'singlecell_marker_subset.h5ad',compression='gzip')
    pd.DataFrame(qc).to_csv(ROOT/'singlecell_QC_summary.csv',index=False);pd.DataFrame(coverage).to_csv(ROOT/'singlecell_program_gene_coverage.csv',index=False)
    print('scRNA:',len(a),'PGAM5 detected',int(a.obs.PGAM5_raw.gt(0).sum()),'Mac candidates',int(a.obs.Mac_RNA_supported.sum()),'PGAM5+ Mac candidates',int(positive.sum()),'T candidates',int(a.obs.T_RNA_supported.sum()),'clusters',a.obs.expression_cluster.nunique(),flush=True)
if __name__=='__main__':
    summary,coverage=spatial_core.run(SAMPLES[0]);pd.DataFrame([summary]).to_csv(ROOT/'sample_QC_summary.csv',index=False);pd.DataFrame(coverage).to_csv(ROOT/'program_gene_coverage.csv',index=False)
    s=SAMPLES[0];a=ad.read_h5ad(ROOT/f'{s["sample"]}_spatial_marker_subset.h5ad');pref=WORK/'inputs'/f'{s["gsm"]}_{s["sample"]}_processed'
    im=np.asarray(Image.open(io.BytesIO(gzip.decompress(Path(str(pref)+'_spatial_tissue_lowres_image.png.gz').read_bytes()))).convert('RGB'))
    a.uns['spatial']={s['sample']:{'images':{'lowres':im},'scalefactors':dict(a.uns['scalefactors']),'metadata':{'source':'GSE203612'}}};a.write_h5ad(ROOT/f'{s["sample"]}_spatial_marker_subset.h5ad',compression='gzip')
    Image.fromarray(im).save(ROOT/'NYU_LIHC1_Vis_tissue_lowres_image.png');singlecell()
