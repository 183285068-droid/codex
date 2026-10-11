"""Expression-only analysis. Spatial coordinates are unavailable, never inferred."""
from pathlib import Path
import gzip
import json
import os

os.environ.setdefault('NUMBA_CACHE_DIR', '/tmp/GSE283406_numba')
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent
WORK = Path(os.environ.get('GSE283406_WORK', ROOT/'work'))
INPUT = WORK/'inputs'
PROGRAMS = json.loads((ROOT/'programs.json').read_text())
SAMPLES = [('GSM8661789','J259B','NextSeq'), ('GSM8661790','J259X','NextSeq'),
           ('GSM8661791','J259B','AVITI'), ('GSM8661792','J259X','AVITI')]
SEED = 283406

def dense(x):
    return np.asarray(x.toarray() if sp.issparse(x) else x)

def load(gsm, section, platform):
    prefix = INPUT/f'{gsm}_{section}'
    features = pd.read_csv(str(prefix)+'.features.tsv.gz', sep='\t', header=None)
    barcodes = pd.read_csv(str(prefix)+'.barcodes.tsv.gz', header=None)[0].astype(str)
    with gzip.open(str(prefix)+'.matrix.mtx.gz', 'rb') as handle:
        counts = scipy.io.mmread(handle, spmatrix=True).T.tocsr()
    counts.eliminate_zeros()
    assert counts.shape == (len(barcodes), len(features))
    assert barcodes.is_unique and features[0].is_unique
    assert features[1].eq('PGAM5').sum() == 1
    assert np.all(counts.data >= 0) and np.all(counts.data == np.floor(counts.data))
    a = ad.AnnData(counts.astype(np.float32), obs=pd.DataFrame(index=pd.Index(barcodes.values, name='barcode')),
                   var=pd.DataFrame({'gene_id':features[0].values, 'original_gene_symbol':features[1].values}, index=pd.Index(features[1].values, name='gene')))
    # Three duplicated symbols map to distinct Ensembl IDs; retain both features and counts.
    a.var_names_make_unique()
    a.obs_names.name = 'barcode'
    a.obs['gsm'] = gsm
    a.obs['section'] = section
    a.obs['platform'] = platform
    a.obs['tissue_scope'] = 'HCC_tumor_section' if section == 'J259X' else 'HCC_and_adjacent_mixed_section'
    a.obs['total_counts'] = np.asarray(counts.sum(1)).ravel()
    a.obs['n_genes'] = counts.getnnz(axis=1)
    a.obs['pct_MT_targeted_panel'] = dense(a.X[:,a.var_names.str.startswith('MT-')].sum(1)).ravel()/a.obs.total_counts*100
    a.obs['PGAM5_raw'] = dense(a[:,['PGAM5']].X).ravel().astype(int)
    # FFPE probe panel: mitochondrial proportion is reported, not used as a poly(A) QC cutoff.
    a.obs['pass_qc'] = a.obs.total_counts.ge(500) & a.obs.n_genes.ge(200)
    a.obs['pass_qc_relaxed'] = a.obs.total_counts.ge(250) & a.obs.n_genes.ge(100)
    return a

def score(a, coverage):
    a.obs['PGAM5_log1p_CP10k'] = dense(a[:,['PGAM5']].X).ravel()
    for program, genes in PROGRAMS.items():
        present = [g for g in genes if g in a.var_names]
        coverage.append({'gsm':a.obs.gsm.iloc[0], 'program':program, 'requested':len(genes),
                         'present':len(present), 'present_genes':';'.join(present),
                         'missing_genes':';'.join(g for g in genes if g not in a.var_names)})
        assert present, program
        a.obs[program] = dense(a[:,present].X).mean(1)
        a.obs[program+'_n_detected'] = (dense(a[:,present].layers['counts']) > 0).sum(1)
    a.obs['Mac_RNA_enriched'] = (a.obs.Macrophage_core.ge(a.obs.Macrophage_core.quantile(.75)) &
        a.obs.Macrophage_core_n_detected.ge(3) & a.obs.Pan_myeloid.ge(a.obs.Pan_myeloid.median()))
    a.obs['T_RNA_enriched'] = a.obs.T_lineage.ge(a.obs.T_lineage.quantile(.75)) & a.obs.T_lineage_n_detected.ge(2)
    for k in ['Activation_early','Exhaustion_core','Effector_cytotoxic']:
        gate = a.obs.T_RNA_enriched
        cutoff = a.obs.loc[gate,k].quantile(.75)
        a.obs[k+'_supported'] = gate & a.obs[k].ge(cutoff) & a.obs[k].gt(0) & a.obs[k+'_n_detected'].ge(2)
    gate = a.obs.Mac_RNA_enriched
    positive = gate & a.obs.PGAM5_raw.gt(0)
    threshold = float(a.obs.loc[positive,'PGAM5_log1p_CP10k'].quantile(.75))
    a.obs['PGAM5_mac_RNA_group'] = 'Outside_mac_RNA_gate'
    a.obs.loc[gate, 'PGAM5_mac_RNA_group'] = 'Undetected'
    a.obs.loc[positive, 'PGAM5_mac_RNA_group'] = 'Other_detected'
    a.obs.loc[positive & a.obs.PGAM5_log1p_CP10k.ge(threshold),'PGAM5_mac_RNA_group'] = 'PGAM5_high'
    a.uns['PGAM5_high_threshold_log1p_CP10k'] = threshold

if __name__ == '__main__':
    cache, coverage, summary, group_rows, correlations = {}, [], [], [], []
    for gsm, section, platform in SAMPLES:
        original = load(gsm,section,platform)
        original.obs.to_csv(ROOT/f'{gsm}_all_spot_QC.csv')
        a = original[original.obs.pass_qc].copy()
        a.layers['counts'] = a.X.copy()
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        score(a,coverage)
        cache[(section,platform)] = a
        row = dict(gsm=gsm,section=section,platform=platform,input_spots=len(original),qc_spots=len(a),
                   relaxed_qc_spots=int(original.obs.pass_qc_relaxed.sum()),genes_in_panel=a.n_vars,
                   median_UMI=a.obs.total_counts.median(),median_genes=a.obs.n_genes.median(),
                   PGAM5_detected=int(a.obs.PGAM5_raw.gt(0).sum()),PGAM5_detection_pct=100*a.obs.PGAM5_raw.gt(0).mean(),
                   PGAM5_mean_log1p_CP10k=a.obs.PGAM5_log1p_CP10k.mean(),
                   Mac_RNA_enriched=int(a.obs.Mac_RNA_enriched.sum()),T_RNA_enriched=int(a.obs.T_RNA_enriched.sum()),
                   T_activation_supported=int(a.obs.Activation_early_supported.sum()),
                   T_exhaustion_supported=int(a.obs.Exhaustion_core_supported.sum()),
                   PGAM5_high_threshold=a.uns['PGAM5_high_threshold_log1p_CP10k'])
        summary.append(row)
        for name,mask in [('All_QC_spots',np.ones(len(a),dtype=bool)),('Mac_RNA_enriched',a.obs.Mac_RNA_enriched),
                          ('T_RNA_enriched',a.obs.T_RNA_enriched),('Activation_supported',a.obs.Activation_early_supported),
                          ('Exhaustion_supported',a.obs.Exhaustion_core_supported)]:
            d=a.obs.loc[mask]
            group_rows.append(dict(gsm=gsm,section=section,platform=platform,group=name,n_spots=len(d),
                PGAM5_detected=int(d.PGAM5_raw.gt(0).sum()),PGAM5_detection_pct=100*d.PGAM5_raw.gt(0).mean(),
                PGAM5_mean_log1p_CP10k=d.PGAM5_log1p_CP10k.mean(),PGAM5_median_log1p_CP10k=d.PGAM5_log1p_CP10k.median()))
        print('QC',row,flush=True)
        if platform == 'AVITI':
            h=a.copy()
            sc.pp.filter_genes(h,min_cells=10)
            h=h[:,~(h.var_names.str.startswith('MT-')|h.var_names.str.match(r'^RP[SL]\d'))].copy()
            sc.pp.highly_variable_genes(h,n_top_genes=2000,flavor='seurat')
            h=h[:,h.var.highly_variable].copy()
            h.var.to_csv(ROOT/f'{section}_AVITI_HVGs.csv',index_label='gene')
            sc.pp.scale(h,max_value=10)
            sc.tl.pca(h,n_comps=40,random_state=SEED)
            sc.pp.neighbors(h,n_neighbors=15,n_pcs=30,random_state=SEED)
            sc.tl.leiden(h,resolution=.6,key_added='cluster',random_state=SEED,flavor='igraph',n_iterations=2,directed=False)
            sc.tl.umap(h,random_state=SEED)
            a.obs['expression_cluster'] = pd.Categorical(['E'+str(x) for x in h.obs.cluster])
            a.obsm['X_umap'] = h.obsm['X_umap']
            a.obsm['X_pca'] = h.obsm['X_pca']
            sc.tl.rank_genes_groups(a,groupby='expression_cluster',method='wilcoxon',tie_correct=True,use_raw=False,n_genes=100,pts=True)
            markers=sc.get.rank_genes_groups_df(a,group=None)
            markers.to_csv(ROOT/f'{section}_AVITI_top100_markers.csv',index=False)
            profiles=[]
            for cluster in a.obs.expression_cluster.cat.categories:
                d=a.obs[a.obs.expression_cluster.eq(cluster)]
                profiles.append(dict(section=section,expression_cluster=cluster,n_spots=len(d),
                    PGAM5_detection_pct=100*d.PGAM5_raw.gt(0).mean(),PGAM5_mean_log1p_CP10k=d.PGAM5_log1p_CP10k.mean(),
                    **{k:d[k].mean() for k in PROGRAMS}))
            pd.DataFrame(profiles).to_csv(ROOT/f'{section}_AVITI_cluster_profiles.csv',index=False)
            # Portable marker subset preserves full-panel normalization and original counts.
            extras=['PGAM5','MYO19','SPP1','MMP12','LPL','STAB1','MARCO','CD163','AFP','GPC3','GLUL','ALDOB','RBP4','FGA','FGB','FGG','MALAT1','HMGB1','MTRNR2L8','CXCL8','CCL2','CCL8','CD274','PVR','NECTIN2','CD44','LGALS9','HAVCR2','TGFB1']
            top=markers.groupby('group',observed=True).head(12)['names'].tolist()
            panel=list(dict.fromkeys(extras+[g for genes in PROGRAMS.values() for g in genes]+top))
            portable=a[:,[g for g in panel if g in a.var_names]].copy()
            portable.uns.pop('rank_genes_groups',None)
            portable.uns['coordinates_available'] = False
            portable.uns['analysis_kind'] = 'expression_only_not_spatial_map'
            portable.write_h5ad(ROOT/f'{section}_AVITI_expression_marker_subset.h5ad',compression='gzip')
            a.write_h5ad(WORK/f'{section}_AVITI_full_processed.h5ad',compression='gzip')
            row['expression_clusters']=a.obs.expression_cluster.nunique()
        a.obs.to_csv(ROOT/f'{gsm}_QC_expression_annotations.csv',float_format='%.17g')

    technical=[]
    for section in ['J259B','J259X']:
        nextseq,aviti=cache[(section,'NextSeq')],cache[(section,'AVITI')]
        shared=nextseq.obs_names.intersection(aviti.obs_names)
        nx=nextseq.obs.loc[shared];av=aviti.obs.loc[shared]
        paired=pd.DataFrame({'barcode':shared, 'section':section,
            'NextSeq_total_UMI':nx.total_counts.values,'AVITI_total_UMI':av.total_counts.values,
            'NextSeq_PGAM5_raw':nx.PGAM5_raw.values,'AVITI_PGAM5_raw':av.PGAM5_raw.values,
            'NextSeq_PGAM5_log1p_CP10k':nx.PGAM5_log1p_CP10k.values,
            'AVITI_PGAM5_log1p_CP10k':av.PGAM5_log1p_CP10k.values})
        paired.to_csv(ROOT/f'{section}_paired_platform_expression.csv',index=False)
        npos=nx.PGAM5_raw.gt(0).to_numpy();apos=av.PGAM5_raw.gt(0).to_numpy()
        technical.append(dict(section=section,shared_QC_spots=len(shared),
            PGAM5_normalized_spearman=spearmanr(nx.PGAM5_log1p_CP10k,av.PGAM5_log1p_CP10k).statistic,
            total_UMI_spearman=spearmanr(nx.total_counts,av.total_counts).statistic,
            PGAM5_both_detected=int((npos & apos).sum()),PGAM5_only_AVITI=int((~npos & apos).sum()),
            PGAM5_only_NextSeq=int((npos & ~apos).sum()),PGAM5_neither_detected=int((~npos & ~apos).sum())))
    pd.DataFrame(summary).to_csv(ROOT/'sample_QC_summary.csv',index=False)
    pd.DataFrame(coverage).to_csv(ROOT/'program_gene_coverage.csv',index=False)
    pd.DataFrame(group_rows).to_csv(ROOT/'PGAM5_lineage_program_summary.csv',index=False)
    pd.DataFrame(technical).to_csv(ROOT/'paired_platform_summary.csv',index=False)
    print('Complete: expression-only analysis; no physical coordinates or adjacency inferred.',flush=True)
