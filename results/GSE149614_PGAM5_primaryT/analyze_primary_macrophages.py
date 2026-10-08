from pathlib import Path
import json
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
R=Path('/workspace/scratch/GSE149614')
markers=pd.read_csv(R/'original_myeloid_cluster_markers.csv');frac=markers.pivot(index='cluster',columns='gene',values='fraction')
# Conservative mature macrophage annotation, independent of PGAM5 grouping.
keep=frac.index[frac[['C1QA','C1QB','C1QC']].min(axis=1).ge(.8)&frac.CD68.ge(.75)].tolist()
expected=[5,6,21,23,38,39,41,44,46];assert keep==expected
labels={10:'DC-enriched/mixed antigen-presenting myeloid',16:'FCN1/S100A8/A9 monocyte-enriched',26:'FCN1/S100A8/A9 monocyte-enriched',52:'GZMB/JCHAIN pDC-enriched',53:'TPSAB1/TPSB2 mast-cell-enriched'}
annotation=pd.DataFrame({'cluster':frac.index,'included_macrophage_cluster':frac.index.isin(keep),'annotation':[('C1Q/CD68 macrophage-enriched' if c in keep else labels.get(c,'excluded myeloid')) for c in frac.index]})
annotation.to_csv(R/'myeloid_cluster_annotation.csv',index=False)
a=ad.read_h5ad(R/'all_sites_myeloid_counts.h5ad')
selected=a.obs['sample'].str.endswith('T').to_numpy()&a.obs['site'].eq('Tumor').to_numpy()&a.obs['res.3'].isin(keep).to_numpy()
a=a[selected].copy();raw_n=len(a);raw_pos=int(a.obs.PGAM5_detected.sum())
a.obs.to_csv(R/'primary_T_macrophage_metadata_before_QC.csv.gz')
a=a[a.obs.n_genes.ge(500).to_numpy()&a.obs.pct_mt.lt(20).to_numpy()].copy()
a.write_h5ad(R/'primary_T_macrophage_counts_QC.h5ad',compression='gzip')
pos=a.obs.PGAM5_counts.to_numpy()>0
assert pos.any() and (~pos).any()
a.obs['PGAM5_group']=pd.Categorical(np.where(pos,'detected','undetected'))
a.obs.to_csv(R/'primary_T_macrophage_analysis_metadata.csv.gz')
a.obs.groupby('sample',observed=True).agg(macrophages=('Cell','size'),PGAM5_detected=('PGAM5_detected','sum'),PGAM5_detection_rate=('PGAM5_detected','mean'),median_counts=('total_counts','median')).to_csv(R/'primary_T_macrophages_by_sample.csv')
fp=np.asarray((a.X[pos]>0).mean(0)).ravel();fn=np.asarray((a.X[~pos]>0).mean(0)).ravel()
sc.pp.normalize_total(a,target_sum=1e4);mp=np.asarray(a.X[pos].mean(0)).ravel();mn=np.asarray(a.X[~pos].mean(0)).ravel();sc.pp.log1p(a)
sc.tl.rank_genes_groups(a,groupby='PGAM5_group',groups=['detected'],reference='undetected',method='wilcoxon',tie_correct=True,use_raw=False,corr_method='benjamini-hochberg')
r=sc.get.rank_genes_groups_df(a,group='detected').set_index('names').rename(columns={'scores':'Wilcoxon_score','logfoldchanges':'log2FC','pvals':'p_value','pvals_adj':'FDR'})
r['gene']=a.var.loc[r.index,'gene_symbol'].values
aux=pd.DataFrame({'fraction_PGAM5_positive':fp,'fraction_PGAM5_undetected':fn,'mean_10k_detected':mp,'mean_10k_undetected':mn},index=a.var_names)
r=r.join(aux);r['defining_gene']=r.gene.eq('PGAM5');r['direction']=np.where(r.log2FC>0,'up','down');r.to_csv(R/'GSE149614_primaryT_PGAM5_full_DE.csv')
passed=r.FDR.lt(.05)&r.log2FC.abs().ge(1)&r[['fraction_PGAM5_positive','fraction_PGAM5_undetected']].max(axis=1).ge(.1)
passing=r[passed].sort_values(['FDR','log2FC'],ascending=[True,False]);passing.to_csv(R/'GSE149614_primaryT_PGAM5_DE_including_PGAM5.csv')
candidates=passing[~passing.defining_gene];up=candidates[candidates.log2FC.ge(1)];down=candidates[candidates.log2FC.le(-1)].sort_values(['FDR','log2FC'])
up.to_csv(R/'GSE149614_primaryT_PGAM5_upregulated.csv');down.to_csv(R/'GSE149614_primaryT_PGAM5_downregulated.csv');up.head(20).to_csv(R/'GSE149614_primaryT_PGAM5_top20_upregulated.csv')
summary={'primary_T_samples':10,'primary_T_cells':34414,'primary_T_myeloid_cells':8209,'macrophage_clusters':keep,'raw_primary_macrophages':raw_n,'raw_PGAM5_detected_macrophages':raw_pos,'QC_primary_macrophages':len(a),'QC_samples_with_macrophages':int(a.obs['sample'].nunique()),'PGAM5_positive_cells':int(pos.sum()),'PGAM5_undetected_cells':int((~pos).sum()),'genes_tested':len(r),'upregulated_excluding_PGAM5':len(up),'downregulated_excluding_PGAM5':len(down),'top20_up':up.gene.head(20).tolist(),'top10_down':down.gene.head(10).tolist(),'criteria':'FDR <0.05, |Scanpy approximate log2FC| >=1, >=10% detection in either group; PGAM5 excluded from candidate lists','method':'pooled-cell Wilcoxon with tie correction and BH FDR; no patient pairing, depth matching or sample/depth covariates','QC':'n_genes >=500 and pct_mt <20; total-count normalize to 10000 then log1p','annotation_rule':'Each C1QA/B/C detected in >=80%, and CD68 in >=75%, of cells in original global myeloid cluster; monocyte/DC/mast clusters excluded'}
(R/'analysis_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
print(up[['gene','log2FC','FDR','fraction_PGAM5_positive','fraction_PGAM5_undetected']].head(20).to_string(),flush=True)
# Primary-T-only marker sanity checks for the retained annotation, not DE evidence.
records=[]
for cluster in keep:
 mask=a.obs['res.3'].eq(cluster).to_numpy()
 for g in ['C1QA','C1QB','C1QC','CD68','CSF1R','CD163','SPP1','MMP9','FCN1','S100A8','CD1C','FCER1A','GZMB','TPSAB1']:
  idx=np.flatnonzero(a.var.gene_symbol.eq(g));v=np.asarray(a.X[mask][:,idx].sum(1)).ravel()
  records.append({'cluster':cluster,'gene':g,'cells':int(mask.sum()),'fraction':float((v>0).mean()),'mean_log1p_10k':float(v.mean())})
pd.DataFrame(records).to_csv(R/'primary_T_macrophage_cluster_markers.csv',index=False)
