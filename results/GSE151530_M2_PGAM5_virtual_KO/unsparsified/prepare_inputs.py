from pathlib import Path
import json,hashlib,sys
import anndata as ad,scanpy as sc,numpy as np,pandas as pd
R=Path(__file__).resolve().parent;I=R/'inputs';I.mkdir(exist_ok=True)
if not (I/'M2_raw.h5ad').exists():
 source=Path(sys.argv[1]) if len(sys.argv)>1 else R.parent/'GSE151530_HCC_hierarchical_clustering/macrophages_annotated_full.h5ad'
 b=ad.read_h5ad(source);a=b[b.obs.mac_cluster.astype(str).eq('M2')].copy();a.X=a.layers['counts'].copy();a.layers.clear();a.obsp.clear();a.uns.clear();a.obsm.clear();a.write_h5ad(I/'M2_raw.h5ad',compression='gzip')
 context=b.obs[['mac_cluster','subtype_EN','PGAM5_counts','total_counts']].copy();context[['UMAP1','UMAP2']]=b.obsm['X_umap'];context.to_csv(I/'macrophage_context.csv.gz')
else:a=ad.read_h5ad(I/'M2_raw.h5ad')
assert len(a)==171 and a.obs.mac_cluster.astype(str).eq('M2').all()
assert np.array_equal(np.asarray(a.X.sum(1)).ravel(),a.obs.total_counts)
assert np.array_equal(a.X[:,a.var_names.get_loc('PGAM5')].toarray().ravel(),a.obs.PGAM5_counts)
a.obs.to_csv(I/'target_cell_metadata.csv.gz');z=a.copy();sc.pp.normalize_total(z,target_sum=1e4);sc.pp.log1p(z);sc.pp.highly_variable_genes(z,n_top_genes=1500,flavor='seurat')
det=np.asarray((a.X>0).mean(0)).ravel();eligible=(det>=.05)&~a.var_names.str.match(r'^(MT-|RPL|RPS)');hv=z.var.loc[eligible].sort_values('dispersions_norm',ascending=False).head(1500).index;genes=sorted(set(hv)|{'PGAM5'});x=pd.DataFrame(a[:,genes].X.toarray().T,index=genes,columns=a.obs_names);x.to_csv(I/'network_raw_counts.csv.gz');pd.DataFrame(dict(gene=a.var_names,detected_fraction=det,eligible=eligible,network_selected=a.var_names.isin(genes))).to_csv(R/'gene_selection.csv',index=False)
a.obs.groupby(['Sample','patient_proxy'],observed=True).agg(cells=('PGAM5_counts','size'),PGAM5_detected=('PGAM5_counts',lambda v:(v>0).sum())).reset_index().to_csv(R/'sample_summary.csv',index=False)
print('Inputs checked:',len(a),'cells;',len(genes),'genes;',a.obs.PGAM5_counts.gt(0).sum(),'PGAM5 detected')
