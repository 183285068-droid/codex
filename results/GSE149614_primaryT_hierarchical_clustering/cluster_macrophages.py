from pathlib import Path
import scanpy as sc,numpy as np,pandas as pd,json,harmonypy
from threadpoolctl import threadpool_limits
from sklearn.metrics import adjusted_rand_score
R=Path(__file__).resolve().parent
if (R/'HCC_allcell_clustered.h5ad').exists():
 b=sc.read_h5ad(R/'HCC_allcell_clustered.h5ad',backed='r')
 defs=pd.read_csv(R/'allcell_cluster_annotations.csv',dtype={'cluster':str});keep=defs.loc[defs.cell_type_EN.eq('Macrophages'),'cluster'].tolist()
 gate_genes=['C1QA','C1QB','C1QC','CD68'];gate_counts=b.layers['counts'][:,b.var_names.get_indexer(gate_genes)].toarray();records=[]
 for c in sorted(b.obs.loc[b.obs.author_type.eq('Myeloid'),'res.3'].unique()):
  ix=(b.obs.author_type.eq('Myeloid')&b.obs['res.3'].eq(c)).to_numpy();det=(gate_counts[ix]>0).mean(0);records.append(dict(original_cluster=int(c),cells=int(ix.sum()),**dict(zip(gate_genes,det)),included=bool((det[:3]>=.8).all() and det[3]>=.75)))
 gate=pd.DataFrame(records);gate.to_csv(R/'primaryT_original_myeloid_marker_gate.csv',index=False);legacy=gate.loc[gate.included,'original_cluster'].tolist()
 selected=b.obs.author_type.eq('Myeloid')&b.obs['res.3'].isin(legacy)&b.obs.cluster.astype(str).isin(keep)
 excluded=b.obs.loc[b.obs.author_type.eq('Myeloid')&~selected].copy();excluded['original_C1Q_CD68_gate_pass']=excluded['res.3'].isin(legacy);excluded['allcell_macrophage_group']=excluded.cluster.astype(str).isin(keep);excluded['exclusion_reason']=np.where(excluded.original_C1Q_CD68_gate_pass,'Not in newly confirmed all-cell macrophage groups','T-sample original-group C1Q/CD68 gate failed');excluded.to_csv(R/'excluded_myeloid_cells.csv.gz')
 a=b[selected.to_numpy(),:].to_memory();b.file.close();a.obs['parent_cluster']=a.obs.cluster.astype(str)
else:
 if (R/'inputs/core_macrophages_raw.h5ad').exists():
  a=sc.read_h5ad(R/'inputs/core_macrophages_raw.h5ad');a.layers['counts']=a.X.copy()
 else:
  a=sc.read_h5ad(R/'macrophages_annotated_full.h5ad');a.X=a.layers['counts'].copy()
a.obs['Sample']=a.obs.Sample.astype('category').cat.remove_unused_categories();assert len(a)>500;a.X=a.layers['counts'].copy();a.write_h5ad(R/'macrophages_raw.h5ad',compression='gzip');a.obsm.clear();a.obsp.clear();a.uns.clear();sc.settings.n_jobs=4
sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a);sc.pp.highly_variable_genes(a,n_top_genes=2000,flavor='seurat',batch_key='Sample');a.var.loc['PGAM5','highly_variable']=False;v=a[:,a.var.highly_variable].copy();v.layers.clear(keep_x=True);sc.pp.scale(v,max_value=10);sc.tl.pca(v,n_comps=40,random_state=149614);a.obsm['X_pca']=v.obsm['X_pca'].copy();del v
sc.pp.neighbors(a,n_neighbors=20,n_pcs=30,key_added='unintegrated',random_state=149614);sc.tl.umap(a,neighbors_key='unintegrated',random_state=149614);a.obsm['X_umap_unintegrated']=a.obsm['X_umap'].copy()
with threadpool_limits(limits=4):h=harmonypy.run_harmony(a.obsm['X_pca'][:,:30],a.obs,'Sample',random_state=149614,max_iter_harmony=20)
z=np.asarray(h.Z_corr);a.obsm['X_pca_harmony']=z.T if z.shape[0]==30 else z;sc.pp.neighbors(a,n_neighbors=20,use_rep='X_pca_harmony',random_state=149614);sc.tl.umap(a,random_state=149614)
for r in [.4,.6,.8]:sc.tl.leiden(a,resolution=r,key_added='leiden_'+str(r),random_state=149614,flavor='igraph',n_iterations=2,directed=False)
a.obs['mac_cluster']=a.obs['leiden_0.6'].map(lambda c:'M'+str(c)).astype('category');sc.tl.leiden(a,resolution=.6,key_added='mac_seed42',random_state=42,flavor='igraph',n_iterations=2,directed=False);print('Mac clusters',a.obs.mac_cluster.value_counts().to_dict(),flush=True)
sc.tl.rank_genes_groups(a,groupby='mac_cluster',method='wilcoxon',tie_correct=True,n_genes=50,use_raw=False);sc.get.rank_genes_groups_df(a,group=None).to_csv(R/'macrophage_top50_markers.csv.gz',index=False)
a.write_h5ad(R/'macrophages_clustered.h5ad',compression='gzip');a.obs.to_csv(R/'macrophage_metadata.csv.gz');pd.DataFrame(a.obsm['X_umap'],index=a.obs_names,columns=['UMAP1','UMAP2']).to_csv(R/'macrophage_UMAP_coordinates.csv.gz');a.obs.groupby(['mac_cluster','Sample'],observed=True).size().unstack(fill_value=0).to_csv(R/'macrophage_cluster_sample_counts.csv')
genes='C1QA C1QB C1QC CSF1R CD68 MERTK FOLR2 CD163 LYVE1 MRC1 TREM2 APOE SPP1 CD5L MARCO HLA-DRA CD74 IL1B CXCL8 FCN1 VCAN S100A8 S100A9 MKI67 TOP2A ISG15 IFIT1 CXCL9 CXCL10 MT1G MT2A HMOX1 NKG7 CD3D CD1C FCER1A CLEC10A ALB KRT19 PECAM1 PGAM5'.split();genes=[g for g in genes if g in a.var_names];rows=[]
for c in a.obs.mac_cluster.cat.categories:
 ix=a.obs.mac_cluster.eq(c).to_numpy()
 for g in genes:
  j=a.var_names.get_loc(g);rows.append(dict(cluster=c,gene=g,mean_log1p_CP10k=float(a.X[ix,j].mean()),detection_percent=float((a.layers['counts'][ix,j]>0).mean()*100)))
pd.DataFrame(rows).to_csv(R/'macrophage_marker_evidence.csv',index=False);summary=dict(cells=len(a),selection='Author Myeloid, original res.3 clusters with each C1QA/B/C detected >=80% and CD68 >=75%, recomputed only in the 10 primary T samples, intersect newly confirmed macrophage all-cell clusters; competing myeloid excluded',samples=int(a.obs.Sample.nunique()),patient_proxies=int(a.obs.patient_proxy.nunique()),clusters=int(a.obs.mac_cluster.nunique()),resolution=.6,resolution_sensitivity={str(r):int(a.obs['leiden_'+str(r)].nunique()) for r in [.4,.6,.8]},seed_ARI=adjusted_rand_score(a.obs.mac_cluster,a.obs.mac_seed42),PGAM5_detected=int(a.obs.PGAM5_detected.sum()),HVG='2000 batch-aware Seurat, PGAM5 not used in geometry');(R/'macrophage_summary.json').write_text(json.dumps(summary,indent=2));print('DONE',summary,flush=True)
