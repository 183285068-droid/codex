from pathlib import Path
import json,numpy as np,pandas as pd,anndata as ad,scanpy as sc
R=Path(__file__).resolve().parent
q=pd.read_csv(R/'global_cluster_marker_QA.csv',dtype={'cluster':str})
rule=q[['fraction_C1QA','fraction_C1QB','fraction_C1QC']].min(axis=1).ge(.75)&q.fraction_CSF1R.ge(.60)&q.fraction_TYROBP.ge(.90)&q.fraction_CD1C.lt(.40)&q.fraction_FCN1.lt(.50)
q['included_as_macrophage']=rule;q.to_csv(R/'cluster_annotation_decisions.csv',index=False)
clusters=q.loc[rule,'cluster'].tolist();assert clusters
o=pd.read_csv(R/'all_tumor_QC_cluster_metadata.csv.gz',index_col=0,dtype={'leiden':str})
a=ad.read_h5ad(R/'all_tumor_counts_QC.h5ad');assert list(a.obs_names)==list(o.index)
a.obs['leiden']=o.leiden.to_numpy()
t=a[a.obs.leiden.isin(clusters).to_numpy()].copy();del a
t.obs['annotation']='marker_inferred_macrophage'
t.obs.to_csv(R/'macrophage_analysis_metadata.csv.gz');t.write_h5ad(R/'macrophage_counts_QC.h5ad',compression='gzip')
t.obs.groupby('Sample',observed=True).agg(macrophages=('barcode','size'),PGAM5_positive=('PGAM5_detected','sum'),rate=('PGAM5_detected','mean')).to_csv(R/'PGAM5_by_sample.csv')
pos=t.obs.PGAM5_counts.to_numpy()>0;assert min(pos.sum(),(~pos).sum())>=2
print('Macrophage clusters',clusters,'cells',len(t),'PGAM5+',pos.sum(),'PGAM5 zero',(~pos).sum(),flush=True)
fp=np.asarray((t.X[pos]>0).mean(0)).ravel();fn=np.asarray((t.X[~pos]>0).mean(0)).ravel()
t.obs['PGAM5_group']=pd.Categorical(np.where(pos,'detected','undetected'))
sc.pp.normalize_total(t,target_sum=1e4);mp=np.asarray(t.X[pos].mean(0)).ravel();mn=np.asarray(t.X[~pos].mean(0)).ravel();sc.pp.log1p(t)
sc.tl.rank_genes_groups(t,groupby='PGAM5_group',groups=['detected'],reference='undetected',method='wilcoxon',tie_correct=True,corr_method='benjamini-hochberg',use_raw=False)
r=sc.get.rank_genes_groups_df(t,group='detected').set_index('names').rename(columns={'scores':'Wilcoxon_score','logfoldchanges':'log2FC','pvals':'p_value','pvals_adj':'FDR'})
r['gene']=t.var.loc[r.index,'gene'].values
aux=pd.DataFrame({'fraction_PGAM5_positive':fp,'fraction_PGAM5_undetected':fn,'mean_10k_detected':mp,'mean_10k_undetected':mn},index=t.var_names);r=r.join(aux);r['defining_gene']=r.gene.eq('PGAM5')
r.to_csv(R/'GSE242889_HCC_tumor_PGAM5_full_DE.csv')
mask=r.FDR.lt(.05)&r.log2FC.abs().ge(1)&r[['fraction_PGAM5_positive','fraction_PGAM5_undetected']].max(axis=1).ge(.1)
passed=r[mask].sort_values(['FDR','log2FC'],ascending=[True,False]);passed.to_csv(R/'GSE242889_HCC_tumor_PGAM5_DE_including_PGAM5.csv')
cand=passed[~passed.defining_gene];up=cand[cand.log2FC.ge(1)];down=cand[cand.log2FC.le(-1)].sort_values(['FDR','log2FC'])
up.to_csv(R/'GSE242889_HCC_tumor_PGAM5_upregulated.csv');down.to_csv(R/'GSE242889_HCC_tumor_PGAM5_downregulated.csv');up.head(20).to_csv(R/'GSE242889_HCC_tumor_PGAM5_top20_upregulated.csv')
qa=[]
for gene in ['C1QA','C1QB','C1QC','CD68','CSF1R','LST1','TYROBP','SPP1','ALB','APOA1','TTR','EPCAM','KRT19','CD3D','TRAC','MS4A1','FCN1','S100A8']:
 ix=np.flatnonzero(t.var.gene.eq(gene))
 if not len(ix):continue
 v=t.X[:,ix].toarray().sum(1)
 for label,sel in [('PGAM5_detected',pos),('PGAM5_undetected',~pos)]:qa.append({'gene':gene,'group':label,'cells':int(sel.sum()),'fraction':float((v[sel]>0).mean()),'mean_log1p_10k':float(v[sel].mean())})
pd.DataFrame(qa).to_csv(R/'TAM_marker_QA.csv',index=False)
all_o=pd.read_csv(R/'all_tumor_cell_metadata_before_QC.csv.gz',index_col=0)
s={'HCC_tumor_samples':int(all_o.Sample.nunique()),'HCC_tumor_patients':int(all_o.patient.nunique()),'source_tumor_cells':len(all_o),'QC_tumor_cells':int(all_o.passes_QC.sum()),'QC_inferred_macrophages':len(t),'PGAM5_positive_cells':int(pos.sum()),'PGAM5_undetected_cells':int((~pos).sum()),'samples_with_macrophages':int(t.obs.Sample.nunique()),'samples_with_positive_macrophages':int(t.obs.loc[pos,'Sample'].nunique()),'tested_genes':len(r),'macrophage_clusters':clusters,'upregulated_excluding_PGAM5':len(up),'downregulated_excluding_PGAM5':len(down),'top20_up':up.gene.head(20).tolist(),'top10_down':down.gene.head(10).tolist(),'scope':'All five HCC tumor samples 1T-5T, MVI+ and MVI- included without MVI filtering; adjacent NT excluded following preceding tumor-only scope','annotation':'Exploratory marker-inferred C1Q-rich macrophage clusters; no author cell annotation supplied by GEO. Not guaranteed to capture C1Q-low macrophages.','annotation_rule':'Global Leiden cluster C1QA/B/C each >=75% detection, CSF1R>=60%, TYROBP>=90%, CD1C<40%, FCN1<50%; target-agnostic marker review, PGAM5 excluded from HVG features','clustering':'QC all tumor cells, normalize to 10000, log1p, 2000 Seurat HVGs minus PGAM5, 30 PCs, 15 neighbors, igraph Leiden resolution 1, random seed 0','QC':'n_genes>=500, pct_mt<20','method':'Pooled-cell Wilcoxon with tie correction, BH across all genes; no patient pairing, depth matching or patient/sample/MVI/depth covariates','criteria':'FDR<0.05, |Scanpy approximate log2FC|>=1, detection>=10% in either group; PGAM5 excluded from candidates','gene_merge':'Outer union of Ensembl IDs; missing rows set to zero; symbol from first source annotation. Input feature lists differ.','MVI_note':'All five tumors retained. GEO sample characteristics do not provide per-sample MVI labels; no unverified sample-level MVI assignment was inferred.'}
(R/'analysis_summary.json').write_text(json.dumps(s,indent=2));print(json.dumps(s,indent=2),flush=True)
