from pathlib import Path
import pandas as pd,numpy as np,json,anndata as ad
from scipy import sparse
R=Path(__file__).resolve().parent
# Locked marker-based curation of the saved 21-cluster solution. No reclustering or coordinate changes.
rows=[
(0,'Macrophages','巨噬细胞','FOLR2/CD163 macrophages','FOLR2/CD163巨噬细胞','supported_RNA','C1QA,C1QB,C1QC,CSF1R,CD68','C1Q genes almost universal; CSF1R 85%; FOLR2/CD163 are top enriched genes; ALB alone is not used for lineage.',False),
(1,'CD1C DC-like myeloid','CD1C树突样髓系细胞','CD1C/FCER1A DC-like myeloid','CD1C/FCER1A树突样髓系细胞','qualified_RNA','CD1C,FCER1A,CLEC10A,LST1,TYROBP','CD1C/FCER1A/CLEC10A/FLT3 enriched; concurrent C1Q/CSF1R prevents claiming a pure cDC2 population.',True),
(2,'Macrophages','巨噬细胞','C1QA/C1QB macrophages','C1QA/C1QB巨噬细胞','supported_RNA','C1QA,C1QB,C1QC,CD68,LST1','C1QA 90%; C1Q/APOC1/HLA-II dominant; CSF1R detection is lower (28%); macrophage RNA identity retained.',False),
(3,'Monocytes','单核细胞','FCN1/VCAN monocytes','FCN1/VCAN单核细胞','supported_RNA','FCN1,CD14,S100A8,S100A9,LST1','FCN1 82%; FCN1/VCAN/CD300E/S100A8/S100A9 enrichment distinguishes the monocyte-like program.',False),
(4,'T cells','T细胞','FOXP3/CTLA4 Treg','FOXP3/CTLA4调节性T细胞','supported_RNA','CD3D,CD3E,TRAC,FOXP3,CTLA4,IL2RA','CD3/TRAC-positive; FOXP3 44%, CTLA4 58%, IL2RA 34%; FOXP3/CTLA4/IKZF2/TNFRSF18 enriched.',False),
(5,'T cells','T细胞','GZMH/NKG7 cytotoxic T','GZMH/NKG7细胞毒性T细胞','supported_RNA','CD3D,CD3E,CD3G,TRAC,NKG7,CD8A,CD8B','Multiple CD3 genes 64-69%; NKG7 91%; GZMH/CCL5/CD8A/CD8B enriched. NK-like cytotoxic genes do not alone establish NK identity or a pure CD8 subset.',False),
(6,'Hepatocyte-like cells','肝细胞样细胞','AOX1/TSPAN8 hepatocyte-like','AOX1/TSPAN8肝细胞样细胞','supported_RNA','ALB,APOA1,TTR,EPCAM,KRT8,KRT18','Coherent hepatic/epithelial program; AOX1/TSPAN8/SLCO1B3/MAL2 enriched; no malignancy inferred from RNA alone.',False),
(7,'Hepatocyte-like cells','肝细胞样细胞','APOC3/APOA2 hepatocyte-like','APOC3/APOA2肝细胞样细胞','supported_RNA','ALB,APOA1,TTR,CPS1,ASGR1','APOC3/APOA2/FABP1/ALB/ALDOB coherent hepatic program; C1Q RNA alone not sufficient to relabel as macrophages.',False),
(8,'Hepatocyte-like cells','肝细胞样细胞','CYP3A5/PCK1 hepatocyte-like','CYP3A5/PCK1肝细胞样细胞','supported_RNA','ALB,APOA1,TTR,CPS1,ASGR1','CYP3A5/ACSM2B/PCK1/GPAM/ONECUT1 metabolic hepatic program; no malignant-versus-normal claim.',False),
(9,'NK cells','NK细胞','FGFBP2/GNLY NK','FGFBP2/GNLY NK细胞','supported_RNA','NKG7,GNLY,KLRD1,PRF1','GNLY 92%; FGFBP2/GNLY/GZMB/KLRD1/PRF1 enriched; T lineage scores much lower.',False),
(10,'Cycling mixed lineage','增殖混合谱系群','MKI67/TOP2A cycling mixed','MKI67/TOP2A增殖混合谱系群','mixed_RNA','MKI67,TOP2A,CENPF,C1QA,CD3D,TRAC','Strong cycling RNA; C1QA 65% and CD3/TRAC about26%; contains more than one lineage and is not assigned as pure cycling macrophages.',True),
(11,'T cells','T细胞','IL7R/LTB T cells','IL7R/LTB T细胞','supported_RNA','CD3D,CD3E,TRAC,IL7R,TCF7','CD3D 63%, TRAC 51%; IL7R/LTB/CD2/CD40LG/BCL11B enriched. No pure CD4-only claim.',False),
(12,'Endothelial cells','内皮细胞','KDR/PLVAP endothelial','KDR/PLVAP内皮细胞','supported_RNA','PECAM1,VWF,KDR,EMCN','PECAM1 78%; KDR/PLVAP/EMCN/CDH5/ADGRL4 enrichment supports endothelial lineage.',False),
(13,'Endothelial cells','内皮细胞','SEMA3G/GJA5 endothelial','SEMA3G/GJA5内皮细胞','supported_RNA','PECAM1,VWF,KDR,EMCN','PECAM1 86%; SEMA3G/GJA5/CLDN5/PTPRB enrichment supports an arterial-like endothelial RNA state.',False),
(14,'Mural cells','血管壁细胞','RGS5/ACTA2 mural cells','RGS5/ACTA2血管壁细胞','supported_RNA','RGS5,ACTA2,TAGLN,MYH11','RGS5/ACTA2/TAGLN each about80%; ECM program coexists. Pericyte/smooth-muscle-like, not claimed as a pure fibroblast subtype.',False),
(15,'B/plasma lineage','B/浆细胞谱系群','MS4A1/CD79A B-plasma lineage','MS4A1/CD79A B/浆细胞谱系群','qualified_RNA','MS4A1,CD79A,CD79B,MZB1,JCHAIN','MS4A1 60%, CD79A 64%; MZB1/JCHAIN about30%; B lineage is supported, but this resolution does not isolate a pure plasma population.',True),
(16,'Macrophages','巨噬细胞','MT1G/TREM2 macrophages','MT1G/TREM2巨噬细胞','supported_RNA','C1QA,C1QB,C1QC,CSF1R,CD68','C1QA 95%, CSF1R 70%; MT1-family/TREM2/GPNMB/HMOX1 enriched; descriptive metal/stress RNA state.',False),
(17,'NK cells','NK细胞','XCL1/GZMK NK','XCL1/GZMK NK细胞','supported_RNA','NKG7,KLRD1,PRF1','NKG7 93%, KLRD1 82%; XCL1/XCL2/KLRC1/GZMK enriched; CD3D only3%.',False),
(18,'Neutrophils','中性粒细胞','FCGR3B/CXCR2 neutrophils','FCGR3B/CXCR2中性粒细胞','supported_RNA','FCGR3B,CSF3R,S100A8,S100A9','FCGR3B 78%; FCGR3B/CXCR2/FPR1/FPR2/S100A12 enriched; distinct from FCN1 monocytes.',False),
(19,'Mac-T/NK mixed RNA','巨噬/T-NK混合RNA群','C1QC/NKG7 Mac-T/NK mixed','C1QC/NKG7巨噬/T-NK混合RNA群','mixed_RNA','C1QA,C1QB,C1QC,CSF1R,CD3D,NKG7','C1Q about95-98%, CSF1R77%, CD3D47%, NKG765%; macrophage and T/NK signals coexist; possible multiplets/admixture, not confirmed doublets.',True),
(20,'Hepatocyte-like cells','肝细胞样细胞','MAGEA1/GLS2 hepatocyte-like','MAGEA1/GLS2肝细胞样细胞','patient_dominated_RNA','ALB,APOA1,TTR,CPS1,ASGR1','Hepatic metabolic markers with MAGEA1/GLS2; 98% from one patient; no generalizable subtype or malignancy claim.',False)]
d=pd.DataFrame(rows,columns=['cluster','cell_type_EN','cell_type_CN','subtype_EN','subtype_CN','annotation_confidence','anchor_genes','rationale','boundary_or_mixed_flag']);d.cluster=d.cluster.astype(str);profiles=pd.read_csv(R/'cluster_profiles.csv',dtype={'cluster':str});d=d.merge(profiles[['cluster','cells','patients','largest_patient_fraction']],on='cluster',validate='one_to_one');d['cluster_label']='C'+d.cluster+' '+d.subtype_EN;d.to_csv(R/'cluster_annotations.csv',index=False)
a=pd.read_csv(R/'cell_metadata_clustered.csv.gz',index_col=0,dtype={'cluster':str});mapping=d.set_index('cluster');assert set(a.cluster)==set(mapping.index)
for col in ['cell_type_EN','cell_type_CN','subtype_EN','subtype_CN','annotation_confidence','boundary_or_mixed_flag','cluster_label']:a[col]=a.cluster.map(mapping[col]);assert a[col].notna().all()
old=pd.read_csv(R/'previous_macrophage_clusters.csv.gz',index_col=0,usecols=['cell_id_index','primary_cluster'],dtype={'primary_cluster':str});a['previous_macrophage_cluster']=old.primary_cluster.reindex(a.index).fillna('');a.index.name='cell_id';a.to_csv(R/'cell_annotations.csv.gz')
# Compact portable AnnData: marker genes only, original full-library normalization denominator.
raw=sparse.load_npz(R/'marker_raw_counts.npz');genes=pd.read_csv(R/'marker_genes.csv').gene.tolist();z=raw.astype(np.float32).multiply(10000/a.total_counts.to_numpy()[:,None]).tocsr();z.data=np.log1p(z.data);compact=ad.AnnData(z,obs=a.copy(),var=pd.DataFrame(index=pd.Index(genes,name='gene')));compact.layers['counts']=raw;u=pd.read_csv(R/'UMAP_coordinates.csv.gz',index_col=0);assert u.index.equals(a.index);compact.obsm['X_umap']=u[['UMAP1','UMAP2']].to_numpy();compact.obsm['X_umap_unintegrated']=u[['unintegrated_UMAP1','unintegrated_UMAP2']].to_numpy();compact.uns['scope']='MARKER GENES ONLY; full-library CP10k denominator; not a full expression matrix';compact.write_h5ad(R/'HCC_annotation_marker_subset.h5ad',compression='gzip')
# Full matrix remains local, while portable annotation is published.
full=ad.read_h5ad(R/'HCC_normalized_clustered.h5ad');assert full.obs_names.equals(a.index)
for col in a.columns:
 if col not in full.obs:full.obs[col]=a[col].to_numpy()
full.write_h5ad(R/'HCC_annotated_full.h5ad',compression='gzip')
print(a.groupby('cell_type_EN').size().to_string(),flush=True)
