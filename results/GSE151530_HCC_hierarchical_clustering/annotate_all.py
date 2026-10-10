from pathlib import Path
import pandas as pd
R=Path(__file__).resolve().parent
names={
0:('T cells','CCL5/CD8A cytotoxic T','CCL5/CD8A细胞毒性T细胞','supported'),
1:('T cells','IL7R/LTB T cells','IL7R/LTB T细胞','supported'),
2:('T cells','FOXP3/CTLA4 Treg','FOXP3/CTLA4调节性T细胞','supported'),
3:('Mixed RNA, unresolved','KIT/ALAS2 mixed RNA','KIT/ALAS2混合RNA群（含T/肥大/红系信号，待定）','mixed'),
4:('NK/T lineage, mixed','FGFBP2/GNLY NK-T lineage','FGFBP2/GNLY NK/T谱系混合群','mixed'),
5:('NK/T lineage, mixed','KLRD1/TRDC NK-T lineage','KLRD1/TRDC NK/γδT样混合群','mixed'),
6:('Tumor/hepatocyte-like (author)','GSTA1/KRT8 epithelial-tumor-like','GSTA1/KRT8上皮/肿瘤样细胞（作者恶性标签）','author-supported'),
7:('B cells','MS4A1/CD79A B cells','MS4A1/CD79A B细胞','supported'),
8:('Monocyte/DC-like myeloid','FCN1/LGALS2 mono-DC-like myeloid','FCN1/LGALS2单核/DC样髓系细胞','qualified'),
9:('Macrophage-like','MSR1/C1QC macrophage-like','MSR1/C1QC巨噬细胞样群','supported'),
10:('T cells','MKI67/TOP2A cycling T','MKI67/TOP2A增殖T细胞','supported'),
11:('Macrophage-like','C1QA/HLA-DRA macrophage-like','C1QA/HLA-DRA巨噬细胞样群','supported'),
12:('Tumor/hepatocyte-like (author)','ACSM2B/CPS1 hepatocyte-tumor-like','ACSM2B/CPS1肝细胞/肿瘤样细胞（作者恶性标签）','author-supported'),
13:('Epithelial/myeloid mixed RNA','PGA5/CHGA epithelial-mixed RNA','PGA5/CHGA上皮/髓系混合RNA群','mixed'),
14:('Endothelial cells','PLVAP/CLEC14A endothelial','PLVAP/CLEC14A内皮细胞','supported'),
15:('Tumor/hepatocyte-like (author)','BIRC5/TOP2A cycling tumor-like','BIRC5/TOP2A增殖肿瘤样细胞（作者恶性标签）','author-supported'),
16:('Tumor/hepatocyte-like (author)','VNN1/MUC13 epithelial-tumor-like','VNN1/MUC13上皮/肿瘤样细胞（作者恶性标签）','author-supported'),
17:('Tumor/hepatocyte-like (author)','CA12/IGF2BP3 tumor-like','CA12/IGF2BP3肿瘤样细胞（作者恶性标签）','author-supported'),
18:('Stromal/mural cells','COL1A2/RGS5 stromal-mural','COL1A2/RGS5基质/血管壁细胞','qualified'),
19:('Endothelial cells','CLEC4G/STAB2 sinusoidal endothelial','CLEC4G/STAB2肝窦内皮样细胞','supported'),
20:('pDC-like cells','CLEC4C/LILRA4 pDC-like','CLEC4C/LILRA4浆细胞样DC','supported'),
21:('Plasma cells','MZB1/JCHAIN plasma cells','MZB1/JCHAIN浆细胞','supported'),
22:('Epithelial-like cells, qualified','EPCAM/KRT19 epithelial-like','EPCAM/KRT19胆管/上皮样细胞','qualified'),
23:('T cells','GIMAP5/SLFN12L T cells','GIMAP5/SLFN12L T细胞','supported')}
c=pd.read_csv(R/'allcell_cluster_author_counts.csv',index_col=0);samples=pd.read_csv(R/'allcell_cluster_sample_counts.csv',index_col=0);rows=[]
for group,(typ,en,cn,status) in names.items():
 v=c.loc[group];s=samples.loc[group];rows.append(dict(cluster=str(group),cell_type_EN=typ,subtype_EN=en,subtype_CN=cn,status=status,cells=int(v.sum()),author_majority=v.idxmax(),author_majority_fraction=float(v.max()/v.sum()),samples=int(s.gt(0).sum()),largest_sample_fraction=float(s.max()/s.sum())))
pd.DataFrame(rows).to_csv(R/'allcell_cluster_annotations.csv',index=False)
