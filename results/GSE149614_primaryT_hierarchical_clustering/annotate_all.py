from pathlib import Path
import pandas as pd
R=Path(__file__).resolve().parent
names={
0:('T cells','IL7R/CD3E T cells','IL7R/CD3E T细胞','supported'),
1:('NK/gamma-delta T-like','KLRD1/TRDC NK-gamma-delta T-like','KLRD1/TRDC NK/γδT样细胞','qualified'),
2:('Monocyte/DC-like myeloid','FCN1/CD1C mono-DC-like myeloid','FCN1/CD1C单核/DC样髓系细胞','qualified'),
3:('Macrophages','MSR1/C1QC macrophages','MSR1/C1QC巨噬细胞','supported'),
4:('DC-like cells','CLEC9A/XCR1 DC-like, pDC admixture','CLEC9A/XCR1 DC样群（伴pDC信号）','qualified'),
5:('Macrophages','TIMD4/CD5L resident-like macrophages','TIMD4/CD5L驻留样巨噬细胞','qualified'),
6:('Endothelial cells','CLEC14A/EMCN endothelial','CLEC14A/EMCN内皮细胞','supported'),
7:('Fibroblast/mural cells','COL1A2/PDGFRB stromal-mural','COL1A2/PDGFRB基质/血管壁细胞','qualified'),
8:('Hepatocyte/epithelial-like','SLCO1B3/CYP2A6 hepatocyte-like','SLCO1B3/CYP2A6肝细胞样群','supported'),
9:('T cells','MKI67/TOP2A cycling T-like','MKI67/TOP2A增殖T细胞样群','qualified'),
10:('Mast cells','TPSAB1/CPA3 mast cells','TPSAB1/CPA3肥大细胞','supported'),
11:('Hepatocyte/epithelial-like','ASGR1/CPS1 hepatocyte-like','ASGR1/CPS1肝细胞样群','supported'),
12:('Hepatocyte/epithelial-like','PRSS3/AGR2 epithelial-like','PRSS3/AGR2上皮样群','qualified'),
13:('Hepatocyte/epithelial-like','KRT19/EPCAM epithelial-like','KRT19/EPCAM上皮样群','supported'),
14:('T cells','FOXP3/CTLA4 Treg-like','FOXP3/CTLA4调节性T细胞样群','supported'),
15:('T/hepatic mixed RNA','CD8A/APOC3 T-hepatic mixed RNA','CD8A/APOC3 T/肝源混合RNA群','mixed'),
16:('Macrophages','MKI67/TOP2A cycling macrophage-like','MKI67/TOP2A增殖巨噬细胞样群','qualified'),
17:('Hepatocyte/epithelial-like','TOP2A/BIRC5 cycling epithelial-like','TOP2A/BIRC5增殖上皮样群','qualified'),
18:('Hepatocyte/epithelial-like','GSTA1/KRT18 hepatocyte-like','GSTA1/KRT18肝细胞样群','supported'),
19:('B cells','MS4A1/CD79A B cells','MS4A1/CD79A B细胞','supported'),
20:('Plasma/B lineage','MZB1/JCHAIN plasma-B lineage','MZB1/JCHAIN浆细胞/B谱系群','qualified'),
21:('Mac/hepatic mixed RNA','C1QB/CYP2D6 Mac-hepatic mixed RNA','C1QB/CYP2D6巨噬/肝源混合RNA群','mixed')}
c=pd.read_csv(R/'allcell_cluster_author_counts.csv',index_col=0);samples=pd.read_csv(R/'allcell_cluster_sample_counts.csv',index_col=0);assert set(c.index)==set(names);rows=[]
for group,(typ,en,cn,status) in names.items():
 v=c.loc[group];s=samples.loc[group];rows.append(dict(cluster=str(group),cell_type_EN=typ,subtype_EN=en,subtype_CN=cn,status=status,cells=int(v.sum()),author_majority=v.idxmax(),author_majority_fraction=float(v.max()/v.sum()),samples=int(s.gt(0).sum()),largest_sample_fraction=float(s.max()/s.sum())))
pd.DataFrame(rows).to_csv(R/'allcell_cluster_annotations.csv',index=False)
