# GSE202642巨噬细胞重新降维聚类

从此前全细胞C0–C20注释中提取C0、C2、C16，共7例HCC肿瘤的11,536个细胞。原全细胞C10增殖混合群和C19巨噬/T-NK混合群未纳入。沿用之前QC（n_genes≥500，线粒体计数比例<20%），没有再按PGAM5筛选细胞。

## 结果

重新计算高变基因、PCA、Harmony、邻接图和UMAP，Leiden分辨率0.6得到9个探索性群，使用M0–M8编号。各群均覆盖7个患者标签。PGAM5 RNA在712个细胞中检出，分布于多个群；M8增殖群检出率最高（15.29%）。不能凭此认定PGAM5构成独立亚群。

|新亚群|RNA标志命名|细胞数|PGAM5检出率|
|---|---|---:|---:|
|M0|FABP5/TREM2巨噬细胞|2282|8.76%|
|M1|MRC1/CD163巨噬细胞|3209|4.74%|
|M2|C1QA/APOC1巨噬细胞（肝源RNA富集）|1869|2.03%|
|M3|NPR3/MYCT1巨噬富集群（待定）|51|0.00%|
|M4|CXCL9/CXCL10干扰素相关巨噬细胞|350|6.57%|
|M5|HSPA1A/DNAJB1应激相关巨噬细胞|2977|7.49%|
|M6|MT1G/HMOX1巨噬细胞|508|8.46%|
|M7|C1QA/CD3D混合RNA群|133|6.77%|
|M8|MCM4/PCLAF增殖巨噬细胞|157|15.29%|

## 方法与稳定性

- 原始全转录组计数→每细胞10,000计数归一化→log1p。
- 按sample_name批次筛选2,500个Seurat高变基因；PGAM5不作为驱动聚类的高变特征，仍保留其表达用于展示和统计。
- 高变基因缩放（截断10），40个PCA，前30个PC输入Harmony（sample_name，种子202642，最多20轮）。
- 20近邻，UMAP及Leiden种子202642；Leiden使用igraph无向图、2轮。
- 分辨率0.4、0.6、0.8分别得到7、9、11群。固定分辨率0.6换种子42，与主划分ARI=0.6206，一致性中等，边界有可变性；9群是探索性RNA状态划分。
- 归一化表达矩阵用于Wilcoxon标志基因排序（tie correction，BH调整），输出每群前50基因；不以PGAM5标签强制合并或移动坐标。

## 注释限制

名称描述RNA程序，不直接断言抗肿瘤/促肿瘤功能。M5应激程序可能受解离过程影响。M2存在ALB/APOA等肝源RNA富集，巨噬身份仍有C1Q/CD68支持，尚未进行环境RNA校正。M3只有51个细胞，其NPR3/MYCT1与血管相关表达需进一步核查，保留待定注释。M7同时具有C1QA和CD3D/CD3E/TRAC信号，保留混合RNA标签，不将其作为明确的新巨噬亚型或已证实双细胞。M6约73.4%来自同一患者。未进行正式双细胞过滤。未检出PGAM5不解释为蛋白阴性。

## 文件与复现

UMAP_named_macrophages.png/pdf为命名聚类图；UMAP_PGAM5.png/pdf为PGAM5逐细胞归一化表达和RNA检出图；UMAP_sample_and_integration.png/pdf比较样本组成及未整合UMAP。marker_dotplot颜色按基因在各群间缩放至0–1，点大小为检出比例；完整未缩放证据见marker_evidence.csv。

cell_annotations.csv.gz保留逐细胞新标签、原C群、3种分辨率标签及第二随机种子结果。cluster_summary.csv为群大小及PGAM5统计，cluster_sample_counts.csv为患者组成，top50_markers.csv.gz为标志基因排序。

inputs/macrophages_raw_full.h5ad为完整36,591基因原始稀疏计数，可从头复现；macrophages_marker_subset.h5ad是仅含标志基因的轻量注释对象，X为使用全库大小计算的log1p(CP10k)，counts层为原始计数，不宜按该子集总和重新归一化。其完整表达对象可运行脚本自行生成，未将超过100 MiB的完整输出上传GitHub。

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python analyze.py
python plot.py
```

检查包含：细胞名单及顺序与原C0/C2/C16完全一致、原始全库总计数一致、PGAM5计数与原注释逐细胞一致、群计数总和正确、降维坐标有限、PNG人工检查。无需额外下载原矩阵。
