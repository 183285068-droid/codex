# GSE151530 HCC：全细胞与巨噬细胞分层聚类、PGAM5表达可视化

本轮沿用此前指定的**全部HCC肿瘤样本**范围，排除ICC。32个HCC样本共有50,023个细胞；以检测基因数≥500、线粒体计数比例<20%质控后，保留44,776个细胞，覆盖全部32个样本（25个患者编号代理）。先独立完成全细胞聚类，再核对作者标签与标志基因、提取巨噬细胞重新降维聚类。样本数不等于独立患者数。

## 第一步：全细胞降维、聚类与注释

得到**24个群C0–C23**，涵盖T、NK/T、B、浆细胞、巨噬、单核/DC样、pDC、内皮、基质、肿瘤/肝细胞样及待定混合RNA群。注释结合作者Type及本轮RNA标志；肿瘤身份沿用作者支持，本轮未另做CNV推断。

- [细胞类型UMAP](01_allcell_cell_types_UMAP.png)
- [24群标志命名UMAP](01_allcell_named_clusters_UMAP.png)
- [作者标签UMAP](01_allcell_author_annotation_UMAP.png)
- [标志基因点图](01_allcell_marker_dotplot.png)
- [质控图](00_HCC_QC.png)、[整合前后对照](01_allcell_integration_comparison.png)、[样本着色UMAP](01_allcell_samples_UMAP.png)
- [逐细胞注释](allcell_cell_annotations.csv.gz)、[群名称及证据状态](allcell_cluster_annotations.csv)、[标志表达证据](allcell_marker_evidence.csv)、[每群前40个标志基因](allcell_top40_markers.csv.gz)

## 第二步：巨噬细胞重新降维聚类

质控后作者TAMs标签有4,494个细胞。全细胞C8包含明显FCN1/S100A8和CD1C/FCER1A竞争髓系信号，不能将全部TAMs直接当作纯巨噬细胞。本轮主分析取**作者TAMs标签与C1Q富集的全细胞C9、C11相交的3,095个细胞**，覆盖28个样本、21个患者编号代理；其余1,399个作者TAMs保存在[排除表](excluded_author_TAMs.csv.gz)，可追溯。

重新聚类得到**16个探索性RNA群M0–M15**。本轮编号与此前数据集或旧聚类版本的编号没有直接对应关系。

- [标志命名巨噬细胞UMAP](02_macrophage_named_clusters_UMAP.png)
- [巨噬及竞争谱系标志点图](02_macrophage_marker_dotplot.png)
- [整合前后对照](02_macrophage_integration_comparison.png)、[样本构成UMAP](02_macrophage_samples_UMAP.png)
- [逐细胞命名表](macrophage_cell_annotations.csv.gz)、[每群前50个标志基因](macrophage_top50_markers.csv.gz)、[群样本构成](macrophage_cluster_sample_counts.csv)

M14/M15显示巨噬与T/NK混合RNA，M6为C1Q富集髓系边界群，M3/M12等仍待复核；表中`qualified`表示注释有限定条件，`mixed`表示混合RNA。因此不能将16群全部视为已确认的纯巨噬亚型。未进行额外双细胞检测或环境RNA校正；混合信号也不直接等同于已证实双细胞。

标志点图颜色按每个基因在群间缩放到0–1，点大小为检出比例；颜色用于该基因的群间比较，不能比较不同基因的绝对表达。

## 第三步：PGAM5表达比较

**保留所有零值**。PGAM5不参与高变基因几何构建，以减少按目标基因分群的循环解释，但保留在原始计数、归一化矩阵、标志表及全部表达比较中。

- [逐细胞PGAM5表达及RNA检出UMAP](03_PGAM5_expression_UMAP.png)：左侧log1p(CP10k)，右侧原始计数>0，灰色为未检出。
- [各群平均表达与检出率UMAP](03_PGAM5_cluster_mean_and_detection_UMAP.png)：每个点使用所在群的汇总值着色，不是该细胞的实际表达值。
- [各群PGAM5比较条形图](03_PGAM5_cluster_comparison.png)：左侧全群平均CP10k，右侧检出率及检出/总细胞数。
- [精确统计表](PGAM5_by_macrophage_cluster.csv)、[群×样本统计](PGAM5_by_macrophage_cluster_sample.csv)

CP10k = PGAM5原始计数 / 原始完整转录组总计数 × 10,000；UMAP表达色阶使用自然对数log1p。所有3,095个细胞参与组内均值，153个检出（4.94%）；各组中位数均为0。均值与检出率回答不同问题，例如M5仅2个检出细胞，但较高单细胞值仍影响组均值。

M3均值最高（0.1271 CP10k，5/26检出），但仅26个细胞、84.6%来自同一样本且功能待定，不能据此确立稳定PGAM5亚型。M2 MKI67/TOP2A增殖型均值0.1010 CP10k，22/171检出（12.87%），为本轮有较明确程序支持的群中均值最高。M8 SPP1/GPNMB均值0.0834 CP10k，19/264检出（7.20%）。这些是RNA表达描述，不是细胞浸润丰度或功能因果证据。

| 群 | 暂定RNA名称 | 细胞数 | 平均PGAM5 CP10k（含零） | 检出率 | 证据状态 |
|---|---|---:|---:|---:|---|
| M0 | FCGR3A/PSAP巨噬细胞 | 730 | 0.0570 | 7.26% | supported |
| M1 | RPL12/RPS27核糖体表达偏高巨噬细胞 | 686 | 0.0465 | 1.02% | qualified |
| M2 | MKI67/TOP2A增殖型巨噬细胞 | 171 | 0.1010 | 12.87% | supported |
| M3 | LINC01419/NOVA1巨噬富集群（待定） | 26 | 0.1271 | 19.23% | qualified |
| M4 | ALOX5AP/COTL1巨噬细胞 | 256 | 0.0579 | 8.59% | supported |
| M5 | CCL3/CCL4趋化因子相关巨噬细胞（伴肝源RNA） | 278 | 0.0523 | 0.72% | qualified |
| M6 | CLEC10A/FCN1 C1Q富集髓系边界群 | 170 | 0.0011 | 0.59% | qualified |
| M7 | FOLR2/SEPP1驻留样巨噬细胞 | 235 | 0.0370 | 2.98% | supported |
| M8 | SPP1/GPNMB巨噬细胞 | 264 | 0.0834 | 7.20% | supported |
| M9 | TIMD4/LYVE1驻留样巨噬细胞 | 39 | 0.0542 | 2.56% | qualified |
| M10 | CXCL9/CXCL10干扰素相关巨噬细胞 | 60 | 0.0520 | 8.33% | supported |
| M11 | CHIT1/CCL18巨噬细胞 | 44 | 0.0749 | 4.55% | supported |
| M12 | HOXA13/UBE2C增殖混合RNA群 | 23 | 0.0000 | 0.00% | qualified |
| M13 | HLA-A/HLA-C巨噬细胞（患者相关） | 32 | 0.0469 | 3.12% | qualified |
| M14 | C1QA/CD3D巨噬/T细胞混合RNA群 | 50 | 0.0331 | 12.00% | mixed |
| M15 | C1QA/EOMES巨噬/NK混合RNA群 | 31 | 0.0000 | 0.00% | mixed |

本轮完成表达描述比较，未进行亚群表达差异的推断性检验，图中无p值。标志筛选使用逐细胞Wilcoxon辅助注释，其p值不代表独立患者验证。本轮未沿用此前signature筛选的FDR/log2FC阈值，因为这里的主要任务是聚类注释及可视化。

## 方法与核查

原始基因名重复项求和，18,667特征合并为18,661个唯一基因。全细胞采用按Sample的3,000个Seurat高变基因、截断scale=10、50PC；40PC构建20邻居图。按Sample进行Harmony，Leiden分辨率0.8。巨噬细胞重新选择2,000个高变基因、40PC，30PC构建20邻居图，Harmony后Leiden分辨率0.6。两级UMAP及主聚类随机种子151530，igraph无向图、2次迭代。原始归一化RNA用于表达比较，不用Harmony坐标计算表达。

固定图上改用种子42：全细胞ARI=0.908，巨噬细胞ARI=0.777。巨噬分辨率0.4/0.6/0.8分别得到13/16/17群，提示精细边界仍有敏感性，名称为探索性注释。极少细胞样本仍纳入批次高变特征选择，可能产生小批次相关警告，未排除这些样本。

完整源数据独立核查**30项通过**：细胞选择、计数与总计数、完整文库归一化、坐标导出、保存图重建聚类ARI=1以及16群表达汇总。见[完整数据审计](full_data_audit_summary.json)、[逐项检查](full_data_independent_audit_checks.csv)。下载包内可运行29项核查；完整全细胞原始计数和矩阵因体积未打包，重跑全流程后可执行第30项。包内原始巨噬输入已在独立目录重新运行完整巨噬聚类流程，逐细胞标签、归一化RNA和UMAP坐标均与本轮结果完全一致（ARI=1），见[独立复跑核查](portable_reclustering_verification.json)。图形已人工查看，PNG与PDF均提供。

## 下载及复现

[完整结果ZIP](GSE151530_HCC_hierarchical_clustering_results.zip)（约92 MiB）。在GitHub打开ZIP后点击 **Download raw file**。完整全细胞原始/归一化h5ad和约291 MiB原始矩阵未纳入ZIP；包内有坐标、逐细胞标签、标志矩阵、完整巨噬归一化+计数矩阵及巨噬原始计数，支持重新绘图和巨噬独立重聚类。

Python 3.12；建议新建虚拟环境，按requirements安装。在解压目录内运行：

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 NUMBA_NUM_THREADS=4
export MPLCONFIGDIR="$PWD/.matplotlib"
python audit.py
python visualize.py both
```

`allcell_marker_subset.h5ad`只含标志子集，其归一化分母仍为完整文库，勿用子集总计数重新归一化。为节省体积，h5ad删去冗余距离及未整合图，保留主connectivities和两套UMAP坐标，仍可重建主分区。`macrophages_annotated_full.h5ad`保存完整18,661基因计数和归一化表达。

只重新聚类当前已选择的巨噬细胞，可直接：

```bash
python cluster_macrophages.py
python annotate_macrophages.py
python visualize.py mac
```

从官方矩阵重新计算全部流程（建议≥32 GiB内存，额外预留约2 GiB磁盘）：

```bash
python download_matrix.py
python cluster_all.py
python annotate_all.py
python cluster_macrophages.py
python annotate_macrophages.py
python visualize.py both
python audit.py
```

`annotate_*.py`记录本轮人工标志注释对应关系。若版本、随机算法或输入改变导致聚类编号变化，必须重新查看标志并核对名称，不能机械套用旧编号。当前按28个样本处理巨噬批次，不假设28个独立患者。

原始数据：[GEO GSE151530](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151530)。矩阵官方下载地址见download_matrix.py；当前复用已上传/缓存矩阵，SHA256核实为`50dad60e779b3332350e2552ed4260369ccba8db4f73918bbe2a9f6549db6403`。基因、条码、作者信息及GEO诊断映射已打包于inputs；元数据保留作者细胞编号与原始矩阵列编号，可回填Seurat/AnnData。
