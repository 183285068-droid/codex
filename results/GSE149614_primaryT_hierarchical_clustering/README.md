# GSE149614原发HCC肿瘤T样本：全细胞与巨噬细胞分层聚类、PGAM5表达

限定HCC01T–HCC10T这**10个原发HCC肿瘤样本**，对应10个作者患者编号，排除末尾N、P、L的癌旁、门静脉癌栓及淋巴结样本。原始输入为GEO提供的计数矩阵及updated元数据；34,414个肿瘤样本细胞均满足本轮检测基因数≥500、线粒体计数比例<20%的条件。本轮读取GEO提供的处理后计数矩阵，这些阈值没有额外删除细胞。分析包含25,712个唯一基因。源文件SHA256见[source_summary.json](source_summary.json)。

## 1. 全细胞降维聚类、RNA类型注释

得到**22个群C0–C21**，涵盖T/NK、B/浆细胞、巨噬、单核/DC、肥大、内皮、基质/血管壁、肝细胞/上皮样群及混合RNA群。注释使用作者celltype、当前群标志及谱系标志表达证据。C16为C1Q富集的MKI67/TOP2A增殖巨噬样群，作者标签存在分歧，名称带有限定。C15与C21明确标为混合RNA群。肝细胞/上皮样名称不是本轮独立确认的恶性身份，未另做CNV分析。

- [细胞类型UMAP](01_allcell_cell_types_UMAP.png) · [标志命名的22群UMAP](01_allcell_named_clusters_UMAP.png)
- [标志基因点图](01_allcell_marker_dotplot.png) · [作者标签](01_allcell_author_annotation_UMAP.png)
- [质控](00_HCC_QC.png) · [样本UMAP](01_allcell_samples_UMAP.png) · [样本类型构成](01_allcell_sample_composition.png)
- [整合前后对照](01_allcell_integration_comparison.png)
- [逐细胞类型表](allcell_cell_annotations.csv.gz) · [群名称和注释证据状态](allcell_cluster_annotations.csv)
- [每群前40标志](allcell_top40_markers.csv.gz) · [谱系标志表达证据](allcell_marker_evidence.csv)

## 2. 巨噬细胞独立重新降维聚类

采用保守谱系筛选：先在这10个T样本的作者Myeloid细胞中，按原作者res.3群重算各C1QA/B/C及CD68检出比例。各C1Q基因均≥80%、CD68≥75%的群为5、6、21、23、38、39、41、44、46，共6,151个细胞。再与本轮全细胞注释确认的巨噬群C3、C5、C16相交，提取**5,605个巨噬富集细胞**，覆盖全部10个T样本。其余546个旧规则候选分散于单核/DC、混合或其他RNA群，未纳入本轮主分析；因此数目与此前6,151个版本不同。

[原作者群的T样本谱系筛选证据](primaryT_original_myeloid_marker_gate.csv)和[排除Myeloid表](excluded_myeloid_cells.csv.gz)支持追溯。此筛选倾向于C1Q富集的成熟巨噬细胞，不能代表所有单核来源的髓系过渡状态。

独立重新选择高变基因、PCA、Harmony、邻居图及UMAP，得到**9个探索性RNA群M0–M8**，名称采用标志基因及相关RNA程序。

- [巨噬亚群命名UMAP](02_macrophage_named_clusters_UMAP.png)
- [巨噬及竞争谱系标志点图](02_macrophage_marker_dotplot.png)
- [巨噬样本UMAP](02_macrophage_samples_UMAP.png) · [整合前后对照](02_macrophage_integration_comparison.png)
- [逐细胞命名表](macrophage_cell_annotations.csv.gz) · [群名称](macrophage_cluster_annotations.csv)
- [前50标志](macrophage_top50_markers.csv.gz) · [群×样本细胞数](macrophage_cluster_sample_counts.csv)

M6显示C1Q及CD3D/CD3E/TRAC共同信号，标为巨噬/T混合RNA群，不能当作已确认的纯巨噬亚型。M3有91.6%细胞来自HCC01T；M4有71.3%来自HCC05T；M8仅42个细胞、功能待定。FOS/JUN及热应激表达也可能受到取材或解离过程影响。未额外进行双细胞或环境RNA校正，混合信号不等同于已证实双细胞。细胞群名称表示RNA标志程序，不能直接等同于经过实验确认的功能。

## 3. 各巨噬RNA群的PGAM5表达比较

**保留全部零值**。CP10k = PGAM5原始计数 / 原始完整转录组文库总计数 × 10,000。RNA检出为原始计数>0；未检出不解释为蛋白阴性。表达统计使用未经Harmony改写的RNA矩阵。PGAM5只从几何构建的高变基因名单中移除，原始计数、归一化表达、标志点图和下列表格均保留该基因。

- [逐细胞PGAM5表达及检出UMAP](03_PGAM5_expression_UMAP.png)：左侧自然对数log1p(CP10k)，右侧RNA检出；灰色为未检出。
- [各群PGAM5平均表达及检出率UMAP](03_PGAM5_cluster_mean_and_detection_UMAP.png)：每个点按其群汇总值着色，不是该细胞的实际表达。
- [平均表达和检出率比较条形图](03_PGAM5_cluster_comparison.png)
- [各群×样本的PGAM5均值图](03_PGAM5_cluster_sample_support.png)：空格表示缺少该群，点大小随细胞数增加。
- [精确群统计](PGAM5_by_macrophage_cluster.csv) · [群×样本统计](PGAM5_by_macrophage_cluster_sample.csv)

共136/5,605个细胞检出（2.43%），各群中位数均为0。**M1 MKI67/TOP2A增殖型均值最高，0.1070 CP10k；19/208检出（9.13%），覆盖全部10个T样本，最大单一样本占25.5%。** PGAM5检出并不限于增殖群，各群均有检出细胞。M8的较高均值仅对应1个检出细胞；M3受单一样本主导，不能据其排序确立稳定亚型。

| 群 | 暂定RNA名称 | 细胞数 | 平均PGAM5 CP10k（含零） | 检出率 | 注释状态 |
|---|---|---:|---:|---:|---|
| M0 | FOS/JUN即时早期反应巨噬细胞 | 1428 | 0.0297 | 1.75% (25/1428) | qualified |
| M1 | MKI67/TOP2A增殖型巨噬细胞 | 208 | 0.1070 | 9.13% (19/208) | supported |
| M2 | C1QC/HLA-DRA抗原呈递相关巨噬细胞 | 1287 | 0.0311 | 2.33% (30/1287) | supported |
| M3 | JUND/HSPA1A应激相关巨噬细胞（样本主导） | 523 | 0.0546 | 3.25% (17/523) | qualified |
| M4 | TIMD4/CD5L驻留样巨噬细胞 | 356 | 0.0291 | 1.97% (7/356) | qualified |
| M5 | SPP1/LGALS1重塑相关巨噬细胞 | 1477 | 0.0370 | 2.03% (30/1477) | supported |
| M6 | C1QA/CD3D巨噬/T混合RNA群 | 149 | 0.0435 | 3.36% (5/149) | mixed |
| M7 | MT1G/MT2A金属响应相关巨噬细胞 | 135 | 0.0103 | 1.48% (2/135) | supported |
| M8 | GRM4/AIM1L巨噬富集群（待定） | 42 | 0.0645 | 2.38% (1/42) | qualified |

本轮进行描述性表达比较，未进行亚群间推断性显著性检验，图中无p值。标志筛选的逐细胞Wilcoxon检验仅辅助注释，不能当作独立患者验证。聚类注释及可视化未机械沿用此前signature的FDR/log2FC阈值。

## 方法与复现核查

全细胞：Sample分层Seurat高变基因3,000个（随后移除PGAM5）、scale截断10、50PC，40PC构建20邻居图；Harmony按Sample整合，最多20轮，本轮14轮收敛；UMAP及Leiden分辨率0.8，种子149614。巨噬细胞独立重复流程，高变基因2,000、40PC、30PC构建20邻居图，Harmony本轮5轮收敛，Leiden分辨率0.6。主聚类采用igraph无向图、2次迭代。

固定邻居图改用种子42：全细胞ARI=0.857，巨噬ARI=0.853。巨噬分辨率0.4/0.6/0.8分别得到6/9/13群，细分边界存在敏感性。本轮编号与其他数据集或旧版本编号没有直接对应关系。标志点图颜色按每个基因在群间缩放到0–1，只能比较同一基因的群间表达；点大小为检出率。

[完整计数审计](full_data_audit_summary.json)、[完整逐项核查](full_data_independent_audit_checks.csv)、[下载矩阵核查](audit_summary.json)、[独立巨噬复跑验证](portable_reclustering_verification.json)检查细胞选择、原始计数、完整文库归一化、保存坐标、固定图重建及逐群/样本统计。图形已人工检查。PNG与PDF均提供。

## 下载与运行

[完整结果ZIP](GSE149614_primaryT_hierarchical_clustering_results.zip)：包括图、表、源信息、代码、全细胞标志子集矩阵、巨噬标志矩阵及完整25,712基因的巨噬计数+归一化矩阵。GitHub打开ZIP后点击 **Download raw file**。体积较大的官方计数文本、全细胞完整计数/归一化h5ad未打包；全细胞全流程重算需下载官方计数文件。

Python3.12，建议≥16 GiB内存并预留至少3 GiB磁盘。解压后：

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 NUMBA_NUM_THREADS=4
export MPLCONFIGDIR="$PWD/.matplotlib"
export XDG_CACHE_HOME="$PWD/.cache"
python audit.py
python visualize.py both
python visualize_sample_support.py
```

包内标志子集矩阵保留完整文库归一化分母，勿按子集总计数重新归一化。为减小体积，删除冗余距离和未整合邻居图，保留主connectivities及两套UMAP坐标；可重建主分区。完整巨噬矩阵为`macrophages_annotated_full.h5ad`，其中`layers['counts']`是原始计数、X是log1p(CP10k)。从它可直接重新运行巨噬聚类：

```bash
python cluster_macrophages.py
python annotate_macrophages.py
python visualize.py mac
python visualize_sample_support.py
```

重新计算全部流程：

```bash
python download_inputs.py
python prepare_counts.py
python cluster_all.py
python annotate_all.py
python cluster_macrophages.py
python annotate_macrophages.py
python visualize.py both
python visualize_sample_support.py
python audit.py
```

`annotate_*.py`记录本轮人工标志命名映射。改变输入、参数或软件版本后若编号变化，必须重新查看标志并注释，不能直接套用旧名称。源数据：[GEO GSE149614](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614)；官方URL及SHA256校验由download_inputs.py和source_summary.json固定。metadata中的原始Cell、sample、patient及res.3均保留，CSV以cell_id作为索引，便于回填分析对象。
