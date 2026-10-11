# GSE281759：SPP1⁺巨噬细胞相关RNA位点与T细胞的空间交互

本轮复用此前**4位HCC患者、16张肿瘤切片**（15张原发、1张复发）的主QC、真实Visium空间坐标、巨噬/T细胞RNA门控。没有把位点或切片当作独立患者，没有改动前轮PGAM5结果。

**可观察到空间邻近及候选SPP1–CD44共检出，但当前没有通过多重校正的增强邻接或T耗竭相关证据。**2,291个SPP1检出巨噬RNA位点中，1,500个有T-RNA直接邻居，300个有耗竭支持T-RNA邻居，81个有早期活化支持T-RNA邻居。主All-T邻接检验的最小BH q=0.2624；全部176项计划检验的最小联合q=0.8096，均未达到0.05。无显著结果不等于没有生物作用。

## 图表入口

- [SPP1检出与T邻接比例总览](figures/SPP1_detection_T_neighborhood_overview.png)
- [空间邻接差异及q值](figures/SPP1_T_adjacency_statistics.png)：主检测状态、耗竭、活化、排除T阳性源、SPP1高表达敏感性。
- [SPP1与邻近T程序的校正相关及q值](figures/SPP1_neighbor_T_program_association.png)
- [4位患者的描述性效应汇总](figures/patient_descriptive_effect_summary.png)
- [候选配体受体RNA共检出](figures/SPP1_mac_T_candidate_LR_codetection.png)
- 原组织图叠加：[A009](figures/A009_SPP1_mac_T_spatial_neighborhoods.png) · [B013](figures/B013_SPP1_mac_T_spatial_neighborhoods.png) · [C004原发/复发](figures/C004_SPP1_mac_T_spatial_neighborhoods.png) · [C008](figures/C008_SPP1_mac_T_spatial_neighborhoods.png)。每片四列：全位点SPP1、巨噬RNA门内SPP1、T-RNA状态、SPP1检出巨噬RNA源到T邻居。

全部9组图均有PNG和PDF。空间连线是**直接位点邻接，无自环**，不是已证实的细胞接触、信号方向或配体传递。共享SPP1色标；跨片Q99最大值用于显示上限，超限只绘图截断，统计使用原值。灰色检验格为不可检验，不能读作阴性。图上q为联合176项BH；精确p和预先单列的主16项BH见CSV。

## 来源、范围和质控

控制数据源为[NCBI GEO GSE281759](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281759)的官方计数及Visium附件。复用已发表到本仓库的[固定版本基础分析](https://github.com/183285068-droid/codex/tree/6bb037700d74331b25f5ab5c21d167f178927fb3/results/GSE281759_spatial_HCC_PGAM5)：17个必要输入（16个marker AnnData和范围表）的字节数、SHA256和固定提交下载URL保存在[reference_manifest.json](reference_manifest.json)。本轮再次核验全部引用文件。

排除A009N邻近非肿瘤切片。GEO的B0013T1–T4与B013T5统一患者标签B013，保留原样本名。C004T1P为原发、C004T1R为复发，另提供只含15张原发的患者描述与主BH15，最小q=0.246，仍不显著。所有切片状态见[sample_scope.json](sample_scope.json)。

沿用主QC：in_tissue=1、UMI≥500、检出基因≥200、MT计数比例≤25%；25,510个输入位点保留25,354。原论文MT≤5%规则会显著减少某些切片覆盖；本轮没有重新完成MT5%下的所有空间检验，因此不能声称结果通过该QC敏感性验证。真实组织位置、缩放系数、H&E和图结构均来自原数据。SPP1归一化使用原**36,601个完整feature的总计数**作CP10k，再log1p；不能从便携marker子集重新计算分母。

## RNA定义和空间邻接

- 巨噬RNA门：C1QA/B/C、CD68、CSF1R、MERTK、MSR1均值≥切片Q75，≥3核心基因检出；LST1/TYROBP/FCER1G/AIF1/LYZ/CTSS均值≥切片中位数。
- **SPP1⁺相关源：上述巨噬RNA门内SPP1原始计数>0。**原始计数0称“未检出”，不等于RNA完全不存在或蛋白阴性。巨噬门不含SPP1，避免用被检验基因决定谱系门。
- SPP1-high：上述检出源内，归一化SPP1≥切片检出源Q75（边界并列保留）。共578个；阈值逐切片固定，不是跨片统一绝对丰度。主分析是所有SPP1检出源，高表达是敏感性分析。
- T-RNA门：CD3D/E/G、TRAC、TRBC1/2均值≥切片Q75，≥2谱系基因检出。4,004个位点满足此门。
- 耗竭支持：T门内PDCD1/LAG3/TIGIT/CTLA4/TOX/CXCL13均值≥T门Q75、>0、≥2程序基因检出；465个位点。早期活化用CD69/IL2RA/TNFRSF9/ICOS/CD40LG/IL2同类门；74个位点。两者可以重叠，不是互斥的纯T细胞群。
- 去CTLA4/TIGIT耗竭敏感性程序为PDCD1/LAG3/TOX/CXCL13/HAVCR2/ENTPD1；另有扩展耗竭及细胞毒效应程序。固定基因列表见[programs.json](programs.json)和[analysis_config.json](analysis_config.json)；程序均值包含零值。去除两基因也不能保证彻底区分Treg与耗竭T细胞。
- 邻接用Visium六边格第一圈：(0,±2)、(±1,±1)，主QC后的真实直接邻居，每源最多6个，没有自环。主邻接检验和比例图排除空间度数0的源；全位点描述计数仍保留它们。邻近T程序均值仅取T门邻居，含表达零值；无T邻居为缺失，不能置0。不是UMAP距离。
- nonT敏感性排除同时满足T-RNA门的巨噬RNA源，减轻同一混合位点T成分的影响；不能清除所有混合RNA或确定单一细胞来源。

## 患者和切片数量

以下是位点数，不是细胞数或经过反卷积校准的浸润丰度。源可同时有活化、耗竭邻居，不能相加作互斥总数。

| patient | sections | mac_RNA_spots | SPP1_mac_detected | SPP1_detected_with_T_neighbor | SPP1_detected_with_Tex_neighbor | SPP1_detected_with_Tact_neighbor |
| --- | --- | --- | --- | --- | --- | --- |
| A009 | 5 | 1360 | 615 | 476 | 189 | 14 |
| B013 | 5 | 1384 | 1112 | 715 | 69 | 54 |
| C004 | 2 | 223 | 34 | 23 | 0 | 0 |
| C008 | 4 | 1124 | 530 | 286 | 42 | 13 |

| sample | tumor_status | mac_RNA_spots | SPP1_mac_detected | SPP1_mac_high | SPP1_detected_with_T_neighbor | SPP1_detected_with_Tex_neighbor | SPP1_detected_with_Tact_neighbor |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A009T1 | primary | 301 | 139 | 35 | 114 | 54 | 1 |
| A009T2 | primary | 128 | 90 | 23 | 48 | 10 | 0 |
| A009T3 | primary | 93 | 60 | 15 | 42 | 3 | 0 |
| A009T4 | primary | 494 | 217 | 55 | 189 | 90 | 7 |
| A009T5 | primary | 344 | 109 | 28 | 83 | 32 | 6 |
| B0013T1 | primary | 199 | 162 | 41 | 110 | 1 | 7 |
| B0013T2 | primary | 295 | 295 | 74 | 178 | 37 | 18 |
| B0013T3 | primary | 182 | 150 | 38 | 74 | 1 | 0 |
| B0013T4 | primary | 350 | 150 | 38 | 113 | 16 | 5 |
| B013T5 | primary | 358 | 355 | 89 | 240 | 14 | 24 |
| C004T1P | primary | 146 | 22 | 6 | 14 | 0 | 0 |
| C004T1R | recurrent | 77 | 12 | 3 | 9 | 0 | 0 |
| C008T1 | primary | 476 | 40 | 10 | 34 | 21 | 10 |
| C008T2 | primary | 216 | 119 | 30 | 91 | 7 | 0 |
| C008T3 | primary | 144 | 103 | 26 | 13 | 0 | 0 |
| C008T4 | primary | 288 | 268 | 67 | 148 | 14 | 3 |

SPP1在所有主QC位点中有12,711个检出，其中10,420个位于巨噬RNA门之外（约82%）；不能把组织SPP1都归因于巨噬细胞。门外也可能有掉落、低巨噬RNA或混合位点，不能反向确认其来自某一其他谱系。

## 邻接统计：SPP1检出与未检出源

统计量是两组中“至少有1个直接T-RNA邻居”的源比例差，单位**百分点**（SPP1检出减未检出），不是细胞丰度差。下表原始比例与检验分母一致，排除空间孤立源；p来自条件置换。

| sample | n_positive | n_comparator | positive_neighbor_pct | comparator_neighbor_pct | statistic | permutation_p | primary_BH_q_16_sections | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A009T1 | 139 | 162 | 82.0144 | 85.1852 | -3.1708 | 0.4456 | 1.0 | tested |
| A009T2 | 90 | 38 | 53.3333 | 36.8421 | 16.4912 | 0.9244 | 1.0 | tested |
| A009T3 | 60 | 33 | 70.0 | 69.697 | 0.303 | 1.0 | 1.0 | tested |
| A009T4 | 217 | 277 | 87.0968 | 91.3357 | -4.239 | 0.5748 | 1.0 | tested |
| A009T5 | 109 | 235 | 76.1468 | 75.7447 | 0.4021 | 1.0 | 1.0 | tested |
| B0013T1 | 162 | 37 | 67.9012 | 78.3784 | -10.4771 | 0.2604 | 1.0 | tested |
| B0013T2 | 294 | 0 | 60.5442 | — | — | — | 1.0 | not_testable |
| B0013T3 | 148 | 32 | 50.0 | 81.25 | -31.25 | 0.0328 | 0.2624 | tested |
| B0013T4 | 150 | 200 | 75.3333 | 63.0 | 12.3333 | 0.0308 | 0.2624 | tested |
| B013T5 | 354 | 3 | 67.7966 | 33.3333 | 34.4633 | — | 1.0 | not_testable |
| C004T1P | 22 | 124 | 63.6364 | 60.4839 | 3.1525 | 1.0 | 1.0 | tested |
| C004T1R | 12 | 65 | 75.0 | 66.1538 | 8.8462 | 0.7164 | 1.0 | tested |
| C008T1 | 40 | 436 | 85.0 | 86.9266 | -1.9266 | 1.0 | 1.0 | tested |
| C008T2 | 117 | 97 | 77.7778 | 80.4124 | -2.6346 | 0.9528 | 1.0 | tested |
| C008T3 | 103 | 41 | 12.6214 | 12.1951 | 0.4262 | 1.0 | 1.0 | tested |
| C008T4 | 268 | 19 | 55.2239 | 47.3684 | 7.8555 | 0.4432 | 1.0 | tested |

B0013T2巨噬RNA源全部检出SPP1，没有未检出对照；B013T5只有3个未检出源，未达到至少10个对照，明确不可比较。其他恒定T状态、缺乏分层内可交换标签或稀疏群亦标为不可检验，未伪造p值。14/16张切片可检验主All-T邻接，11/16可检验主耗竭邻接，8/16可检验主活化邻接。

## 归一化SPP1与邻近T耗竭程序

在巨噬RNA源且至少有1个T-RNA邻居的位点中保留SPP1零值，偏Spearman调整9项秩协变量：源logUMI、MT%、巨噬核心、肝细胞、上皮、T谱系；T邻居平均T谱系、logUMI、MT%。效应是条件相关，不是表达改变的因果量。

| sample | n_sources | statistic | permutation_p | BH_q_176_planned_tests | status |
| --- | --- | --- | --- | --- | --- |
| A009T1 | 252 | -0.0487 | 0.2936 | 1.0 | tested |
| A009T2 | 62 | -0.0249 | 0.91 | 1.0 | tested |
| A009T3 | 65 | -0.0541 | 0.652 | 1.0 | tested |
| A009T4 | 442 | -0.003 | 0.2112 | 1.0 | tested |
| A009T5 | 261 | -0.0151 | 0.7332 | 1.0 | tested |
| B0013T1 | 139 | 0.0109 | 0.6112 | 1.0 | tested |
| B0013T2 | 178 | -0.208 | 0.0092 | 0.8096 | tested |
| B0013T3 | 100 | 0.018 | 0.8788 | 1.0 | tested |
| B0013T4 | 239 | 0.0176 | 0.8092 | 1.0 | tested |
| B013T5 | 241 | -0.029 | 0.2876 | 1.0 | tested |
| C004T1P | 89 | 0.1265 | 0.0444 | 1.0 | tested |
| C004T1R | 52 | 0.0404 | 0.4872 | 1.0 | tested |
| C008T1 | 413 | 0.0586 | 0.3468 | 1.0 | tested |
| C008T2 | 169 | -0.1348 | 0.19 | 1.0 | tested |
| C008T3 | 18 | — | — | 1.0 | not_testable |
| C008T4 | 157 | 0.1406 | 0.1432 | 1.0 | tested |

B0013T2的耗竭核心关联为rho=-0.2080，名义p=0.0092，但联合q=0.8096；不能据此声称SPP1抑制T耗竭。各切片有正有负，患者汇总只取可检验切片效应的等权中位数；不把位点合并扩大患者样本量，不报告4患者的伪精确总体p或置信区间。

## 置换、敏感性和多重检验

每项4,999次置换，固定原目标门及邻居结果，只在**全切片坐标3×3区块×当前分析源UMI二分位×巨噬核心二分位**内打乱SPP1表达/标签。两侧等尾p带+1校正，允许条件零分布均值非0；没有假定零分布以0为中心。最少30源；检测与未检出各≥10，高表达与其余检出各≥5/10；连续关联剩余自由度≥15。退化零分布不报告显著性。

计划176项：每切片5个邻接检验（全T、耗竭、活化、nonT耗竭、高表达vs其余检出耗竭）和6个连续关联检验（5套T程序、nonT耗竭核心）。145项可检验，31项不可检验；后者p保留缺失，仅在BH时按1计入计划家族。主All-T邻接另列16项BH，以及15张原发主检验BH。阈值/程序/家族在本轮计算前写入配置，没有据结果择优调整。

分层置换依赖层内近似可交换性；粗区块不能完全消除空间自相关、共享T目标、组织结构及残余混合RNA。9协变量用于连续关联，邻接检验仅使用上述分层。缺少病理ROI，未单独识别肿瘤细胞区域或MVI区域。当前结果支持病例内RNA邻域探索，不是跨患者已验证的通信推断。

## SPP1–CD44及候选免疫调节通路

固定8个候选来自前轮LIANA 1.10.0资源快照：SPP1–CD44、CD274–PDCD1、PVR–TIGIT、LGALS9–HAVCR2、TGFB1–TGFBR1/TGFBR2、IL10–IL10RA/IL10RB、CXCL9–CXCR3、CD80–CD28。仅统计巨噬RNA源到直接T-RNA邻居的原始表达共检出；多亚基受体须全部检出。几何CP10k分数是自定义描述，**没有运行CellChat/LIANA通信算法，没有配体受体显著性或结合/抑制证据**。

SPP1检出源到所有T邻居共有2,949条有向边，其中1,429条同时检出目标CD44；耗竭支持目标有381条边，其中186条检出CD44。重复源/目标不能当作独立观察，不能用这些边数作患者样本量。因为源已按SPP1检出选择，SPP1–CD44共检出比例实际主要反映目标CD44检出；CD44并非T细胞专属。下表是耗竭支持目标子集，百分比分母为eligible_edges；无边为缺失，不是0。

| sample | eligible_edges | unique_sources | unique_targets | supported_edges | supported_sources | supported_targets | edge_codetection_pct |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A009T1 | 72 | 54 | 40 | 24 | 22 | 12 | 33.33 |
| A009T2 | 10 | 10 | 8 | 0 | 0 | 0 | 0.0 |
| A009T3 | 3 | 3 | 3 | 0 | 0 | 0 | 0.0 |
| A009T4 | 119 | 90 | 70 | 38 | 34 | 23 | 31.93 |
| A009T5 | 43 | 32 | 35 | 12 | 12 | 7 | 27.91 |
| B0013T1 | 1 | 1 | 1 | 0 | 0 | 0 | 0.0 |
| B0013T2 | 38 | 37 | 10 | 37 | 36 | 9 | 97.37 |
| B0013T3 | 1 | 1 | 1 | 0 | 0 | 0 | 0.0 |
| B0013T4 | 16 | 16 | 9 | 11 | 11 | 5 | 68.75 |
| B013T5 | 18 | 14 | 7 | 18 | 14 | 7 | 100.0 |
| C004T1P | 0 | 0 | 0 | 0 | 0 | 0 | — |
| C004T1R | 0 | 0 | 0 | 0 | 0 | 0 | — |
| C008T1 | 38 | 21 | 30 | 33 | 17 | 27 | 86.84 |
| C008T2 | 8 | 7 | 3 | 7 | 6 | 2 | 87.5 |
| C008T3 | 0 | 0 | 0 | 0 | 0 | 0 | — |
| C008T4 | 14 | 14 | 4 | 6 | 6 | 2 | 42.86 |

checkpoint受体也与T状态门共享基因，相关共检出不是独立的耗竭验证。[candidate_LR_codetection_summary.csv](candidate_LR_codetection_summary.csv)提供全部源组/目标状态/通路及唯一源、唯一目标分母；各片边表可按barcode回填。未将候选通路RNA共检出解读为细胞发送/接收功能。

## 验证及复现

本轮独立审核**10,240项通过**：[audit_results.json](audit_results.json)。包括17个引用文件字节/SHA、SPP1全库归一化和全部RNA程序、原巨噬/T门、Q75和并列、独立KD-tree重建六边邻接、逐源显式邻居均值、所有统计分母、独立OLS秩残差相关及首3次置换、等尾p、三种BH、患者计数/等权中位数、全部候选复合受体共检出和9对PNG/PDF。人工复核统计热图及真实H&E空间图。

结果包：[GSE281759_SPP1_macrophage_T_interaction_results.zip](GSE281759_SPP1_macrophage_T_interaction_results.zip)。原输入AnnData含图像、marker原始计数及全库归一化表达，依据固定引用清单下载；包内不重复所有原对象。

```bash
python -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
python fetch_reference.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python analyze.py
OPENBLAS_NUM_THREADS=1 python plot.py
OPENBLAS_NUM_THREADS=1 python audit.py
python build_report.py
python package_results.py
```

脚本目录旁将建立`GSE281759_spatial_HCC_PGAM5`引用目录；已有文件只校验、不覆盖。结果CSV、注释及审核数组由脚本重建。固定版本基准为4患者/16片，不需要重新降维或按SPP1标签重排空间位置。
