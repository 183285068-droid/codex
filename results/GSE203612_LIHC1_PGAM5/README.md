# GSE203612：NYU_LIHC1原发HCC的单细胞与空间RNA分析

已完成**1位人类原发HCC患者**的inDrop单细胞和配对Visium分析。6,000个原始单细胞条码中保留1,783个，分为14个暂定RNA群；1,661个Visium位点中保留1,634个，分为3个混合RNA区域。

**PGAM5结论：**所有单细胞中66个检出PGAM5（3.70%）；C3/C6巨噬细胞相关群共325个细胞、8个检出。进一步联合巨噬谱系RNA门控后，246个候选细胞中仅6个检出（2.44%），**均为1个原始计数**。这些检出标签不能证明稳定PGAM5⁺巨噬亚群。66个检出细胞中39个来自C0/C1肝细胞样群，说明PGAM5并非巨噬细胞谱系专属。

空间数据有71/1,634个位点检出PGAM5（4.35%）。以163个有T-RNA邻居的巨噬RNA源位点作主关联检验，保留PGAM5零值并调整9项协变量：偏Spearman rho=0.00393，空间置换p=0.7804，7项联合BH q=1。去除CTLA4/TIGIT的耗竭程序也未显著。**当前病例没有显示PGAM5巨噬RNA邻域与T耗竭程序的统计可靠正关联；不显著也不等于没有作用。**

## 图表

- [单细胞14群、PGAM5表达及3文库UMAP](figures/singlecell_UMAP_PGAM5_libraries.png) · [PDF](figures/singlecell_UMAP_PGAM5_libraries.pdf)
- [单细胞标志基因dotplot](figures/singlecell_marker_dotplot.png) · [PDF](figures/singlecell_marker_dotplot.pdf)
- [各群PGAM5检出与巨噬候选定位](figures/PGAM5_detection_and_macrophage_candidates.png) · [PDF](figures/PGAM5_detection_and_macrophage_candidates.pdf)
- [Visium真实组织坐标九面板图](figures/Visium_spatial_RNA_maps.png) · [PDF](figures/Visium_spatial_RNA_maps.pdf)：区域、PGAM5、巨噬、T、耗竭、去CTLA4/TIGIT耗竭、早期活化、Treg及IFN程序。
- [Visium表达UMAP及dotplot](figures/Visium_UMAP_marker_dotplot.png) · [PDF](figures/Visium_UMAP_marker_dotplot.pdf)
- [PGAM5与邻近T程序关联及p/q](figures/Visium_PGAM5_T_exhaustion_association.png) · [PDF](figures/Visium_PGAM5_T_exhaustion_association.pdf)
- [PGAM5巨噬RNA源与耗竭T-RNA空间邻域](figures/Visium_PGAM5_Tex_spatial_neighborhoods.png) · [PDF](figures/Visium_PGAM5_Tex_spatial_neighborhoods.pdf)

Dotplot点面积为原始计数>0的比例，颜色是**包含零值的实际平均log1p(CP10k)**，不是按基因z分数。单细胞图显示每细胞实际表达；空间程序图显示实际程序均值，不是把均值回填成每细胞值。空间颜色各面板独立Q99显示截断，统计未截断。UMAP是表达嵌入，不是组织坐标。

## 纳入范围和来源

| GEO编号 | 标识 | 类型 | 数量 |
| --- | --- | --- | --- |
| GSM6177611 | NYU_LIHC1 | inDrop单细胞 | 同一病例3个文库，各2,000个原始条码 |
| GSM6177612 | NYU_LIHC1_Vis | 10x Visium | 同一病例1张空间切片 |

GSE203612是泛癌研究，其他人类癌种、黑色素瘤PDX及小鼠模型全部排除；小鼠`Exp3_Liver`是胰腺癌肝脏植入/转移相关模型，**不是人类原发肝癌病例**。两个纳入记录和论文均指向LIHC NYU1，文库/细胞/空间位点不能替代独立患者数。本文没有跨病例重复验证、生存分析或临床分期推断。

13个官方附件完整取得，共22.59MiB，长度与GEO官方文件清单一致，下载SHA256已记录：[source_manifest.json](source_manifest.json)。包含3组单细胞矩阵+基因表，以及Visium filtered H5、tissue_positions、scalefactors、低/高分辨图像及对齐图。

单细胞MatrixMarket注释行含**原始基因顺序和原始细胞ID**，已恢复；不是无barcode文件而人为生成未知对应的细胞。不同文库同名条码加lib前缀，例如`lib1|cell_134-323`，保留`original_cell_id`可回填。单细胞19,936个feature主要为原处理基因集，Visium33,538个feature；两平台各自使用**所提供完整基因矩阵总计数**作CP10k，再log1p。由于平台和基因覆盖差异，不能直接把跨平台均值差解释为表达增减。重复symbol保留唯一feature ID，PGAM5每个平台均单一feature。

## 质控和方法差异

主QC固定为UMI≥500、检出基因≥200、MT计数比例≤25%、核糖体基因比例≤30%；空间另要求in_tissue=1。MT使用`MT-`，核糖体使用`^RP[SL]\d`。原论文采用UMI≥500、MT/核糖体比例≤30%，使用Seurat SCTransform及其他参考流程；本轮增加基因200和更严MT25，并使用CP10k/Scanpy，**不是原论文的精确复现**。

| library | input_cells | qc_cells | qc_author_rule_cells | median_input_UMI |
| --- | --- | --- | --- | --- |
| lib1 | 2000 | 691 | 706 | 339.5 |
| lib2 | 2000 | 540 | 554 | 292.0 |
| lib3 | 2000 | 552 | 559 | 208.0 |

空间主QC保留1,634；原论文规则计数为1,635。单细胞原论文规则共1,819；主QC保留1,783。原始条码大部分低于500 UMI，不能通过放宽到极低测序深度来扩大PGAM5候选数并假定身份更可靠。逐条码QC表保留，未删除输入原始文件。

分群每平台独立：2,000高变基因；单细胞至少5细胞、空间至少10位点检出，HVG排除MT和典型核糖体基因；标准化clip10、PCA40、30PC/15邻居、Leiden0.6、UMAP，种子203612。没有依据PGAM5标签强制形成岛，也没有使用批次校正人为生成群。3个文库按同病例合并、以文库UMAP检查分布；未将其作为3位患者。Wilcoxon/top100标志用于RNA群解释，属于病例内探索，不能作为跨患者重复亚群检验。

## 单细胞RNA群注释

| cluster | RNA_annotation | n_cells | PGAM5_detected | PGAM5_detection_pct | PGAM5_mean_log1p_CP10k |
| --- | --- | --- | --- | --- | --- |
| C0 | C0: HRG/CYP3A5 &#124; Hepatocyte-like | 214 | 15 | 7.0093 | 0.1058 |
| C1 | C1: CYP2E1/CYP3A4 &#124; Hepatocyte-like | 316 | 24 | 7.5949 | 0.0987 |
| C2 | C2: MZB1/JCHAIN &#124; Plasma-cell-related | 103 | 4 | 3.8835 | 0.0889 |
| C3 | C3: C1QA/CD163 &#124; Macrophage-like | 221 | 4 | 1.81 | 0.0439 |
| C4 | C4: MKI67/TOP2A &#124; Cycling; lineage unresolved | 45 | 1 | 2.2222 | 0.053 |
| C5 | C5: FCN1/VCAN &#124; Monocyte/DC-related | 147 | 2 | 1.3605 | 0.0219 |
| C6 | C6: CD5L/MARCO &#124; Macrophage-like | 104 | 4 | 3.8462 | 0.0449 |
| C7 | C7: VWF/KDR &#124; Endothelial-like | 112 | 3 | 2.6786 | 0.0492 |
| C8 | C8: IL7R/CCL5 &#124; T/NK-related | 375 | 7 | 1.8667 | 0.0534 |
| C9 | C9: FCGR3B/S100A8 &#124; Granulocyte-like | 43 | 0 | 0.0 | 0.0 |
| C10 | C10: COL1A2/NOTCH3 &#124; Stromal/perivascular-related | 16 | 0 | 0.0 | 0.0 |
| C11 | C11: FOXP3/IKZF2 &#124; Treg-related | 46 | 1 | 2.1739 | 0.0648 |
| C12 | C12: MS4A1/CD79A &#124; B-cell-related | 31 | 1 | 3.2258 | 0.0772 |
| C13 | C13: CPA3/TPSB2 &#124; Mast-cell-related | 10 | 0 | 0.0 | 0.0 |

C3：C1QA/CD163，C6：CD5L/MARCO，支持巨噬细胞相关身份；C5为FCN1/VCAN单核/DC连续状态。C11具有FOXP3/IKZF2/TIGIT标志，暂定Treg相关；C8为IL7R/CCL5的T/NK相关群，未将该混合群强行全部判为耗竭T细胞。C4虽然MKI67/TOP2A增殖信号明显，但混合谱系信号不能确定来源，**未直接命名为增殖型巨噬细胞**。

C0/C1是肝细胞样RNA群，不因来自肿瘤组织就把每个细胞判为恶性。未运行CNV/SingleR、没有官方逐细胞恶性标签；10–16细胞的小群注释亦仅暂定。名称依据公开标志和程序复核，映射表为[RNA_group_annotations.csv](RNA_group_annotations.csv)，规则为[annotation_rules.json](annotation_rules.json)。

保守逐细胞巨噬门：至少2个C1QA/B/C、CD68、CSF1R、MERTK、MSR1检出，至少2个泛髓系基因检出，并限制T/B/DC和强肝细胞信号。初始门288个RNA候选，进一步与C3/C6巨噬相关群交集为246个，6个PGAM5检出。C3/C6全部325细胞中另2个PGAM5检出未通过更严门；结果同时报告，而不是偷偷替换分母。门控仍可能受dropout/ambient RNA影响，不是经过独立验证的细胞百分比。

## 空间RNA区域

| cluster | RNA_annotation | n_spots | PGAM5_detected | PGAM5_detection_pct | PGAM5_mean_log1p_CP10k |
| --- | --- | --- | --- | --- | --- |
| S0 | S0: CYP3A4/TF &#124; Hepatocyte-RNA-rich mixed region | 720 | 13 | 1.8056 | 0.0353 |
| S1 | S1: DCXR/GLUL &#124; Metabolic hepatocyte-RNA mixed region | 399 | 15 | 3.7594 | 0.0633 |
| S2 | S2: CD74/HLA-DRB1 &#124; Immune/stromal-RNA mixed region | 515 | 43 | 8.3495 | 0.0935 |

S0/S1以肝细胞/代谢RNA为主，S2显示较强CD74/HLA及免疫/间质混合RNA。是真实组织坐标上的表达区域，不能解释为纯细胞群、克隆、确定恶性区域或MVI区域。没有独立病理ROI标签。

空间巨噬RNA门：核心分数≥切片Q75、至少3核心基因检出，泛髓系分数≥切片中位数；共188个位点，其中21个检出PGAM5，6个属于检出位点表达最高25%（分位包含边界并保留原门阈值）。T-RNA门为谱系分数≥Q75且≥2谱系基因检出，共223个位点。

活化/耗竭支持状态须在T门内，分数≥T门Q75、>0且≥2程序基因检出，允许重叠。耗竭支持14、活化支持9。43个巨噬RNA源有耗竭支持T-RNA直接邻居；**本轮6个PGAM5高源中没有这样的直接邻居**。高源有任意T邻居的仅4个，因此不能把零共定位作为PGAM5抑制耗竭的证据。完整邻域和候选配体受体共检出表公开；21个固定LIANA资源候选仅为自定义RNA共检出描述，未运行CellChat/LIANA通信推断或声称互作显著性。

## PGAM5空间模式和邻近T程序统计

全主QC空间位点PGAM5 Moran统计：

| method | Moran_I | conditional_expected_I | permutation_p | BH_q_2_tests |
| --- | --- | --- | --- | --- |
| unrestricted | 0.02606 | -0.00061 | 0.0632 | 0.1264 |
| UMI_MT_stratified | 0.02606 | 0.0065 | 0.1798 | 0.1798 |

每类4,999次置换，分别为全片无约束及UMI五分位×MT三分位层内置换；双侧p围绕各模型精确期望值计算，再2项BH。PGAM5全位点空间自相关未通过校正；该检验不是巨噬细胞专属PGAM5富集检验。

关联统计：

| program | cohort | contrast | n_sources | statistic | permutation_p | BH_q_7_planned_tests | status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Exhaustion_core | all_mac_including_zero | partial_Spearman | 163 | 0.00393 | 0.7804 | 1.0 | tested |
| Exhaustion_expanded | all_mac_including_zero | partial_Spearman | 163 | 0.02034 | 0.952 | 1.0 | tested |
| Activation_early | all_mac_including_zero | partial_Spearman | 163 | 0.08327 | 0.1712 | 1.0 | tested |
| Effector_cytotoxic | all_mac_including_zero | partial_Spearman | 163 | 0.04407 | 0.6796 | 1.0 | tested |
| Exhaustion_no_CTLA4_TIGIT | all_mac_including_zero | partial_Spearman | 163 | 0.03624 | 0.8296 | 1.0 | tested |
| Exhaustion_core | PGAM5_detected_only | partial_Spearman | 19 | — | — | 1.0 | not_testable |
| Exhaustion_core | PGAM5_detected_only | adjusted_high_minus_other | 19 | — | — | 1.0 | not_testable |

暴露为巨噬RNA源PGAM5归一化表达，主分析**保留零值**；结局为该源所有直接T-RNA邻居的程序均值，保留T位点程序零值，没有T邻居的源记缺失，不填0。每源一条观察，但多个源可共享T靶；没有把边当独立患者。

偏Spearman移除9项秩协变量：源log1p(UMI)、MT%、巨噬核心、肝细胞、上皮、T谱系，以及T邻居的T谱系、log1p(UMI)、MT%均值。3×3原始array坐标网格×分析源UMI二分位中置换PGAM5，4,999次；双侧等尾随机化p，避免强行假定零模型中心0。粗空间分块只近似保留局部结构，层内可交换性及剩余混杂未完全解决；协变量校正不意味着识别了巨噬细胞内PGAM5或因果效果。

耗竭核心为PDCD1/LAG3/TIGIT/CTLA4/TOX/CXCL13；扩展加HAVCR2/ENTPD1。由于单细胞发现明确Treg相关群，新增**去CTLA4/TIGIT**的PDCD1/LAG3/TOX/CXCL13/HAVCR2/ENTPD1敏感性程序，rho=0.03624、p=0.8296，仍不显著。去两基因不能保证程序专属于耗竭T细胞，剩余基因亦有状态重叠。

预设7项（五套含零值程序关联、仅检出耗竭关联、高vs其他检出耗竭），实际可检验5项；联合BH全q=1。主检验单独p=0.7804。仅检出范围19源（4高、15其他）小于固定源数30，未进行检出位点/高低组推断；缺失p保留`not_testable`，只在BH计算临时设1。未以低样本量强行给出p、普通独立位点t/Wilcoxon检验、伪置信区间或患者层面预后结论。

## 生物学含义和限制

这个病例可以支持“巨噬细胞相关RNA群及PGAM5检出状态”的描述，但6个严格候选均1 UMI、同病例、PGAM5广泛见于肝细胞样和其他群，不能定义稳定/特异PGAM5⁺巨噬细胞亚群或可靠TCGA丰度参考。单细胞PGAM5归一化高值也受稀疏检出和文库深度影响。

空间RNA没有显示可靠的PGAM5–邻域T耗竭正关联；不能据此说促进或抑制耗竭。混合Visium位点不定位RNA到其中某种细胞；Treg/活化与耗竭标志重叠，程序不等于功能衰竭。IFN程序是本轮固定RNA面板，不是原论文癌细胞IFN模块的精确复现。未使用原论文MNN/NNLS细胞比例、CODEX蛋白或临床证据，未把失败的旧PGAM5细胞丰度参考应用于此病例。

## 复现与独立核验

通过2,008项检查：[audit_results.json](audit_results.json)。核验13原文件长度/SHA、MatrixMarket基因/细胞注释顺序、独立原始矩阵QC、原始计数与全库CP10k、26套程序、RNA名称映射、真实坐标/scale/组织图、KDTree独立六邻接、Moran期望/p/BH、21候选受体亚基，以及逐源邻近结局、独立最小二乘偏相关、置换和BH7/图像校验。

```bash
cd GSE203612_LIHC1_PGAM5
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
export GSE203612_WORK="$PWD/work"
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
export MPLCONFIGDIR="$GSE203612_WORK/matplotlib"
export XDG_CACHE_HOME="$GSE203612_WORK/cache"
mkdir -p "$GSE203612_WORK"
python download_inputs.py
python preprocess.py
python annotate.py
python spatial_statistics.py
python association.py
python plot.py
python audit.py
python build_report.py
python package_results.py
```

Python3.12；[requirements.txt](requirements.txt)固定库版本。全部便携对象和逐条码数据公开，h5ad保留原计数、全原矩阵总量归一化表达、UMAP及真实空间信息，但**仅为marker子集**，不得用子集总计数重新归一化或代替全基因HVG/差异分析。全基因对象保存在scratch/work，原始数据可由脚本自动重取；复现分群需该步骤。

## 下载和文献

[下载全部图表、便携对象、逐条码注释、统计和代码](GSE203612_LIHC1_PGAM5_results.zip)，打开后点击 **Download raw file**。[包大小与SHA256](package_manifest.json)。不含原始22.59MiB附件，`download_inputs.py`可自动下载并核验。

- [官方GEO GSE203612](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203612)；[GSM6177611](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6177611)；[GSM6177612](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6177612)。SOFT和官方文件清单在[provenance](provenance)。
- 原论文 *Cancer cell states recur across tumor types and form specific interactions with the tumor microenvironment*，Nature Genetics2022；[PMID35931863](https://pubmed.ncbi.nlm.nih.gov/35931863/)，[DOI10.1038/s41588-022-01141-9](https://doi.org/10.1038/s41588-022-01141-9)，PMC9886402。官方摘要/全文XML保留，本轮与作者分析方法的差异如上。
