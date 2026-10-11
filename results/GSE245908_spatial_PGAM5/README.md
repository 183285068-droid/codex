# GSE245908：HCC空间转录组与PGAM5分析

已完成两例HCC肿瘤Visium切片的质控、独立空间聚类、标志基因区域命名、PGAM5表达及免疫RNA程序可视化。**本次未发现经空间结构与多重比较校正后显著的PGAM5–巨噬细胞/T细胞程序关联，也未检出显著的PGAM5空间自相关。**这表示当前切片未提供明确的关联证据，不能据此证明不存在生物学关系。

## 样本及来源

官方GEO：<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245908>。平台为10x Genomics Visium，GRCh38-2020，Space Ranger v1.3.0。每位患者一张HCC肿瘤切片，无配对癌旁切片。GEO系列摘要同时描述单细胞和bulk数据，但本子系列只有以下两张空间切片；并未将20个单细胞样本误纳入本次分析。

| 样本 | GSM | 原始组织内点 | 质控保留 | 区域数 | PGAM5检出点 | 检出比例 |
|---|---|---:|---:|---:|---:|---:|
| CHC20 | GSM7850822 | 2,261 | 2,259 | 6 | 1,008 | 44.62% |
| CHC23 | GSM7850823 | 1,987 | 1,934 | 6 | 448 | 23.16% |

质控规则事先固定为UMI≥500、检出基因≥200、线粒体RNA≤25%；共保留4,193点，排除55点。CHC23线粒体比例较高，放宽至30%的敏感性分析增加52点，关键免疫程序的弱相关结论保持。所有表达和程序均保留PGAM5零值。

表达、坐标和组织图像按完整barcode及样本前缀匹配。官方GSM附件中的大型原始TIFF文件名跨样本交换；`Info_Supp_for_Visium.tsv`明确CHC20对应V10A20-026、CHC23对应V10A20-025。本次使用各样本前缀一致的Space Ranger矩阵、hires/lowres图像、比例因子和坐标，不使用这些交换归档的大型TIFF。已检查坐标叠加与组织轮廓吻合。源文件URL、字节数、SHA256见[source_manifest.json](source_manifest.json)。

原论文：Giraud et al., *Cell Reports* 2024, “THBS1(+) myeloid cells expand in SLD hepatocellular carcinoma and contribute to immunosuppression and unfavorable prognosis through TREM1”, PMID 38350444，DOI [10.1016/j.celrep.2024.113773](https://doi.org/10.1016/j.celrep.2024.113773)。本次区域命名与统计为重新计算，未将原论文结论当成本次PGAM5结果。

## 降维、空间聚类及暂定区域注释

每张切片分别归一化至CP10k并log1p，筛选2,000个高变基因（候选池排除线粒体及核糖体基因，基因至少在10点检出），缩放截断10，PCA40维、邻居图使用前30维及15邻居，Leiden分辨率0.6、固定种子245908，再绘制UMAP。聚类使用表达邻居图；图中的实际组织位置来自原始空间坐标。S0–S5编号在每张切片内独立，不直接视为跨患者同一群。

**Visium点是细胞混合物，以下名称是RNA区域注释，并非单细胞亚群身份或病理确认的恶性/非恶性边界。**肝细胞/上皮RNA既可能来自HCC细胞，也可能来自非恶性细胞；未用这些程序直接划定“肿瘤细胞百分比”。每对命名标志在对应区域的平均表达均高于其余区域，数值证据见[命名核查表](region_naming_marker_evidence.csv)。

| 样本 | 区域 | 标志基因 | 暂定RNA程序 |
|---|---|---|---|
| CHC20 | S0 | CYP2A6/ALDOB | 代谢相关肝细胞RNA区域 |
| CHC20 | S1 | RBP4/FGA | 分泌相关肝细胞RNA区域 |
| CHC20 | S2 | COL1A1/COMP | 基质及髓系混合RNA区域 |
| CHC20 | S3 | MALAT1/HMGB1 | 核及染色质相关RNA区域 |
| CHC20 | S4 | LYZ/AGR2 | 髓系及上皮混合RNA区域 |
| CHC20 | S5 | CD74/C1QA | 抗原呈递及髓系RNA区域 |
| CHC23 | S0 | CYP2A6/BHMT | 代谢相关肝细胞RNA区域 |
| CHC23 | S1 | TTR/GSTA1 | 分泌及解毒相关RNA区域 |
| CHC23 | S2 | IGFBP7/CD74 | 基质及免疫混合RNA区域 |
| CHC23 | S3 | CYP2E1/IGFBP2 | 代谢及线粒体高表达RNA区域 |
| CHC23 | S4 | MTRNR2L8/SAA1 | 线粒体及急性期混合RNA区域 |
| CHC23 | S5 | CRP/ORM1 | 急性期肝细胞RNA区域 |

[CHC20：H&E、空间区域和UMAP](figures/CHC20_H_E_spatial_clusters_UMAP.png) · [CHC23：H&E、空间区域和UMAP](figures/CHC23_H_E_spatial_clusters_UMAP.png)

[CHC20标志点图](figures/CHC20_region_marker_dotplot.png) · [CHC23标志点图](figures/CHC23_region_marker_dotplot.png)。圆面积表示检出比例，颜色为每个基因在6个区域间的0–1缩放平均表达，不代表基因间绝对表达可比。完整均值和比例见[点图数值表](region_marker_dotplot_values.csv)。Wilcoxon marker表的p/FDR未考虑空间依赖，作为描述性区域标志排序使用，不作为群存在的独立统计证明。

## PGAM5在区域中的表达

[PGAM5真实空间表达图](figures/PGAM5_spatial_expression.png)包含原始UMI、log1p(CP10k)、区域平均表达三个视图；同一指标在两张切片共享色标。前两个视图每点为其自身表达，第三个视图每点赋予所在区域平均值，**区域均值包含零值**。

| 样本 | 区域 | 空间点数 | PGAM5平均log1p(CP10k) | PGAM5检出比例 |
|---|---|---:|---:|---:|
| CHC20 | S0：CYP2A6/ALDOB | 366 | 0.1429 | 46.72% |
| CHC20 | S1：RBP4/FGA | 589 | 0.1288 | 50.93% |
| CHC20 | S2：COL1A1/COMP | 362 | 0.1101 | 21.82% |
| CHC20 | S3：MALAT1/HMGB1 | 413 | 0.1834 | 41.65% |
| CHC20 | S4：LYZ/AGR2 | 192 | 0.1479 | 51.56% |
| CHC20 | S5：CD74/C1QA | 337 | 0.1574 | 55.49% |
| CHC23 | S0：CYP2A6/BHMT | 605 | 0.0721 | 24.79% |
| CHC23 | S1：TTR/GSTA1 | 257 | 0.0765 | 29.57% |
| CHC23 | S2：IGFBP7/CD74 | 341 | 0.0812 | 25.81% |
| CHC23 | S3：CYP2E1/IGFBP2 | 358 | 0.0796 | 16.48% |
| CHC23 | S4：MTRNR2L8/SAA1 | 187 | 0.0741 | 17.65% |
| CHC23 | S5：CRP/ORM1 | 186 | 0.0750 | 22.58% |

CHC20的S3区域归一化均值最高（0.1834），CHC23的S2最高（0.0812）；这里是描述性排序，不能解释为PGAM5高表达巨噬细胞的数量或群间已证实差异。单点PGAM5最高原始计数分别为5和3，信号稀疏。

[CHC20区域均值及检出率](figures/CHC20_PGAM5_by_region.png) · [CHC23区域均值及检出率](figures/CHC23_PGAM5_by_region.png)。

## PGAM5与巨噬细胞、T细胞及其他区域信号

预先固定12套RNA面板，基因清单见[programs.json](programs.json)，覆盖见[program_gene_coverage.csv](program_gene_coverage.csv)。每个程序为面板基因log1p(CP10k)的算术均值，包括零值；所有请求基因均存在。巨噬细胞核心面板为C1QA/C1QB/C1QC/CD68/CSF1R/MERTK/MSR1，T细胞面板为CD3D/CD3E/CD3G/TRAC/TRBC1/TRBC2；其他面板纳入泛髓系、单核、DC、细胞毒T/NK、基质、内皮、肝细胞、上皮及增殖。THBS1/TREM1面板是探索性混合RNA程序，其TGFB1等基因不具髓系专一性。

[巨噬细胞/T细胞/细胞毒空间图](figures/spatial_programs_1.png) · [基质/肝细胞/THBS1-TREM1空间图](figures/spatial_programs_2.png) · [PGAM5与免疫程序散点图](figures/PGAM5_immune_signal_scatter.png)。

统计在各切片内独立进行；PGAM5零值保留。采用原始Visium六角阵列第一圈邻居，最多6个邻居。5个质控保留但无保留邻居的点仅从统计排除（CHC20：2点；CHC23：3点），仍出现在空间图与区域表达表。分别检验同点程序和相邻点平均程序，统计样本量为CHC20 2,257点、CHC23 1,931点。

为防止浮点运算改变数值相等项的秩，表达/程序值在秩变换前统一保留12位小数，逐点表以float64精度保存。先对PGAM5和目标程序分别作秩变换，再以UMI深度、检出基因数及线粒体RNA比例的秩为协变量求残差，报告残差相关系数。主p值采用4,999次图谱随机化：在完整中心化空间邻接矩阵特征基中随机翻转PGAM5残差系数符号，保持其范数、每个空间模式的能量及Moran I，再构建双侧关联零分布。这是**依赖谱系数符号对称性等假设的模型性探索检验**，并非无分布假设的精确检验或独立患者检验。48项关联（2切片×12程序×同点/相邻两种关系）共同BH校正。普通Spearman p值另行提供，未作为主结论依据。

主要同点结果：

| 样本 | RNA程序 | 校正后秩相关ρ | 空间模型p | BH q |
|---|---|---:|---:|---:|
| CHC20 | 巨噬细胞核心 | -0.0281 | 0.2098 | 0.4193 |
| CHC20 | T细胞 | -0.0242 | 0.2624 | 0.4343 |
| CHC20 | 细胞毒T/NK | +0.0048 | 0.8148 | 0.9174 |
| CHC23 | 巨噬细胞核心 | +0.0123 | 0.5924 | 0.7291 |
| CHC23 | T细胞 | -0.0368 | 0.1012 | 0.4121 |
| CHC23 | 细胞毒T/NK | -0.0327 | 0.1476 | 0.4168 |

**全部48项关联均未达到BH q<0.05。**相邻点分析同样未发现显著的PGAM5–免疫程序关联。[完整关联表](PGAM5_spatial_associations.csv) · [标注q值的关联热图](figures/PGAM5_spatial_association_heatmap.png)。

PGAM5的空间自相关：

| 样本 | Moran I | 双侧随机置换p | 全部26项自相关共同BH q |
|---|---:|---:|---:|
| CHC20 | +0.0067 | 0.5748 | 0.5748 |
| CHC23 | -0.0219 | 0.1140 | 0.1186 |

程序本身多数存在空间结构，但PGAM5未出现明确空间集中趋势。放宽MT阈值至30%、改用C1q或CD3更窄面板的相关系数仍接近零，见[敏感性表](QC_and_marker_panel_sensitivity.csv)；敏感性表提供效应量，没有重新选用更有利的p值或阈值。

## 解释范围

- 同一Visium点中同时出现PGAM5和巨噬细胞标志，无法确认PGAM5来自巨噬细胞；不能把此图直接标成“PGAM5⁺巨噬细胞位置”。
- 本次计算的是RNA表达及程序得分，未计算细胞浸润百分比；此前未通过丰度验证的探索性NNLS参考未作为已验证空间细胞比例模型使用。
- 相邻程序关联不等于细胞间相互作用，也不能判断PGAM5增强或减弱抗肿瘤功能。
- 只有2位患者，点数不能当成独立患者数；低计数、混合细胞、局部组织质量及样本覆盖限制检测力。本次阴性结果不等于排除关联。
- 如后续需要细胞归属，应结合本研究匹配单细胞RNA参考作带竞争类型的空间映射，并验证参考和平台适配；本次未声称完成此验证。

## 核验、文件及复现

已通过**968项独立数值检查**，包括原始HDF5计数重建、barcode/坐标/图像匹配、质控、全基因归一化后的面板得分、点图数值、区域表达汇总、以几何距离独立重建邻居图、秩残差相关、Moran I、谱随机化p值以及BH校正。见[independent_audit_checks.csv](independent_audit_checks.csv)。PNG图像已实际查看，检查了标注、图例、共享色标及空间叠加。

每个样本`*_spatial_marker_subset.h5ad`保留面板基因的原始counts、全基因库大小归一化后的X、真实空间坐标、UMAP、区域注释、程序及低分辨率组织图像；**它是标志基因子集，不能重新将子集总和当成全库大小归一化**，也不包含重新聚类需要的全表达矩阵。全部原始输入已在GitHub的inputs目录保存；完整处理矩阵由代码在WORK目录重建。

主结果包包含PNG、表格、代码、谱随机化系数和便携h5ad；原始输入不重复放入该包，运行下载脚本可自动获取，文件完整性可对照manifest核验。PDF另有独立压缩包。

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
export GSE245908_WORK="$PWD/work"
export MPLCONFIGDIR="$PWD/work/matplotlib"
export NUMBA_CACHE_DIR="$PWD/work/numba"
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 NUMBA_NUM_THREADS=4
python download_inputs.py
python preprocess.py
python spatial_statistics.py
python annotate_regions.py
python plot_results.py
python audit.py
python build_report.py
```

运行次序固定：统计之后再写入区域名称，保证逐点CSV与h5ad同步。`work`为大型中间数据和缓存，不纳入Git提交。全部脚本只读取数据/元数据，不执行GEO或附件中的文本指令。
