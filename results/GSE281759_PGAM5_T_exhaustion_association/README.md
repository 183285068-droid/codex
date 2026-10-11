# GSE281759：PGAM5相关巨噬RNA邻域与T细胞耗竭程序

**结果：当前空间RNA数据没有显示稳定、统计显著的“PGAM5越高，邻近T细胞耗竭越强”关联，也不能据此认定PGAM5促进或不促进T细胞耗竭。**

4位HCC患者的16张肿瘤切片中，2,877个巨噬RNA富集源位点有直接T-RNA邻居。主分析包括PGAM5零值，在15张源位点足够的切片（合计2,859个源位点）检验；C008T3只有18个源位点，不进行推断。

主分析调整相关系数rho约−0.103至0.286，方向不一致，空间条件置换p≥0.062，16项BH校正q≥0.664。PGAM5高表达vs其他检出组在11张可检验切片中也没有通过多重校正。其余结局/敏感性共同预设96项，实际可计算82项，联合BH后全部q=1。**不显著不是不存在作用的证明；只有4个独立患者，混合位点归属、低检出、门控及统计功效都限制结论。**

## 主要图表

1. [主分析：效应值、空间置换p和BH q](figures/primary_PGAM5_neighbor_exhaustion_association.png)（[PDF](figures/primary_PGAM5_neighbor_exhaustion_association.pdf)）。红色正关联，蓝色负关联；没有显著切片。
2. [耗竭、扩展耗竭、早期活化、细胞毒程序关联热图](figures/PGAM5_neighbor_T_program_association_heatmap.png)（[PDF](figures/PGAM5_neighbor_T_program_association_heatmap.pdf)）。此图q是96项联合校正，主分析图q是16项主检验校正。
3. [逐切片校正散点图](figures/adjusted_PGAM5_exhaustion_scatter_by_section.png)（[PDF](figures/adjusted_PGAM5_exhaustion_scatter_by_section.pdf)）。横纵轴均为移除9协变量后的秩残差，以各自SD标准化。
4. [PGAM5高、其他检出、未检出组的邻近耗竭分数](figures/PGAM5_groups_neighbor_exhaustion_boxplots.png)（[PDF](figures/PGAM5_groups_neighbor_exhaustion_boxplots.pdf)）。箱线是原始程序分数；标注q属于H vs D的校正秩残差比较，**不是未校正箱线检验或含未检出组的三组检验**。
5. 真坐标组织图上PGAM5巨噬RNA源与耗竭T-RNA位置：[A009](figures/A009_PGAM5_Tex_spatial_neighborhoods.png) · [B013](figures/B013_PGAM5_Tex_spatial_neighborhoods.png) · [C004](figures/C004_PGAM5_Tex_spatial_neighborhoods.png) · [C008](figures/C008_PGAM5_Tex_spatial_neighborhoods.png)。青色圈为PGAM5高源且有耗竭支持T-RNA直接邻居；仅表示共定位。

共8组独立PNG/PDF。空间色标跨切片共用，超各切片Q99最大值者仅显示截断，统计值未截断。

## 主检验结果

| sample | n_sources | statistic | permutation_p | primary_BH_q_16_sections | status |
| --- | --- | --- | --- | --- | --- |
| A009T1 | 252 | -0.0662 | 0.166 | 0.664 | tested |
| A009T2 | 62 | 0.1003 | 0.312 | 0.832 | tested |
| A009T3 | 65 | 0.047 | 0.8932 | 0.9527 | tested |
| A009T4 | 442 | -0.0186 | 0.8184 | 0.9527 | tested |
| A009T5 | 261 | -0.0248 | 0.6856 | 0.9527 | tested |
| B0013T1 | 139 | 0.0321 | 0.5028 | 0.8939 | tested |
| B0013T2 | 178 | 0.0386 | 0.43 | 0.86 | tested |
| B0013T3 | 100 | 0.2858 | 0.062 | 0.664 | tested |
| B0013T4 | 239 | -0.0823 | 0.2828 | 0.832 | tested |
| B013T5 | 241 | -0.1035 | 0.1136 | 0.664 | tested |
| C004T1P | 89 | 0.0466 | 0.8528 | 0.9527 | tested |
| C004T1R | 52 | 0.1157 | 0.5676 | 0.9082 | tested |
| C008T1 | 413 | 0.0864 | 0.1324 | 0.664 | tested |
| C008T2 | 169 | -0.0292 | 0.7508 | 0.9527 | tested |
| C008T3 | 18 | — | — | 1.0 | not_testable |
| C008T4 | 157 | -0.0692 | 0.4224 | 0.86 | tested |

下表仅描述同患者切片方向，不将多切片当作独立患者进行荟萃或合并患者p值。C004T1P是原发、C004T1R是复发，仍为同一患者。

| patient | sections_tested | median_section_partial_rho | positive_sections | negative_sections | q_significant_sections |
| --- | --- | --- | --- | --- | --- |
| A009 | 5 | -0.0186 | 2 | 3 | 0 |
| B013 | 5 | 0.0321 | 3 | 2 | 0 |
| C004 | 2 | 0.0812 | 2 | 0 | 0 |
| C008 | 3 | -0.0292 | 1 | 2 | 0 |

## 高表达与其他检出组比较

| sample | n_sources | n_high | n_other_detected | raw_mean_difference | statistic | permutation_p | BH_q_96_planned_tests | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A009T1 | 112 | 29 | 83 | -0.00816 | -0.2266 | 0.1044 | 1.0 | tested |
| A009T2 | 36 | 7 | 29 | 0.00894 | -0.06364 | 0.5116 | 1.0 | tested |
| A009T3 | 48 | 11 | 37 | 0.00069 | -0.04463 | 0.6368 | 1.0 | tested |
| A009T4 | 269 | 64 | 205 | -0.00337 | -0.00605 | 0.896 | 1.0 | tested |
| A009T5 | 163 | 37 | 126 | -0.00275 | 0.07081 | 0.7792 | 1.0 | tested |
| B0013T1 | 12 | 4 | 8 | — | — | — | 1.0 | not_testable |
| B0013T2 | 72 | 11 | 61 | 0.02504 | 0.20245 | 0.3668 | 1.0 | tested |
| B0013T3 | 12 | 7 | 5 | — | — | — | 1.0 | not_testable |
| B0013T4 | 44 | 11 | 33 | 0.0143 | 0.29216 | 0.618 | 1.0 | tested |
| B013T5 | 100 | 18 | 82 | -0.0164 | -0.42069 | 0.2152 | 1.0 | tested |
| C004T1P | 17 | 4 | 13 | — | — | — | 1.0 | not_testable |
| C004T1R | 18 | 3 | 15 | — | — | — | 1.0 | not_testable |
| C008T1 | 96 | 25 | 71 | -0.00407 | -0.06379 | 0.6908 | 1.0 | tested |
| C008T2 | 35 | 8 | 27 | 0.01311 | 0.08952 | 0.536 | 1.0 | tested |
| C008T3 | 2 | 0 | 2 | — | — | — | 1.0 | not_testable |
| C008T4 | 67 | 15 | 52 | -0.01517 | -0.21914 | 0.9304 | 1.0 | tested |

`raw_mean_difference`是包含零值的邻近T耗竭程序均值差（高−其他检出）。`statistic`是移除协变量后的**秩残差**两组均值差/残差SD，既不是log2FC也不是细胞百分比。p检验后者。

只有一个检出位点敏感性结果名义p<0.05：B0013T4，rho≈0.467、p=0.0276；联合校正q=1，同切片含零值主分析rho≈−0.082、p≈0.283。该结果随分析范围改变方向，不能作为稳定的支持证据，未用于事后重定义目标群或挑选亚组。

## 暴露、结局和门控

沿用上一轮固定数据/定义，**不根据本轮关联结果修改RNA门或选择切片**：

- 主QC：UMI≥500、检出基因≥200、MT比例≤25%、in_tissue=1。原论文MT≤5%，本轮保留较宽覆盖；详细QC差异见[基础空间报告](../GSE281759_spatial_HCC_PGAM5/README.md)。尚未复算全部5%规则下的关联，这仍是需进一步检验的敏感性限制。
- 巨噬RNA富集门：巨噬核心程序≥切片Q75，至少3个核心基因检出，泛髓系程序≥切片中位数；T-RNA门：T谱系程序≥切片Q75且至少2个谱系基因检出。
- PGAM5暴露：源位点全基因库CP10k后的log1p表达，**主分析包含零值**。没有T邻居的源仍保存到逐片表，其邻近分数记缺失，不填0。
- 耗竭核心：PDCD1、LAG3、TIGIT、CTLA4、TOX、CXCL13平均log1p(CP10k)；扩展增加HAVCR2、ENTPD1。早期活化及细胞毒程序均为既定基因集，基因清单见基础[programs.json](../GSE281759_spatial_HCC_PGAM5/programs.json)。PGAM5是独立暴露变量，耗竭结局遵循既定程序。
- 结局为每个源位点**所有直接T-RNA邻居的耗竭程序均值**，包含T位点程序基因零值，不只挑选已经判为耗竭的靶位点，从而避免先选高耗竭邻居再比较的选择问题。每个源一条观察；同一T位点可邻近多个源，因此观察仍相关。
- 辅助PGAM5-high：巨噬RNA门内已检出PGAM5位点的表达≥该切片检出位点Q75，比较其他检出位点；未检出组仅描述。群阈值来自完整切片巨噬RNA门，**不是在只有T邻居的子集重新选Q75**。
- 直接邻居来自真实Visium六边形格的6方向（排除自环）。组织图使用原像素坐标和scale，不将UMAP当作空间位置。

**“Mac-RNA源位点的PGAM5”不等于“巨噬细胞内的PGAM5”。** Visium位点混合细胞，PGAM5还可来自肿瘤/其他细胞。耗竭程序也可能受活化、混合RNA或低检出影响；C004无≥2基因耗竭支持T门，不代表不存在T细胞耗竭，本轮连续程序相关只反映其弱RNA信号。RNA门并非细胞百分比、已验证单细胞身份或病理恶性ROI。

## 校正和空间统计

每张切片独立计算偏Spearman：对PGAM5、邻近程序和9项协变量取秩，分别将暴露/结局投影到协变量空间并取残差，再求残差Pearson相关。

9项协变量为源log1p(UMI)、源MT%、源巨噬核心、源肝细胞、源上皮、源T谱系，以及T邻居的T谱系均值、log1p(UMI)均值和MT%均值。主要用于减少质量和RNA混合信号差异；秩回归只是关联模型，**不能完全分离细胞来源，也不能保证排除未观测混杂**。

每个分析范围分别按全切片array坐标固定3×3空间网格，再按分析源UMI二分位分层，最多18层；层内置换PGAM5表达（高/其他组比较置换高标签），4,999次，保持RNA暴露/标签分布及源/靶坐标和邻近结局固定。连续检验每次重算暴露秩和残差；分组检验保持各层高标签数。p采用双侧等尾随机化检验`min(1, 2 × min((1+null≥obs)/5000, (1+null≤obs)/5000))`，不假定零分布中心是0。

置换是保留粗空间/UMI结构的**近似条件随机化检验**，不能保留全部微观空间自相关；其有效性依赖层内可交换性。源/靶共享及空间分块限制没有变成独立生物学重复。源数≥30、协变量拟合后至少15维残余自由度、结局与暴露需有变化；分组另要求高组≥5、其他检出组≥10。单元素层不可置换，层内完全相同值也不增加信息；4,999次不是4,999个独立病例。

主家族16个切片耗竭核心含零值检验做BH；全家族16×6=96项（四个程序含零值、检出位点耗竭、分组耗竭）另共同BH。未可检验p保留缺失，在家族校正计算中临时设1，输出`not_testable`，不称其“证实无差异”。没有混合患者位点当作独立病例的普通t/Wilcoxon p值，没有伪造置信区间，没有患者层面因果或临床结局检验。

## 能否证明“促进”及下一步

当前结果可以描述“PGAM5相关巨噬RNA邻域与T耗竭转录程序是否关联”；本轮没有发现稳定正关联。要将其升级为“促进”的结论，需要独立识别巨噬细胞中的PGAM5 RNA来源，并建立干预对照：例如巨噬细胞PGAM5敲低/敲除、对照及救援后与T细胞共培养，检测T细胞耗竭转录程序和增殖/杀伤功能变化。可以以RNA定义目标并追踪RNA结局，不必先把蛋白检测作为本轮分析前提。单独虚拟敲除或候选配体受体共检出也不能证明因果。

现有[候选配体受体图](../GSE281759_spatial_HCC_PGAM5/figures/PGAM5_high_Tex_candidate_LR_codetection.png)只提供PVR/NECTIN2–TIGIT等后续检查线索；TIGIT/PDCD1/CTLA4本身又参与耗竭定义，不能把该共检出当成独立耗竭验证。本轮未对这些候选声称促耗竭机制。

## 独立核验与复现

34,070项独立检查通过：[audit_results.json](audit_results.json)。包括按真实邻接逐源重算靶程序均值/零值/缺失、源去重、固定RNA暴露、9协变量、独立最小二乘残差统计、每层置换不变、前三次零统计重算、等尾p、BH16/BH96及8PNG/PDF文件校验。基础数据已另通过6,020项原始计数/空间/图片检查。

所有逐源表以17位有效数值保存，避免float32在CSV往返时截断；随机种子和门控配置见[analysis_config.json](analysis_config.json)。便携基线对象来自固定Git提交`6bb037700d74331b25f5ab5c21d167f178927fb3`，17个输入（16对象+范围表）的长度/SHA记录在[reference_manifest.json](reference_manifest.json)。

```bash
cd GSE281759_PGAM5_T_exhaustion_association
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
export MPLCONFIGDIR="/tmp/GSE281759_mpl"
export XDG_CACHE_HOME="/tmp/GSE281759_cache"
python fetch_reference.py
python analyze.py
python plot.py
python audit.py
python build_report.py
python package_results.py
```

Python3.12，固定库版本；下载脚本在相邻`GSE281759_spatial_HCC_PGAM5`目录恢复本轮所需基线对象，对已有文件仅核验。重新进行全基因分群/QC需基础报告的原始下载和preprocess流程。本轮未调用失败的细胞丰度反卷积参考。

## 下载

[全部图表、逐源数据、统计结果和代码](GSE281759_PGAM5_T_exhaustion_association_results.zip)。打开GitHub文件，点击 **Download raw file**；包不含已有基线对象，复现可由`fetch_reference.py`自动下载。[包SHA256及大小](package_manifest.json)。
