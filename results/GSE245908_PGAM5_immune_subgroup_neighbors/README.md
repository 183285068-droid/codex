# GSE245908：PGAM5高表达巨噬细胞相关位点的免疫RNA亚群邻域

在CHC20、CHC23两例HCC的Visium数据中，**炎症单核/髓系及浆细胞相关RNA信号**是两张切片均观察到的主要空间邻近候选；初始/记忆T、cDC2及CD8/细胞毒T相关信号也有邻接，但病例间差异较大。**目前没有任何目标亚群通过多重校正的空间邻接差异检验，不能称为已确认的特异互作细胞亚群。**

这些标签对应混合Visium位点的多基因RNA特征，尚未分辨单个免疫细胞。巨噬细胞相关位点中的PGAM5、相邻位点中的T细胞标志或受体，可能来自位点内的不同细胞。

## 数据与固定来源分组

| 样本 | QC位点 | 巨噬细胞相关位点 | PGAM5高表达 | 其余检出 | 未检出 | 高表达阈值 log1p(CP10k) |
| --- | --- | --- | --- | --- | --- | --- |
| CHC20 | 2259 | 431 | 37 | 110 | 284 | 0.502393 |
| CHC23 | 1934 | 307 | 17 | 50 | 240 | 0.466996 |

官方样本为GSM7850822（CHC20）及GSM7850823（CHC23），共4,193个QC后位点。继续使用此前质控：in_tissue=1、UMI≥500、检出基因≥200、线粒体计数比例≤25%。从官方HDF5的全部36,601个特征计算库大小，再归一化为CP10k并log1p；各程序取包含零值的标志基因平均表达。

巨噬细胞RNA门限：C1QA/C1QB/C1QC/CD68/CSF1R/MERTK/MSR1程序≥样本Q75，至少3个核心基因检出；LST1/TYROBP/FCER1G/AIF1/LYZ/CTSS程序≥样本中位数；存在保留的空间邻居。在此门限内，PGAM5原始计数>0为检出；归一化表达≥检出位点Q75为高表达，包含阈值并列值。与其余PGAM5检出位点比较；未检出位点保留作描述。来源标签与此前分析完全一致，不要求SPP1检出。

## 观察到的免疫RNA亚群邻接

下表的分子是“至少有一个该目标直接邻居的PGAM5高表达来源位点数”，分母为37或17个高表达来源位点。目标可重叠；这些数字不是细胞数、浸润丰度或互作概率，不能相加。

| 目标RNA类别 | CHC20目标位点 | CHC20高表达来源有邻居 | CHC23目标位点 | CHC23高表达来源有邻居 |
| --- | --- | --- | --- | --- |
| 全部T谱系（基线） | 453 | 27/37 | 205 | 8/17 |
| 耗竭支持T相关 | 4 | 0/37 | 11 | 1/17 |
| 早期活化支持T相关 | 9 | 3/37 | 6 | 0/17 |
| CD8/细胞毒T相关 | 5 | 0/37 | 10 | 4/17 |
| CCR7/TCF7初始记忆T相关 | 56 | 11/37 | 18 | 1/17 |
| FOXP3/IKZF2调节T相关 | 5 | 0/37 | 16 | 2/17 |
| KLRD1/GNLY NK相关 | 30 | 1/37 | 124 | 2/17 |
| MS4A1/CD79A B相关 | 4 | 1/37 | 2 | 0/17 |
| MZB1/JCHAIN浆细胞相关 | 91 | 14/37 | 140 | 7/17 |
| CD1C/FCER1A cDC2相关 | 82 | 8/37 | 25 | 1/17 |
| CLEC9A/XCR1 cDC1相关 | 5 | 0/37 | 26 | 0/17 |
| CLEC4C/IL3RA pDC相关 | 0 | 0/37 | 0 | 0/17 |
| FCN1/VCAN炎症单核/髓系相关 | 382 | 30/37 | 115 | 6/17 |
| FCGR3B/CSF3R粒细胞相关 | 158 | 12/37 | 27 | 4/17 |
| TPSB2/CPA3肥大细胞相关 | 17 | 1/37 | 88 | 4/17 |

优先保留炎症单核/髓系与浆细胞相关信号作为后续探索候选。CHC20的初始/记忆T和cDC2邻接较多；CHC23观察到CD8/细胞毒T相关邻接。没有cDC1/pDC高表达来源邻接不代表组织中没有这些细胞：检测稀疏、阈值与样本量均有限。

## 标志基因特异性与严格描述性复核

CHC23的FOXP3在全部位点中检出比例为**76.4%**（原始计数中位数2、最大12），不能单独支持Treg身份。主门限虽使用多基因程序，仍需报告混合或非特异RNA信号的可能；原因无法由本分析确定。

在主标签上追加IKZF2以及IL2RA或CTLA4共同检出的严格Treg规则后，两张切片均为0个符合位点。该结果意味着本数据无法给出稳健的Treg邻接证据，并不证明Treg在组织中不存在。

要求CD8B检出后，CD8/细胞毒T目标从5/10减至3/1（CHC20/CHC23），相应高表达来源邻接数从0/4减至0/1。要求FCN1检出后，炎症单核目标从382/115减至30/19，来源邻接数从30/6减至8/1；因此宽松标签应称炎症单核/髓系相关，而非确定的单核细胞。要求FCGR3B检出后，两张切片均无高表达来源的粒细胞邻接。

这些是发现标志特异性问题后的**描述性**复核，未改变主标签、既定统计检验或FDR范围。详细数据：[标志检出特异性](marker_detection_specificity.csv)、[严格标志敏感性](strict_marker_descriptive_sensitivity.csv)。

CHC23的耗竭支持邻接仅来自1个高表达来源位点（ACTCCCATTCCTAAAG-1），连接2个目标位点（ACCTGCGTGTCATGTT-1、GCAACACACTAGAACT-1）。两目标同时满足CD8/细胞毒和初始/记忆T标签，一个还满足浆细胞/髓系标签，另一个满足宽松Treg标签。不能把这2个位点解释成已分辨的耗竭T细胞或可重复的特异互作群。

## 空间统计

预设90项检验中，75项可检验，15项因来源不足、目标过少、恒定结果或条件零分布退化而不可检验。**无BH校正q<0.05的结果。主分析30项检验的最小q=0.486；联合90项的最小q=0.423。**

使用真实Visium六角阵列的直接相邻关系，排除自身，每个位点最多6个邻居。主结局为每个来源位点的目标邻居数/全部保留邻居数，对比高表达与其余检出来源的平均差。补充结局为是否存在目标邻居；敏感性检验排除自身满足同一目标标签的来源位点。所有来源保留零邻接，不把边或目标位点重复计为独立病例。

在全切片3×3坐标分区×区内候选来源UMI中位数分层内置换高表达标签，保留各层高表达数。标签组合空间≤4,999时完整枚举，否则4,999次随机置换；双侧p检验观测值相对**精确条件期望**的偏差，该期望可以不为0。至少3个高表达及3个比较来源、至少3个目标位点才可检验；恒定或退化零分布不输出p。不可检验p为缺失，仅在预设BH校正族中以1占位。

CHC20 cDC2主分析为负向差异−3.48个百分点，未经校正p=0.0162，主分析q=0.486，不能称显著富集或显著减少。观察邻接较多也不等于PGAM5高表达组选择性增加该邻域，空间组织与测序深度条件化会改变零假设期望。

分层仅近似控制局部组织与深度；共享邻居和残余空间依赖仍可能存在。仅2位患者，位点检验不能代替独立患者层面的重复证据。完整统计：[spatial_conditional_tests.csv](spatial_conditional_tests.csv)。

## 候选配体–受体RNA支持

从固定LIANA 1.10.0 omni_resource快照核验15对候选配体–受体。下表仅描述相邻来源与目标的RNA共同检出：来源配体>0，目标全部受体亚基>0。未运行CellChat/LIANA通信求解器，未进行配体–受体显著性检验。边得分为sqrt(来源配体CP10k×目标最小受体亚基CP10k)，包括零支持的全部合格边，不是互作概率或通路活性。

| 样本 | 目标RNA类别 | 候选配体–受体 | 共同检出边/全部合格边 | 支持来源数 | 支持目标数 |
| --- | --- | --- | --- | --- | --- |
| CHC20 | FCN1/VCAN炎症单核/髓系相关 | SPP1 → CD44 | 50/98 | 25 | 43 |
| CHC23 | FCN1/VCAN炎症单核/髓系相关 | CXCL12 → CXCR4 | 10/14 | 5 | 8 |
| CHC20 | MZB1/JCHAIN浆细胞相关 | CXCL12 → CXCR4 | 3/15 | 3 | 3 |
| CHC23 | MZB1/JCHAIN浆细胞相关 | CXCL12 → CXCR4 | 11/19 | 5 | 9 |
| CHC20 | CCR7/TCF7初始记忆T相关 | SPP1 → CD44 | 8/12 | 7 | 6 |
| CHC23 | CD8/细胞毒T相关 | CXCL12 → CXCR4 | 4/5 | 3 | 4 |
| CHC23 | 耗竭支持T相关 | PVR → TIGIT | 2/2 | 1 | 2 |

SPP1–CD44及CXCL12–CXCR4可作为部分候选邻域的RNA线索。耗竭支持邻域中的PVR–TIGIT只有1个来源和2个目标，尽管2/2边共同检出，也不是重复验证或PGAM5促进耗竭的证据。某些T相关混合位点检出CSF1R等髓系受体，不能将这些受体归属到T细胞。详细结果：[候选配体–受体](candidate_LR_by_immune_subgroup.csv)、[数据库定义](selected_candidate_LR_pairs.csv)。

## 图表

图中所有目标均应读作“RNA相关位点”。空间图使用原始病理图与真实坐标；金色环为PGAM5高表达来源，蓝色为目标RNA门限，紫色线为直接邻接。标志dotplot为门限内检出比例及包含零的平均表达，不能作为独立于门限的细胞身份验证。统计热图标联合90项BH q；配体–受体热图中灰色为无合格边，黄色0表示有边但未共同检出，稀少边的100%需结合原始分母。

- [免疫亚群邻接覆盖比较 PNG](figures/immune_subgroup_neighborhood_coverage.png) · [PDF](figures/immune_subgroup_neighborhood_coverage.pdf)
- [空间条件检验与q值 PNG](figures/immune_subgroup_spatial_statistics.png) · [PDF](figures/immune_subgroup_spatial_statistics.pdf)
- [CHC20全部亚群空间邻接 PNG](figures/CHC20_immune_subgroup_spatial_neighborhoods.png) · [PDF](figures/CHC20_immune_subgroup_spatial_neighborhoods.pdf)
- [CHC23全部亚群空间邻接 PNG](figures/CHC23_immune_subgroup_spatial_neighborhoods.png) · [PDF](figures/CHC23_immune_subgroup_spatial_neighborhoods.pdf)
- [CHC20标志基因dotplot PNG](figures/CHC20_immune_subgroup_marker_dotplot.png) · [PDF](figures/CHC20_immune_subgroup_marker_dotplot.pdf)
- [CHC23标志基因dotplot PNG](figures/CHC23_immune_subgroup_marker_dotplot.png) · [PDF](figures/CHC23_immune_subgroup_marker_dotplot.pdf)
- [候选配体–受体RNA共同检出 PNG](figures/immune_subgroup_candidate_LR_codetection.png) · [PDF](figures/immune_subgroup_candidate_LR_codetection.pdf)

## 数据定义与可复现执行

17套固定程序包含2套来源谱系程序、全部T基线及14个候选目标标签。所有标志均可在原始矩阵中找到。T状态要求全部T门限：CD3D/CD3E/CD3G/TRAC/TRBC1/TRBC2平均≥样本Q75且至少2个检出；状态程序≥T门限内Q75且至少2个状态标志检出，CD8/细胞毒额外要求CD8A或CD8B，Treg额外要求FOXP3。非T目标程序≥全样本Q75且>0，至少2个标志检出（浆细胞、pDC至少3个），并要求相应锚定标志。NK排除全部T门限；炎症单核门限额外要求泛髓系程序≥中位数。完整程序在[programs.json](programs.json)，锚定规则在[analyze.py](analyze.py)，统计契约在[analysis_config.json](analysis_config.json)。

独立核验脚本通过7,596项检查，涵盖原始HDF5计数/特征ID、库大小、门限、真实坐标、独立重建六角邻接、置换分层与标签数、精确条件期望、p/BH重算、配体–受体各亚基及来源/目标分母、严格标志复核、PNG/PDF文件。部分与此前发布文件的同一性检查仅在此前结果目录存在时执行；单独下载后有效检查数可略有不同。图形布局已人工查看。详见[audit_results.json](audit_results.json)。

下载包包含代码、全部7套PNG/PDF、CSV、标志子集AnnData对象（原始标志计数、全库归一化值、注释、邻接及低分辨率病理图）、置换审计数组与固定配体–受体资源。未重复打包全部官方HDF5；prepare_inputs.py从已发布基础分析的固定提交下载缺失依赖并逐项验证SHA256。子集对象不可按子集库大小重新归一化。

解压后在此目录执行（本轮运行Python 3.12.14；完整复算会下载原始依赖）：

```bash
python -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
python prepare_inputs.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python analyze.py
OPENBLAS_NUM_THREADS=1 python marker_specificity.py
OPENBLAS_NUM_THREADS=1 python plot.py
OPENBLAS_NUM_THREADS=1 python audit.py
python build_report.py
python package_results.py
```

已保存当前运行版本和随机种子；统计含Monte Carlo，跨平台数值细节可能有微小差异。绘图不需要重跑统计，可直接读取下载包中的结果与子集对象执行plot.py。

## 来源

- [GSE245908官方GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245908)；GSM7850822、GSM7850823官方表达、坐标、比例因子及组织图。
- 基础分析固定提交：`9cccc00079d26e55eb9d4f1d746cca9394366b76`；[依赖清单](source_inputs/base_requirements.json)、[输入SHA256](input_manifest.json)。
- 关联单细胞GSE245906采用CD45+、panTCRαβ−、CD19−先天免疫分选，不能提供完整T/B竞争参考，因此本轮使用透明的固定RNA门限，不宣称完成单细胞分辨率映射。
- 配体–受体数据库固定快照：LIANA 1.10.0 omni_resource.csv；已保存本地快照与SHA256。

[下载完整结果包](GSE245908_PGAM5_immune_subgroup_neighbors_results.zip)（在GitHub文件页点击Download raw file；或使用聊天中的直接下载链接）。
