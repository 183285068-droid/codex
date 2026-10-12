# GSE238264：Q75 PGAM5高表达巨噬RNA源与T细胞亚群的空间互作候选

**已完成7张切片的CD8/CD4/Treg候选分类、CD8耗竭/活化/效应/初始记忆程序、真实邻接和22条候选配体–受体描述。本轮未发现校正后显著的High组T亚群邻域富集。**HCC4R的High源有7条CD8候选邻接边，其中2条连接同一个耗竭相关目标；Treg仅1条边。稀疏且共享目标的RNA支持不能证明调控机制。

## 范围和分组

沿用原7位患者（HCC1R–HCC4R应答，HCC5NR–HCC7NR未应答）的治疗后切除Visium切片，固定14,946个QC点。每切片1位患者，点/边不能作独立患者重复。本轮不另做应答/未应答总体检验，原标签只用于背景。

PGAM5源门与原分析逐点相同：巨噬核心RNA Q75、至少3核心标志、泛髓系≥中位数且有保留邻居；PGAM5原计数>0的巨噬RNA源内log1p(CP10k)Q75定义High。源数依次8/1/10/23/4/36/2，总84。PGAM5零值未删除。未使用此前Q20/Q25等阈值替代本轮Q75。

从原始完整H5获取所需基因计数，重复symbol合并；全库UMI作CP10k分母。来源、原QC和患者标签见[基线](../GSE238264_PGAM5_spatial_response/README.md)，免疫总体候选见[前一轮](../GSE238264_PGAM5_Q75_immune_neighbors/README.md)。本轮为空间RNA代理，没有单细胞去卷积或确定细胞身份。

## 系统分类

T-RNA富集门沿用CD3D/E/G、TRAC、TRBC1/2的切片Q75、至少2谱系标志及有邻居。

- CD8候选：T门内CD8A检出，CD8A/B平均log1p(CP10k)>CD4表达。CD4候选：T门内CD4检出且CD4表达>CD8A/B均值；其他保留为身份未解析。CD4 RNA可来自混合点内髓系，表达相对比较是保守操作规则，不是纯细胞判定。
- Treg候选：CD4候选内FOXP3检出，并在FOXP3/IKZF2/IL2RA/CTLA4至少2基因检出。其余CD4为CD4_other；CD8、CD4_other、Treg、未解析四身份互斥并覆盖全部T门。
- CD8耗竭：PDCD1/LAG3/TIGIT/CTLA4/TOX/CXCL13。
- CD8早期活化：CD69/IL2RA/TNFRSF9/ICOS/CD40LG/IL2。
- CD8效应：IFNG/PRF1/GZMB/GZMA/CCL5/TBX21。
- CD8初始/记忆：CCR7/SELL/TCF7/LEF1/IL7R。

四状态在每张切片CD8候选内要求程序均值≥Q75、>0、至少2基因检出；状态可重叠，不把同一目标重复计作独立细胞。初始/记忆是合并程序，不能严格区分初始、中央记忆、效应记忆或驻留状态。免疫检查点也可能反映活化，不能等同功能衰竭或终末耗竭。完整基因覆盖和状态阈值公开。

### 每切片候选点数

| 候选亚群 | HCC1R | HCC2R | HCC3R | HCC4R | HCC5NR | HCC6NR | HCC7NR |
|---|---:|---:|---:|---:|---:|---:|---:|
| All_T | 591 | 280 | 151 | 566 | 82 | 12 | 39 |
| All_CD8 | 62 | 30 | 10 | 54 | 0 | 0 | 0 |
| CD8_exhaustion | 2 | 0 | 0 | 11 | 0 | 0 | 0 |
| CD8_activation | 2 | 1 | 1 | 2 | 0 | 0 | 0 |
| CD8_effector | 11 | 7 | 3 | 13 | 0 | 0 | 0 |
| CD8_naive_memory | 12 | 3 | 3 | 14 | 0 | 0 | 0 |
| CD4_other | 105 | 45 | 70 | 216 | 5 | 1 | 8 |
| Treg | 0 | 0 | 1 | 4 | 0 | 0 | 0 |
| T_unresolved | 424 | 205 | 70 | 292 | 77 | 11 | 31 |

### Q75 High源到候选的有向邻接边数

| 候选亚群 | HCC1R | HCC2R | HCC3R | HCC4R | HCC5NR | HCC6NR | HCC7NR |
|---|---:|---:|---:|---:|---:|---:|---:|
| All_T | 8 | 3 | 3 | 83 | 3 | 1 | 0 |
| All_CD8 | 2 | 0 | 1 | 7 | 0 | 0 | 0 |
| CD8_exhaustion | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| CD8_activation | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| CD8_effector | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| CD8_naive_memory | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| CD4_other | 2 | 0 | 0 | 38 | 1 | 0 | 0 |
| Treg | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| T_unresolved | 4 | 3 | 2 | 37 | 2 | 1 | 0 |

“点”与“边”不是单细胞数量；All_T包括所有身份，All_CD8包括全部CD8状态，状态之间可能共享同一目标，因此不能跨行相加。没有合格直接边不代表不存在远距离或功能互作。

## 分类敏感性及证据边界

CD8A-only宽门控的候选点为67/40/27/108/1/0/0；CD8A/CD8B同时检出的严格门控为4/2/4/5/0/0/0，**各切片High源均无严格CD8门直接边**。这表明CD8身份支持稀疏且对门控敏感，不能将未应答切片的主CD8门0点解释为组织中没有CD8细胞，也不能据此证明治疗组差异。

Treg要求≥3标志的严格门在所有切片均为0点，说明本轮Treg候选亦很有限。另在CD8内保留去CTLA4/TIGIT的耗竭敏感性门（PDCD1/LAG3/TOX/CXCL13/HAVCR2/ENTPD1）逐点标签，不新增p值；HAVCR2等仍可能来自混合髓系，因此该门也不保证耗竭特异。

[门控敏感性表](gate_sensitivity.csv) · CD8_state_thresholds.csv · 各切片T_subtypes.csv.gz

## 空间条件统计

用真实Visium六角第一圈邻居，每源计算目标覆盖（目标候选邻居/全部保留邻居）；没有该亚群邻居计0。连续CD8状态结局仅在至少1个CD8候选邻居的源上计算，没有CD8邻居记NA，不填0。每源先汇总，避免把共享边当作独立观察。

High vs Other detected在各切片3×3位置网格及格内UMI中位数分层交换标签，保留分层High数；双侧检验围绕精确条件期望。小配置空间全枚举，大空间4,999次固定种子置换。每组源数≥3且有可交换标签才检验。该探索性检验假设层内可交换，不能完全解决剩余空间相关、细胞混合或临床混杂。

预设91项：7×（9亚群覆盖＋4连续CD8程序），联合BH，不可检验p仅在校正时设1，p/q保留NA。本轮42项不可检验，49项可计算（含18项所有分数0），31项有变异；全部有效q=1，无校正后显著关联。HCC2R/HCC7NR High数不足且标签不能交换。连续CD8程序只有HCC4R可比较（6 High、11 Other源），但条件空间仅2种配置，检测力极低。不能据不显著认定相同，也不能靠放宽标签去追求更小p。

[空间p/q图](figures/spatial_T_tests.png) · [全部检验](spatial_T_subtype_tests.csv) · [置换支持度](permutation_support.csv) · [覆盖比较图](figures/T_subtype_neighborhood_coverage.png)

## 候选配体–受体RNA支持

22条资源支持候选覆盖检查点、共刺激、趋化及其他信号，见selected_LR_pairs.csv。全部逐项核对LIANA 1.10.0 omni_resource快照，未运行LIANA/CellChat通讯算法。本轮只提供描述性候选支持，不额外对稀疏LR亚群边做显著性推断。

源配体和相邻目标受体以CP10k计算；复合受体全部亚基要求检出，最低亚基为代理。边RNA兼容得分=sqrt(源配体×受体代理)，零表达边保留分母，每源先平均再组内平均。NA表示没有合格目标边，与有边但无共检出不同。RNA方向是分析约定，不证明分泌方向或细胞来源，尤其CXCL12可来自源点内成纤维等细胞。

- HCC4R：High→全部CD8有7边（6源、6目标），CXCL12–CXCR4为5/7、SPP1–CD44为2/7、TGFB1–TGFBR1/TGFBR2为2/7共检出。
- HCC4R：耗竭相关CD8有**2源连接同1个目标**的2条边；CXCL12–CXCR4及TGFB1复合受体均2/2，PVR–TIGIT、NECTIN2–TIGIT、LGALS9–HAVCR2各1/2。这不是两个独立耗竭细胞或重复机制证据；耗竭定义含TIGIT等，本身与部分受体支持相关。
- HCC4R：Treg仅1源→1目标边，NECTIN2–TIGIT、SPP1–CD44、CXCL12–CXCR4各1/1。百分比100%来自同一条边，不能解释为强稳定通讯。
- HCC1R：High→CD8为2源→2目标的2边，SPP1–CD44为2/2，CXCL12–CXCR4为1/2。HCC3R仅1条CD8边，CXCL12–CXCR4为1/1。
- 所有切片Q75 High源未直接邻接多标志活化或效应CD8支持点；不能解释为PGAM5抑制其活化/效应。

[HE候选互作连线](figures/HCC4R_HE_candidate_T_LR_edges.png) · [候选LR热图](figures/T_subtype_LR_codetection.png) · [全部亚群LR表](LR_by_T_subtype.csv) · 各切片T_LR_edges.csv.gz（源/目标条形码及22条分值）

## HE融合图

- [HCC1R：HE T亚群空间图](figures/HCC1R_HE_T_subtypes.png)
- [HCC2R：HE T亚群空间图](figures/HCC2R_HE_T_subtypes.png)
- [HCC3R：HE T亚群空间图](figures/HCC3R_HE_T_subtypes.png)
- [HCC4R：HE T亚群空间图](figures/HCC4R_HE_T_subtypes.png)
- [HCC5NR：HE T亚群空间图](figures/HCC5NR_HE_T_subtypes.png)
- [HCC6NR：HE T亚群空间图](figures/HCC6NR_HE_T_subtypes.png)
- [HCC7NR：HE T亚群空间图](figures/HCC7NR_HE_T_subtypes.png)

每切片9面板：All_T、全部CD8、四CD8状态、CD4_other、Treg、未解析。蓝色半透明为目标RNA候选，橙圈为固定Q75 High源。按官方tissue_hires_scalef转换原pixel_row/col，点半径来自spot_diameter_fullres，无手动旋转或平移；所有坐标位于匹配图像范围。HE不能替代细胞身份验证，未自行标注TLS/恶性/病理边界。

## 核验和复现

4873项核验通过：SHA256、源分组冻结、完整库归一化、所有原始面板/门控/身份分区、独立KDTree几何、每源和亚群邻域汇总、原始LR及受体亚基/边数、条件p和91项BH。HE官方配准检查另见HE_registration_checks.json。

需完整仓库相邻基线及GSE245908固定资源快照，原7样本解压目录由GSE238264_INPUT指定。结果ZIP不重复大原始文件；代码、逐点/逐边数据、图表和审计均提供。

```bash
export GSE238264_INPUT=/你的路径/inputs
export MPLCONFIGDIR=/tmp/gse238264_mpl XDG_CACHE_HOME=/tmp/gse238264_cache
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
python analyze.py
python audit.py
python plot.py
python plot_statistics.py
python plot_edges.py
python package_results.py
```

requirements.txt固定环境。本轮结论仍是探索性空间RNA候选关系，不能确立PGAM5高表达巨噬细胞诱导T细胞耗竭、调节Treg或造成治疗抵抗。
