# GSE245908：PGAM5 高表达巨噬 RNA 区域与 T 细胞亚群候选的空间关系

本次在此前活化/耗竭分析基础上补充 CD8、CD4、Treg 和初始/记忆候选分类。**未发现 PGAM5 高表达组特异增强某一 T 亚群邻近覆盖：20 项可计算检验全部 BH q=1；另4项因有效对照不足不检验。**预设比较共24项。候选配体–受体共检出存在，但样本和边数较少，不能解释为巨噬细胞诱导 T 细胞耗竭或活化。

## 数据及观察单位

两张 HCC Visium 切片 CHC20、CHC23（GSM7850822/23），沿用基线全部4,193个质控点和完整基因库 UMI 归一化。每张切片分别比较，不把空间点或邻接边当作患者独立重复；不进行只有2位患者的总体推断。High/Other detected/Undetected源点分别为CHC20：37/110/284；CHC23：17/50/240，与上轮完全一致。源点为多标志支持的巨噬/髓系RNA富集空间区域，不能确认每个点中PGAM5来自哪一个巨噬细胞。

源分组定义及原先21条LR总体T邻居比较见[上轮分析](../GSE245908_PGAM5_Tcell_spatial_interactions/README.md)。本轮不调整源阈值、不重选样本；T门控仍为CD3D/E/G、TRAC、TRBC1/2平均log1p(CP10k)达本切片Q75、至少2标志检出且有保留邻居。CHC20 T富集点453，CHC23为205。

## 分类规则与局限

- CD8候选：T门内CD8A检出，CD8A/B平均归一化表达高于CD4表达。CD4候选：T门内CD4检出且高于CD8A/B均值。其余保留为身份未解析；四个身份类别（CD8、其余CD4、Treg、未解析）互斥并覆盖T门。
- Treg候选：CD4候选内FOXP3检出，FOXP3/IKZF2/IL2RA/CTLA4至少2基因检出。没有仅凭IL2RA或CTLA4判定Treg。
- CD8内的耗竭相关（PDCD1/LAG3/TIGIT/CTLA4/TOX/CXCL13）、早期活化（CD69/IL2RA/TNFRSF9/ICOS/CD40LG/IL2）、效应（IFNG/PRF1/GZMB/GZMA/CCL5/TBX21）、初始/记忆（CCR7/SELL/TCF7/LEF1/IL7R）程序，均以各自CD8候选内Q75、得分>0、至少2基因检出定义。

CD8状态允许重叠，不能相加为总细胞数；初始/记忆仅为联合转录程序，不能凭该面板严格区分初始、中央记忆、效应记忆。检查点也可反映活化，不等于确认终末耗竭。CD4转录本可来自混合点内髓系细胞，该比较是保守、操作性的RNA代理门控，不能替代单细胞注释。CD8身份检出不足时不从细胞毒程序直接推断CD8。标志均存在于完整输入。

分类敏感性揭示明显检测限制：仅要求CD8A的宽门控为30/32点；同时要求CD8A和CD8B的严格门控只有1/2点。Treg至少3标志的严格敏感性门控两张切片均为0点。因此Treg与精细CD8结论尤其不稳定，不将某门控阈值当作生物学真值，也不按结果选择门控。

## 亚群分布与High源点相邻边

| 亚群/状态候选 | CHC20 点数 | CHC23 点数 | CHC20 High→目标边 | CHC23 High→目标边 |
|---|---:|---:|---:|---:|
| 全部 CD8 候选 | 17 | 18 | 3 | 5 |
| CD8 耗竭相关 | 1 | 3 | 0 | 2 |
| CD8 活化 | 0 | 1 | 0 | 0 |
| CD8 细胞毒效应 | 2 | 5 | 0 | 4 |
| CD8 初始/记忆 | 2 | 3 | 0 | 2 |
| 其余 CD4 候选 | 198 | 71 | 34 | 8 |
| Treg 候选 | 1 | 10 | 0 | 1 |
| T 身份未解析 | 237 | 106 | 27 | 6 |

这些为混合空间点与有向邻接边数量，不是单个免疫细胞数量。CD8各状态行包含在“全部CD8”内且彼此可能重叠。CHC20 High→CD8仅3条边；CHC23仅5条。CHC23耗竭与初始/记忆的2条边为同一组目标，不能当作两份独立支持。

[真实空间位置图](figures/spatial_T_subtypes.png) · [邻近覆盖图](figures/subtype_neighbor_coverage.png)

## 空间统计

沿用真实Visium阵列六角第一圈邻接，每个源点先计算邻居覆盖比例（分母为所有质控保留邻居），再比较High与Other detected。无某亚群邻居的源点保留并计0；无CD8邻居的源点不进入连续CD8状态程序比较。

各切片3×3位置格及格内UMI中位数分层交换标签，保持各层High个数；主统计衡量观察差值偏离精确条件零期望的程度。可配置空间≤4,999时全枚举，否则4,999次Monte Carlo，固定种子24590842/43。普通位置/深度分层不能完全解决细空间相关和细胞混合，属探索性检验。

固定24项检验＝2切片×（8亚群覆盖＋4CD8邻居连续程序），仅20项具有可计算p值，统一BH校正；全部q=1。CHC23连续程序有4个High、仅2个Other源点，低于每组至少3个的规则，不检验。CHC20连续程序也只有3个High、13个Other，条件空间仅2种配置，4项p均1，检测力极低。

全部CD8候选覆盖：CHC20 High 1.35% vs Other 2.00%，条件p=0.4812；CHC23为5.39% vs0.67%，条件p=0.7512。CHC23未调整差异较大，但空间条件比较不支持High特异增强。没有合格边/检出也不表示对应细胞或通路生物学不存在。

[全部比较表](spatial_subtype_tests.csv) · [标注p/q的图](figures/subtype_spatial_conditional_comparison.png) · [置换支持度](permutation_support.csv)

## 候选配体–受体

沿用上轮从LIANA 1.10.0固定资源核对的21条候选，未运行CellChat/LIANA通讯算法。源为巨噬RNA富集点，目标为真实相邻T亚群候选；复合受体要求全部亚基检出，以最小亚基CP10k为代理。边得分=sqrt(源配体CP10k×目标受体代理CP10k)，零表达边保留在分母。缺少边标记NA，区别于存在边但没有共检出。

- CHC20 SPP1–CD44：High→CD8为3/3共检出边，仅3个源点和3个目标点；High→CD4_other为18/34。100%来自极少数边，不能理解为强或稳定通讯。
- CHC23 PVR–TIGIT：High→CD8为2/5共检出边，均来自同一个源点、两个目标点；在耗竭候选中为2/2。耗竭门控本身含TIGIT，属于相关定义，并非独立机制验证。这两条也计入效应和初始/记忆重叠状态。
- CHC23 CCL5–CCR5有2条High→CD8共检出边，CXCL9–CXCR3有1条；源点是混合组织区，不能确认CCL5等配体来自巨噬细胞。
- CHC23 High→Treg只有1条邻接边，预设21候选均无共检出；证据不足以判定互作存在或不存在。

本轮亚群LR仅提供完整21条候选的描述性RNA支持，不为稀疏边再做显著性检验或给出通讯概率。上轮总体T门21条LR比较也未发现校正后显著差异。方向表示选定的源配体/目标受体表达，并不证明实际分泌方向或PGAM5的调节作用。

[候选LR亚群图](figures/LR_subtype_codetection.png) · [SPP1–CD44 CD8边空间图](figures/SPP1_CD44_CD8_spatial_edges.png) · [全部21条亚群LR表](LR_by_subtype.csv)

## 核验与复现

1259项核验全部通过：输入SHA256、源分组冻结、原始全库UMI、面板得分、分类分区、几何邻接独立重建、逐边原始表达及复合亚基、所有组/亚群LR聚合、条件p及BH复算。详见audit_checks.csv。原始HDF5直接读取，不依赖scratch缓存。程序提供逐点分类、全部源/目标条形码、状态重叠标记及候选边得分。

在完整仓库内运行，需保留相邻基线与上轮结果目录；ZIP为本轮结果，未重复打包约数十MB原始输入。缺少基线输入时先执行上轮prepare_inputs.py（其固定哈希下载从已发布仓库恢复），然后：

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
export MPLCONFIGDIR=/tmp/gse245908_mpl XDG_CACHE_HOME=/tmp/gse245908_cache
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
python analyze.py
python audit.py
python plot.py
python build_report.py
```

本轮扩展没有确认细胞层面的机制或患者总体差异。需要含T细胞的匹配单细胞参考、更高空间分辨率或功能实验验证；同研究先天免疫分选单细胞数据不构成完整T细胞参考。
