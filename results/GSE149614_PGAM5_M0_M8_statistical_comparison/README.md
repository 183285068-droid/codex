# GSE149614：PGAM5 亚群表达丰度统计比较

## 结果与图

**主要结论：未检出 BH 校正后 q<0.05 的亚群间差异。** 可检验的配对共 12/36 对，最小 q=0.5625。这不等同于证明各亚群表达相同；仅有10个样本，部分群覆盖不足，检验效能有限。

细胞合并描述中，M1（MKI67/TOP2A 增殖型）平均 CP10k=0.1070、检出率9.13%，均为9群最高。
细胞层面探索性检验有 6/36 对 q<0.05。**这些细胞层面的星号不能替代样本层面的证据。**

- [主图：样本表达分布和全部两两比较显著性](outputs/01_sample_aware_expression_significance.png) · [SVG](outputs/01_sample_aware_expression_significance.svg)
- [补充图：合并表达丰度和细胞层面显著性](outputs/02_exploratory_cell_expression_significance.png) · [SVG](outputs/02_exploratory_cell_expression_significance.svg)
- [柱状图：直接标注显著性括号、星号与q值](outputs/03_barplot_with_significance.png) · [SVG](outputs/03_barplot_with_significance.svg)
- [主分析完整36对统计](outputs/primary_sample_paired_comparisons.csv)
- [细胞层面完整36对统计](outputs/exploratory_cell_MannWhitney_comparisons.csv)
- [各样本各亚群统计](outputs/sample_cluster_statistics.csv) · [样本细胞数](outputs/sample_cluster_cell_counts.csv)
- [各群表达汇总](outputs/cluster_expression_summary.csv)
- [敏感性分析](outputs/sensitivity_sample_comparisons.csv) · [检出率次要比较](outputs/secondary_sample_detection_comparisons.csv)

## 数据与统计方法

来源：同仓库 `results/GSE149614_PGAM5_M0_M8_expression_supplement/inputs/cell_expression_coordinates.csv.gz`。
沿用5,605个细胞、10个原发HCC肿瘤T样本、M0–M8原标签；不重新聚类、不插补、不筛掉表达零值。
PGAM5 在136/5605个细胞中检出。M0–M8是聚类编号，M1编号不代表经典M1极化。

1. CP10k = PGAM5原始计数 / 全转录组文库总计数 × 10,000。每个细胞先取 log1p，再对同一样本、同一亚群内所有细胞求均值，得到主要终点。
2. 主分析要求每个样本-亚群至少10个细胞，每对比较至少5个共享合格样本；缺失亚群不当作表达0。每个样本只贡献一个汇总值，且在比较中等权。每对的完整样本名单和均值差写入CSV。
3. 对配对样本差值进行双侧精确符号翻转/配对标签交换检验，枚举全部2^n种标签交换，以绝对均值差为检验量。在零假设下要求配对标签可交换（差值符号对称）；本分析为观察性探索，并非随机试验。
4. 主分析36个预定比较组成一个Benjamini–Hochberg（BH）校正族。不能检验的比较以P=1参与校正以保留整个比较族，但输出P/q留空并标为NT，不能解释为不显著。
5. 图中的 `*`、`**`、`***`、`****` 分别指 q<0.05、0.01、0.001、0.0001；`ns`指 q>=0.05；`NT`指合格配对样本不足。主图同时标明配对样本数。
6. 敏感性分析：最小细胞数1/10/20；主要终点、算术平均CP10k、按计数汇总的pseudobulk CP10k。每个终点和阈值组合分别对36对做BH，不据此挑选显著结果。检出率的样本级比较另列为次要分析，单独校正。
7. 细胞层面为补充探索：双侧Mann–Whitney U（渐近法，校正并列值及连续性）比较表达分布，36对BH校正；列出秩二列相关效应量。另保存Fisher精确检验的检出率比较。合并Kruskal–Wallis H=45.59，P=2.839e-07，亦不调整样本内相关性。
8. 柱状图直接标注M1与其余8群的比较（含不显著结果），q值仍来自全部36对比较的BH校正，没有缩小比较族。其余群之间的两两比较保留在完整热图和CSV中。括号检验的是细胞表达分布，柱高描述算术平均CP10k；星号不应解读为样本层面的均值检验。

## 细胞层面探索性显著比较

| 比较 | 原始P | BH校正q | 标注 |
|---|---:|---:|---|
| M0 vs M1 | 1.044e-09 | 3.759e-08 | **** |
| M1 vs M2 | 3.342e-07 | 4.011e-06 | **** |
| M1 vs M3 | 0.001253 | 0.009019 | ** |
| M1 vs M4 | 0.000106 | 0.0009536 | *** |
| M1 vs M5 | 1.547e-08 | 2.784e-07 | **** |
| M1 vs M7 | 0.003641 | 0.02185 | * |

## 解释限制

PGAM5非常稀疏，所有群的细胞中位数为0，因此仅用中位数箱线图会隐藏差别。主图展示样本均值的分布，补充图保留原始均值量纲。
M3有91.6%细胞来自HCC01T，M6是巨噬/T混合RNA群，M8仅42个细胞且1个检出；群名称和功能均属探索性注释。
样本细胞数阈值10是本次分析约定的稳定性门槛，不是经验证的通用标准；同时交付1和20阈值结果供检查。
M1是已有结果中观察到最高的群，不另行挑选其8个比较组成更小校正族。
这里比较的是标准化RNA表达，不是蛋白丰度、通路活性或因果作用；未调整细胞状态、测序深度等其他混杂因素。

## 复现与核查

```bash
python -m pip install -r requirements.txt
MPLCONFIGDIR=/tmp/pgam5-matplotlib python analyze.py
```

输入保存在inputs目录，并核验SHA-256；9群均值与之前的CSV一致。
精确检验通过与SciPy全枚举置换检验及已知答案核查，BH通过已知例子核查。核查结果见 `outputs/verification.json`。
图为英文标签，中文方法在此说明。未将本次统计结果改写回原分析目录。

## 方法参考

- Zimmerman et al. A practical solution to pseudoreplication bias in single-cell studies. Nat Commun 12, 738 (2021). https://doi.org/10.1038/s41467-021-21038-1
- SciPy Mann–Whitney U: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html
- SciPy BH correction: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html
