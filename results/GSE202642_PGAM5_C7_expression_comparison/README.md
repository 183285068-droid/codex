# GSE202642：PGAM5表达丰度与C7定向比较

基于原始 `PGAM5_cluster_mean_expression_UMAP.pdf` 对应的12,135个细胞、7位患者、C0–C14原始分群；未重新聚类或更改细胞成员。

## 图和统计表

- [主柱状图：C7与其他巨噬细胞相关群](outputs/01_C7_macrophage_barplot_significance.png) · [PDF](outputs/01_C7_macrophage_barplot_significance.pdf) · [SVG](outputs/01_C7_macrophage_barplot_significance.svg)
- [对照柱状图：全部15个原始RNA群，混合群另行标记](outputs/02_C7_all15_RNA_groups_barplot_significance.png) · [PDF](outputs/02_C7_all15_RNA_groups_barplot_significance.pdf) · [SVG](outputs/02_C7_all15_RNA_groups_barplot_significance.svg)
- [患者层面图与比较结果](outputs/03_C7_patient_aware_comparison.png) · [PDF](outputs/03_C7_patient_aware_comparison.pdf) · [SVG](outputs/03_C7_patient_aware_comparison.svg)
- [各群表达汇总](outputs/cluster_expression_summary.csv)
- [C7定向细胞层面比较](outputs/C7_vs_others_cell_exploratory.csv) · [患者层面主分析](outputs/C7_vs_others_patient_primary.csv)
- [患者层面敏感性分析](outputs/C7_patient_sensitivity.csv) · [检出率次要分析](outputs/C7_patient_detection_secondary.csv) · [探索性Fisher检验](outputs/C7_detection_Fisher_exploratory.csv)

## 主要结果

C7共403个细胞、79个PGAM5检出细胞；平均CP10k=0.161737，平均log1p(CP10k)=0.108005，检出率=19.60%。
在巨噬细胞相关群的9个定向比较中，细胞层面有8个q<0.05；C4与C7不显著。全部15群中有12/14个定向比较显著，C4和C14不显著。
**这些是细胞层面的探索性显著性；同一患者内细胞不是独立生物学重复。**
患者层面的14个比较中，只有2个满足预设覆盖门槛；均未达到BH q<0.05。未通过门槛的比较标NT，不能解释成表达相同或不显著。

| 定向比较 | 原始P | BH q（14个比较） | 柱状图标注 |
|---|---:|---:|---|
| C7 vs C0 | 1.585e-17 | 7.398e-17 | **** |
| C7 vs C2 | 2.956e-14 | 1.035e-13 | **** |
| C7 vs C3 | 1.531e-26 | 2.143e-25 | **** |
| C7 vs C4 | 0.2482 | 0.2673 | ns |
| C7 vs C6 | 1.854e-18 | 1.297e-17 | **** |
| C7 vs C9 | 2.656e-13 | 7.438e-13 | **** |
| C7 vs C11 | 1.233e-07 | 2.157e-07 | **** |
| C7 vs C12 | 6.322e-12 | 1.475e-11 | **** |
| C7 vs C13 | 1.343e-08 | 2.686e-08 | **** |

## 表达量与统计定义

1. CP10k = PGAM5原始计数 / 原始全转录组文库总计数 × 10,000。**柱高为每个细胞先取log1p(CP10k)后求群内算术平均，保留全部零值，与原UMAP均值一致。** 同时导出原始量纲CP10k均值与检出率；两种均值不混用。C7的平均CP10k为15群最高，但平均log1p(CP10k)低于混合RNA群C14，不声称C7在所有量纲均最高。
2. 用户指定C7作为参照，事先定义C7与其他14个原始群的双侧Mann–Whitney U检验（渐近法、并列值和连续性校正）。BH对整个14比较族校正。主图只展示其中9个巨噬细胞相关群对比，**没有重新缩小校正族或重算更有利的q值**。U检验比较表达分布，星号不应当作对算术均值的独立生物学重复检验。
3. 主图保留C0、C2、C3、C4、C6、C7、C9、C11、C12、C13。依据原注释，C1为低复杂度混合RNA，C5为内皮RNA混合，C8为T/NK RNA混合，C10为粒细胞RNA混合，C14为肝细胞RNA混合，均待身份复核；它们不作为已确认巨噬细胞亚群命名，保留在15群对照图并以灰色斜线和问号标明。C0/C2仍保留原巨噬细胞相关注释和肝源/肿瘤RNA限定。此为基于原注释的展示分层，不删除原始细胞。
4. 患者层面主分析：每患者-亚群至少10个细胞，每对至少5位共享合格患者；每位患者贡献一个含零值的平均log1p(CP10k)，等权配对比较。不存在的亚群不补零。用双侧精确配对标签交换检验，枚举2^n种符号翻转并比较绝对均值差；零假设要求配对标签可交换/差值符号对称。
5. 患者比较同样按14个预定比较做BH；不能检验的比较以P=1参与校正，输出P/q留空并标NT。次要检出率和敏感性分析分别组成各自14比较校正族。敏感性涵盖细胞数阈值1/10/20、平均log1p(CP10k)、平均CP10k、汇总计数pseudobulk CP10k，不据此挑选显著结论。
6. 7位患者的独立性沿用原分析中用户已确认的样本映射，详见源目录原有HCC_patient_sample_mapping.csv。这里没有重新识别或额外验证患者身份。10个细胞、5位患者为本次操作门槛，不是通用生物学标准；原分析候选群方向一致性使用的20细胞规则另在敏感性表中保留。

星号：`* q<0.05`、`** q<0.01`、`*** q<0.001`、`**** q<0.0001`；`ns`表示q>=0.05；`NT`表示配对患者不足。

## 解释限制和复现

许多群由单一患者主导；细胞层面极小P值不能取代跨患者重复。原数据尚未校正环境RNA/双细胞，未控制测序深度、细胞周期等混杂。本结果衡量RNA表达，不是蛋白丰度、功能活性或因果作用，也不自动验证TCGA反卷积用途。

```bash
python -m pip install -r requirements.txt
MPLCONFIGDIR=/tmp/gse202642-mpl python analyze_C7.py
```

本目录仅发布统计汇总、图和代码，不重复发布逐细胞数据、患者映射或含患者标识的明细。输入读取同仓库相邻的 `GSE202642_PGAM5_subpopulation_assessment` 原目录，source_provenance.json记录源提交和SHA-256；请克隆完整仓库后运行。也可通过环境变量 `GSE202642_SOURCE_DIR` 指定原数据目录。运行将重生成汇总outputs和此README。
原始计数/全库分母已重算，15群均值、检出数与原表及UMAP映射逐项一致；精确置换和BH实现经过独立数值核对。outputs/verification.json保存核查结果。

方法参考：[SciPy Mann–Whitney U](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html)；[BH校正](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html)；[单细胞伪重复问题](https://doi.org/10.1038/s41467-021-21038-1)。
