# GSE149614 M1：PGAM5虚拟敲除与可视化

固定GSE149614的10个原发HCC肿瘤T样本、最近9群版本中的**M1：MKI67/TOP2A增殖型巨噬细胞**。共208个细胞，19个检出PGAM5，189个零值全部保留；覆盖10个样本/患者编号，最大样本贡献25.48%。不重新聚类、不按PGAM5检出筛除细胞，也不将原始观测UMAP称为敲除后UMAP。

## 模型结果

完成5次稀疏主模型（网络保留阈值q=0.95）、5次关闭全局稀疏化的敏感性模型（网络阈值q=0），以及10次相同WT网络的无敲除数值对照。网络阈值q与统计学BH q是不同参数。

| 模型 | 各次模型BH q<0.05的下游基因数 | 至少4/5次显著 | 5/5次显著 | 预设技术稳定性 |
|---|---|---:|---:|---|
| 稀疏主模型 | 28, 30, 28, 29, 26 | 24 | 20 | 通过 |
| 关闭稀疏化 | 26, 27, 27, 26, 23 | 19 | 15 | 通过 |

主模型有效干预5/5次；最低/中位两两位移排序rho=0.817/0.859，前50 Jaccard最低/中位=0.639/0.770。敏感性模型有效干预5/5次；最低/中位rho=0.888/0.908，Jaccard最低/中位=0.538/0.695。

两方案中位位移的排序rho=0.848，TOP20交集10/20；[两方案均至少4/5次显著的共同候选](cross_parameter_recurrent_candidates.csv)共12个。技术种子重复和网络参数敏感性不能替代独立患者或真实敲除实验验证。通过本轮技术稳定标准也不证明真实因果调控、PGAM5蛋白作用或抗肿瘤功能变化。

## 图形与统计

- [M1定位及敲除前实际PGAM5表达UMAP](01_M1_target_and_observed_PGAM5_UMAP.png)、[样本构成](02_M1_sample_representation.png)。
- **[主模型TOP20条形图](03_PGAM5_virtual_KO_TOP20.png)**、[主模型TOP20重复热图](05_TOP20_cross_seed_heatmap.png)、[主模型稳定性和无敲除对照](04_repetition_stability_and_controls.png)。
- [关闭稀疏化TOP20](unsparsified/03_PGAM5_virtual_KO_TOP20.png)、[重复热图](unsparsified/05_TOP20_cross_seed_heatmap.png)、[稳定性及对照](unsparsified/04_repetition_stability_and_controls.png)。
- [网络参数敏感性对比](06_model_configuration_comparison.png)。所有PNG均有同名PDF。
- [主模型完整基因汇总](all_genes_perturbation_summary.csv)、[TOP20精确数值](TOP20_downstream_genes.csv)、[重复候选](recurrent_candidates.csv)、[跨参数逐基因对比](cross_parameter_gene_comparison.csv)。

TOP20条形表示5次模型的中位**无方向网络位移**，误差线为Q1–Q3；蓝色表示5/5次模型BH q<0.05，橙色4/5，灰色≤3/5。Q1–Q3不是患者置信区间。位移与模型FC均不是实际RNA表达fold change，不能据此将基因判为上调或下调。若有效信号少于20个，只展示通过数值噪声过滤的基因。

| 排名 | 下游基因 | 中位无方向网络位移 | 模型q<0.05次数 |
|---:|---|---:|---:|
| 1 | S100A6 | 3.96226e-06 | 5/5 |
| 2 | S100A4 | 3.48411e-06 | 5/5 |
| 3 | ANXA1 | 3.09428e-06 | 5/5 |
| 4 | LGALS1 | 2.90533e-06 | 5/5 |
| 5 | HMGB2 | 2.82689e-06 | 5/5 |
| 6 | MNDA | 2.65817e-06 | 5/5 |
| 7 | CRIP1 | 2.62242e-06 | 5/5 |
| 8 | HSP90AA1 | 2.54606e-06 | 5/5 |
| 9 | C15orf48 | 2.26595e-06 | 5/5 |
| 10 | HSPH1 | 2.04601e-06 | 5/5 |
| 11 | DNAJB1 | 2.03766e-06 | 4/5 |
| 12 | HSPA1B | 1.99633e-06 | 4/5 |
| 13 | HSPE1 | 1.73321e-06 | 5/5 |
| 14 | UBE2C | 1.51345e-06 | 5/5 |
| 15 | FCGR2B | 1.44607e-06 | 5/5 |
| 16 | LINC01272 | 1.43592e-06 | 5/5 |
| 17 | TUBA1A | 1.43091e-06 | 5/5 |
| 18 | ZFAND2A | 1.38885e-06 | 3/5 |
| 19 | GAPDH | 1.33367e-06 | 5/5 |
| 20 | HK2 | 1.31244e-06 | 5/5 |

PGAM5始终保留在输入、完整统计和[被强制干预靶基因表](PGAM5_forced_target_summary.csv)中；由于其被直接干预，不计入独立下游TOP20或稳定候选。

## 方法与预设标准

scTenifoldpy 0.5.1：PC回归构建30个有放回重采样网络（每个208个细胞、3成分），CP张量秩5、最大1,000次、容差1e-5，网络归一化后保留6位小数并转置、去自环，再将PGAM5对应调控行设为0；其他行与WT完全不变。30维流形对齐计算无方向位移，模型卡方p及BH q描述模型内相对扰动，不是实际敲除实验差异表达检验。

原始矩阵包含25,712基因。按检出率≥5%排除MT-/RPL/RPS，取log1p(CP10k)的Seurat标准化离散度前1,500基因，再确保加入PGAM5，共1,501基因。网络使用CPM，分母始终为完整25,712基因原始文库总计数。没有额外回归患者、测序深度或细胞周期；M1仅19个PGAM5检出细胞，存在稀疏表达和网络推断限制。

两模型使用相同细胞、基因、种子149614/42/101/202/303及其余参数。预设标准：全部张量收敛、全部WT PGAM5行非零且位移高于数值噪声、全部无敲除对照无显著、所有两两rho≥0.8、所有前50 Jaccard≥0.5。模型BH q<0.05在至少4/5次出现作为重复候选；重复频率及中位q不是合并p值。

独立审计：主模型86/86项、敏感性模型86/86项通过。包括全部208个细胞和原始计数、完整文库分母、19检出/189零值、实际KO矩阵SHA256、仅PGAM5行置零、对齐位移、独立卡方p和BH q重算、数值噪声过滤、10次无敲除对照及TOP20排序。审计证明实现正确，不自动证明生物学结论。

## 下载与复现

- [主模型完整包](GSE149614_M1_PGAM5_virtual_KO_results.zip)：输入、5次网络及数值结果、图、代码、运行依赖源码快照。
- [关闭稀疏化敏感性包](GSE149614_M1_PGAM5_virtual_KO_unsparsified_results.zip)：独立输入、5次网络、图、统计和代码。两包独立解压，无须二进制合并。

GitHub打开ZIP点击Download raw file。原始计数不改变；输入M1_raw.h5ad为完整25,712基因，network_raw_counts.csv.gz仅为选中1,501基因，不可用其行和重新归一化。原分析[完整巨噬细胞矩阵](../GSE149614_primaryT_hierarchical_clustering/macrophages_annotated_full.h5ad)的哈希记录在input_summary.json。

Python3.12，建议4线程及≥8 GiB内存。依赖版本见requirements.txt。vendor/保存本轮实际执行的scTenifold源码及许可证、SHA256，run.py优先加载该快照，避免未来同名版本变化。

```bash
pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
export MPLCONFIGDIR="$PWD/.matplotlib"
python visualize.py
python audit.py
```

从头重算请在新目录复制inputs/、vendor/、脚本、sample_summary.csv和requirements.txt后运行prepare_inputs.py、run.py、visualize.py、audit.py。run.py跳过已有metrics.json的完成种子；重新运行全部模型需使用无models/的新目录。manifest.json与archive_verification.json提供完整性记录。
