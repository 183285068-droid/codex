# 原15群版本C7：PGAM5虚拟敲除

## 结论

已在GSE202642原15群版本的C7（MKI67/TOP2A增殖型巨噬细胞）完成scTenifoldKnk虚拟敲除。403个细胞来自7个患者标签，其中79个检出PGAM5 RNA；包含所有PGAM5零值细胞。

种子202642得到26个非靶基因的模型FDR<0.05候选；种子42得到0个，因此没有在两次运行中均显著的下游候选。不能据此建立稳定的PGAM5因果调控结论，也不能判定敲除会增强或减弱抗肿瘤免疫或增殖。

主运行的领先候选包括HSP90AA1、DNAJB1、HSPA1B、PLAUR、NFKBIA、FOSB、PPP1R15A、HSPH1、NR4A1、DUSP1等。预定义热休克/蛋白稳态程序的探索性富集BH q=4.19e-15，即时反应程序q=8.87e-8；二者在第二种子均无显著基因支持。这些富集仅描述主运行候选，不能作为重复验证的通路效应。

## 方法

完整36,591基因原始计数提取自已有GSE202642全细胞原始对象，逐细胞与原15群C7名单、PGAM5计数和全部2,102个便携源特征计数核对一致。选择在≥5%目标细胞中检出的基因，排除MT-/RPL/RPS，按log1p(CP10k)Seurat归一化离散度选前1,500基因，再加入PGAM5，共1,501基因。网络范围因此受筛选限制，未建模全部转录组。

scTenifoldpy 0.5.1默认scTenifoldKnk行置零方法：每次10个PC网络，每网络有放回抽取300个细胞、PC回归3成分、保留绝对边权高于95%分位的边；CP张量分解rank=5、tol=1e-5、最多1,000轮，采用包默认1位小数精度；流形对齐30维。输入网络的CPM用原始全转录组总计数作为分母（明确覆盖包QC阶段按筛选后库大小归一化的默认结果），没有把PGAM5表达直接删零后重新计算DE。

WT网络中的PGAM5出边在种子202642为293条、种子42为28条，网络估计显著依赖重采样和张量初始化；默认张量精度、稀疏化也可能影响弱边。敲除将PGAM5所在行连接置零，其他行完全一致。相同WT网络对齐作为无敲除对照，0个显著基因，最大位移约7.17e-15。两次运行非靶基因位移排序Spearman rho=0.5795。

“Distance”是无方向的流形位移；包输出“FC”是位移平方相对于背景期望的比值，**不是RNA表达fold change**。输出p值基于包的卡方模型假设、BH校正，是网络模型统计量；没有真实KO细胞或独立患者级实验对照，不是敲除后表达差异的实验p值。PGAM5本身受强制干预，表中保留且标记为目标，不计入下游候选或富集证据。

程序使用之前免疫程序定义加两套明确列出的热休克及即时应答基因集，均保留PGAM5成员；检验下游证据时单独去除这个自动受干预的靶点，背景为其他1,500个实际建模基因，单侧Fisher、14个程序BH校正。不是完整GO/KEGG富集。没有额外控制患者、测序深度或细胞周期；患者标签组成见sample_summary.csv。

## 结果文件

- primary_candidate_genes.png/pdf：主运行候选位移图（标注第二种子未复现）。
- seed_sensitivity.png/pdf、seed_comparison.csv：重复种子比较。
- differential_regulation_seed*.csv：两次完整模型统计；primary_exploratory_non_target_hits.csv：主运行26个候选。
- exploratory_program_enrichment.csv、program_definitions.json：探索性程序分析。
- WT_network_seed*.csv.gz、KO_network_seed*.csv.gz：敲除前后网络。
- no_KO_control.csv、validation_summary.json：无敲除对照和核查。
- C7_raw.h5ad：目标403细胞的完整原始计数；network_raw_counts.csv.gz：1,501建模基因原始计数；gene_selection.csv：筛选记录。

## 复现

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
OPENBLAS_NUM_THREADS=4 python run.py
OPENBLAS_NUM_THREADS=4 python audit_report.py
```

随包已有C7_raw.h5ad，无需重新下载GEO；run.py会跳过已存在的种子统计。要从头重算，请在新目录复制C7_raw.h5ad、run.py、audit_report.py、program_definitions.json及requirements.txt后执行。稀疏特征筛选和重采样固定种子；底层特征求解器及不同平台可能产生微小数值差异。所有代码、参数、原始计数和完整结果保留，正结果及未复现情况均报告。
