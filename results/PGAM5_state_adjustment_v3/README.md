# PGAM5相关巨噬细胞状态与分层模型调整

完成跨患者一致性筛选、竞争谱系特异性筛选、有限细胞周期排除敏感性、分层反卷积及非细胞周期NMF状态探索。原有FDR<0.05、|log2FC|>=1差异表达表未更改；没有患者配对DE或测序深度匹配/协变量。

## 筛选结论

联合训练筛选：稳定候选1个（CENPK），谱系特异候选0个，去除固定细胞周期清单后稳定候选0个。CENPK在14位合格训练患者中10位方向一致，但其最大竞争群为Cycling_T_NK，未通过特异性筛选。

筛选使用训练患者中位表达比、方向一致性及竞争组分表达上限，并非新的统计显著性DE检验。完整阈值见protocol.json。每次患者留出重新筛选，少于5个候选时明确记为不可测试，没有把缺失当作零误差。此前21基因共识仅作为探索性对照；其中20个在本次有限细胞周期清单内，去除后仅剩HCFC1，未构成可测试多基因非增殖评分。该21基因来自此前含GSE189903的分析，不能宣称完全独立的外部特征发现。

## 分层反卷积

第一层按广义谱系同时估计总巨噬细胞与竞争细胞；第二层从混合表达中扣除拟合的非巨噬细胞贡献，再用两类巨噬细胞参考估计条件比例。未检出巨噬细胞参考按患者等权，并按细胞数合并其循环与非循环成员。最终目标=总巨噬细胞比例×条件目标比例。该形式仍然需要验证，不能因公式成立而认定丰度有效。

| 外部等RNA测试 | MAE（百分点） | r | 零目标预测95分位（%） |
|---|---:|---:|---:|
| direct_NNLS | 2.798 | -0.087493259393389 | 7.641 |
| hierarchical_NNLS | 1.700 | 不可计算 | 0.000 |

分层外部等RNA预测全部为0，无法识别非零目标；零目标误判变低不是验证成功。RNA权重模拟、内部留出及单核/DC/Mixed_APC挑战另见validation_summary.csv。只做NNLS，未执行CIBERSORTx、MuSiC或BayesPrism。20组匹配两状态参考的混合准确恢复，独立有界求解器的第二层系数与NNLS一致，MAE逐条重新计算；这不证明注释或参考的生物学准确性。

## 连续状态探索

在全部训练巨噬细胞中，排除PGAM5、固定细胞周期清单及MT/RPL/RPS基因，从训练表达中选取500个变异基因，拟合4个非负表达程序。仅依据训练数据挑选与连续PGAM5表达相关性最高的程序，使用其30个高载荷基因构建连续评分。分别执行151530→149614、149614→151530跨队列验证，再联合训练后复测189903。参数不按外部结果调整。

| 测试 | 合格患者 | 中位AUC（RNA检出） | 中位Spearman相关 | rho>=0.1患者比例 | 支持条件全部满足 |
|---|---:|---:|---:|---:|---|
| cohort_holdout_149_to_151 | 5 | 0.523 | 0.017 | 20.0% | False |
| cohort_holdout_151_to_149 | 9 | 0.437 | -0.032 | 0.0% | False |
| joint_external_retest | 2 | 0.431 | -0.032 | 0.0% | False |

模型基因、训练标准化参数与全NMF载荷已导出。训练相关不等于已验证的PGAM5驱动状态；PGAM5 RNA检出不是蛋白或功能金标准。未控制深度，因此深度相关信号仍可能存在。非巨噬细胞评分高于巨噬细胞的逐患者检查见program_support_checks.csv。NMF未收敛者明确记录，不能视为稳定基因程序。

## 使用限制与结果文件

当前不提供经验证可用于TCGA-LIHC丰度的签名，也不将探索性状态分数包装为细胞丰度。GSE189903仅2位患者有>=5个目标RNA检出巨噬细胞，且无符合规则的纯DC群；本轮复用了已看过的外部数据，不是新增独立验证。状态评分在bulk中的平台/单位校准也未完成。

- training_gene_stability_specificity.csv、per_donor_gene_contrasts.csv.gz：全部基因的训练一致性/特异性证据。
- donor_target_consistency_QA.csv：目标检出率、单计数阳性占比、循环标记占比的患者审计。
- panel_fold_testability.csv、fold_selected_genes.csv：每个折的候选数量与不可测试状态。
- predictions.csv、validation_summary.csv、abundance_gate_checks.csv：直接/分层模型完整结果。
- program_score_parameters.csv、*_NMF_basis.tsv、all_programs_training_association.csv：探索性状态评分参数与程序。
- program_validation_by_patient.csv、program_support_checks.csv、program_lineage_specificity.csv：跨队列和外部复测。
- heldout_cell_state_scores.csv.gz、program_cell_scores.csv.gz：可审计的逐细胞探索性分数。

复现：依赖上一轮三个harmonized_counts_QC.h5ad及其train_joint_reference.py，路径默认/workspace/scratch/PGAM5_myeloid_reference_v2。先python analysis.py，再python discover_programs.py，最后python report.py；OPENBLAS_NUM_THREADS=1。大型矩阵不打包。输入和代码哈希见来源清单。细胞周期清单是记录来源的有限人工转录列表及扩展，非完整GO注释；外部下载403/404未用未验证内容替代，不能把“未列入”理解为确定非增殖。

下一步更有价值的是新增独立、巨噬细胞/单核/DC覆盖充分的HCC队列，并以蛋白或功能证据定位PGAM5相关状态。现有结果不支持继续放宽阈值得到最终TCGA签名。

## 追加：训练谱系筛选后的非增殖程序

初始训练NMF高载荷出现ALB/APOA1/TTR后，追加仅依据训练数据的谱系条件：巨噬细胞患者等权均值不低于每个有支持的非巨噬细胞群，且巨噬细胞检出率>=10%，再执行同样的NMF流程。这是在训练标记审计后的探索性改进，不是预注册独立验证。

| 测试 | 合格患者 | 中位AUC | 中位Spearman相关 | rho>=0.1患者比例 | 全部条件满足 |
|---|---:|---:|---:|---:|---|
| cohort_holdout_149_to_151 | 5 | 0.686 | 0.101 | 60.0% | False |
| cohort_holdout_151_to_149 | 9 | 0.516 | 0.011 | 0.0% | False |
| joint_external_retest | 2 | 0.607 | 0.051 | 0.0% | False |

联合模型的研究用30基因程序含C1QA/B/C、CTSB、CTSD、PSAP、TREM2、SPP1、CXCL9/10等巨噬细胞相关基因，减少了初始程序的肝细胞相关基因，但不能证明环境RNA或双细胞问题已纠正；但这些基因组成不能自行证明PGAM5调控、功能因果或bulk特异性。未做GO/通路富集检验。基因列表和分数参数见lineage_restricted_programs/program_score_parameters.csv。

不把未达到关联/特异性/患者覆盖要求的程序命名为已验证PGAM5状态签名。所有程序和负结果均保留，后续可结合独立蛋白或功能证据检验；此轮未运行TCGA状态评分。复现追加程序：在report.py之前运行python discover_lineage_restricted_programs.py。

实现复核完成：从原始计数独立重算60条细胞评分，最大误差2.45e-15以内；32项AUC通过独立Mann–Whitney方法复核。研究用30基因清单另单独导出为EXPLORATORY_30_gene_macrophage_program_NOT_VALIDATED.csv。
