# HCC 全部样本合并：收紧差异阈值

沿用已完成的细胞层面 Wilcoxon 检验（tie correction、BH FDR），不做患者配对或测序深度匹配。分组、QC 和常规归一化不变：PGAM5 RNA 检出 180 个细胞、未检出 4,314 个细胞。

新阈值为 FDR <0.05 且 |近似 log2FC| >=1；保留此前至少一组检出率 >=10% 的表达过滤。只调整筛选阈值不会改变原始检验的 p 值、FDR 或 log2FC，所以直接使用保留的完整检验结果重新筛选，无需重复运行相同统计检验。

PGAM5 本身保留在 DE_including_PGAM5 表中，因用于定义分组，从上调/下调候选表排除。上调表中保留线粒体和核糖体等类别，不增加额外功能过滤。候选列表是探索性结果，未经独立验证。

详细数量与参数见 summary.json。原始完整检验表为 HCC_all_samples_PGAM5_full_DE.csv，执行脚本为 unadjusted_cell_analysis.py。
