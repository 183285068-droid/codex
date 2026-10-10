# GSE149614：M0–M8的PGAM5表达丰度补充

沿用10个原发HCC肿瘤T样本的5,605个巨噬富集细胞、原始M0–M8标签和UMAP坐标。使用保存矩阵中的原始PGAM5计数及完整文库总计数重新计算统计，与此前9群的结果逐项一致。未重新聚类或改变UMAP坐标。

[标注表达数值的UMAP（PNG）](PGAM5_M0_M8_abundance_labeled_UMAP.png) · [PDF](PGAM5_M0_M8_abundance_labeled_UMAP.pdf) · [精确CSV](PGAM5_M0_M8_expression_summary.csv) · [下载补充结果包](GSE149614_PGAM5_M0_M8_expression_supplement.zip)

## 按GSE202642模板补充

新增[群平均表达UMAP（模板样式）](PGAM5_cluster_mean_expression_UMAP_template_style.png) · [PDF](PGAM5_cluster_mean_expression_UMAP_template_style.pdf)，以及[逐细胞表达UMAP（同样式）](PGAM5_expression_UMAP_template_style.png) · [PDF](PGAM5_expression_UMAP_template_style.pdf)。采用黄色至深红色渐变、两侧的编号/标志基因/功能名称引线和底部水平色标，与用户提供的GSE202642图一致。

模板群均值图使用**每个细胞log1p(CP10k)的算术平均值（保留零值）**，与先对群平均CP10k取log1p是不同的汇总方法。M1群平均log1p(CP10k)=0.0686，为本轮最高；原始平均CP10k=0.1070仍保留在统计表中。色阶覆盖本数据的完整观察范围，两个数据集各自设定色阶，不能仅凭颜色深浅比较跨数据集的绝对表达量。图中每个点的群均值着色表示所在群的汇总值；逐细胞图则按该细胞表达值着色，未检出为灰色。

运行`python plot_template_style.py`可重绘这两张模板样式图。坐标及9群名称沿用此前GSE149614结果；核查记录见[template_style_verification.json](template_style_verification.json)。

## 群表达数值

两幅UMAP均标注M0–M8编号、标志基因及功能程序名称。左图按每群平均PGAM5 CP10k着色并标注均值；右图按每群RNA检出率着色并标注百分比。图中每个点的颜色代表所在群的汇总值。平均表达保留所有零值，CP10k = PGAM5原始计数 / 完整转录组文库总计数 × 10,000。RNA检出定义为PGAM5原始计数>0。

| 群 | 暂定RNA名称 | 细胞数 | 平均PGAM5 CP10k（含零） | 检出率（检出/总数） |
|---|---|---:|---:|---:|
| M0 | FOS/JUN即时早期反应巨噬细胞 | 1428 | 0.0297 | 1.75% (25/1428) |
| M1 | MKI67/TOP2A增殖型巨噬细胞 | 208 | 0.1070 | 9.13% (19/208) |
| M2 | C1QC/HLA-DRA抗原呈递相关巨噬细胞 | 1287 | 0.0311 | 2.33% (30/1287) |
| M3 | JUND/HSPA1A应激相关巨噬细胞（样本主导） | 523 | 0.0546 | 3.25% (17/523) |
| M4 | TIMD4/CD5L驻留样巨噬细胞 | 356 | 0.0291 | 1.97% (7/356) |
| M5 | SPP1/LGALS1重塑相关巨噬细胞 | 1477 | 0.0370 | 2.03% (30/1477) |
| M6 | C1QA/CD3D巨噬/T混合RNA群 | 149 | 0.0435 | 3.36% (5/149) |
| M7 | MT1G/MT2A金属响应相关巨噬细胞 | 135 | 0.0103 | 1.48% (2/135) |
| M8 | GRM4/AIM1L巨噬富集群（待定） | 42 | 0.0645 | 2.38% (1/42) |

M1 MKI67/TOP2A增殖型的平均表达及检出率最高。M8均值对应仅1个检出细胞；M3此前有91.6%细胞来自HCC01T；M6为巨噬/T混合RNA群。群名称是探索性RNA注释。本轮补充描述性表达统计及可视化，未新增推断性显著性检验。

## 复现

ZIP包含用于重绘的逐细胞计数、文库总计数、UMAP坐标、群注释和代码。沿用[原分析](../GSE149614_primaryT_hierarchical_clustering/README.md)的环境（Python3.12、numpy2.3.5、pandas3.0.6、matplotlib3.10.8），解压后运行；或先用pip安装requirements.txt。有Noto Sans CJK字体时标注中文功能名，否则标注对应英文名称。

```bash
export MPLCONFIGDIR="$PWD/.matplotlib"
export XDG_CACHE_HOME="$PWD/.cache"
python plot_expression.py
```

`prepare_inputs.py`用于从相邻的原分析目录重新提取输入，需要anndata0.13.4；单独重绘无需运行该脚本。来源校验见[source_provenance.json](source_provenance.json)，统计核对见[verification.json](verification.json)。
