# GSE202642：全细胞与巨噬细胞亚群dotplot统一格式

采用与GSE151530、GSE149614相同的Scanpy DotPlot绘图样式：红色连续色阶，基因标签竖排，点大小表示检出比例，右侧同时显示点大小图例和横向颜色图例。分组标志面板显示顶部括号。所有群按编号顺序完整显示，命名版保留“编号＋标志基因＋暂定RNA亚群名称”。

## 图形下载

| 聚类版本 | 图形 | PDF | 群数 |
|---|---|---|---:|
| 全细胞C0–C20，编号版 | [PNG](allcell_C0_C20_marker_dotplot.png) | [PDF](allcell_C0_C20_marker_dotplot.pdf) | 21 |
| 全细胞C0–C20，名称版 | [PNG](allcell_C0_C20_named_dotplot.png) | [PDF](allcell_C0_C20_named_dotplot.pdf) | 21 |
| 最近重聚类巨噬细胞M0–M8，名称版 | [PNG](macrophage_M0_M8_named_dotplot.png) | [PDF](macrophage_M0_M8_named_dotplot.pdf) | 9 |
| 此前功能命名巨噬细胞C0–C14，名称版 | [PNG](macrophage_C0_C14_legacy_named_dotplot.png) | [PDF](macrophage_C0_C14_legacy_named_dotplot.pdf) | 15 |

[下载全部图、数值及代码](GSE202642_dotplot_style_unified.zip)。GitHub页面点击Download raw file。

两个巨噬细胞版本采用不同的细胞筛选和聚类，不能互换编号或合并当作同一分群。C0–C20属于全细胞分析，与巨噬细胞内部C0–C14也不对应。原有混合RNA、患者主导和待定标签全部保留。

## 数值定义与核验

保留原图全部标志基因：全细胞43个，最近巨噬重聚类24个，此前功能命名巨噬36个；全细胞图按生物学标志面板组织列。数据范围、原群成员、编号、名称和RNA表达统计不重新计算聚类。

颜色表示逐细胞自然log1p(CP10k)的群平均值，均值包含所有零值，再对**每个基因在各群间进行0–1最小–最大缩放**，与其他数据集的`standard_scale='var'`一致。所有群恒定的基因颜色为0。颜色用于同一基因的群间相对比较，不能比较不同基因的绝对表达量。原全细胞图仅除以基因最大值、此前15群图使用绝对均值色阶；此次统一了显示缩放，原始均值仍完整导出。

点大小来自原始RNA计数>0的细胞比例，固定0–100%范围，采用Scanpy默认大小映射。零值不排除，也不插补或平滑。CP10k归一化使用原完整转录组总计数，而非标志物子矩阵行和。

`*_plot_values.csv`包含每个点的细胞数、检出比例、未缩放均值及绘图颜色值。`independent_value_audit.json`记录逐点与原始计数矩阵/完整文库归一化的独立核对；原全细胞均值以float32计算，因此复核允许2e-5数值容差。`verification.json`检查所有群标签完整显示、实际点图数值正确、大小与颜色图例存在。图形已检查名称和基因标签无裁切。

## 复现

```bash
python -m pip install -r requirements.txt
python plot_dotplots.py
```

包内`inputs`含已核对的群×基因汇总和命名信息，足以离线重绘。Scanpy的布局载体只有每群一行；实际点大小与颜色从真实汇总表显式传入，不从该布局载体估计，布局载体不用于生物学分析。需要重新提取/核对时，在相邻原分析目录存在的情况下运行`prepare_inputs.py`和`audit_plot_values.py`。

源分析：[全细胞21群](../GSE202642_HCC_allcell_annotation/README.md)、[最近巨噬9群](../GSE202642_macrophage_reclustering/README.md)、[此前巨噬15群](../GSE202642_PGAM5_subpopulation_assessment/README.md)。补充版本集中存放在本目录，原分析下载包保留历史图形。
