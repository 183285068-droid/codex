from pathlib import Path
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
import pandas as pd

R = Path(__file__).resolve().parent
d = pd.read_csv(R / 'inputs/cell_expression_coordinates.csv.gz', index_col=0)
defs = pd.read_csv(R / 'inputs/cluster_annotations.csv').set_index('cluster')
previous = pd.read_csv(R / 'inputs/previous_cluster_statistics.csv').set_index('cluster')
assert d.index.is_unique and len(d) == 5605
assert set(d.mac_cluster) == {f'M{i}' for i in range(9)}
assert (d.total_counts > 0).all() and (d.PGAM5_raw_counts >= 0).all()
d['PGAM5_CP10k'] = d.PGAM5_raw_counts / d.total_counts * 10000
d['PGAM5_log1p_CP10k'] = np.log1p(d.PGAM5_CP10k)
rows = []
for c in [f'M{i}' for i in range(9)]:
    v = d[d.mac_cluster.eq(c)]
    rows.append(dict(cluster=c, subtype_EN=defs.loc[c, 'subtype_EN'], subtype_CN=defs.loc[c, 'subtype_CN'], status=defs.loc[c, 'status'], cells=len(v), PGAM5_detected=int(v.PGAM5_raw_counts.gt(0).sum()), PGAM5_detection_percent=float(v.PGAM5_raw_counts.gt(0).mean()*100), PGAM5_mean_CP10k_all_cells=float(v.PGAM5_CP10k.mean()), PGAM5_mean_log1p_CP10k=float(v.PGAM5_log1p_CP10k.mean()), PGAM5_median_CP10k_all_cells=float(v.PGAM5_CP10k.median())))
s = pd.DataFrame(rows).set_index('cluster')
fields = ['cells', 'PGAM5_detected', 'PGAM5_detection_percent', 'PGAM5_mean_CP10k_all_cells', 'PGAM5_mean_log1p_CP10k', 'PGAM5_median_CP10k_all_cells']
np.testing.assert_allclose(s[fields], previous.loc[s.index, fields], rtol=1e-12, atol=1e-12)
assert s.cells.sum() == 5605 and s.PGAM5_detected.sum() == 136
s.to_csv(R / 'PGAM5_M0_M8_expression_summary.csv')
xy = d[['UMAP1', 'UMAP2']].to_numpy()
assert np.isfinite(xy).all()
order = np.random.default_rng(149614).permutation(len(d))
cjk = Path('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
font = FontProperties(fname=str(cjk)) if cjk.exists() else None
label_names = {
    'M0': ('FOS/JUN', '即时早期反应', 'immediate-early'),
    'M1': ('MKI67/TOP2A', '增殖型', 'cycling'),
    'M2': ('C1QC/HLA-DRA', '抗原呈递相关', 'antigen-presentation'),
    'M3': ('JUND/HSPA1A', '应激相关（样本主导）', 'stress, sample-dominated'),
    'M4': ('TIMD4/CD5L', '驻留样', 'resident-like'),
    'M5': ('SPP1/LGALS1', '重塑相关', 'remodeling-associated'),
    'M6': ('C1QA/CD3D', '巨噬/T混合RNA', 'Mac-T mixed RNA'),
    'M7': ('MT1G/MT2A', '金属响应相关', 'metal-response'),
    'M8': ('GRM4/AIM1L', '巨噬富集（待定）', 'macrophage-enriched, unresolved'),
}
fig, axs = plt.subplots(1, 2, figsize=(18, 8))
settings = [('PGAM5_mean_CP10k_all_cells', 'Mean PGAM5 CP10k (all zeros included)', 'Group mean CP10k', '.4f'), ('PGAM5_detection_percent', 'PGAM5 RNA detection by group', 'Detected cells (%)', '.2f')]
for ax, (field, title, barlabel, fmt) in zip(axs, settings):
    values = d.mac_cluster.map(s[field]).to_numpy()
    points = ax.scatter(*xy[order].T, c=values[order], cmap='viridis', s=5, vmin=0, vmax=s[field].max(), linewidths=0, rasterized=True)
    fig.colorbar(points, ax=ax, label=barlabel, shrink=.85)
    for c in s.index:
        center = np.median(xy[d.mac_cluster.eq(c)], axis=0)
        val = format(s.loc[c, field], fmt) + ('%' if field.endswith('percent') else '')
        gene, name_cn, name_en = label_names[c]
        label = c+' '+gene+'\n'+(name_cn if font else name_en)+'\n'+val
        ax.text(*center, label, fontproperties=font, ha='center', va='center', fontsize=9, linespacing=1.1, bbox=dict(facecolor='white', edgecolor='#666666', linewidth=.4, alpha=.85, boxstyle='round,pad=.25'))
    ax.set_title(title, fontsize=12)
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    ax.set_xticks([])
    ax.set_yticks([])
fig.suptitle('GSE149614 | 10 primary HCC tumor samples | M0-M8', fontsize=15)
fig.text(.5, .035, 'Colors show group summaries at the original cell coordinates. Labels give the value for each group.\nCP10k uses the full-transcriptome library total; RNA detection means PGAM5 raw count > 0.', ha='center', fontsize=10)
fig.tight_layout(rect=(0,.09,1,.95))
for ext in ['png', 'pdf']:
    fig.savefig(R / ('PGAM5_M0_M8_abundance_labeled_UMAP.'+ext), dpi=200)
plt.close(fig)
(R / 'verification.json').write_text(json.dumps({'cells': len(d), 'PGAM5_detected': int(s.PGAM5_detected.sum()), 'clusters_verified': 9, 'previous_statistics_match_rtol': 1e-12, 'all_zeros_retained': True, 'coordinates_not_recomputed': True, 'finite_coordinates': True, 'all_nine_group_names_labeled': True, 'label_language': 'Chinese' if font else 'English'}, indent=2))
table = '\n'.join(f'| {x.Index} | {x.subtype_CN} | {x.cells} | {x.PGAM5_mean_CP10k_all_cells:.4f} | {x.PGAM5_detection_percent:.2f}% ({x.PGAM5_detected}/{x.cells}) |' for x in s.itertuples())
(R / 'README.md').write_text('''# GSE149614：M0–M8的PGAM5表达丰度补充

沿用10个原发HCC肿瘤T样本的5,605个巨噬富集细胞、原始M0–M8标签和UMAP坐标。使用保存矩阵中的原始PGAM5计数及完整文库总计数重新计算统计，与此前9群的结果逐项一致。未重新聚类或改变UMAP坐标。

[标注表达数值的UMAP（PNG）](PGAM5_M0_M8_abundance_labeled_UMAP.png) · [PDF](PGAM5_M0_M8_abundance_labeled_UMAP.pdf) · [精确CSV](PGAM5_M0_M8_expression_summary.csv) · [下载补充结果包](GSE149614_PGAM5_M0_M8_expression_supplement.zip)

两幅UMAP均标注M0–M8编号、标志基因及功能程序名称。左图按每群平均PGAM5 CP10k着色并标注均值；右图按每群RNA检出率着色并标注百分比。图中每个点的颜色代表所在群的汇总值。平均表达保留所有零值，CP10k = PGAM5原始计数 / 完整转录组文库总计数 × 10,000。RNA检出定义为PGAM5原始计数>0。

| 群 | 暂定RNA名称 | 细胞数 | 平均PGAM5 CP10k（含零） | 检出率（检出/总数） |
|---|---|---:|---:|---:|
'''+table+'''

M1 MKI67/TOP2A增殖型的平均表达及检出率最高。M8均值对应仅1个检出细胞；M3此前有91.6%细胞来自HCC01T；M6为巨噬/T混合RNA群。群名称是探索性RNA注释。本轮补充描述性表达统计及可视化，未新增推断性显著性检验。

## 复现

ZIP包含用于重绘的逐细胞计数、文库总计数、UMAP坐标、群注释和代码。沿用[原分析](../GSE149614_primaryT_hierarchical_clustering/README.md)的环境（Python3.12、numpy2.3.5、pandas3.0.6、matplotlib3.10.8），解压后运行；或先用pip安装requirements.txt。有Noto Sans CJK字体时标注中文功能名，否则标注对应英文名称。

```bash
export MPLCONFIGDIR="$PWD/.matplotlib"
export XDG_CACHE_HOME="$PWD/.cache"
python plot_expression.py
```

`prepare_inputs.py`用于从相邻的原分析目录重新提取输入，需要anndata0.13.4；单独重绘无需运行该脚本。来源校验见[source_provenance.json](source_provenance.json)，统计核对见[verification.json](verification.json)。
''')
print(s[['cells', 'PGAM5_mean_CP10k_all_cells', 'PGAM5_detection_percent']].to_string())
