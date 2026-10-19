"""Reproduce sample-aware PGAM5 comparisons and static figures from saved counts."""
from pathlib import Path
from itertools import combinations, product
import hashlib
import json
import platform
import importlib.metadata as metadata

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'outputs'
OUT.mkdir(exist_ok=True)
GROUPS = [f'M{i}' for i in range(9)]
NAMES = ['FOS/JUN', 'MKI67/TOP2A', 'C1QC/HLA-DRA', 'JUND/HSPA1A',
         'TIMD4/CD5L', 'SPP1/LGALS1', 'C1QA/CD3D', 'MT1G/MT2A', 'GRM4/AIM1L']
MIN_CELLS = 10
MIN_PAIRS = 5

def stars(q):
    if not np.isfinite(q): return 'NT'
    return '****' if q < .0001 else '***' if q < .001 else '**' if q < .01 else '*' if q < .05 else 'ns'

def exact_signflip(x):
    """Two-sided paired label-swap test of the mean difference; enumerate all 2^n swaps."""
    x = np.asarray(x, dtype=float)
    observed = abs(x.mean())
    null = np.asarray(list(product([-1., 1.], repeat=len(x)))) @ x / len(x)
    return float(np.mean(np.abs(null) >= observed - 1e-12))

def adjust(table, p='p_value'):
    # All 36 planned comparisons stay in the family. Untestable p values enter
    # the adjustment as 1, but are displayed as NT with q left unavailable.
    values = table[p].fillna(1).to_numpy()
    table['q_BH'] = stats.false_discovery_control(values, method='bh')
    table.loc[table[p].isna(), 'q_BH'] = np.nan
    table['significance'] = table.q_BH.map(stars)
    return table

def paired_table(aggregated, endpoint, minimum_cells):
    usable = aggregated[aggregated.cells >= minimum_cells]
    wide = usable.pivot(index='Sample', columns='cluster', values=endpoint).reindex(columns=GROUPS)
    rows = []
    for a, b in combinations(GROUPS, 2):
        pair = wide[[a, b]].dropna()
        delta = pair[a].to_numpy() - pair[b].to_numpy()
        n = len(pair)
        rows.append(dict(group_a=a, group_b=b, endpoint=endpoint,
                         min_cells_per_sample_cluster=minimum_cells,
                         paired_samples=n, mean_difference_a_minus_b=float(delta.mean()) if n else np.nan,
                         median_difference_a_minus_b=float(np.median(delta)) if n else np.nan,
                         samples=';'.join(pair.index),
                         p_value=exact_signflip(delta) if n >= MIN_PAIRS else np.nan,
                         status='tested' if n >= MIN_PAIRS else 'insufficient_paired_samples'))
    return adjust(pd.DataFrame(rows))

def save(fig, stem):
    for ext in ['png', 'svg']:
        fig.savefig(OUT / f'{stem}.{ext}', dpi=210, facecolor='white')
    plt.close(fig)

def significance_matrix(ax, table, sample_level):
    ax.set_xlim(-.5, 8.5); ax.set_ylim(8.5, -.5); ax.set_aspect('equal')
    ax.set_xticks(range(9), GROUPS); ax.set_yticks(range(9), GROUPS)
    ax.tick_params(length=0, labelsize=10)
    norm = Normalize(0, 4)
    cmap = plt.get_cmap('Blues')
    for i, a in enumerate(GROUPS):
        for j, b in enumerate(GROUPS):
            if i <= j:
                ax.add_patch(plt.Rectangle((j-.5, i-.5), 1, 1, facecolor='#f5f5f5', edgecolor='white'))
                continue
            r = table[(table.group_a == b) & (table.group_b == a)].iloc[0]
            q = r.q_BH
            color = '#e5e7eb' if not np.isfinite(q) else cmap(norm(min(4, -np.log10(max(q, 1e-300)))))
            ax.add_patch(plt.Rectangle((j-.5, i-.5), 1, 1, facecolor=color, edgecolor='white', lw=1))
            qtext = f'{q:.2g}' if np.isfinite(q) else 'NT'
            label = f'{qtext}\n{stars(q)}' if np.isfinite(q) else 'NT'
            if sample_level: label += f' | n={r.paired_samples}' if np.isfinite(q) else f'\nn={r.paired_samples}'
            ax.text(j, i, label, ha='center', va='center', fontsize=7.5 if sample_level else 9,
                    color='white' if np.isfinite(q) and q < .01 else '#17212b', linespacing=1.5)
    for spine in ax.spines.values(): spine.set_visible(False)
    cb = ax.figure.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, fraction=.045, pad=.04, shrink=.75)
    cb.set_label('-log10(BH-adjusted P)', fontsize=9)
    ax.set_title('All 36 pairwise comparisons\nBH-adjusted P and significance', fontsize=12, pad=12)

def main():
    source = ROOT / 'inputs/cell_expression_coordinates.csv.gz'
    d = pd.read_csv(source)
    expected_sha = '068e7ce434dc5b9cbec3fd9282eacd1c58747f2c3e0698ed92c26eeccba95a2c'
    assert hashlib.sha256(source.read_bytes()).hexdigest() == expected_sha
    assert len(d) == 5605 and d.cell_id.is_unique and set(d.mac_cluster) == set(GROUPS)
    assert not d[['cell_id', 'Sample', 'mac_cluster', 'PGAM5_raw_counts', 'total_counts']].isna().any().any()
    assert (d.total_counts > 0).all() and (d.PGAM5_raw_counts >= 0).all()
    assert (d.PGAM5_raw_counts <= d.total_counts).all()
    assert d.Sample.nunique() == 10
    d['CP10k'] = d.PGAM5_raw_counts / d.total_counts * 10000
    d['log1p_CP10k'] = np.log1p(d.CP10k)
    d['detected'] = d.PGAM5_raw_counts > 0
    assert int(d.detected.sum()) == 136
    agg = d.groupby(['Sample', 'mac_cluster']).agg(
        cells=('cell_id', 'size'), detected=('detected', 'sum'),
        mean_CP10k=('CP10k', 'mean'), mean_log1p_CP10k=('log1p_CP10k', 'mean'),
        PGAM5_counts=('PGAM5_raw_counts', 'sum'), library_counts=('total_counts', 'sum')).reset_index()
    agg = agg.rename(columns={'mac_cluster': 'cluster'})
    agg['detection_percent'] = agg.detected / agg.cells * 100
    agg['pseudobulk_CP10k'] = agg.PGAM5_counts / agg.library_counts * 10000
    agg['eligible_primary'] = agg.cells >= MIN_CELLS
    agg.to_csv(OUT / 'sample_cluster_statistics.csv', index=False)
    agg.pivot(index='Sample', columns='cluster', values='cells').reindex(columns=GROUPS).fillna(0).astype(int).to_csv(OUT/'sample_cluster_cell_counts.csv')
    summary = d.groupby('mac_cluster').agg(cells=('cell_id', 'size'), detected=('detected', 'sum'),
        samples=('Sample', 'nunique'), mean_CP10k=('CP10k', 'mean'),
        mean_log1p_CP10k=('log1p_CP10k', 'mean'), median_CP10k=('CP10k', 'median')).reindex(GROUPS)
    summary.index.name = 'cluster'
    summary['detection_percent'] = summary.detected / summary.cells * 100
    summary['marker_label'] = NAMES
    summary['eligible_samples_primary'] = agg[agg.eligible_primary].groupby('cluster').size().reindex(GROUPS, fill_value=0)
    summary.to_csv(OUT / 'cluster_expression_summary.csv')
    # Reconcile all cluster means/counts with the previously published statistics.
    prior = pd.read_csv(ROOT/'inputs/previous_expression_summary.csv').set_index('cluster').loc[GROUPS]
    np.testing.assert_array_equal(summary.cells, prior.cells)
    np.testing.assert_array_equal(summary.detected, prior.PGAM5_detected)
    np.testing.assert_allclose(summary.mean_CP10k, prior.PGAM5_mean_CP10k_all_cells, rtol=1e-12)
    np.testing.assert_allclose(summary.mean_log1p_CP10k, prior.PGAM5_mean_log1p_CP10k, rtol=1e-12)
    primary = paired_table(agg, 'mean_log1p_CP10k', MIN_CELLS)
    primary.to_csv(OUT/'primary_sample_paired_comparisons.csv', index=False)
    secondary = paired_table(agg, 'detection_percent', MIN_CELLS)
    secondary.to_csv(OUT/'secondary_sample_detection_comparisons.csv', index=False)
    sensitivity = pd.concat([paired_table(agg, endpoint, threshold)
        for endpoint in ['mean_log1p_CP10k', 'mean_CP10k', 'pseudobulk_CP10k']
        for threshold in [1, 10, 20]], ignore_index=True)
    sensitivity.to_csv(OUT/'sensitivity_sample_comparisons.csv', index=False)
    cell_rows, detection_rows = [], []
    for a, b in combinations(GROUPS, 2):
        x, y = [d[d.mac_cluster.eq(c)] for c in [a, b]]
        u = stats.mannwhitneyu(x.log1p_CP10k, y.log1p_CP10k, alternative='two-sided', method='asymptotic', use_continuity=True)
        cell_rows.append(dict(group_a=a, group_b=b, cells_a=len(x), cells_b=len(y), U=float(u.statistic),
            rank_biserial_a_minus_b=float(2*u.statistic/(len(x)*len(y))-1),
            mean_difference_log1p_CP10k=float(x.log1p_CP10k.mean()-y.log1p_CP10k.mean()), p_value=float(u.pvalue)))
        fisher = stats.fisher_exact([[x.detected.sum(), (~x.detected).sum()], [y.detected.sum(), (~y.detected).sum()]])
        detection_rows.append(dict(group_a=a, group_b=b, odds_ratio=float(fisher.statistic), p_value=float(fisher.pvalue)))
    cell = adjust(pd.DataFrame(cell_rows)); cell.to_csv(OUT/'exploratory_cell_MannWhitney_comparisons.csv', index=False)
    adjust(pd.DataFrame(detection_rows)).to_csv(OUT/'exploratory_cell_detection_Fisher_comparisons.csv', index=False)
    kw = stats.kruskal(*[d[d.mac_cluster.eq(c)].log1p_CP10k for c in GROUPS])
    # Cross-check the exact permutation implementation against SciPy, on real
    # matched observations and a known all-positive toy case.
    wide = agg[agg.cells >= MIN_CELLS].pivot(index='Sample', columns='cluster', values='mean_log1p_CP10k')
    pair = wide[['M0', 'M1']].dropna()
    perm = stats.permutation_test((pair.M0.to_numpy(), pair.M1.to_numpy()), lambda a,b: np.mean(a-b),
        permutation_type='samples', n_resamples=np.inf, alternative='two-sided')
    np.testing.assert_allclose(exact_signflip((pair.M0-pair.M1).to_numpy()), perm.pvalue)
    assert np.isclose(exact_signflip(np.ones(5)), 2/32)
    assert exact_signflip(np.zeros(5)) == 1
    np.testing.assert_allclose(stats.false_discovery_control([.01,.04,.03]), [.03,.04,.04])

    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':10, 'axes.spines.top':False,
        'axes.spines.right':False, 'axes.titleweight':'bold', 'svg.fonttype':'none'})
    fig, axs = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={'width_ratios':[1.15,1]})
    rng = np.random.default_rng(149614)
    for i, c in enumerate(GROUPS):
        v = agg[agg.cluster.eq(c)]
        good = v[v.eligible_primary]; low = v[~v.eligible_primary]
        color = '#c28732' if c == 'M1' else '#386b8d'
        if len(good):
            bp = axs[0].boxplot([good.mean_log1p_CP10k], positions=[i], widths=.52, patch_artist=True,
                showfliers=False, manage_ticks=False, medianprops={'color':'#17212b','lw':1.8})
            bp['boxes'][0].set(facecolor=color, alpha=.25, edgecolor=color)
            axs[0].scatter(i+rng.uniform(-.15,.15,len(good)),good.mean_log1p_CP10k,s=32,c=color,zorder=3,edgecolor='white',lw=.4)
        if len(low): axs[0].scatter(i+rng.uniform(-.15,.15,len(low)),low.mean_log1p_CP10k,s=30,facecolor='none',edgecolor='#a0a0a0',lw=.9,zorder=2)
    axs[0].set_xticks(range(9), [f'{c}\n{summary.loc[c,"eligible_samples_primary"]}/10' for c in GROUPS])
    axs[0].set_xlabel('Cluster | eligible samples / total samples')
    axs[0].set_ylabel('Sample mean PGAM5 log1p(CP10k), zeros included')
    axs[0].set_ylim(bottom=0); axs[0].grid(axis='y', color='#e8e8e8'); axs[0].set_axisbelow(True)
    axs[0].set_title('A  Sample-level expression\nEach point = one sample-cluster', loc='left', fontsize=12)
    significance_matrix(axs[1], primary, True)
    fig.suptitle('PGAM5 across macrophage-enriched clusters | sample-aware primary analysis', fontsize=16, y=.98)
    fig.text(.07,.06,'Filled points: >=10 cells per sample-cluster; hollow grey points: <10 cells (descriptive only).\n'
        'Paired two-sided exact label-swap test; >=5 matched samples required. BH family: all 36 pairs.\n'
        '* q<0.05; ** q<0.01; *** q<0.001; **** q<0.0001; ns: q>=0.05; NT: insufficient matched samples.',fontsize=10)
    fig.tight_layout(rect=[.025,.17,.99,.94]); save(fig,'01_sample_aware_expression_significance')

    fig, axs = plt.subplots(1,2,figsize=(16,8),gridspec_kw={'width_ratios':[1.15,1]})
    colors = ['#c28732' if c=='M1' else '#386b8d' for c in GROUPS]
    axs[0].bar(range(9), summary.mean_CP10k,color=colors,width=.65)
    for i,c in enumerate(GROUPS):
        axs[0].text(i,summary.loc[c,'mean_CP10k']+.002,f'{summary.loc[c,"mean_CP10k"]:.4f}',ha='center',fontsize=9)
    axs[0].set_xticks(range(9),[f'{c}\nn={summary.loc[c,"cells"]}' for c in GROUPS],fontsize=9)
    axs[0].set_ylim(0,summary.mean_CP10k.max()*1.22);axs[0].set_ylabel('Pooled cell mean PGAM5 CP10k, zeros included')
    axs[0].set_xlabel('Cluster | cell count');axs[0].grid(axis='y',color='#e8e8e8');axs[0].set_axisbelow(True)
    axs[0].set_title('A  Pooled descriptive abundance\nBars are cell-weighted; no biological-replicate error bars',loc='left',fontsize=12)
    significance_matrix(axs[1],cell,False)
    fig.suptitle('PGAM5 | exploratory cell-level comparisons',fontsize=16,y=.98)
    fig.text(.07,.06,f'Two-sided Mann-Whitney U tests with tie/continuity correction; BH over 36 pairs. Kruskal-Wallis P={kw.pvalue:.2g}.\n'
        'Heatmap tests cell-level distributions of log1p(CP10k); the bar panel shows arithmetic mean CP10k.\n'
        'Cells from the same sample are correlated. These stars are exploratory and do not establish sample-level significance.',fontsize=10)
    fig.tight_layout(rect=[.025,.17,.99,.94]);save(fig,'02_exploratory_cell_expression_significance')

    tested = primary[primary.status.eq('tested')]
    verification = {'input_sha256':expected_sha,'cells':len(d),'samples':d.Sample.nunique(),
        'PGAM5_detected_cells':int(d.detected.sum()),'primary_planned_pairs':36,
        'primary_testable_pairs':len(tested),'primary_significant_pairs_q_lt_0_05':int((tested.q_BH<.05).sum()),
        'primary_smallest_q':float(tested.q_BH.min()),'exploratory_cell_significant_pairs':int((cell.q_BH<.05).sum()),
        'exploratory_Kruskal_Wallis_H':float(kw.statistic),'exploratory_Kruskal_Wallis_P':float(kw.pvalue),
        'previous_statistics_reconciled':True,'exact_test_checked_against_scipy':True,
        'missing_sample_clusters_not_filled_as_zero':True,'all_cell_zeros_retained':True,
        'python':platform.python_version(),'packages':{p:metadata.version(p) for p in ['numpy','pandas','scipy','matplotlib']}}
    (OUT/'verification.json').write_text(json.dumps(verification,indent=2)+'\n')
    sig = cell[cell.q_BH < .05]
    sig_lines = '\n'.join(f'| {r.group_a} vs {r.group_b} | {r.p_value:.4g} | {r.q_BH:.4g} | {r.significance} |' for r in sig.itertuples())
    main_text = '未检出 BH 校正后 q<0.05 的亚群间差异。' if not (tested.q_BH<.05).any() else '存在 BH 校正后 q<0.05 的亚群间差异，详见完整统计表。'
    readme = f'''# GSE149614：PGAM5 亚群表达丰度统计比较

## 结果与图

**主要结论：{main_text}** 可检验的配对共 {len(tested)}/36 对，最小 q={tested.q_BH.min():.4g}。这不等同于证明各亚群表达相同；仅有10个样本，部分群覆盖不足，检验效能有限。

细胞合并描述中，M1（MKI67/TOP2A 增殖型）平均 CP10k=0.1070、检出率9.13%，均为9群最高。
细胞层面探索性检验有 {len(sig)}/36 对 q<0.05。**这些细胞层面的星号不能替代样本层面的证据。**

- [主图：样本表达分布和全部两两比较显著性](outputs/01_sample_aware_expression_significance.png) · [SVG](outputs/01_sample_aware_expression_significance.svg)
- [补充图：合并表达丰度和细胞层面显著性](outputs/02_exploratory_cell_expression_significance.png) · [SVG](outputs/02_exploratory_cell_expression_significance.svg)
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
7. 细胞层面为补充探索：双侧Mann–Whitney U（渐近法，校正并列值及连续性）比较表达分布，36对BH校正；列出秩二列相关效应量。另保存Fisher精确检验的检出率比较。合并Kruskal–Wallis H={kw.statistic:.4g}，P={kw.pvalue:.4g}，亦不调整样本内相关性。

## 细胞层面探索性显著比较

| 比较 | 原始P | BH校正q | 标注 |
|---|---:|---:|---|
{sig_lines}

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
'''
    (ROOT/'README.md').write_text(readme)
    print(json.dumps(verification,indent=2))
    print(sig[['group_a','group_b','q_BH','significance']].to_string(index=False))

if __name__ == '__main__':
    main()
