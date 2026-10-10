"""Render all-cell normalized and raw-count panels without dropping zeros."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def render(root):
    root = Path(root)
    d = pd.read_csv(root / 'macrophage_expression.csv.gz')
    s = pd.read_csv(root / 'correlation_statistics.csv').set_index('analysis')
    counts = json.loads((root / 'detection_summary.json').read_text())
    both = d.PGAM5_detected & d.MYO19_detected
    titles = ['Primary: log1p(CP10k)', 'Sensitivity: raw integer counts']
    suffixes = ['log1p_CP10k', 'raw_count']
    fig, axs = plt.subplots(1, 2, figsize=(12.6, 7.6))
    fig.subplots_adjust(left=.07, right=.98, bottom=.41, top=.82, wspace=.25)
    for i, (ax, suffix) in enumerate(zip(axs, suffixes)):
        x, y = d['PGAM5_' + suffix], d['MYO19_' + suffix]
        ax.scatter(x, y, s=19, alpha=.25, color='#538aa2', edgecolors='none',
                   label='All macrophages (zeros retained)')
        ax.scatter(x[both], y[both], s=35, alpha=.85, color='#a36837', edgecolors='none',
                   label=f'Both detected (n={int(both.sum())})')
        ax.set(xlabel=f'PGAM5 {titles[i].split(": ", 1)[1]}',
               ylabel=f'MYO19 {titles[i].split(": ", 1)[1]}', title=titles[i])
        statistic_suffix = 'raw_counts' if suffix == 'raw_count' else suffix
        sp = s.loc['Spearman_all_cells_' + statistic_suffix]
        pe = s.loc['Pearson_all_cells_' + statistic_suffix]
        label = (f'Spearman rho = {sp.correlation:.4f}; two-sided asymptotic p = {sp.asymptotic_p:.3g}\n')
        if i == 0:
            label += (f'Monte Carlo permutation p = {sp.permutation_p:.3g}\n'
                      f'({int(sp.extreme_permutations)} extreme / 200,000 shuffles; plus-one estimate at resolution limit)\n\n')
        else:
            label += f'Spearman secondary BH q = {sp.secondary_BH_q:.3g}\n\n'
        label += (f'Pearson r = {pe.correlation:.4f}; two-sided p = {pe.asymptotic_p:.3g}\n'
                  f'Pearson secondary BH q = {pe.secondary_BH_q:.3g}')
        fig.text(ax.get_position().x0, .325, label, ha='left', va='top', fontsize=9)
        ax.set_xlim(-max(float(x.max()) * .035, .04), float(x.max()) * 1.06)
        ax.set_ylim(-max(float(y.max()) * .035, .04), float(y.max()) * 1.06)
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(alpha=.12)
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(.52, .12),
               ncol=2, frameon=False, fontsize=9)
    fig.suptitle('GSE149614: PGAM5 and MYO19 RNA expression in HCC macrophages', fontsize=14, y=.97)
    fig.text(.5, .91, 'All 7,448 macrophages from 10 primary tumor T samples; all zeros retained',
             ha='center', fontsize=11)
    fig.text(.07, .04,
             f'Detected: PGAM5 {counts["PGAM5_detected"]}; MYO19 {counts["MYO19_detected"]}; '
             f'both {counts["both_detected"]}; neither {counts["neither_detected"]:,}. Identical coordinates overlap.\n'
             'All tests use 7,448 cells. Cell-level p-values do not establish an independent patient-level association.',
             fontsize=9, ha='left')
    fig.savefig(root / 'PGAM5_MYO19_all_macrophages.png', dpi=220)
    fig.savefig(root / 'PGAM5_MYO19_all_macrophages.pdf')
    plt.close(fig)


if __name__ == '__main__':
    render(Path(__file__).resolve().parent)
