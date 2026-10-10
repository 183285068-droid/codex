"""Standalone figures from saved evidence; no feature or model selection."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

R = Path(__file__).resolve().parent
COHORTS = json.loads((R/'config.json').read_text())['primary_cohorts']
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                     'axes.spines.right': False, 'pdf.fonttype': 42})


def save(fig, name):
    fig.savefig(R/f'{name}.png', dpi=180, bbox_inches='tight', facecolor='white')
    fig.savefig(R/f'{name}.pdf', bbox_inches='tight', facecolor='white')
    plt.close(fig)


def main():
    candidates = pd.read_csv(R/'cross_cohort_recurrent_candidates.csv').set_index('gene')
    genes = ['TOP2A', 'CENPK', 'TYMS', 'MKI67', 'NUSAP1', 'DHFR']
    values = np.array([[candidates.loc[g, c+'_log2FC'] for c in COHORTS] for g in genes])
    fig, ax = plt.subplots(figsize=(9.5, 5.4))
    im = ax.imshow(values, aspect='auto', cmap='RdBu_r', norm=TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2.5))
    for i, g in enumerate(genes):
        for j, c in enumerate(COHORTS):
            star = '*' if candidates.loc[g, c+'_strict_up'] else ''
            ax.text(j, i, f'{values[i,j]:.2f}{star}', ha='center', va='center',
                    color='white' if abs(values[i,j]) >= 1.6 else 'black')
    ax.set_xticks(range(5), COHORTS, rotation=20, ha='right')
    ax.set_yticks(range(6), genes)
    ax.set_title('Independent recurrent candidates: exact linear-mean log2 fold change', pad=16)
    fig.colorbar(im, ax=ax, fraction=.035, pad=.03, label='log2FC: RNA detected / undetected macrophages')
    fig.subplots_adjust(bottom=.26)
    fig.text(.08, .035, '* FDR <0.05, log2FC >=1, target detection >=10%. Recurrence >=3/5 cohorts.\n'
             'MKI67, NUSAP1, DHFR pass the >=70% donor-direction filter; PGAM5 is a definition anchor.\n'
             'These are exploratory candidates, not a validated PGAM5-specific signature.', fontsize=9)
    save(fig, 'candidate_evidence')

    pred = pd.read_csv(R/'nested_outer_predictions.csv.gz')
    summary = pd.read_csv(R/'nested_validation_summary.csv').set_index(['dataset','unit'])
    units = ['equalized_cell_fraction', 'library_RNA_contribution']
    labels = ['Idealized equal-cell normalization', 'Observed library RNA contribution']
    fig, axes = plt.subplots(2, 5, figsize=(18, 8.7))
    for i, unit in enumerate(units):
        d = pred[pred.unit.eq(unit) & pred.scenario.eq('standard')]
        limit = max(d.truth.max(), d.predicted.max())*100*1.1
        for j, c in enumerate(COHORTS):
            ax = axes[i,j]
            sub = d[d.dataset.eq(c)]
            ax.scatter(sub.truth*100, sub.predicted*100, s=17, alpha=.48, color='#2369a1', edgecolors='none')
            ax.plot([0, limit], [0, limit], linestyle='--', color='#777777', lw=1)
            ax.set(xlim=(-.2, limit), ylim=(-.2, limit), title=c,
                   xlabel='Known target mixture fraction (%)')
            s = summary.loc[c, unit]
            ax.text(.04, .96, f'r = {s.standard_Pearson_r:.3f}\nMAE = {s.standard_MAE_pp:.2f} pp\n'
                    f'Zero-target P95 = {s.all_zero_target_p95_percent:.2f}%', transform=ax.transAxes,
                    va='top', fontsize=9, bbox={'facecolor':'white', 'alpha':.8, 'edgecolor':'none'})
            ax.grid(alpha=.13)
            if j == 0:
                ax.set_ylabel(labels[i]+'\nEstimated target coefficient x 100')
    fig.suptitle('Nested whole-cohort holdouts: target recovery fails validation', fontsize=16)
    fig.subplots_adjust(top=.90, bottom=.14, hspace=.50, wspace=.36)
    fig.text(.05, .035, 'Dots: standard mixtures from held-out donors, with repeats; models selected inside each training fold.\n'
             'Zero-target P95 includes all target-absent stress scenarios, not only the plotted standard mixtures.\n'
             'Observed RNA contribution and idealized equal-cell mixtures are separate units; neither calibrates TCGA cell abundance.', fontsize=10)
    save(fig, 'nested_recovery')

    fig, ax = plt.subplots(figsize=(10, 5.8))
    x = np.arange(5)
    for i, (unit, color) in enumerate(zip(units, ['#2369a1', '#d07a31'])):
        y = [summary.loc[c, unit].all_zero_target_p95_percent for c in COHORTS]
        bars = ax.bar(x+(i-.5)*.36, y, width=.36, label=labels[i], color=color)
        ax.bar_label(bars, fmt='%.2f', padding=3, fontsize=9)
    ax.axhline(.5, linestyle='--', color='#b82d33', lw=1.6, label='Required upper bound: 0.5%')
    ax.set_xticks(x, COHORTS)
    ax.set(ylabel='95th percentile of estimated coefficient x 100', ylim=(0, 23),
           title='False target signal when no PGAM5 RNA-detected macrophages are present')
    ax.legend(loc='upper left', frameon=False, fontsize=9)
    fig.subplots_adjust(bottom=.22)
    fig.text(.10, .035, 'All target-absent cases: standard zero-dose, macrophage-dose challenges, PGAM5-detected\n'
             'nonmacrophages, cycling nonmacrophages and pure competitors. Mixtures are resampled, not independent patients.\n'
             'Every cohort/unit exceeds the fixed 0.5% validation gate.', fontsize=9)
    save(fig, 'zero_target_false_signal')

    donor = pd.read_csv(R/'state_discrimination_by_donor.csv')
    ss = pd.read_csv(R/'state_discrimination_summary.csv').set_index('dataset')
    fig, ax = plt.subplots(figsize=(10, 5.5))
    for j, c in enumerate(COHORTS):
        d = donor[donor.dataset.eq(c) & donor.descriptive_support_n1_ge5]
        y = d.RNA_detection_AUC.dropna().to_numpy()
        if len(y):
            ax.scatter(j+np.linspace(-.08,.08,len(y)), y, color='#2369a1', s=34)
            ax.hlines(np.median(y), j-.18, j+.18, color='black', lw=2)
            ax.text(j, .88, f'{len(y)} donor labels\n{int(ss.loc[c,"selected_genes"])} genes', ha='center', fontsize=9)
        else:
            ax.text(j, .48, 'NA\n0 selected genes', ha='center', color='#777777')
    ax.axhline(.5, color='#777777', linestyle='--', label='Chance discrimination')
    ax.axhline(.7, color='#b82d33', linestyle=':', label='Descriptive AUC benchmark: 0.70')
    ax.set_xticks(range(5), COHORTS)
    ax.set(ylim=(0,1), ylabel='AUC: PGAM5 RNA detected versus undetected',
           title='Held-cohort independent candidate scores provide little discrimination')
    ax.legend(loc='lower right', frameon=False, fontsize=9)
    fig.subplots_adjust(bottom=.22)
    fig.text(.10, .035, 'Features selected using other cohorts only. PGAM5 is excluded from this independent diagnostic score.\n'
             'Donors shown require >=5 detected and >=20 undetected macrophages. All panels have <5 genes;\n'
             'these descriptive scores do not meet the minimum panel-support requirement.', fontsize=9)
    save(fig, 'heldout_state_discrimination')
    pd.DataFrame({'figure':['candidate_evidence','nested_recovery','zero_target_false_signal','heldout_state_discrimination'],
                  'source':['cross_cohort_recurrent_candidates.csv','nested_outer_predictions.csv.gz; nested_validation_summary.csv',
                            'nested_validation_summary.csv','state_discrimination_by_donor.csv; state_discrimination_summary.csv']}).to_csv(R/'figure_sources.csv',index=False)


if __name__ == '__main__':
    main()
