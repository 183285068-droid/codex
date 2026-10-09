"""Render the stage figure from saved results without rerunning statistics."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt


def render_plot(results_dir):
    root = Path(results_dir)
    source = 'primary_Toil_clinical'
    stages = ['I', 'II', 'III-IV']
    patients = pd.read_csv(root/'primary_Toil_clinical_analysis_patients.csv')
    summary = pd.read_csv(root/'stage_descriptive_statistics.csv')
    summary = summary[summary.source.eq(source)].set_index('stage').loc[stages]
    tests = pd.read_csv(root/'stage_overall_and_secondary_tests.csv')
    kw = tests[tests.source.eq(source) & tests.test.eq('Kruskal_Wallis_all_stages')].iloc[0]
    pairs = pd.read_csv(root/'stage_pairwise_tests.csv')
    pairs = pairs[pairs.source.eq(source)].copy()
    assert len(pairs) == 3
    assert set(zip(pairs.stage_A, pairs.stage_B)) == {('I','II'), ('I','III-IV'), ('II','III-IV')}

    fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.6), constrained_layout=True)
    colors = ['#377da3', '#4d9875', '#b67c3d']
    positions = np.arange(1, 4)
    values = [patients.loc[patients.stage.eq(stage), 'exploratory_target_fraction'].to_numpy()*100 for stage in stages]
    boxes = axs[0].boxplot(values, positions=positions, widths=.5, showfliers=False,
                           patch_artist=True, medianprops={'color':'black', 'linewidth':1.4})
    rng = np.random.default_rng(502)
    for i, y in enumerate(values):
        boxes['boxes'][i].set_facecolor(colors[i])
        boxes['boxes'][i].set_alpha(.4)
        axs[0].scatter(np.full(len(y), i+1)+rng.uniform(-.17,.17,len(y)), y,
                       s=13, color=colors[i], alpha=.5, edgecolors='none')
    axs[0].set(xticks=positions,
               xticklabels=[f'{stage}\n(n={len(y)})' for stage, y in zip(stages, values)],
               xlabel='Pathological AJCC stage',
               ylabel='Exploratory relative mixture estimate (%)')
    axs[0].set_title(f'All estimates, including zeros\nKruskal-Wallis p={kw.asymptotic_p:.4g}; permutation p={kw.permutation_p:.4g}\nBrackets: permutation p and BH-adjusted q', fontsize=10)

    # Brackets apply to the coefficient distributions on the left. The right
    # panel describes nonzero patient proportions, a different endpoint.
    highest = max(y.max() for y in values)
    spacing = max(4., highest*.1)
    baseline = highest+spacing*.75
    height = spacing*.18
    annotations = []
    for row in pairs.itertuples():
        x1, x2 = stages.index(row.stage_A)+1, stages.index(row.stage_B)+1
        y = baseline + (spacing*1.9 if x2-x1 == 2 else 0)
        label = f'p={row.permutation_p:.5f}\nq(BH)={row.pairwise_BH_q_3:.5f}'
        axs[0].plot([x1,x1,x2,x2], [y,y+height,y+height,y], color='black', linewidth=1.05)
        axs[0].text((x1+x2)/2, y+height+spacing*.10, label,
                    ha='center', va='bottom', fontsize=9)
        annotations.append({'stage_A':row.stage_A, 'stage_B':row.stage_B,
                            'permutation_p':row.permutation_p,
                            'BH_q_3':row.pairwise_BH_q_3,
                            'figure_label':label,
                            'test':'Two-sided Mann-Whitney U, 20000 permutations; BH across three pairs'})
    axs[0].set_ylim(-highest*.05, baseline+spacing*3.3)

    proportion = summary.nonzero_percent.to_numpy()
    lower = summary.nonzero_Wilson95CI_lower_percent.to_numpy()
    upper = summary.nonzero_Wilson95CI_upper_percent.to_numpy()
    axs[1].bar(positions, proportion, color=colors, alpha=.75)
    axs[1].errorbar(positions, proportion, yerr=np.vstack([proportion-lower, upper-proportion]),
                    fmt='none', ecolor='black', capsize=4)
    for i, (_, row) in enumerate(summary.iterrows()):
        axs[1].text(i+1, min(101,row.nonzero_Wilson95CI_upper_percent+3),
                    f'{int(row.nonzero_estimates)}/{int(row.patients)}', ha='center', fontsize=9)
    axs[1].set(xticks=positions, xticklabels=stages, xlabel='Pathological AJCC stage',
               ylabel='Patients with nonzero estimate (%)', ylim=(0,105),
               title='Nonzero estimate proportion (Wilson 95% CI)')
    for ax in axs:
        ax.spines[['top','right']].set_visible(False)
        ax.grid(axis='y', alpha=.15)
    fig.suptitle('TCGA-LIHC: GSE151530-only exploratory estimates by stage', fontsize=12)
    fig.savefig(root/'stage_comparison.png', dpi=200)
    fig.savefig(root/'stage_comparison.pdf')
    plt.close(fig)
    pd.DataFrame(annotations).to_csv(root/'figure_annotation_values.csv', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-dir', type=Path, default=Path(__file__).resolve().parent)
    render_plot(parser.parse_args().results_dir)
