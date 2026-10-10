"""Show source-backed identity and naming-marker expression for named RNA groups."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parent


def main():
    names = pd.read_csv(ROOT / 'functional_subgroup_names.csv',
                        dtype={'primary_cluster': str})
    obs = pd.read_csv(ROOT / 'cell_annotations.csv.gz', index_col=0,
                      dtype={'primary_cluster': str})
    genes = pd.read_csv(ROOT / 'portable_feature_genes.csv').gene.to_numpy()
    raw = sparse.load_npz(ROOT / 'portable_raw_feature_counts.npz').tocsr()
    norm = raw.multiply((1e4 / obs.total_counts.to_numpy())[:, None]).tocsr()
    log = norm.copy()
    log.data = np.log1p(log.data)
    selected = ['C1QA', 'C1QB', 'C1QC', 'CD68', 'CSF1R', 'MSR1', 'PGAM5']
    for row in names.itertuples():
        for gene in [row.marker_gene_1, row.marker_gene_2]:
            if gene not in selected:
                selected.append(gene)
    points = []
    for y, row in enumerate(names.itertuples()):
        mask = obs.primary_cluster.eq(row.primary_cluster).to_numpy()
        assert mask.sum() == row.cells
        for x, gene in enumerate(selected):
            indices = np.flatnonzero(genes == gene)
            assert len(indices) == 1
            j = int(indices[0])
            points.append({'original_cluster_id': row.original_cluster_id,
                           'functional_name_EN': row.functional_name_EN,
                           'gene': gene, 'x': x, 'y': y,
                           'cells': int(mask.sum()),
                           'detection_fraction': float((raw[mask, j] > 0).mean()),
                           'mean_CP10k': float(norm[mask, j].mean()),
                           'mean_log1p_CP10k': float(log[mask, j].mean())})
    dots = pd.DataFrame(points)
    assert len(dots) == len(names) * len(selected)
    dots.to_csv(ROOT / 'functional_subgroup_marker_dotplot_values.csv', index=False)
    # Check the complete plot's naming-marker entries against the independent
    # 30-marker evidence export, not a second fitted model or test.
    evidence = pd.read_csv(ROOT / 'functional_subgroup_marker_evidence.csv')
    for row in evidence.itertuples():
        value = dots[dots.original_cluster_id.eq(row.original_cluster_id) & dots.gene.eq(row.gene)].iloc[0]
        assert np.isclose(value.mean_CP10k, row.cluster_mean_CP10k)
        assert np.isclose(value.detection_fraction, row.cluster_detection_fraction)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'pdf.fonttype': 42})
    fig, ax = plt.subplots(figsize=(19, 9))
    scatter = ax.scatter(dots.x, dots.y, s=110 * dots.detection_fraction,
                         c=dots.mean_log1p_CP10k, cmap='viridis', vmin=0,
                         vmax=dots.mean_log1p_CP10k.max(), edgecolors='none')
    ax.set(xticks=np.arange(len(selected)), xticklabels=selected,
           yticks=np.arange(len(names)),
           yticklabels=[n.replace(' — ', ' | ') for n in names.functional_name_EN],
           ylim=(len(names) - .4, -.6), xlim=(-.7, len(selected) - .3))
    ax.tick_params(axis='x', rotation=65, labelsize=9)
    ax.tick_params(axis='y', labelsize=9)
    ax.axvline(6.5, color='#AAB3BA', lw=1)
    for row in range(len(names)):
        ax.axhline(row, color='#ECEFF1', lw=.5, zorder=-1)
    cax = fig.add_axes([.87, .38, .012, .37])
    colorbar = fig.colorbar(scatter, cax=cax)
    colorbar.set_label('Mean log1p(CP10k), including zeros', fontsize=10)
    for fraction in [.1, .5, 1.]:
        ax.scatter([], [], s=110 * fraction, c='#777777', label=f'{fraction:.0%} detected')
    ax.legend(frameon=False, bbox_to_anchor=(1.08, 1), loc='upper left', fontsize=9)
    ax.spines[['top', 'right', 'left']].set_visible(False)
    fig.text(.03, .955, 'GSE202642: macrophage identity and functional-name marker evidence', fontsize=16)
    fig.text(.03, .918, 'Naming markers describe RNA patterns; mixed-RNA groups remain provisional.', fontsize=11)
    fig.text(.03, .035, 'Dot area: raw RNA detection fraction. Color: mean log1p(full-library CP10k), including all zeros.\n'
             'The first seven genes show macrophage identity and PGAM5; the remaining genes support the descriptive names.', fontsize=10)
    fig.subplots_adjust(left=.32, right=.85, bottom=.22, top=.87)
    for extension in ['png', 'pdf']:
        fig.savefig(ROOT / ('functional_subgroup_marker_dotplot.' + extension), dpi=240, facecolor='white')
    plt.close(fig)
    print(json.dumps({'named_clusters': len(names), 'genes': len(selected),
                      'plot_values': len(dots), 'naming_marker_evidence_agreement': True}))


if __name__ == '__main__':
    main()
