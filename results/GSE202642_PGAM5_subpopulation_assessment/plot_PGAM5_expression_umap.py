"""Overlay PGAM5 RNA intensity on the saved, PGAM5-free macrophage UMAP.

No embedding, clustering, imputation, smoothing or expression clipping is run.
Normalization uses each cell's full-library total, not the partial matrix sum.
"""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parent
ZERO_COLOR = '#C9CDD2'


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def label_clusters(ax, frame):
    for cluster, cells in frame.groupby('primary_cluster', sort=False):
        x, y = cells[['UMAP1', 'UMAP2']].median()
        ax.text(x, y, 'C' + cluster, ha='center', va='center', fontsize=9,
                color='#18232C', weight='bold',
                path_effects=[pe.withStroke(linewidth=3, foreground='white')])


def finish_map(fig, ax, scatter, label, caption, name, legend=False):
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_aspect('equal', adjustable='box')
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(scatter, ax=ax, fraction=.035, pad=.025, shrink=.80)
    cbar.set_label(label, fontsize=11)
    cbar.ax.tick_params(labelsize=9)
    if legend:
        ax.legend(handles=[Line2D([0], [0], linestyle='', marker='o',
                                 color=ZERO_COLOR, markersize=6,
                                 label='PGAM5 not detected (raw count = 0)')],
                  loc='lower left', bbox_to_anchor=(0, -.07), frameon=False,
                  fontsize=9)
    fig.text(.09, .055, caption, ha='left', va='bottom', fontsize=9,
             color='#47525C', linespacing=1.6)
    fig.subplots_adjust(left=.07, right=.90, bottom=.16, top=.89)
    fig.savefig(ROOT / (name + '.png'), dpi=240, facecolor='white')
    fig.savefig(ROOT / (name + '.pdf'), facecolor='white')
    plt.close(fig)


def main():
    coord_path = ROOT / 'UMAP_coordinates.csv.gz'
    coord_hash = sha256(coord_path)
    obs = pd.read_csv(ROOT / 'cell_annotations.csv.gz', index_col=0,
                      dtype={'primary_cluster': str})
    coords = pd.read_csv(coord_path)
    assert obs.cell_id.is_unique and coords.cell_id.is_unique
    assert set(obs.cell_id) == set(coords.cell_id)
    assert obs.total_counts.gt(0).all()
    genes = pd.read_csv(ROOT / 'portable_feature_genes.csv').gene
    gene_index = np.flatnonzero(genes.eq('PGAM5').to_numpy())
    assert len(gene_index) == 1
    counts = sparse.load_npz(ROOT / 'portable_raw_feature_counts.npz')
    assert counts.shape == (len(obs), len(genes))
    pg_raw = counts[:, gene_index[0]].toarray().ravel().astype(np.float64)
    assert np.array_equal(pg_raw, obs.PGAM5_counts.to_numpy())
    cp10k = pg_raw / obs.total_counts.to_numpy(np.float64) * 1e4
    norm_error = float(np.max(np.abs(cp10k - obs.PGAM5_CP10k.to_numpy())))
    assert np.allclose(cp10k, obs.PGAM5_CP10k, rtol=1e-6, atol=1e-7)
    obs['PGAM5_CP10k_recomputed'] = cp10k
    obs['PGAM5_log1p_CP10k'] = np.log1p(cp10k)
    frame = obs[['cell_id', 'patient_id', 'sample_name', 'primary_cluster',
                 'PGAM5_counts', 'total_counts', 'PGAM5_CP10k_recomputed',
                 'PGAM5_log1p_CP10k']].merge(coords, on='cell_id', how='left',
                                           validate='one_to_one', sort=False)
    assert frame.cell_id.tolist() == obs.cell_id.tolist()
    assert np.isfinite(frame[['UMAP1', 'UMAP2', 'PGAM5_log1p_CP10k']]).all().all()
    frame['PGAM5_detected'] = frame.PGAM5_counts.gt(0)
    assert np.array_equal(frame.PGAM5_detected, obs.PGAM5_detected)
    grouped = frame.groupby('primary_cluster').agg(
        cells=('cell_id', 'size'),
        PGAM5_detected_cells=('PGAM5_detected', 'sum'),
        PGAM5_detection_fraction=('PGAM5_detected', 'mean'),
        mean_PGAM5_CP10k=('PGAM5_CP10k_recomputed', 'mean'),
        mean_PGAM5_log1p_CP10k=('PGAM5_log1p_CP10k', 'mean'),
        median_PGAM5_CP10k=('PGAM5_CP10k_recomputed', 'median'))
    grouped = grouped.loc[sorted(grouped.index, key=int)]
    assert grouped.cells.sum() == len(frame)
    assert grouped.PGAM5_detected_cells.sum() == frame.PGAM5_detected.sum()
    weighted_mean = np.average(grouped.mean_PGAM5_log1p_CP10k,
                               weights=grouped.cells)
    assert np.isclose(weighted_mean, frame.PGAM5_log1p_CP10k.mean())
    original = pd.read_csv(ROOT / 'cluster_PGAM5_association.csv',
                           dtype={'cluster': str})
    original = original[original.geometry.eq('primary')].set_index('cluster')
    assert np.allclose(grouped.mean_PGAM5_CP10k,
                       original.loc[grouped.index].cluster_PGAM5_mean_CP10k,
                       rtol=1e-6, atol=1e-8)
    frame['cluster_mean_PGAM5_log1p_CP10k'] = frame.primary_cluster.map(
        grouped.mean_PGAM5_log1p_CP10k)
    frame.to_csv(ROOT / 'PGAM5_UMAP_expression_values.csv.gz', index=False)
    grouped.rename_axis('cluster').reset_index().to_csv(
        ROOT / 'PGAM5_expression_by_cluster.csv', index=False)

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    # Trim only the colormap's near-white portion; expression values use the
    # complete observed range, without quantile clipping or saturation.
    cmap = LinearSegmentedColormap.from_list(
        'PGAM5_YlOrRd', plt.colormaps['YlOrRd'](np.linspace(.20, 1, 256)))
    fig, ax = plt.subplots(figsize=(10, 8))
    zero = frame[~frame.PGAM5_detected]
    positive = frame[frame.PGAM5_detected].sort_values('PGAM5_log1p_CP10k')
    ax.scatter(zero.UMAP1, zero.UMAP2, s=3, color=ZERO_COLOR,
               edgecolors='none', rasterized=True)
    vmax = float(frame.PGAM5_log1p_CP10k.max())
    sc = ax.scatter(positive.UMAP1, positive.UMAP2, s=13,
                    c=positive.PGAM5_log1p_CP10k, cmap=cmap,
                    norm=Normalize(0, vmax), edgecolors='none', rasterized=True)
    label_clusters(ax, frame)
    ax.set_title('GSE202642 macrophages: PGAM5 RNA expression',
                 fontsize=15, loc='left', pad=16)
    finish_map(fig, ax, sc, 'Per-cell PGAM5: log1p(CP10k)',
               '12,135 macrophages | 868 PGAM5-detected cells | 7 HCC patients\n'
               'Saved PGAM5-free UMAP; C0–C14 are RNA clusters. '
               'All zero values retained; no smoothing or imputation.',
               'PGAM5_expression_UMAP', legend=True)

    fig, ax = plt.subplots(figsize=(10, 8))
    has_expression = frame.cluster_mean_PGAM5_log1p_CP10k.gt(0)
    ax.scatter(frame.loc[~has_expression, 'UMAP1'],
               frame.loc[~has_expression, 'UMAP2'], s=5, color=ZERO_COLOR,
               edgecolors='none', rasterized=True)
    average_max = float(grouped.mean_PGAM5_log1p_CP10k.max())
    sc = ax.scatter(frame.loc[has_expression, 'UMAP1'],
                    frame.loc[has_expression, 'UMAP2'], s=5,
                    c=frame.loc[has_expression, 'cluster_mean_PGAM5_log1p_CP10k'],
                    cmap=cmap, norm=Normalize(0, average_max),
                    edgecolors='none', rasterized=True)
    label_clusters(ax, frame)
    ax.set_title('GSE202642 macrophages: PGAM5 mean expression by cluster',
                 fontsize=14, loc='left', pad=16)
    finish_map(fig, ax, sc, 'Cluster mean of log1p(CP10k)',
               'Each point receives its cluster mean, including zero-expression cells.\n'
               'This panel shows group averages, not each cell’s RNA level. '
               'C1 has no PGAM5 detection and is grey.',
               'PGAM5_cluster_mean_expression_UMAP')
    assert sha256(coord_path) == coord_hash
    result = {
        'dataset': 'GSE202642', 'macrophages': len(frame),
        'PGAM5_detected': int(frame.PGAM5_detected.sum()),
        'PGAM5_undetected': int((~frame.PGAM5_detected).sum()),
        'RNA_clusters': len(grouped), 'patient_count': int(frame.patient_id.nunique()),
        'UMAP_coordinate_SHA256_before_and_after': coord_hash,
        'UMAP_recomputed': False, 'one_to_one_cell_ID_join_passed': True,
        'raw_PGAM5_matches_portable_matrix': True,
        'normalization': 'raw PGAM5 / original full-library total_counts * 10000',
        'plotted_single_cell_value': 'natural log1p(CP10k)',
        'maximum_abs_difference_vs_stored_CP10k': norm_error,
        'single_cell_color_range': [0.0, vmax],
        'cluster_mean_color_range': [0.0, average_max],
        'cluster_mean': 'arithmetic mean of per-cell log1p(CP10k), including all zeros',
        'single_cell_zero_color': ZERO_COLOR,
        'positive_cells_drawn_after_zero_cells_in_ascending_expression_order': True,
        'expression_clipping': False, 'imputation': False, 'smoothing': False,
        'existing_cluster_mean_CP10k_agreement_passed': True,
        'cluster_weighted_mean_agreement_passed': True,
        'highest_mean_CP10k_cluster': 'C' + grouped.mean_PGAM5_CP10k.idxmax(),
        'highest_mean_log1p_CP10k_cluster': 'C' + grouped.mean_PGAM5_log1p_CP10k.idxmax(),
        'checks_passed': True,
    }
    (ROOT / 'PGAM5_expression_UMAP_validation.json').write_text(
        json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
