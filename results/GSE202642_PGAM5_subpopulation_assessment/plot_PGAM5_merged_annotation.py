"""Merge PGAM5 RNA-detected cells into one display label on the saved UMAP.

This is an annotation merge, not a new unsupervised cluster or new embedding.
Original functional names, mixed-RNA flags and all cell coordinates are retained.
"""
from pathlib import Path
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parent
TARGET = 'Macrophage_PGAM5_RNA_detected'
TARGET_COLOR = '#C62828'
ZERO_COLOR = '#CCD1D5'


def configure(ax):
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_aspect('equal', adjustable='box')
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def marker(color, text):
    return Line2D([0], [0], linestyle='', marker='o', color=color,
                  markersize=6, label=text)


def save(fig, name):
    for extension in ['png', 'pdf']:
        fig.savefig(ROOT / (name + '.' + extension), dpi=240, facecolor='white')
    plt.close(fig)


def main():
    path = ROOT / 'UMAP_coordinates.csv.gz'
    coordinate_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    cells = pd.read_csv(ROOT / 'functional_subgroup_cell_annotations.csv.gz')
    coords = pd.read_csv(path)
    assert cells.cell_id.is_unique and coords.cell_id.is_unique
    assert set(cells.cell_id) == set(coords.cell_id)
    frame = cells.merge(coords, on='cell_id', validate='one_to_one', sort=False)
    obs = pd.read_csv(ROOT / 'cell_annotations.csv.gz', index_col=0).set_index('cell_id')
    assert frame.cell_id.tolist() == cells.cell_id.tolist()
    assert np.isfinite(frame[['UMAP1', 'UMAP2']]).all().all()
    # Check detection against portable raw counts in their saved row order.
    genes = pd.read_csv(ROOT / 'portable_feature_genes.csv').gene
    indices = np.flatnonzero(genes.eq('PGAM5').to_numpy())
    assert len(indices) == 1
    raw = sparse.load_npz(ROOT / 'portable_raw_feature_counts.npz')
    source_pg = pd.Series(raw[:, int(indices[0])].toarray().ravel(), index=obs.index)
    assert np.array_equal(source_pg.loc[frame.cell_id], frame.PGAM5_counts)
    positive = frame.PGAM5_counts.gt(0)
    assert np.array_equal(positive, frame.PGAM5_detected)
    frame['merged_PGAM5_RNA_group'] = np.where(positive, TARGET,
                                               'Macrophage_PGAM5_RNA_undetected')
    frame['merged_display_annotation'] = np.where(positive, TARGET, frame.functional_name_EN)
    frame['annotation_merge_not_unsupervised_cluster'] = True
    assert frame.loc[positive, 'merged_display_annotation'].nunique() == 1
    assert frame.loc[~positive, 'merged_display_annotation'].equals(
        frame.loc[~positive, 'functional_name_EN'])
    frame.to_csv(ROOT / 'PGAM5_merged_group_cell_annotations.csv.gz', index=False)
    origins = frame.groupby(['PGAM5_free_cluster', 'functional_name_EN', 'functional_name_CN'],
                            sort=False).agg(cells=('cell_id', 'size'),
                                            PGAM5_detected_cells=('PGAM5_detected', 'sum')).reset_index()
    origins = origins.sort_values('PGAM5_free_cluster', key=lambda s: s.str[1:].astype(int))
    origins['PGAM5_undetected_cells'] = origins.cells - origins.PGAM5_detected_cells
    origins['fraction_of_merged_PGAM5_detected_group'] = origins.PGAM5_detected_cells / positive.sum()
    assert origins.PGAM5_detected_cells.sum() == positive.sum()
    origins.to_csv(ROOT / 'PGAM5_merged_group_origin_summary.csv', index=False)

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'pdf.fonttype': 42})
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_axes([.08, .20, .84, .66])
    ax.scatter(frame.loc[~positive, 'UMAP1'], frame.loc[~positive, 'UMAP2'],
               s=4, color=ZERO_COLOR, edgecolors='none', rasterized=True)
    ax.scatter(frame.loc[positive, 'UMAP1'], frame.loc[positive, 'UMAP2'],
               s=13, color=TARGET_COLOR, edgecolors='none', rasterized=True)
    configure(ax)
    fig.text(.06, .955, 'GSE202642: one PGAM5 RNA-detected annotation group', fontsize=15)
    fig.legend(handles=[marker(TARGET_COLOR, f'PGAM5 RNA detected: one label (n={positive.sum():,})'),
                        marker(ZERO_COLOR, f'PGAM5 RNA undetected (n={(~positive).sum():,})')],
               loc='lower left', bbox_to_anchor=(.06, .08), frameon=False, fontsize=10)
    fig.text(.06, .027, 'Raw PGAM5 > 0 defines RNA detection; all zero values retained.\n'
             'Annotation merge only; original PGAM5-free UMAP. A shared label does not establish one independent cluster.',
             fontsize=9, color='#47525C')
    save(fig, 'PGAM5_RNA_detected_merged_group_UMAP')

    # One target color, with undetected cells retaining their functional names.
    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_axes([.06, .18, .53, .68])
    handles = [marker(TARGET_COLOR, f'PGAM5 RNA detected | n={positive.sum():,}')]
    palette = plt.colormaps['tab20']
    for i, row in enumerate(origins.itertuples()):
        mask = (~positive) & frame.PGAM5_free_cluster.eq(row.PGAM5_free_cluster)
        color = tuple(.5 * np.array(palette(i)[:3]) + .5)
        ax.scatter(frame.loc[mask, 'UMAP1'], frame.loc[mask, 'UMAP2'],
                   s=4, color=color, edgecolors='none', rasterized=True)
        handles.append(marker(color, row.functional_name_EN.replace(' — ', ' | ') +
                              f' | n={mask.sum():,}'))
    ax.scatter(frame.loc[positive, 'UMAP1'], frame.loc[positive, 'UMAP2'],
               s=13, color=TARGET_COLOR, edgecolors='none', rasterized=True)
    configure(ax)
    fig.legend(handles=handles, loc='center left', bbox_to_anchor=(.61, .53),
               frameon=False, fontsize=8.5, labelspacing=.9, handletextpad=.7)
    fig.text(.06, .955, 'GSE202642: unified PGAM5 label with functional RNA backgrounds', fontsize=15)
    fig.text(.06, .91, 'All PGAM5-detected cells share the red label; other cells retain their provisional functional names.', fontsize=10)
    fig.text(.06, .035, 'Same saved UMAP and original RNA clusters; no cells removed or coordinates changed.\n'
             'Mixed-RNA flags retained. Merged PGAM5 annotation is an expression status, not a validated independent subtype.',
             fontsize=9, color='#47525C')
    save(fig, 'PGAM5_merged_group_functional_background_UMAP')
    assert hashlib.sha256(path.read_bytes()).hexdigest() == coordinate_hash
    result = {'dataset': 'GSE202642', 'cells': len(frame),
              'merged_PGAM5_RNA_detected_group_cells': int(positive.sum()),
              'PGAM5_RNA_undetected_cells': int((~positive).sum()),
              'origin_RNA_clusters_with_detected_cells': int(origins.PGAM5_detected_cells.gt(0).sum()),
              'annotation_rule': 'raw PGAM5 counts > 0; original inferred macrophage set',
              'all_detected_cells_share_one_display_label': True,
              'raw_count_detection_and_one_to_one_cell_join_passed': True,
              'original_functional_names_and_mixed_RNA_flags_retained': True,
              'new_unsupervised_clustering_performed': False,
              'UMAP_recomputed': False, 'cells_removed': 0,
              'coordinate_SHA256_before_and_after': coordinate_hash,
              'independent_subtype_established_by_label_merge': False,
              'checks_passed': True}
    (ROOT / 'PGAM5_merged_annotation_validation.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
