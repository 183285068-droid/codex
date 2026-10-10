"""Illustrate a forced PGAM5-detection island, explicitly NOT a new UMAP.

Left: saved expression-derived UMAP. Right: arbitrary label-imposed positions.
Positive positions form a sunflower disk using cell ID order only. They contain
no transcriptome, distance, clustering or mechanistic information.
"""
from pathlib import Path
import argparse
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'GSE202642_PGAM5_subpopulation_assessment' / 'PGAM5_merged_group_cell_annotations.csv.gz'
RED = '#C62828'
GREY = '#CDD2D6'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=SOURCE)
    args = parser.parse_args()
    source_hash = hashlib.sha256(args.source.read_bytes()).hexdigest()
    cells = pd.read_csv(args.source).sort_values('cell_id').reset_index(drop=True)
    assert cells.cell_id.is_unique
    assert np.isfinite(cells[['PGAM5_counts', 'UMAP1', 'UMAP2']]).all().all()
    assert cells.PGAM5_counts.ge(0).all()
    pos = cells.PGAM5_counts.gt(0).to_numpy()
    assert np.array_equal(pos, cells.PGAM5_detected)
    n_positive, n_zero = int(pos.sum()), int((~pos).sum())
    assert n_positive > 0 and n_zero > 0
    n_origins = cells.loc[pos, 'PGAM5_free_cluster'].nunique()

    display = np.zeros((len(cells), 2), dtype=float)
    order = np.arange(n_positive)
    radius = .8 * np.sqrt((order + .5) / n_positive)
    angle = order * np.pi * (3 - np.sqrt(5))
    display[pos, 0] = -1.8 + radius * np.cos(angle)
    display[pos, 1] = radius * np.sin(angle)
    # Undetected cells preserve their relative original coordinates, rescaled
    # into the other illustration region. Axes are explicitly arbitrary.
    xy = cells.loc[~pos, ['UMAP1', 'UMAP2']].to_numpy()
    midpoint = (xy.min(axis=0) + xy.max(axis=0)) / 2
    scale = max(np.ptp(xy, axis=0)) / 2
    assert scale > 0
    display[~pos] = (xy - midpoint) / scale + np.array([1.8, 0])
    assert display[pos, 0].max() < display[~pos, 0].min()

    values = cells[['cell_id', 'sample_name', 'PGAM5_counts', 'PGAM5_detected',
                    'PGAM5_CP10k', 'PGAM5_free_cluster', 'functional_name_EN',
                    'mixed_RNA_flag', 'UMAP1', 'UMAP2']].copy()
    values['illustrative_x_NOT_UMAP'] = display[:, 0]
    values['illustrative_y_NOT_UMAP'] = display[:, 1]
    values['label_imposed_schematic_not_a_fitted_embedding'] = True
    values.to_csv(ROOT / 'schematic_coordinates_NOT_UMAP.csv.gz', index=False)

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 2, figsize=(14, 7.5))
    for ax, coord in zip(axes, [cells[['UMAP1', 'UMAP2']].to_numpy(), display]):
        ax.scatter(coord[~pos, 0], coord[~pos, 1], s=4, color=GREY,
                   edgecolors='none', rasterized=True)
        ax.scatter(coord[pos, 0], coord[pos, 1], s=13, color=RED,
                   edgecolors='none', rasterized=True)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    axes[0].set(title='A. Saved expression-derived UMAP', xlabel='UMAP 1', ylabel='UMAP 2')
    axes[1].set(title='B. Label-imposed schematic — NOT UMAP',
                xlabel='Illustrative x (arbitrary)', ylabel='Illustrative y (arbitrary)',
                xlim=(-2.95, 3.05), ylim=(-1.45, 1.45))
    axes[1].text(-1.8, 1.05, 'PGAM5 detected\nArtificial island',
                 ha='center', va='bottom', color=RED, fontsize=11)
    axes[1].text(1.8, 1.05, 'PGAM5 undetected', ha='center', va='bottom',
                 color='#56636D', fontsize=11)
    fig.text(.055, .945, 'GSE202642: original UMAP and a forced PGAM5-detection island', fontsize=16)
    fig.text(.055, .901, f'Same {len(cells):,} cells in both panels; RNA detection = raw PGAM5 count > 0.',
             fontsize=11, color='#47525C')
    handles = [Line2D([], [], linestyle='', marker='o', markersize=6, color=color, label=label)
               for color, label in [(RED, f'PGAM5 RNA detected (n={n_positive:,})'),
                                    (GREY, f'PGAM5 RNA undetected (n={n_zero:,})')]]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .125),
               frameon=False, ncol=2, fontsize=10)
    fig.text(.055, .033,
             'Right-panel separation is imposed by the detection label; it is not learned from expression similarity.\n'
             'The artificial island provides no evidence of an independent subtype or a valid deconvolution signature.',
             fontsize=10, color='#47525C')
    fig.subplots_adjust(left=.06, right=.97, top=.80, bottom=.23, wspace=.23)
    for ext in ['png', 'pdf']:
        fig.savefig(ROOT / ('original_UMAP_vs_label_imposed_schematic.' + ext),
                    dpi=240, facecolor='white')
    plt.close(fig)
    assert hashlib.sha256(args.source.read_bytes()).hexdigest() == source_hash
    provenance = {
        'dataset': 'GSE202642', 'source_file': args.source.name,
        'source_SHA256': source_hash, 'cells': len(cells),
        'PGAM5_RNA_detected_cells': n_positive, 'PGAM5_RNA_undetected_cells': n_zero,
        'original_RNA_clusters_with_PGAM5_detection': int(n_origins),
        'left_panel': 'Unchanged saved expression-derived UMAP coordinates',
        'right_panel': 'Arbitrary label-imposed schematic, NOT UMAP or supervised UMAP',
        'positive_position_rule': 'Sunflower disk in cell-ID order; no expression features or distances',
        'undetected_position_rule': 'Saved UMAP uniformly scaled and translated to a separate illustration region',
        'no_new_UMAP_or_clustering_model_fitted': True,
        'all_source_cells_and_original_annotations_retained': True,
        'source_file_unchanged': True,
        'independent_subtype_or_deconvolution_validity_established': False}
    (ROOT / 'illustration_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps(provenance, indent=2))


if __name__ == '__main__':
    main()
