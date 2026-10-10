"""Reproduce the supplied GSE202642 figure style on saved GSE151530 data.

Run with Python 3.12, numpy, pandas, matplotlib and (input extraction only) anndata.
Portable inputs are sufficient to rerender without the original h5ad.
"""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'GSE151530_HCC_hierarchical_clustering'
ZERO_COLOR = '#C9CDD2'


def sha256(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def prepare_inputs():
    import anndata as ad
    dest = ROOT / 'inputs'
    dest.mkdir(exist_ok=True)
    a = ad.read_h5ad(SOURCE / 'macrophage_marker_subset.h5ad')
    pg = a.layers['counts'][:, a.var_names.get_loc('PGAM5')].toarray().ravel()
    np.testing.assert_array_equal(pg, a.obs.PGAM5_counts)
    coordinates = pd.read_csv(SOURCE / 'macrophage_UMAP_coordinates.csv.gz', index_col=0)
    assert a.obs_names.equals(coordinates.index)
    np.testing.assert_array_equal(a.obsm['X_umap'], coordinates.to_numpy().astype(a.obsm['X_umap'].dtype))
    d = a.obs[['mac_cluster', 'Sample', 'patient_proxy', 'total_counts']].copy()
    d.index.name = 'cell_id'
    d['PGAM5_raw_counts'] = pg
    d[['UMAP1', 'UMAP2']] = a.obsm['X_umap']
    d.to_csv(dest / 'cell_expression_coordinates.csv.gz')
    for name in ['macrophage_cluster_annotations.csv', 'PGAM5_by_macrophage_cluster.csv']:
        pd.read_csv(SOURCE / name).to_csv(dest / name, index=False)
    names = ['macrophage_marker_subset.h5ad', 'macrophage_UMAP_coordinates.csv.gz',
             'macrophage_cluster_annotations.csv', 'PGAM5_by_macrophage_cluster.csv']
    (ROOT / 'source_provenance.json').write_text(json.dumps({
        'dataset': 'GSE151530', 'source_directory': SOURCE.name,
        'scope': 'Saved HCC macrophage-enriched subset: author TAMs intersect C1Q-rich all-cell C9/C11',
        'source_sha256': {n: sha256(SOURCE / n) for n in names},
        'raw_PGAM5_verified_against_matrix': True,
        'UMAP_exact_agreement_with_source_coordinate_file': True,
    }, indent=2) + '\n')


LABELS = {
    'M0': 'FCGR3A/PSAP\nMacrophage RNA program',
    'M1': 'RPL12/RPS27\nRibosome-high',
    'M2': 'MKI67/TOP2A\nCycling',
    'M3': 'LINC01419/NOVA1\nMacrophage-enriched / unresolved',
    'M4': 'ALOX5AP/COTL1\nMacrophage RNA program',
    'M5': 'CCL3/CCL4\nChemokine-associated / hepatic RNA',
    'M6': 'CLEC10A/FCN1\nC1Q-rich myeloid boundary',
    'M7': 'FOLR2/SEPP1\nResident-like',
    'M8': 'SPP1/GPNMB\nMacrophage RNA program',
    'M9': 'TIMD4/LYVE1\nResident-like',
    'M10': 'CXCL9/CXCL10\nIFN-associated',
    'M11': 'CHIT1/CCL18\nMacrophage RNA program',
    'M12': 'HOXA13/UBE2C\nCycling / mixed RNA',
    'M13': 'HLA-A/HLA-C\nPatient-associated',
    'M14': 'C1QA/CD3D\nMac/T-RNA mixed',
    'M15': 'C1QA/EOMES\nMac/NK-RNA mixed',
}


def label_clusters(ax, frame):
    centers = frame.groupby('mac_cluster')[['UMAP1', 'UMAP2']].median()
    ordered = centers.sort_values('UMAP1').index.tolist()
    for side, clusters in [('left', ordered[:8]), ('right', ordered[8:])]:
        rows = centers.loc[clusters].sort_values('UMAP2', ascending=False)
        for ypos, (cluster, center) in zip(np.linspace(.96, .04, len(rows)), rows.iterrows()):
            ax.annotate(LABELS[cluster], xy=(center.UMAP1, center.UMAP2), xycoords='data',
                        xytext=(-.04 if side == 'left' else 1.04, ypos),
                        textcoords='axes fraction', ha='right' if side == 'left' else 'left',
                        va='center', fontsize=9.5, color='#18232C', linespacing=1.35,
                        annotation_clip=False,
                        bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': 1},
                        arrowprops={'arrowstyle': '-|>', 'color': '#87929B',
                                    'mutation_scale': 7, 'lw': .65, 'shrinkA': 4, 'shrinkB': 1})


def finish_map(fig, ax, scatter, label, caption, name, legend=False):
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2', rotation=0, ha='left')
    ax.yaxis.set_label_coords(0, 1.025)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    cax = fig.add_axes([.37, .125, .28, .018])
    cbar = fig.colorbar(scatter, cax=cax, orientation='horizontal')
    cbar.set_label(label, fontsize=11)
    cbar.ax.tick_params(labelsize=9)
    if legend:
        fig.legend(handles=[Line2D([0], [0], linestyle='', marker='o', color=ZERO_COLOR,
                                 markersize=6, label='PGAM5 not detected (raw count = 0)')],
                   loc='lower left', bbox_to_anchor=(.07, .16), frameon=False, fontsize=9)
    fig.text(.07, .025, caption, ha='left', va='bottom', fontsize=9,
             color='#47525C', linespacing=1.6)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    # Check that all 16 names and connecting arrows stay on the exported canvas.
    for annotation in ax.texts:
        box = annotation.get_window_extent(renderer)
        assert box.x0 >= 0 and box.x1 <= fig.bbox.width, annotation.get_text()
        assert box.y0 >= 0 and box.y1 <= fig.bbox.height, annotation.get_text()
    fig.savefig(ROOT / (name + '.png'), dpi=240, facecolor='white')
    fig.savefig(ROOT / (name + '.pdf'), facecolor='white')
    plt.close(fig)


def main():
    input_path = ROOT / 'inputs/cell_expression_coordinates.csv.gz'
    if not input_path.exists():
        prepare_inputs()
    input_hash = sha256(input_path)
    frame = pd.read_csv(input_path, index_col=0)
    assert frame.index.is_unique and len(frame) == 3095
    assert frame.total_counts.gt(0).all() and frame.PGAM5_raw_counts.ge(0).all()
    assert np.isfinite(frame[['UMAP1', 'UMAP2', 'total_counts', 'PGAM5_raw_counts']]).all().all()
    frame['PGAM5_CP10k'] = frame.PGAM5_raw_counts / frame.total_counts * 10000
    frame['PGAM5_log1p_CP10k'] = np.log1p(frame.PGAM5_CP10k)
    frame['PGAM5_detected'] = frame.PGAM5_raw_counts.gt(0)
    grouped = frame.groupby('mac_cluster').agg(
        cells=('Sample', 'size'), PGAM5_detected=('PGAM5_detected', 'sum'),
        PGAM5_detection_percent=('PGAM5_detected', lambda x: x.mean() * 100),
        PGAM5_mean_CP10k_all_cells=('PGAM5_CP10k', 'mean'),
        PGAM5_median_CP10k_all_cells=('PGAM5_CP10k', 'median'),
        PGAM5_mean_log1p_CP10k=('PGAM5_log1p_CP10k', 'mean'),
        samples=('Sample', 'nunique'))
    grouped = grouped.loc[sorted(grouped.index, key=lambda x: int(x[1:]))]
    assert set(grouped.index) == set(LABELS)
    previous = pd.read_csv(ROOT / 'inputs/PGAM5_by_macrophage_cluster.csv').set_index('cluster')
    for col in grouped:
        np.testing.assert_allclose(grouped[col], previous.loc[grouped.index, col], rtol=1e-12, atol=1e-12)
    assert grouped.cells.sum() == len(frame) and grouped.PGAM5_detected.sum() == 153
    assert np.isclose(np.average(grouped.PGAM5_mean_log1p_CP10k, weights=grouped.cells),
                      frame.PGAM5_log1p_CP10k.mean())
    names = pd.read_csv(ROOT / 'inputs/macrophage_cluster_annotations.csv').set_index('cluster')
    grouped = names.loc[grouped.index].join(grouped)
    grouped['plot_label'] = pd.Series(LABELS)
    grouped.rename_axis('cluster').reset_index().to_csv(ROOT / 'PGAM5_expression_by_cluster.csv', index=False)
    frame['cluster_mean_PGAM5_log1p_CP10k'] = frame.mac_cluster.map(grouped.PGAM5_mean_log1p_CP10k)
    frame.to_csv(ROOT / 'PGAM5_UMAP_expression_values.csv.gz')

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    cmap = LinearSegmentedColormap.from_list(
        'PGAM5_YlOrRd', plt.colormaps['YlOrRd'](np.linspace(.20, 1, 256)))
    for mode in ['group_means', 'per_cell']:
        fig = plt.figure(figsize=(14, 10.5))
        ax = fig.add_axes([.22, .23, .54, .64])
        if mode == 'group_means':
            column = 'cluster_mean_PGAM5_log1p_CP10k'
            positive = frame[frame[column].gt(0)]
            zero = frame[frame[column].eq(0)]
            ax.scatter(zero.UMAP1, zero.UMAP2, s=5, color=ZERO_COLOR,
                       edgecolors='none', rasterized=True)
            sc = ax.scatter(positive.UMAP1, positive.UMAP2, s=5, c=positive[column],
                            cmap=cmap, norm=Normalize(0, frame[column].max()),
                            edgecolors='none', rasterized=True)
            title = 'GSE151530 macrophage RNA groups: PGAM5 group means'
            caption = ('Each point receives its cluster mean, including zero-expression cells.\n'
                       'Group averages, not each cell’s RNA level. HOXA13/UBE2C and C1QA/EOMES: no PGAM5 detection (grey).')
            scale = 'Cluster mean of log1p(CP10k)'
            filename = 'PGAM5_cluster_mean_expression_UMAP'
        else:
            column = 'PGAM5_log1p_CP10k'
            zero = frame[~frame.PGAM5_detected]
            positive = frame[frame.PGAM5_detected].sort_values(column)
            ax.scatter(zero.UMAP1, zero.UMAP2, s=3, color=ZERO_COLOR,
                       edgecolors='none', rasterized=True)
            sc = ax.scatter(positive.UMAP1, positive.UMAP2, s=13, c=positive[column],
                            cmap=cmap, norm=Normalize(0, frame[column].max()),
                            edgecolors='none', rasterized=True)
            title = 'GSE151530 macrophage RNA groups: PGAM5 expression'
            caption = ('3,095 macrophage-enriched cells | 153 PGAM5-detected cells | 28 HCC samples\n'
                       'Saved PGAM5-free UMAP. Mixed = other-lineage RNA signal. All zero values retained; no smoothing or imputation.')
            scale = 'Per-cell PGAM5: log1p(CP10k)'
            filename = 'PGAM5_expression_UMAP'
        label_clusters(ax, frame)
        fig.text(.07, .96, title, fontsize=16)
        fig.text(.07, .925, 'Labels: marker genes + observed RNA program (provisional annotations)',
                 fontsize=11, color='#47525C')
        finish_map(fig, ax, sc, scale, caption, filename, legend=(mode == 'per_cell'))
    assert sha256(input_path) == input_hash
    result = {'dataset': 'GSE151530', 'cells': len(frame), 'clusters_labeled': len(grouped),
              'PGAM5_detected': int(frame.PGAM5_detected.sum()),
              'samples': int(frame.Sample.nunique()), 'patient_proxies': int(frame.patient_proxy.nunique()),
              'template': 'GSE202642 supplied figure and original plot_PGAM5_expression_umap.py',
              'labels': 'Marker genes + RNA program, without cluster-number prefix as in template; ID mapping in CSV',
              'zero_expression_groups': grouped.index[grouped.PGAM5_detected.eq(0)].tolist(),
              'input_sha256_before_and_after': input_hash,
              'UMAP_recomputed': False, 'all_zeros_retained': True,
              'statistics_match_original': True, 'all_label_bounds_pass': True,
              'normalization': 'PGAM5 raw counts / full original library total * 10000',
              'group_mean': 'Arithmetic mean of per-cell natural log1p(CP10k), including zeros',
              'group_mean_color_range': [0, float(grouped.PGAM5_mean_log1p_CP10k.max())],
              'per_cell_color_range': [0, float(frame.PGAM5_log1p_CP10k.max())],
              'expression_clipping': False, 'imputation': False, 'smoothing': False}
    (ROOT / 'verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
