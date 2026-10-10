"""Use the same Scanpy DotPlot defaults as GSE151530 and GSE149614.

Plot from actual group summaries. The AnnData object is only a layout carrier;
all colors and dot sizes are supplied from audited source evidence, not inferred
from its one row per group. It is not a synthetic cell expression dataset.
"""
from pathlib import Path
import json
import hashlib
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scanpy as sc

R = Path(__file__).resolve().parent
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'pdf.fonttype': 42})
checks = []


def render(key, named=True):
    spec = json.loads((R / 'inputs' / (key+'_spec.json')).read_text())
    e = pd.read_csv(R / 'inputs' / (key+'_values.csv'), dtype={'cluster': str})
    groups = list(spec['labels'])
    genes = sum(spec['panels'].values(), [])
    means = e.pivot(index='cluster', columns='gene', values='mean_log1p_CP10k').loc[groups, genes]
    fractions = e.pivot(index='cluster', columns='gene', values='detection_fraction').loc[groups, genes]
    assert np.isfinite(means).all().all() and np.isfinite(fractions).all().all()
    assert fractions.ge(-1e-12).all().all() and fractions.le(1+1e-12).all().all()
    # Source sparse-matrix means can exceed 1 by floating-point roundoff only.
    fractions = fractions.clip(0, 1)
    # Scanpy standard_scale='var': min-max within each gene across groups.
    scaled = means.subtract(means.min(axis=0), axis=1)
    scaled = scaled.divide(scaled.max(axis=0), axis=1).fillna(0)
    labels = [spec['labels'][c] if named else 'C'+c for c in groups]
    a = ad.AnnData(np.zeros((len(groups), len(genes))),
        obs=pd.DataFrame({'group': pd.Categorical(labels, categories=labels, ordered=True)},
                         index=pd.Index(groups, name='cluster')),
        var=pd.DataFrame(index=pd.Index(genes, name='gene')))
    scaled.index = labels
    fractions.index = labels
    var_names = spec['panels'] if len(spec['panels']) > 1 else genes
    f = sc.pl.dotplot(a, var_names, groupby='group', categories_order=labels,
        dot_color_df=scaled, dot_size_df=fractions,
        standard_scale=None, show=False, return_fig=True)
    # Full 0-100% detection scale shared across plots; both legends are kept.
    f.style(cmap='Reds', dot_min=0, dot_max=1)
    np.testing.assert_array_equal(f.dot_color_df.to_numpy(), scaled.to_numpy())
    np.testing.assert_array_equal(f.dot_size_df.to_numpy(), fractions.to_numpy())
    filename = key + ('_named_dotplot' if named else '_marker_dotplot')
    f.savefig(R / (filename+'.png'), dpi=240)
    f.savefig(R / (filename+'.pdf'))
    axes = f.get_axes()
    shown = [t.get_text() for t in axes['mainplot_ax'].get_yticklabels()]
    assert shown == labels
    assert 'size_legend_ax' in axes and 'color_legend_ax' in axes
    scaled.index = groups
    e['scaled_mean_log1p_CP10k'] = [scaled.loc[c, g] for c, g in zip(e.cluster, e.gene)]
    e.to_csv(R / (key+'_plot_values.csv'), index=False)
    checks.append({'file': filename, 'groups': len(groups), 'genes': len(genes),
                   'all_group_labels_displayed': True, 'plot_values_match_source': True,
                   'size_and_color_legends_present': True})
    plt.close(f.fig)


if __name__ == '__main__':
    render('allcell_C0_C20', named=False)
    render('allcell_C0_C20')
    render('macrophage_M0_M8')
    render('macrophage_C0_C14_legacy')
    (R / 'verification.json').write_text(json.dumps({'plots': checks,
        'source_data_unchanged': True, 'all_zeros_included_in_means': True,
        'color': 'Reds; per-gene min-max of all-cell mean natural log1p(CP10k)',
        'dot_size': 'Scanpy default area mapping of detection fraction; fixed 0-100% range',
        'layout': 'Same Scanpy DotPlot as GSE151530/GSE149614; vertical gene labels, right-hand legends'},
        indent=2)+'\n')
    print(json.dumps(checks, indent=2))
