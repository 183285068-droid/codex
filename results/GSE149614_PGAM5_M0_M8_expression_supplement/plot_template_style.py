"""Render GSE149614 using the supplied GSE202642 UMAP figure layout."""
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

R = Path(__file__).resolve().parent
input_path = R / 'inputs/cell_expression_coordinates.csv.gz'
input_hash = hashlib.sha256(input_path.read_bytes()).hexdigest()
d = pd.read_csv(input_path, index_col=0)
assert d.index.is_unique and len(d) == 5605
assert d.total_counts.gt(0).all() and d.PGAM5_raw_counts.ge(0).all()
d['PGAM5_CP10k'] = d.PGAM5_raw_counts / d.total_counts * 10000
d['PGAM5_log1p_CP10k'] = np.log1p(d.PGAM5_CP10k)
d['detected'] = d.PGAM5_raw_counts.gt(0)
groups = d.groupby('mac_cluster').agg(cells=('Sample', 'size'), detected=('detected', 'sum'), mean_CP10k=('PGAM5_CP10k', 'mean'), mean_log1p_CP10k=('PGAM5_log1p_CP10k', 'mean'))
assert set(groups.index) == {f'M{i}' for i in range(9)}
previous = pd.read_csv(R / 'PGAM5_M0_M8_expression_summary.csv').set_index('cluster')
np.testing.assert_allclose(groups.mean_log1p_CP10k, previous.loc[groups.index, 'PGAM5_mean_log1p_CP10k'], rtol=1e-12, atol=1e-12)
assert groups.cells.sum() == 5605 and groups.detected.sum() == 136
assert np.isclose(np.average(groups.mean_log1p_CP10k, weights=groups.cells), d.PGAM5_log1p_CP10k.mean())
d['group_mean_log1p_CP10k'] = d.mac_cluster.map(groups.mean_log1p_CP10k)
labels = {
    'M0': 'M0  FOS/JUN\nImmediate-early response',
    'M1': 'M1  MKI67/TOP2A\nCycling',
    'M2': 'M2  C1QC/HLA-DRA\nAntigen presentation',
    'M3': 'M3  JUND/HSPA1A\nStress-associated / sample-dominated',
    'M4': 'M4  TIMD4/CD5L\nResident-like',
    'M5': 'M5  SPP1/LGALS1\nRemodeling-associated',
    'M6': 'M6  C1QA/CD3D\nMac/T-RNA mixed',
    'M7': 'M7  MT1G/MT2A\nMetal response',
    'M8': 'M8  GRM4/AIM1L\nMacrophage-enriched / unresolved',
}
centers = d.groupby('mac_cluster')[['UMAP1', 'UMAP2']].median()
ordered = centers.sort_values('UMAP1').index.tolist()

def label_groups(ax):
    for side, clusters in [('left', ordered[:5]), ('right', ordered[5:])]:
        rows = centers.loc[clusters].sort_values('UMAP2', ascending=False)
        for ypos, (c, center) in zip(np.linspace(.96, .04, len(rows)), rows.iterrows()):
            ax.annotate(labels[c], xy=(center.UMAP1, center.UMAP2), xycoords='data', xytext=(-.04 if side=='left' else 1.04, ypos), textcoords='axes fraction', ha='right' if side=='left' else 'left', va='center', fontsize=9.5, color='#18232C', linespacing=1.35, annotation_clip=False, bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': 1}, arrowprops={'arrowstyle': '-|>', 'color': '#87929B', 'mutation_scale': 7, 'lw': .65, 'shrinkA': 4, 'shrinkB': 1})

def finish(fig, ax, scatter, scale_label, caption, filename, legend=False):
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
    cbar.set_label(scale_label, fontsize=11)
    cbar.ax.tick_params(labelsize=9)
    if legend:
        fig.legend(handles=[Line2D([0],[0],linestyle='',marker='o',color='#C9CDD2',markersize=6,label='PGAM5 not detected (raw count = 0)')], loc='lower left', bbox_to_anchor=(.07,.16), frameon=False, fontsize=9)
    fig.text(.07,.025,caption,ha='left',va='bottom',fontsize=9,color='#47525C',linespacing=1.6)
    for ext in ['png','pdf']:
        fig.savefig(R/(filename+'.'+ext),dpi=240,facecolor='white')
    plt.close(fig)

plt.rcParams.update({'font.family': 'DejaVu Sans','font.size':11,'pdf.fonttype':42,'ps.fonttype':42})
cmap = LinearSegmentedColormap.from_list('PGAM5_YlOrRd', plt.colormaps['YlOrRd'](np.linspace(.20,1,256)))
fig = plt.figure(figsize=(14,10.5))
ax = fig.add_axes([.22,.23,.54,.64])
scatter = ax.scatter(d.UMAP1,d.UMAP2,s=5,c=d.group_mean_log1p_CP10k,cmap=cmap,norm=Normalize(0,float(groups.mean_log1p_CP10k.max())),edgecolors='none',rasterized=True)
label_groups(ax)
fig.text(.07,.96,'GSE149614 macrophage RNA groups: PGAM5 group means',fontsize=16)
fig.text(.07,.925,'Labels: marker genes + observed RNA program (provisional annotations)',fontsize=11,color='#47525C')
finish(fig,ax,scatter,'Cluster mean of log1p(CP10k)','Each point receives its group mean, including zero-expression cells.\n5,605 macrophage-enriched cells | 10 primary HCC tumor T samples | M0-M8.', 'PGAM5_cluster_mean_expression_UMAP_template_style')

fig = plt.figure(figsize=(14,10.5))
ax = fig.add_axes([.22,.23,.54,.64])
zero = d[~d.detected]
positive = d[d.detected].sort_values('PGAM5_log1p_CP10k')
ax.scatter(zero.UMAP1,zero.UMAP2,s=3,color='#C9CDD2',edgecolors='none',rasterized=True)
scatter = ax.scatter(positive.UMAP1,positive.UMAP2,s=13,c=positive.PGAM5_log1p_CP10k,cmap=cmap,norm=Normalize(0,float(d.PGAM5_log1p_CP10k.max())),edgecolors='none',rasterized=True)
label_groups(ax)
fig.text(.07,.96,'GSE149614 macrophage RNA groups: PGAM5 expression',fontsize=16)
fig.text(.07,.925,'Labels: marker genes + observed RNA program (provisional annotations)',fontsize=11,color='#47525C')
finish(fig,ax,scatter,'Per-cell PGAM5: log1p(CP10k)','5,605 macrophage-enriched cells | 136 PGAM5-detected cells | 10 primary HCC tumor T samples\nAll zero values retained. M6 is Mac/T-RNA mixed; subgroup names are provisional RNA annotations.', 'PGAM5_expression_UMAP_template_style',legend=True)
assert hashlib.sha256(input_path.read_bytes()).hexdigest() == input_hash
(R/'template_style_verification.json').write_text(json.dumps({'dataset':'GSE149614','style_reference':'GSE202642 PGAM5 group means supplied by user','cells':len(d),'PGAM5_detected':int(d.detected.sum()),'clusters_labeled':9,'UMAP_recomputed':False,'input_sha256_before_and_after':input_hash,'group_mean':'Arithmetic mean of per-cell natural log1p(CP10k), including all zeros','group_means_match_previous_statistics':True,'weighted_mean_matches_all_cells':True,'group_mean_color_range':[0,float(groups.mean_log1p_CP10k.max())],'per_cell_color_range':[0,float(d.PGAM5_log1p_CP10k.max())],'colormap':'YlOrRd, near-white portion trimmed (.20-1); full expression range','expression_clipping':False},indent=2))
print('Rendered template-style group-mean and per-cell UMAPs; 9 group means verified.')
