"""Describe PGAM5 RNA in source spots; no cell-specific abundance inference."""
from pathlib import Path
import hashlib
import json
import os

os.environ.setdefault('MPLCONFIGDIR', '/tmp/GSE245908_expression_mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
BASE = ROOT.parent
rows, spot_rows, manifest = [], [], {}
checks = 0
for sample in ['CHC20', 'CHC23']:
    paths = [BASE / f'{sample}_spot_states.csv',
             BASE / f'{sample}_directed_LR_edges.csv.gz']
    for path in paths:
        manifest[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    d = pd.read_csv(paths[0]).set_index('barcode')
    e = pd.read_csv(paths[1])
    pairs = [c for c in e.columns if '--' in c]
    mac = d.Mac_enriched
    nearby = mac & d.T_neighbor_count.gt(0)
    assert d.index.is_unique
    assert set(e.source_barcode) == set(d.index[nearby])
    assert d.loc[e.source_barcode.unique(), 'Mac_enriched'].all()
    assert d.loc[e.target_barcode.unique(), 'T_enriched'].all()
    reconstructed = np.log1p(1e4 * d.PGAM5_raw / d.total_counts)
    assert np.allclose(reconstructed, d.PGAM5_log1p_CP10k, atol=1e-7)
    assert d.PGAM5_log1p_CP10k.eq(0).equals(d.PGAM5_raw.eq(0))
    counts = e.groupby('source_barcode').size()
    assert counts.eq(d.loc[counts.index, 'T_neighbor_count']).all()
    checks += 7

    def source_mask(edge_mask):
        return d.index.isin(e.loc[edge_mask, 'source_barcode'].unique())

    masks = {
        'All_mac_RNA_spots': mac,
        'Adjacent_T_RNA_spots': nearby,
        'Without_T_RNA_neighbor': mac & d.T_neighbor_count.eq(0),
        'Any_candidate_LR_supported': source_mask(e[pairs].gt(0).any(axis=1)),
        'SPP1_CD44_supported': source_mask(e['SPP1--CD44'].gt(0)),
        'PVR_TIGIT_supported': source_mask(e['PVR--TIGIT'].gt(0)),
        'Adjacent_activation_supported_T': source_mask(e.target_T_state.isin(['Activation_supported', 'Both_programs'])),
        'Adjacent_exhaustion_supported_T': source_mask(e.target_T_state.isin(['Exhaustion_supported', 'Both_programs'])),
        'PGAM5_high_adjacent_T': nearby & d.source_group.eq('PGAM5_high'),
        'Other_detected_adjacent_T': nearby & d.source_group.eq('Other_detected'),
        'Undetected_adjacent_T': nearby & d.source_group.eq('Undetected'),
    }
    for region, mask in masks.items():
        a = d.loc[mask]
        assert a.Mac_enriched.all()
        checks += 1
        x = a.PGAM5_log1p_CP10k
        pos = a.PGAM5_raw.gt(0)
        rows.append(dict(sample=sample, region=region, n_unique_source_spots=len(a),
                         n_PGAM5_detected=int(pos.sum()),
                         detection_percent=100 * pos.mean(),
                         mean_log1p_CP10k_including_zeros=x.mean(),
                         median_log1p_CP10k_including_zeros=x.median(),
                         mean_log1p_CP10k_detected_only=x[pos].mean(),
                         mean_CP10k_including_zeros=(1e4*a.PGAM5_raw/a.total_counts).mean(),
                         mean_raw_UMI_including_zeros=a.PGAM5_raw.mean(),
                         min_raw_UMI=a.PGAM5_raw.min(), max_raw_UMI=a.PGAM5_raw.max()))
        z = a[['sample', 'PGAM5_raw', 'PGAM5_log1p_CP10k', 'total_counts',
               'source_group', 'T_neighbor_count']].reset_index()
        z.insert(1, 'region', region)
        spot_rows.append(z)
    assert masks['Adjacent_T_RNA_spots'].sum() + masks['Without_T_RNA_neighbor'].sum() == mac.sum()
    checks += 1

summary = pd.DataFrame(rows)
summary.to_csv(ROOT/'region_expression_summary.csv', index=False, float_format='%.10g')
pd.concat(spot_rows, ignore_index=True).to_csv(ROOT/'region_source_spot_expression.csv.gz', index=False)
(ROOT/'input_sha256.json').write_text(json.dumps(manifest, indent=2)+'\n')
(ROOT/'validation.json').write_text(json.dumps({'passed_checks': checks,
    'unit': 'unique mixed-cell source spot within each region',
    'checks': 'unique IDs, original expression formula, zero retention, lineage gates, adjacency counts, source deduplication and partition sizes'}, indent=2)+'\n')

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'pdf.fonttype': 42})
regions = ['All_mac_RNA_spots', 'Adjacent_T_RNA_spots', 'Without_T_RNA_neighbor', 'Any_candidate_LR_supported']
labels = ['All mac-RNA\nspots', 'Adjacent to\nT-RNA spots', 'No T-RNA\nneighbor', 'Candidate LR\nco-detection']
colors = ['#9DA9B5', '#D96370', '#77A4BF', '#C99155']
fig, axs = plt.subplots(2, 2, figsize=(12.8, 9))
for i, sample in enumerate(['CHC20', 'CHC23']):
    a = summary[summary['sample'].eq(sample)].set_index('region').loc[regions]
    for j, (metric, title, ymax) in enumerate([
        ('mean_log1p_CP10k_including_zeros', 'PGAM5 mean log1p(CP10k)', .20),
        ('detection_percent', 'PGAM5 RNA detection (%)', 55)]):
        ax = axs[i,j]
        bars = ax.bar(range(4), a[metric], color=colors, width=.65)
        for k, bar in enumerate(bars):
            r = a.iloc[k]
            value = f'{r[metric]:.4f}' if j == 0 else f'{r[metric]:.2f}%'
            ax.text(k, bar.get_height()+ymax*.025, value, ha='center', fontsize=10)
        ax.set_xticks(range(4), [f'{label}\nn={int(n)}' for label,n in zip(labels, a.n_unique_source_spots)])
        ax.set_ylim(0, ymax)
        ax.set_ylabel(title)
        ax.set_title(sample)
        ax.spines[['top', 'right']].set_visible(False)
fig.suptitle('GSE245908: PGAM5 RNA in macrophage-enriched source regions', fontsize=16)
fig.text(.04, .018, 'Mixed-cell Visium spots; PGAM5 cannot be assigned specifically to macrophages.\nZeros included. Source spots counted once per category; categories overlap. Descriptive results, no significance claim.', fontsize=10)
fig.tight_layout(rect=[0,.085,1,.95])
fig.savefig(ROOT/'PGAM5_expression_T_neighbor_regions.png', dpi=220, bbox_inches='tight')
fig.savefig(ROOT/'PGAM5_expression_T_neighbor_regions.pdf', bbox_inches='tight')
plt.close(fig)
print(summary.to_string(index=False))
print(f'Passed {checks} checks.')
