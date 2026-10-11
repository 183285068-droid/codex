"""Reproduce the intersection figures from the three previously reported TOP20 lists."""
from pathlib import Path
import csv
import hashlib
import itertools
import json
import shutil
import zipfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np

OUT = Path(__file__).resolve().parent
BASE = OUT.parent
SOURCES = {
    'GSE202642': 'GSE202642_C7_PGAM5_virtual_KO_TOP20/TOP20_downstream_genes.csv',
    'GSE149614': 'GSE149614_M1_PGAM5_virtual_KO/TOP20_downstream_genes.csv',
    'GSE151530': 'GSE151530_M2_PGAM5_virtual_KO/unsparsified/TOP20_downstream_genes.csv',
}
LABELS = ['GSE202642\nC7 | v2 ensemble', 'GSE149614\nM1 | sparse main model',
          'GSE151530\nM2 | unfiltered sensitivity model']


def main():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 12,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    (OUT / 'inputs').mkdir(exist_ok=True)
    ranks, provenance = {}, {}
    for dataset, relative in SOURCES.items():
        source = BASE / relative
        local = OUT / 'inputs' / f'{dataset}_TOP20.csv'
        if source.exists():
            shutil.copyfile(source, local)
        with local.open() as f:
            rows = list(csv.DictReader(f))
        ranks[dataset] = {r['gene'].strip(): int(r['rank']) for r in rows}
        assert len(rows) == len(ranks[dataset]) == 20
        assert sorted(ranks[dataset].values()) == list(range(1, 21))
        assert 'PGAM5' not in ranks[dataset]
        provenance[dataset] = {'source': relative, 'sha256': hashlib.sha256(local.read_bytes()).hexdigest()}
    names = list(ranks)
    sets = [set(ranks[d]) for d in names]
    regions = {}
    for mask in itertools.product([0, 1], repeat=3):
        if not any(mask):
            continue
        included = [sets[i] for i, bit in enumerate(mask) if bit]
        excluded = set().union(*(sets[i] for i, bit in enumerate(mask) if not bit))
        regions[''.join(map(str, mask))] = sorted(set.intersection(*included) - excluded)
    assert sum(map(len, regions.values())) == len(set.union(*sets)) == 56
    assert [sum(len(v) for m, v in regions.items() if m[i] == '1') for i in range(3)] == [20, 20, 20]
    assert regions['110'] == ['DNAJB1', 'HSP90AA1', 'HSPH1']
    assert regions['011'] == ['ANXA1'] and not regions['111'] and not regions['101']
    shared = sorted(g for g in set.union(*sets) if sum(g in s for s in sets) >= 2)
    with (OUT / 'shared_genes_rank_table.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['gene', 'datasets_in_TOP20', *[d + '_rank' for d in names]])
        for gene in shared:
            writer.writerow([gene, sum(gene in s for s in sets), *[ranks[d].get(gene, '') for d in names]])
    with (OUT / 'exclusive_intersections.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['membership_GSE202642_GSE149614_GSE151530', 'count', 'genes'])
        for mask, genes in regions.items():
            writer.writerow([mask, len(genes), ';'.join(genes)])
    summary = {'datasets': names, 'sources': provenance, 'list_sizes': [20, 20, 20],
               'union_count': 56, 'exclusive_regions': regions, 'triple_intersection': [],
               'checks': 'Source list sizes, unique genes, ranks, set partition and all expected overlaps verified.',
               'caveat': 'Mixed model configurations. GSE202642 and GSE151530 did not meet overall stability criteria. Network displacement has no expression direction.'}
    (OUT / 'intersection_summary.json').write_text(json.dumps(summary, indent=2) + '\n')

    fig, ax = plt.subplots(figsize=(11, 9))
    fig.subplots_adjust(left=.05, right=.95, top=.86, bottom=.24)
    centers = [(0, .65), (.95, .65), (.475, -.13)]
    colors = ['#327CA8', '#E6A23C', '#B66788']
    for center, color in zip(centers, colors):
        ax.add_patch(Circle(center, .85, facecolor=color, edgecolor=color, linewidth=1.8, alpha=.23))
    locations = {'100': (-.4, .86), '010': (1.35, .86), '001': (.475, -.65),
                 '110': (.475, 1.01), '101': (0, .02), '011': (.95, .02), '111': (.475, .40)}
    for mask, position in locations.items():
        ax.text(*position, str(len(regions[mask])), ha='center', va='center', fontsize=25,
                fontweight='bold', color='#25313D')
    ax.text(-.43, 1.68, LABELS[0] + '\n20 genes', ha='center', va='center', color=colors[0])
    ax.text(1.38, 1.68, LABELS[1] + '\n20 genes', ha='center', va='center', color='#9E6812')
    ax.text(.475, -1.12, LABELS[2] + '\n20 genes', ha='center', va='center', color='#8C4163')
    ax.set(xlim=(-1.08, 2.03), ylim=(-1.32, 1.90))
    ax.set_aspect('equal')
    ax.axis('off')
    fig.suptitle('PGAM5 virtual knockout: overlap of downstream TOP20 lists', fontsize=18, y=.96)
    fig.text(.5, .905, 'GSE202642 C7 | GSE149614 M1 | GSE151530 M2', ha='center', fontsize=13)
    fig.text(.09, .175, 'GSE202642 + GSE149614:  DNAJB1, HSP90AA1, HSPH1\n'
             'GSE149614 + GSE151530:  ANXA1\n'
             'GSE202642 + GSE151530:  none     |     All three:  none', fontsize=12, linespacing=1.6)
    fig.text(.09, .045, 'Numbers are gene counts; circle areas are schematic, not proportional.\n'
             'Mixed model configurations; GSE202642 and GSE151530 failed overall stability criteria.\n'
             'Exploratory network perturbations; no RNA up/down direction is implied.',
             fontsize=10, color='#56616D', linespacing=1.5)
    for suffix in ['png', 'pdf']:
        fig.savefig(OUT / f'01_TOP20_overlap_venn.{suffix}', dpi=220, facecolor='white')
    plt.close(fig)

    data = np.array([[ranks[d].get(g, np.nan) for d in names] for g in shared])
    fig, ax = plt.subplots(figsize=(9.5, 6.3))
    fig.subplots_adjust(left=.19, right=.85, bottom=.28, top=.80)
    cmap = plt.get_cmap('Blues_r').copy()
    cmap.set_bad('#EEEEEE')
    im = ax.imshow(np.ma.masked_invalid(data), cmap=cmap, vmin=1, vmax=20, aspect='auto')
    ax.set_yticks(range(len(shared)), labels=shared)
    ax.set_xticks(range(3), labels=['GSE202642\nC7 | v2 ensemble',
                                   'GSE149614\nM1 | sparse main',
                                   'GSE151530\nM2 | unfiltered\nsensitivity'])
    ax.tick_params(axis='both', length=0, pad=10)
    for i in range(len(shared)):
        for j in range(3):
            rank = data[i, j]
            ax.text(j, i, f'#{int(rank)}' if np.isfinite(rank) else 'Outside\nTOP20',
                    ha='center', va='center', fontsize=12,
                    color='white' if np.isfinite(rank) and rank <= 8 else '#27333F')
    ax.set_xticks(np.arange(-.5, 3, 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(shared), 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=2)
    ax.tick_params(which='minor', bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cb = fig.colorbar(im, ax=ax, pad=.04, ticks=[1, 5, 10, 15, 20])
    cb.set_label('Within-dataset TOP20 rank\n1 = largest median displacement')
    fig.suptitle('Shared TOP20 genes: within-dataset ranks', fontsize=18, y=.95)
    fig.text(.5, .87, 'Four genes occur in two lists; no gene occurs in all three.', ha='center', fontsize=12)
    fig.text(.08, .07, 'Grey: absent from that TOP20 list, not evidence of zero expression or no perturbation.\n'
             'Ranks compare relative positions only; model displacement magnitudes are not pooled.\n'
             'GSE202642 and GSE151530 did not meet overall stability criteria.', fontsize=10, color='#56616D')
    for suffix in ['png', 'pdf']:
        fig.savefig(OUT / f'02_shared_genes_rank_matrix.{suffix}', dpi=220, facecolor='white')
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 6.3))
    fig.subplots_adjust(left=.025, right=.975, top=.78, bottom=.26, wspace=.16)
    for ax, (i, j) in zip(axes, itertools.combinations(range(3), 2)):
        common = sorted(sets[i] & sets[j])
        for center, color in zip([(-.48, 0), (.48, 0)], [colors[i], colors[j]]):
            ax.add_patch(Circle(center, .85, facecolor=color, edgecolor=color, linewidth=1.8, alpha=.23))
        ax.text(-.86, 0, str(len(sets[i] - sets[j])), ha='center', va='center', fontsize=25, fontweight='bold')
        ax.text(.86, 0, str(len(sets[j] - sets[i])), ha='center', va='center', fontsize=25, fontweight='bold')
        ax.text(0, 0, str(len(common)), ha='center', va='center', fontsize=27, fontweight='bold')
        ax.text(-.64, .99, names[i] + '\n20 genes', ha='center', va='center', fontsize=12)
        ax.text(.64, .99, names[j] + '\n20 genes', ha='center', va='center', fontsize=12)
        ax.text(0, -1.12, f'Shared: {len(common)} genes\n' + ('\n'.join(common) if common else 'None'),
                ha='center', va='top', fontsize=12, linespacing=1.4)
        ax.set(xlim=(-1.47, 1.47), ylim=(-1.85, 1.28))
        ax.set_aspect('equal')
        ax.axis('off')
    fig.suptitle('PGAM5 virtual knockout TOP20: pairwise intersections', fontsize=19, y=.96)
    fig.text(.5, .87, 'GSE202642 C7: v2 ensemble  |  GSE149614 M1: sparse main  |  GSE151530 M2: unfiltered sensitivity',
             ha='center', fontsize=11)
    fig.text(.06, .07, 'Numbers are exact gene counts; circle areas are schematic, not proportional.\n'
             'Exploratory candidates from mixed model configurations; GSE202642 and GSE151530 failed overall stability criteria.\n'
             'No RNA up/down direction is implied.', fontsize=10, color='#56616D', linespacing=1.5)
    for suffix in ['png', 'pdf']:
        fig.savefig(OUT / f'03_pairwise_TOP20_intersections.{suffix}', dpi=220, facecolor='white')
    plt.close(fig)
    files = sorted(p for p in OUT.rglob('*') if p.is_file() and p.suffix != '.zip' and p.name != 'SHA256SUMS')
    (OUT / 'SHA256SUMS').write_text(''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(OUT)}\n' for p in files))
    with zipfile.ZipFile(OUT / 'PGAM5_virtual_KO_three_dataset_TOP20_overlap.zip', 'w', zipfile.ZIP_DEFLATED) as z:
        for path in files + [OUT / 'SHA256SUMS']:
            z.write(path, path.relative_to(OUT))
    with zipfile.ZipFile(OUT / 'PGAM5_virtual_KO_three_dataset_TOP20_overlap.zip') as z:
        assert z.testzip() is None
    print(json.dumps({'exclusive_counts': {m: len(g) for m, g in regions.items()},
                      'shared_genes': shared, 'union_count': 56, 'validation': 'passed'}, indent=2))


if __name__ == '__main__':
    main()
