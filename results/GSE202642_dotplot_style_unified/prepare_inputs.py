"""Extract existing expression evidence without changing any cell labels."""
from pathlib import Path
import hashlib
import json
import anndata as ad
import numpy as np
import pandas as pd

R = Path(__file__).resolve().parent
I = R / 'inputs'
I.mkdir(exist_ok=True)
sources = {}


def record(p):
    with p.open('rb') as f:
        sources[str(p.relative_to(R.parent))] = hashlib.file_digest(f, 'sha256').hexdigest()


def export(key, evidence, labels, panels):
    genes = sum(panels.values(), [])
    assert len(genes) == len(set(genes))
    evidence = evidence[evidence.gene.isin(genes)].copy()
    assert not evidence.duplicated(['cluster', 'gene']).any()
    assert len(evidence) == len(labels) * len(genes)
    evidence.to_csv(I / (key + '_values.csv'), index=False)
    (I / (key + '_spec.json')).write_text(json.dumps({
        'labels': labels, 'panels': panels, 'groups': len(labels), 'genes': len(genes)}, indent=2) + '\n')


P = R.parent / 'GSE202642_HCC_allcell_annotation'
defs = pd.read_csv(P / 'cluster_annotations.csv', dtype={'cluster': str})
e = pd.read_csv(P / 'cluster_marker_evidence.csv', dtype={'cluster': str})
e = e.rename(columns={'fraction': 'detection_fraction'})
panels = {
    'T_cell': ['CD3D', 'CD3E', 'TRAC', 'IL7R', 'CD8A', 'FOXP3'],
    'NK': ['NKG7', 'GNLY'], 'B_cell': ['MS4A1', 'CD79A'], 'Plasma': ['MZB1', 'JCHAIN'],
    'Macrophage': ['C1QA', 'C1QB', 'CSF1R', 'CD68', 'LST1'],
    'Monocyte': ['FCN1', 'S100A8', 'CD14'], 'cDC': ['CD1C', 'FCER1A', 'CLEC9A'],
    'Migratory_DC': ['LAMP3'], 'pDC': ['IL3RA', 'TCF4'],
    'Endothelial': ['PECAM1', 'VWF'], 'Fibroblast': ['COL1A1', 'DCN'],
    'Mural': ['RGS5', 'ACTA2'], 'Hepatocyte_like': ['ALB', 'APOA1', 'TTR'],
    'Epithelial': ['EPCAM', 'KRT19'], 'Mast': ['TPSAB1', 'KIT'],
    'Neutrophil': ['FCGR3B', 'CSF3R'], 'Cycling': ['MKI67', 'TOP2A']}
labels = {str(x.cluster): x.cluster_label for x in defs.itertuples()}
assert list(labels) == [str(i) for i in range(21)]
export('allcell_C0_C20', e, labels, panels)
for n in ['cluster_annotations.csv', 'cluster_marker_evidence.csv']:
    record(P / n)

P = R.parent / 'GSE202642_macrophage_reclustering'
a = ad.read_h5ad(P / 'macrophages_marker_subset.h5ad')
names = json.loads((P / 'names.json').read_text())
genes = 'C1QA CSF1R FOLR2 CD163 LYVE1 MRC1 TREM2 SPP1 APOE MARCO CD5L HLA-DRA CD74 IL1B CXCL8 FCN1 S100A8 MKI67 TOP2A ISG15 IFIT1 MT1G MT2A PGAM5'.split()
counts = a.layers['counts'][:, a.var_names.get_indexer(genes)]
assert a.var_names.get_indexer(genes).min() >= 0
log = counts.multiply((10000 / a.obs.total_counts.to_numpy())[:, None]).tocsr()
log.data = np.log1p(log.data)
rows = []
for c in sorted(names, key=lambda x: int(x[1:])):
    mask = a.obs.mac_cluster.astype(str).eq(c).to_numpy()
    frac = np.asarray((counts[mask] > 0).mean(axis=0)).ravel()
    mean = np.asarray(log[mask].mean(axis=0)).ravel()
    for j, g in enumerate(genes):
        rows.append({'cluster': c, 'gene': g, 'cells': int(mask.sum()),
                     'detection_fraction': float(frac[j]), 'mean_log1p_CP10k': float(mean[j])})
export('macrophage_M0_M8', pd.DataFrame(rows), {c: c+' '+names[c] for c in names}, {'Markers': genes})
record(P / 'macrophages_marker_subset.h5ad')
record(P / 'names.json')

P = R.parent / 'GSE202642_PGAM5_subpopulation_assessment'
names = pd.read_csv(P / 'functional_subgroup_names.csv', dtype={'primary_cluster': str})
e = pd.read_csv(P / 'functional_subgroup_marker_dotplot_values.csv')
e['cluster'] = e.original_cluster_id.str.removeprefix('C')
genes = e.sort_values('x').gene.drop_duplicates().tolist()
labels = {x.primary_cluster: x.original_cluster_id+' '+x.functional_name_EN.replace(' — ', ' | ')
          for x in names.itertuples()}
export('macrophage_C0_C14_legacy', e[['cluster', 'gene', 'cells', 'detection_fraction', 'mean_log1p_CP10k']],
       labels, {'Macrophage_identity': genes[:6], 'PGAM5': genes[6:7], 'Naming_markers': genes[7:]})
for n in ['functional_subgroup_names.csv', 'functional_subgroup_marker_dotplot_values.csv']:
    record(P / n)
(R / 'source_provenance.json').write_text(json.dumps({'source_sha256': sources,
    'reclustering': False, 'all_original_gene_sets_retained': True}, indent=2) + '\n')
print('Prepared 21 all-cell groups, 9 re-clustered macrophage groups and 15 legacy macrophage groups.')
