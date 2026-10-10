"""Give saved RNA clusters descriptive marker/program names, with ambiguity flags.

Names describe observed RNA patterns; they do not assert a validated lineage,
mechanism, pure population, or functional assay result. No cells are removed.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parent

# Each pair occurs among that cluster's top 25 descriptive markers.
NAMES = [
    ('0', 'SPP1', 'HLA-DRB5', 'MHC-II / SPP1-high',
     'SPP1/HLA-DRB5高表达、抗原呈递相关巨噬细胞（伴肝源/肿瘤RNA）',
     'MHC-II-associated and SPP1-high RNA pattern',
     'HLA genotype/patient effects and hepatic/tumor RNA background may contribute.'),
    ('1', 'NUPR1', 'AMBP', 'Low-complexity / mixed',
     'NUPR1/AMBP高表达、低复杂度混合RNA群（待复核）',
     'NUPR1-associated stress RNA with hepatic RNA and low library complexity',
     'Not a confirmed functional macrophage subtype; hepatic/tumor RNA background.'),
    ('2', 'DNASE1L3', 'SPIC', 'Resident-like / mixed',
     'DNASE1L3/SPIC高表达、驻留样巨噬细胞（伴肝源RNA）',
     'Resident/scavenging-associated RNA; resident identity remains tentative',
     'FABP1/ALDOB/FGL1 hepatic RNA background; resident-like is an expression analogy.'),
    ('3', 'CD5L', 'HLA-DQA2', 'Antigen presentation',
     'CD5L/HLA-DQA2高表达、抗原呈递相关巨噬细胞',
     'CD5L and MHC-II-associated RNA',
     'Patient/HLA effects possible; gene expression does not measure antigen presentation activity.'),
    ('4', 'SPP1', 'MMP12', 'Matrix remodeling',
     'SPP1/MMP12高表达、基质重塑相关巨噬细胞',
     'Matrix-remodeling-associated RNA with MMP9/OLR1/IL1RN',
     'Patient-dominated cluster; remodeling is inferred from RNA markers.'),
    ('5', 'COL4A1', 'CDH5', 'Endothelial-RNA mixed',
     'COL4A1/CDH5高表达、内皮RNA混合群（待复核）',
     'Endothelial/vascular RNA including FLT1/KDR/NOTCH4',
     'Mixed-lineage RNA; doublet, ambient RNA or cell identity require review, not inferred here.'),
    ('6', 'LPL', 'STAB1', 'Lipid handling',
     'LPL/STAB1高表达、脂质处理相关巨噬细胞',
     'Lipid-handling/scavenging-associated RNA with FCGBP/VSIG4',
     'Patient-dominated; does not establish lipid flux or an exclusive subtype.'),
    ('7', 'MKI67', 'TOP2A', 'Cycling',
     'MKI67/TOP2A高表达、增殖型巨噬细胞',
     'Cell-cycle RNA including CENPF/PCLAF/CDK1/TPX2',
     'PGAM5-enriched RNA candidate; only a subset has detectable PGAM5.'),
    ('8', 'NKG7', 'TRBC2', 'T/NK-RNA mixed',
     'NKG7/TRBC2高表达、T/NK RNA混合群（待复核）',
     'T/NK-associated RNA including CD2/TRBC1/CD3G/CCL5',
     'Not named a cytotoxic macrophage subtype; lymphocyte admixture/doublets not resolved.'),
    ('9', 'CCL2', 'CCL8', 'Chemokine-associated',
     'CCL2/CCL8高表达、趋化因子相关巨噬细胞',
     'Chemokine-associated RNA with GAS6/KCNJ5',
     'Patient-dominated; chemokine secretion/recruitment not directly measured.'),
    ('10', 'S100A8', 'FCGR3B', 'Granulocyte-RNA mixed',
     'S100A8/FCGR3B高表达、粒细胞RNA混合群（待复核）',
     'Granulocyte-associated RNA including S100A9/CSF3R/CXCR2',
     'Only 34 cells; do not equate this with a confirmed inflammatory macrophage subtype.'),
    ('11', 'MT1G', 'MT2A', 'Metal response',
     'MT1G/MT2A高表达、金属应答相关巨噬细胞',
     'Metallothionein/metal-response RNA with MT1H/MT1M/MT1X',
     'Stress/handling and patient effects possible; no metal-response assay performed.'),
    ('12', 'CXCL2', 'NR4A2', 'Early inflammatory response',
     'CXCL2/NR4A2高表达、早期炎症应答相关巨噬细胞',
     'Inflammatory/immediate-response RNA with DUSP2/TNFAIP3/ATF3',
     'Dissociation/immediate stress may contribute; not evidence of a permanently stable lineage.'),
    ('13', 'IL1B', 'CXCL8', 'Inflammatory',
     'IL1B/CXCL8高表达、炎症相关巨噬细胞',
     'Inflammatory RNA with NFKBIZ/CLEC4E/EGR1',
     'Patient-dominated and early-response RNA present; not an M1 designation.'),
    ('14', 'PCK1', 'CYP3A5', 'Hepatocyte-RNA mixed',
     'PCK1/CYP3A5高表达、肝细胞RNA混合群（待复核）',
     'Hepatocyte-associated RNA with ACSM2A/ACSM2B/ONECUT1',
     'Only 38 cells; do not interpret hepatic metabolic genes as confirmed macrophage function.'),
]


def main():
    obs = pd.read_csv(ROOT / 'cell_annotations.csv.gz', index_col=0,
                      dtype={'primary_cluster': str})
    markers = pd.read_csv(ROOT / 'cluster_descriptive_markers.csv',
                          dtype={'cluster': str})
    markers = markers[markers.geometry.eq('primary')]
    genes = pd.read_csv(ROOT / 'portable_feature_genes.csv').gene.to_numpy()
    raw = sparse.load_npz(ROOT / 'portable_raw_feature_counts.npz').tocsr()
    assert raw.shape == (len(obs), len(genes))
    norm = raw.multiply((1e4 / obs.total_counts.to_numpy())[:, None]).tocsr()
    log = norm.copy()
    log.data = np.log1p(log.data)
    rows, evidence = [], []
    for cluster, g1, g2, function, name_cn, program, caveat in NAMES:
        sel = obs.primary_cluster.eq(cluster).to_numpy()
        assert sel.any()
        patient = obs.loc[sel, 'sample_name'].value_counts()
        row = {'primary_cluster': cluster, 'original_cluster_id': 'C' + cluster,
               'marker_gene_1': g1, 'marker_gene_2': g2,
               'functional_name_EN': g1 + '/' + g2 + ' — ' + function,
               'functional_name_CN': name_cn, 'plot_label': g1 + '/' + g2 + '\n' + function,
               'observed_RNA_program': program, 'interpretation_note': caveat,
               'annotation_status': 'Provisional descriptive RNA annotation',
               'mixed_RNA_flag': cluster in {'0', '1', '2', '5', '8', '10', '14'},
               'cells': int(sel.sum()), 'patients_with_cells': len(patient),
               'largest_patient_fraction': float(patient.max() / sel.sum()),
               'median_n_genes': float(obs.loc[sel, 'n_genes'].median()),
               'PGAM5_detected_cells': int(obs.loc[sel, 'PGAM5_detected'].sum())}
        for gene in (g1, g2):
            j = int(np.flatnonzero(genes == gene)[0])
            m = markers[markers.cluster.eq(cluster) & markers.gene.eq(gene)].iloc[0]
            mean = float(norm[sel, j].mean())
            other = float(norm[~sel, j].mean())
            detected = float((raw[sel, j] > 0).mean())
            assert m['rank'] <= 25 and m.descriptive_log2FC_CP10k_pseudocount0p05 > 0
            assert np.isclose(mean, m.cluster_mean_CP10k, rtol=1e-6, atol=1e-7)
            assert np.isclose(detected, m.cluster_detection_fraction, atol=1e-12)
            evidence.append({'original_cluster_id': 'C' + cluster,
                             'functional_name_EN': row['functional_name_EN'], 'gene': gene,
                             'descriptive_marker_rank': int(m['rank']),
                             'cluster_mean_CP10k': mean, 'other_mean_CP10k': other,
                             'cluster_detection_fraction': detected,
                             'mean_log1p_CP10k': float(log[sel, j].mean()),
                             'descriptive_log2FC_CP10k_pseudocount0p05': float(np.log2((mean + .05) / (other + .05)))})
        rows.append(row)
    names = pd.DataFrame(rows)
    assert names.primary_cluster.is_unique and names.functional_name_EN.is_unique
    assert set(names.primary_cluster) == set(obs.primary_cluster)
    assert names.cells.sum() == len(obs)
    names.to_csv(ROOT / 'functional_subgroup_names.csv', index=False)
    pd.DataFrame(evidence).to_csv(ROOT / 'functional_subgroup_marker_evidence.csv', index=False)

    annotation = pd.read_csv(ROOT / 'RNA_subgroup_cell_annotation.csv.gz')
    mapping = names.set_index('original_cluster_id')
    for column in ['functional_name_EN', 'functional_name_CN', 'annotation_status', 'mixed_RNA_flag']:
        annotation[column] = annotation.PGAM5_free_cluster.map(mapping[column])
        assert annotation[column].notna().all()
    assert annotation.cell_id.is_unique and len(annotation) == len(obs)
    annotation.to_csv(ROOT / 'functional_subgroup_cell_annotations.csv.gz', index=False)
    result = {'named_RNA_clusters': len(names), 'annotated_cells': len(annotation),
              'marker_evidence_rows': len(evidence),
              'all_naming_markers_present_and_in_top25': True,
              'raw_marker_means_and_detection_match_existing_marker_table': True,
              'all_cells_have_one_provisional_name': True, 'cells_removed': 0,
              'mixed_RNA_flagged_clusters': names.loc[names.mixed_RNA_flag, 'original_cluster_id'].tolist(),
              'functional_subgroup_names_SHA256': hashlib.sha256((ROOT / 'functional_subgroup_names.csv').read_bytes()).hexdigest(),
              'interpretation': 'Within-cohort RNA descriptions. Mixed RNA is observed, but its source is not established; no functional assay, new clustering, or abundance validation.',
              'checks_passed': True}
    (ROOT / 'functional_subgroup_naming_validation.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
