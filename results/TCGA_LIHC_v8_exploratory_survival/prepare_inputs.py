"""One-time source preparation; source paths are overrideable and not needed to refit."""
from pathlib import Path
import argparse, hashlib, json, shutil
import numpy as np
import pandas as pd

R = Path(__file__).resolve().parent
UNITS = ['library_RNA_contribution', 'equalized_cell_fraction']

def digest(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for c in iter(lambda: f.read(1024*1024), b''): h.update(c)
    return {'bytes': p.stat().st_size, 'sha256': h.hexdigest()}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference-root', type=Path, default=Path('/workspace/scratch/PGAM5_RNA_identity_v8'))
    parser.add_argument('--tcga-root', type=Path, default=Path('/workspace/scratch/PGAM5_signature_optimization_v5'))
    parser.add_argument('--bulk-original', type=Path, default=Path('/workspace/scratch/TCGA_LIHC_candidate30_survival/TCGA24_tumor_featurecounts.txt.gz'))
    args = parser.parse_args()
    assert not (R/'protocol.json').exists(), 'Preserve frozen inputs'
    refs = args.reference_root/'POOLED_ALL_TRAINING/continuous_rank'
    sources = {}
    for unit in UNITS:
        out = R/'reference'/unit
        out.mkdir(parents=True, exist_ok=True)
        for name in [f'{unit}_model.npz', f'EXPLORATORY_reference_{unit}.tsv',
                     f'{unit}_row_scales.csv', f'{unit}_definition.json',
                     f'{unit}_candidate_state_genes.csv', f'{unit}_group_support.csv']:
            p = refs/name
            sources[str(p)] = digest(p)
            shutil.copy2(p, out/name)
    for name in ['guard_solver.py', 'validation_result.json', 'protocol.json']:
        p = args.reference_root/name
        dest = R/('guard_solver.py' if name == 'guard_solver.py' else 'v8_source_'+name)
        shutil.copy2(p, dest)
        sources[str(p)] = digest(p)
    cachepath = args.tcga_root/'TCGA_primary_common_gene_counts.npz'
    cache = np.load(cachepath)
    counts = pd.DataFrame(cache['counts'], index=cache['genes'], columns=cache['patients'])
    assert counts.index.is_unique and counts.columns.is_unique and counts.shape[1] == 371
    assert np.isfinite(counts.to_numpy()).all() and (counts.to_numpy() >= 0).all()
    requested = sorted(set().union(*(set(np.load(refs/f'{u}_model.npz')['genes']) for u in UNITS)))
    present = [g for g in requested if g in counts.index]
    counts.loc[present].to_csv(R/'TCGA_reference_union_raw_counts.csv.gz')
    meta = pd.read_csv(args.tcga_root/'TCGA_primary_aliquot_selection.csv').set_index('patient').loc[counts.columns]
    meta.index.name = 'patient'
    assert meta.index.is_unique and meta['sample'].str[13:15].eq('01').all()
    assert (meta.full_library_total_counts > 0).all()
    meta.to_csv(R/'TCGA_primary_aliquot_selection.csv')
    cp = counts.loc[present].div(meta.full_library_total_counts, axis=1)*10000
    cp.to_csv(R/'TCGA_reference_union_CP10k.csv.gz')
    sources[str(cachepath)] = digest(cachepath)
    metadata = args.tcga_root/'TCGA_input_provenance.json'
    shutil.copy2(metadata, R/'bulk_source_provenance.json')
    for name in ['TCGA_CDR_mirror_survival.csv.gz', 'source_provenance.json']:
        p = args.tcga_root/'survival_exploratory'/name
        shutil.copy2(p, R/(name if name != 'source_provenance.json' else 'clinical_source_provenance.json'))
        sources[str(p)] = digest(p)
    protocol = {
        'reference': 'v8 POOLED_ALL_TRAINING/continuous_rank; no new feature or outcome-driven model selection',
        'primary_unit': UNITS[0], 'sensitivity_unit': UNITS[1],
        'unit_rationale': 'Library-RNA contribution is the primary comparator to bulk RNA; equalized cell-RNA reference is an uncalibrated sensitivity analysis, not a validated cell-fraction model.',
        'reference_size': '395 genes x 20 components in each original model; gene sets differ slightly between units',
        'predictor': 'sum of coefficients at the saved v8 target_mask (TAM_PGAM5_detected); total macrophages is not the primary predictor',
        'solver': 'Unchanged v8 robust nonnegative sum-to-one FCLS, ridge 0.001, 3 reweighting steps, target-only PGAM5 RNA upper budget; NNLS initialization, not a pure NNLS fit.',
        'query': 'GSE62944 FeatureCounts read counts divided by full measured-library total and multiplied by 10000; counts not TPM; no log transform or selected-gene denominator.',
        'missing_genes': 'Use exact measured gene-symbol intersection and retain original row scales. No zero filling and no ambiguous alias conversion. This measured-row projection modifies the available features and does not inherit validation.',
        'patient_selection': 'Reuse outcome-blind primary tumor type 01 selection: highest full-matrix assigned count total, then lexical sample barcode. One patient, one aliquot.',
        'grouping': 'For each unit separately, High > global median over all 371 estimates; Low <= median. Freeze before survival merging; same cutpoint for OS and PFI. No outcome-optimized cutpoint or removal of coefficient zeros.',
        'endpoints': ['OS', 'PFI'], 'PFS_available': False,
        'endpoint_rule': 'Use exact TCGA-CDR OS/OS.time and PFI/PFI.time. Finite time > 0 days and binary event 0 or 1. No imputation. PFI is not relabeled PFS.',
        'primary_tests': 'Unweighted two-sided logrank for OS and PFI; BH across these 2 primary-unit tests.',
        'sensitivity_tests': 'Both endpoints for the other reference unit; BH over its 2 tests separately. Retain all four tests; optional joint-four-test BH also reported.',
        'KM': 'Product-limit curves; Greenwood log-log 95% CI; censor marks; number at risk at start of each displayed month.',
        'effect_estimates': 'Descriptive unadjusted high/low Cox HR, Efron ties; rank-time proportional-hazards diagnostic. No adjusted/independent-prognosis or causal claim.',
        'authorization': 'Current user explicitly requests exploratory TCGA application, superseding archived v8 no-TCGA release gate. Existing validation failures remain recorded.',
        'interpretation': 'Exploratory reference coefficient-survival association only. Neither a coefficient nor prognostic significance establishes PGAM5 specificity or calibrated cellular abundance.',
        'grouping_frozen_before_survival_analysis': True,
    }
    (R/'protocol.json').write_text(json.dumps(protocol, ensure_ascii=False, indent=2)+'\n')
    sources[str(args.bulk_original)] = digest(args.bulk_original)
    (R/'source_input_hashes.json').write_text(json.dumps(sources, indent=2)+'\n')
    print('Frozen protocol and portable inputs: 371 unique primary tumors,', len(present), 'union reference genes.')

if __name__ == '__main__': main()
