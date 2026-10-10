"""Held-cohort pseudobulks with full non-target pools and two distinct units."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import sparse

R = Path(__file__).resolve().parent
CONFIG = json.loads((R/'config.json').read_text())


def main():
    genes = pd.read_csv(R/'common_genes.csv').gene.to_numpy()
    results, rows, coverage, membership_checks = [], [], [], []
    for cindex, cohort in enumerate(CONFIG['primary_cohorts']):
        obs = pd.read_csv(R/f'{cohort}_metadata.csv.gz')
        raw = sparse.load_npz(R/f'{cohort}_raw_common.npz').astype(float)
        total = obs.total_counts.to_numpy(float)
        norm = raw.multiply(10000/total[:, None]).tocsr()
        rng = np.random.default_rng(CONFIG['seed'] + 1000*cindex)
        for donor, d in obs.groupby('donor', sort=True):
            ix = d.index.to_numpy()
            mac = d.harmonized_group.eq('Macrophage').to_numpy()
            pg = d.PGAM5_counts.gt(0).to_numpy()
            pools = {'positive': ix[mac & pg], 'negative': ix[mac & ~pg],
                'other': ix[~mac], 'PGAM5_nonmac': ix[~mac & pg],
                'cycling_nonmac': ix[~mac & d.cycling_flag.to_numpy()]}
            eligible = len(pools['positive']) >= 5 and len(pools['negative']) >= 20 and len(pools['other']) >= 200
            coverage.append({'dataset': cohort, 'donor': donor, 'macrophages': int(mac.sum()),
                'PGAM5_detected': len(pools['positive']), 'nonmacrophages': len(pools['other']),
                'dose_response_eligible': eligible,
                'unit_kind': 'sample_proxy' if cohort == 'GSE202642' else 'donor_label'})
            cases = []
            if eligible:
                cases += [('standard', f, .25, 'other', rep) for f in CONFIG['mixtures']['target_fractions'] for rep in range(3)]
                cases += [('macrophage_fraction_'+str(m), f, m, 'other', rep)
                          for m in [.1, .5] for f in [0, .02, .05] for rep in range(2)]
            if len(pools['negative']) >= 20:
                for group in ['PGAM5_nonmac', 'cycling_nonmac']:
                    if len(pools[group]) >= 20:
                        cases += [('zero_'+group, 0, .25, group, rep) for rep in range(3)]
            for group, part in d.groupby('fine_group', sort=True):
                j = part.index.to_numpy()
                if group == 'TAM_PGAM5_detected' or len(j) < 20:
                    continue
                if obs.loc[j, 'harmonized_group'].eq('Macrophage').any():
                    assert obs.loc[j, 'PGAM5_counts'].eq(0).all()
                key = 'pure_'+group
                pools[key] = j
                cases += [(key, 0, 0, key, rep) for rep in range(2)]
            for scenario, f, m, other, rep in cases:
                if scenario.startswith('pure_'):
                    chosen = rng.choice(pools[other], 2000, replace=True)
                    selected_positive = np.array([], dtype=int)
                else:
                    npg, nneg = int(round(f*2000)), int(round((m-f)*2000))
                    selected_positive = rng.choice(pools['positive'], npg, replace=True) if npg else np.array([], dtype=int)
                    neg = rng.choice(pools['negative'], nneg, replace=True) if nneg else np.array([], dtype=int)
                    chosen = np.r_[selected_positive, neg, rng.choice(pools[other], 2000-npg-nneg, replace=True)]
                assert len(chosen) == 2000
                is_target = obs.loc[chosen, 'harmonized_group'].eq('Macrophage').to_numpy() & obs.loc[chosen, 'PGAM5_counts'].gt(0).to_numpy()
                cell_truth = is_target.mean()
                rna_truth = total[chosen][is_target].sum()/total[chosen].sum()
                is_macrophage = obs.loc[chosen, 'harmonized_group'].eq('Macrophage').to_numpy()
                assert np.isclose(cell_truth, f)
                if f == 0:
                    assert not is_target.any()
                equal = np.asarray(norm[chosen].mean(0)).ravel()
                pooled = np.asarray(raw[chosen].sum(0)).ravel()/total[chosen].sum()*10000
                membership_checks.append({'dataset': cohort, 'donor': donor, 'scenario': scenario,
                    'replicate': rep, 'sampled_cells': len(chosen), 'target_cells': int(is_target.sum()),
                    'raw_count_denominator': float(total[chosen].sum()), 'truth_cell_fraction': cell_truth,
                    'truth_RNA_contribution': rna_truth})
                for unit, expression, truth in [('equalized_cell_fraction', equal, cell_truth),
                                                 ('library_RNA_contribution', pooled, rna_truth)]:
                    mixture_id = len(rows)
                    rows.append({'mixture_id': mixture_id, 'dataset': cohort, 'donor': donor,
                        'scenario': scenario, 'nominal_target_cell_fraction': f,
                        'macrophage_cell_fraction': float(is_macrophage.mean()),
                        'macrophage_RNA_contribution': float(total[chosen][is_macrophage].sum()/total[chosen].sum()),
                        'replicate': rep, 'unit': unit, 'truth': truth})
                    results.append(expression.astype(np.float32))
            print('Mixtures', cohort, donor, len(cases), flush=True)
    pd.DataFrame(rows).to_csv(R/'mixture_index.csv', index=False)
    np.savez_compressed(R/'all_mixture_expression.npz', values=np.stack(results))
    pd.DataFrame(coverage).to_csv(R/'validation_donor_coverage.csv', index=False)
    pd.DataFrame(membership_checks).to_csv(R/'mixture_membership_checks.csv', index=False)
    print('TOTAL', len(rows), 'mixture/unit records', flush=True)


if __name__ == '__main__':
    main()
