"""Nested cohort holdouts evaluate feature/configuration selection without leakage."""
from pathlib import Path
import itertools
import json
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from reference_core import ReferenceBuilder, fit, T, N, COHORTS, CONFIG

R = Path(__file__).resolve().parent


def summary_frame(table, keys):
    results = []
    for values, d in table.groupby(keys, dropna=False, sort=True):
        values = values if isinstance(values, tuple) else (values,)
        zero = d[d.truth.eq(0)]
        standard = d[d.scenario.eq('standard')]
        r = pearsonr(standard.truth, standard.predicted).statistic if len(standard) and standard.predicted.nunique() > 1 and standard.truth.nunique() > 1 else np.nan
        results.append(dict(zip(keys, values)) | {'donor_labels': d.donor.nunique(), 'mixtures': len(d),
            'standard_MAE_pp': abs(standard.predicted-standard.truth).mean()*100 if len(standard) else np.nan,
            'standard_Pearson_r': r, 'all_zero_target_p95_percent': zero.predicted.quantile(.95)*100 if len(zero) else np.nan,
            'all_zero_target_max_percent': zero.predicted.max()*100 if len(zero) else np.nan,
            'scenario_MAE_pp': abs(d.predicted-d.truth).mean()*100})
    return pd.DataFrame(results)


def main():
    builder = ReferenceBuilder()
    index = pd.read_csv(R/'mixture_index.csv')
    expression = np.load(R/'all_mixture_expression.npz')['values']
    assert len(index) == len(expression)
    options = [(top, solver) for top in CONFIG['configurations']['broad_marker_top']
               for solver in CONFIG['configurations']['solvers']]
    selected, tuning, outer, sensitivities, audit_cases, gene_rows = [], [], [], [], [], []
    inner_cache = {}

    def predict(training, held, top, solver, phase, inner=False, pg=True, no_cycle=False, save_cases=False):
        query = index[index.dataset.eq(held)]
        if inner:
            query = query[query.replicate.eq(0)]
        predictions = []
        for unit, d in query.groupby('unit'):
            model = builder.model(training, top, unit, include_pgam5=pg, remove_cycle=no_cycle)
            target = model['groups'].index(T)
            macrophage_columns = [j for j, g in enumerate(model['groups']) if g in [T, N, 'Cycling_'+N]]
            for _, row in d.iterrows():
                y = expression[int(row.mixture_id)].astype(float)
                coefficients, residual, weights = fit(model, y, solver)
                predictions.append(row.to_dict() | {'phase': phase, 'training_cohorts': ','.join(training),
                    'top_markers': top, 'solver': solver, 'include_PGAM5': pg, 'remove_cycle': no_cycle,
                    'features': len(model['features']), 'state_markers': len(model['state_features']),
                    'reference_components': len(model['groups']), 'predicted': float(coefficients[target]),
                    'predicted_total_macrophage': float(coefficients[macrophage_columns].sum()),
                    'relative_residual': residual})
                if save_cases and row.replicate == 0 and row.scenario == 'standard' and row.nominal_target_cell_fraction in [0, .05]:
                    # A bounded, inspectable sample of actual held-out inputs.
                    case = len(audit_cases)
                    if sum(c['dataset'] == held and c['unit'] == unit for c in audit_cases) < 4:
                        np.savez_compressed(R/f'solver_audit_case_{case}.npz', A=model['matrix'],
                            b=y[model['features']]/model['scale'], coefficients=coefficients, weights=weights,
                            solver=solver, target_index=target, gene_indices=model['features'])
                        audit_cases.append({'case': case, 'mixture_id': int(row.mixture_id),
                            'dataset': held, 'unit': unit, 'solver': solver, 'truth': row.truth})
        return pd.DataFrame(predictions)

    def configuration_scores(training):
        records = []
        for held in training:
            inner_training = tuple(c for c in training if c != held)
            for top, solver in options:
                key = (inner_training, held, top, solver)
                if key not in inner_cache:
                    predictions = predict(inner_training, held, top, solver, 'inner_cohort_selection', inner=True)
                    scores = summary_frame(predictions, ['dataset', 'unit'])
                    score_values = scores.standard_MAE_pp/100 + 2*scores.all_zero_target_p95_percent/100 + .02*(1-scores.standard_Pearson_r.fillna(0).clip(lower=0))
                    inner_cache[key] = float(score_values.mean())
                records.append({'held_inner': held, 'training_inner': ','.join(inner_training),
                                'top': top, 'solver': solver, 'score': inner_cache[key]})
        table = pd.DataFrame(records)
        means = table.groupby(['top', 'solver'], sort=True).score.mean()
        best = means.idxmin()
        return int(best[0]), best[1], table

    for held in COHORTS:
        training = tuple(c for c in COHORTS if c != held)
        top, solver, table = configuration_scores(training)
        tuning.append(table.assign(outer_held=held))
        selected.append({'held_dataset': held, 'training': ','.join(training),
                         'top': top, 'solver': solver})
        print('OUTER CONFIGURATION', held, top, solver, flush=True)
        predictions = predict(training, held, top, solver, 'outer_nested', save_cases=True)
        outer.append(predictions)
        for pg, cycle in [(False, False), (True, True)]:
            sensitivities.append(predict(training, held, top, solver,
                'outer_definition_anchor_sensitivity' if not pg else 'outer_cycle_sensitivity',
                pg=pg, no_cycle=cycle))
        model = builder.model(training, top, 'equalized_cell_fraction')
        gene_rows.extend({'held_dataset': held, 'gene': builder.genes[k], 'role': 'state_candidate'}
                         for k in model['state_features'])
        pd.DataFrame(selected).to_csv(R/'outer_configuration_choices.csv', index=False)
        print('OUTER COMPLETED', held, len(predictions), flush=True)
    pd.concat(tuning, ignore_index=True).to_csv(R/'inner_configuration_scores.csv', index=False)
    final_top, final_solver, final_scores = configuration_scores(tuple(COHORTS))
    final_scores.to_csv(R/'final_training_configuration_scores.csv', index=False)
    print('FINAL LOCKED', final_top, final_solver, flush=True)
    state_evidence, choices, contrasts = builder.state_evidence(tuple(COHORTS))
    state_evidence.to_csv(R/'all_gene_candidate_evidence.csv.gz', index=False)
    state_evidence[state_evidence.support_cohorts.ge(3)].to_csv(R/'cross_cohort_recurrent_candidates.csv', index=False)
    state_evidence[state_evidence.stable_candidate].to_csv(R/'patient_consistent_state_candidates.csv', index=False)
    contrasts.to_csv(R/'per_donor_candidate_contrasts.csv.gz', index=True)
    pd.DataFrame(gene_rows).to_csv(R/'outer_fold_state_candidates.csv', index=False)
    # A role-labelled panel includes PGAM5 explicitly without counting it as independent evidence.
    macrophage_anchors = set('C1QA C1QB C1QC CSF1R CD68 TYROBP FCER1G LST1 AIF1 CD163 MSR1 MRC1 SPP1'.split())
    panel = [{'gene': 'PGAM5', 'role': 'RNA_definition_anchor', 'validated_PGAM5_specific_marker': False}]
    panel += [{'gene': builder.genes[k], 'role': 'patient_consistent_state_candidate',
               'validated_PGAM5_specific_marker': False} for k in choices]
    panel += [{'gene': g, 'role': 'macrophage_lineage_anchor', 'validated_PGAM5_specific_marker': False}
              for g in sorted(macrophage_anchors) if g in builder.genes]
    pd.DataFrame(panel).to_csv(R/'EXPLORATORY_signature_gene_roles.csv', index=False)
    locked = {'training_cohorts': COHORTS, 'top': final_top, 'solver': final_solver,
        'feature_selection': 'All-five training, selected by fixed entire-cohort inner CV; nested outer evaluation separately reported',
        'status': 'VALIDATION_REQUIRED; retrospective data, not new independent cohorts',
        'state_genes': builder.genes[choices].tolist(), 'units': {}}
    for unit in ['equalized_cell_fraction', 'library_RNA_contribution']:
        model = builder.model(tuple(COHORTS), final_top, unit)
        pd.DataFrame(model['ref'][model['features']], index=builder.genes[model['features']],
            columns=model['groups']).rename_axis('gene').to_csv(R/f'EXPLORATORY_reference_{unit}.tsv', sep='\t')
        pd.DataFrame({'gene': builder.genes[model['features']], 'scale': model['scale']}).to_csv(R/f'{unit}_row_scales.csv', index=False)
        pd.DataFrame(model['markers']).to_csv(R/f'{unit}_reference_marker_evidence.csv', index=False)
        pd.DataFrame(model['support']).to_csv(R/f'{unit}_reference_group_support.csv', index=False)
        locked['units'][unit] = {'features': builder.genes[model['features']].tolist(), 'groups': model['groups']}
    (R/'locked_reference.json').write_text(json.dumps(locked, indent=2)+'\n')
    out = pd.concat(outer, ignore_index=True)
    out.to_csv(R/'nested_outer_predictions.csv.gz', index=False)
    pd.concat(sensitivities, ignore_index=True).to_csv(R/'nested_sensitivity_predictions.csv.gz', index=False)
    pd.DataFrame(audit_cases).to_csv(R/'solver_audit_case_manifest.csv', index=False)
    summary = summary_frame(out, ['dataset', 'unit'])
    summary.to_csv(R/'nested_validation_summary.csv', index=False)
    summary_frame(out, ['dataset', 'unit', 'scenario']).to_csv(R/'validation_by_scenario.csv', index=False)
    summary_frame(pd.concat(sensitivities), ['phase', 'dataset', 'unit']).to_csv(R/'sensitivity_validation_summary.csv', index=False)
    coverage = pd.read_csv(R/'validation_donor_coverage.csv')
    gates = []
    for _, row in summary.iterrows():
        scope = {'dataset': row.dataset, 'unit': row.unit}
        for name, value, passed in [
            ('standard_MAE_pp<=1', row.standard_MAE_pp, row.standard_MAE_pp <= 1),
            ('standard_Pearson>=0.7', row.standard_Pearson_r, np.isfinite(row.standard_Pearson_r) and row.standard_Pearson_r >= .7),
            ('zero_target_p95_percent<=0.5', row.all_zero_target_p95_percent, row.all_zero_target_p95_percent <= .5)]:
            gates.append(scope | {'criterion': name, 'value': float(value), 'passed': bool(passed)})
        n = int(coverage[coverage.dataset.eq(row.dataset) & coverage.dose_response_eligible].donor.nunique())
        genuine = row.dataset != 'GSE202642'
        gates.append(scope | {'criterion': 'dose_response_independent_donors>=3', 'value': n,
                              'passed': bool(n >= 3 and genuine)})
    for label in ['Monocyte_enriched', 'DC_enriched']:
        test = out[out.scenario.eq('pure_'+label) & ~out.dataset.eq('GSE202642')]
        n = test[['dataset', 'donor']].drop_duplicates().shape[0]
        for unit in out.unit.unique():
            v = test[test.unit.eq(unit)].predicted.quantile(.95) * 100 if len(test) else np.nan
            gates.append({'dataset': 'all_outer', 'unit': unit, 'criterion': label+'_challenge_p95<=0.5_and_donors>=2',
                          'value': float(v), 'donors': n, 'passed': bool(n >= 2 and np.isfinite(v) and v <= .5)})
    gates.append({'dataset': 'all_training', 'unit': 'state_marker_evidence',
                  'criterion': 'patient_consistent_independent_state_markers>=5',
                  'value': len(choices), 'passed': len(choices) >= 5})
    g = pd.DataFrame(gates)
    g.to_csv(R/'validation_gate_checks.csv', index=False)
    result = {'computational_gates_all_passed': bool(g.passed.all()),
        'valid_for_TCGA_cellular_abundance': False,
        'bulk_platform_and_cell_RNA_calibration_completed': False,
        'final_solver': final_solver, 'final_top_markers': final_top,
        'patient_consistent_state_genes': builder.genes[choices].tolist(),
        'failed_gates': int((~g.passed).sum()), 'total_gates': len(g),
        'interpretation': 'Only an exploratory RNA-detected macrophage reference unless all gates and external platform/cell-unit calibration pass. Do not report TCGA coefficients as cellular infiltration percentages.'}
    (R/'validation_result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(summary.to_string(index=False), flush=True)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
