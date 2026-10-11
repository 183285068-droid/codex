"""Independent source reconstruction and numerical verification of the published analysis."""
from pathlib import Path
import csv, gzip, io, json, hashlib, re
import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear
from scipy.stats import chi2

R = Path(__file__).resolve().parent
I = R / 'inputs'
checks = []


def check(name, value):
    checks.append({'check': name, 'passed': bool(value)})
    assert value, name


def source(accession):
    text = gzip.open(I / f'{accession}_series_matrix.txt.gz', 'rt').read()
    rows = [next(csv.reader([line], delimiter='\t')) for line in text.splitlines()
            if line.startswith('!Sample')]
    gsms = next(row[1:] for row in rows if row[0] == '!Sample_geo_accession')
    titles = next(row[1:] for row in rows if row[0] == '!Sample_title')
    frame = pd.DataFrame({'GSM': gsms, 'source_title': titles})
    frame['patient_id'] = [re.search(r'\((?:bkg_|hcc_)(\d+)\)', t).group(1) for t in titles]
    for row in rows:
        if row[0] == '!Sample_characteristics_ch1':
            field = next(v for v in row[1:] if ':' in v).rsplit(':', 1)[0].strip()
            frame[field] = [v.rsplit(':', 1)[-1].strip() for v in row[1:]]
    matrix = pd.read_csv(io.StringIO(text.split('!series_matrix_table_begin')[1]
                         .split('!series_matrix_table_end')[0].strip()), sep='\t', index_col=0)
    return matrix, frame


def main():
    prov = json.loads((R / 'source_provenance.json').read_text())
    for name, expected in prov['input_sha256'].items():
        check(f'Immutable input SHA256: {name}', hashlib.sha256((I / name).read_bytes()).hexdigest() == expected)
    for name, expected in prov.get('compressed_SOFT_payload_sha256', {}).items():
        check(f'Original source SOFT byte payload: {name}',
              hashlib.sha256(gzip.decompress((I / name).read_bytes())).hexdigest() == expected)
    raw, tumor = source('GSE10141')
    _, background = source('GSE10140')
    patients = pd.read_csv(R / 'patient_coefficients_clinical.csv', dtype={'patient_id': str}).set_index('GSM')
    check('Exactly 80 unique tumor patients', len(patients) == 80 and patients.patient_id.is_unique)
    check('Only tumor expression analyzed', set(patients.index) == set(raw.columns) and
          not set(patients.index).intersection(background.GSM))
    clinical = pd.read_excel(io.BytesIO(gzip.decompress((I / 'clinical.xls.gz').read_bytes()))).set_index('!Sample_geo_accession')
    clinical = clinical.loc[patients.index]
    check('All OS times equal official XLS', np.array_equal(patients.OS_time_days, clinical['Survival time (days)']))
    check('All OS events equal official XLS', np.array_equal(patients.OS_event, clinical.iloc[:, 2]))
    check('Patient codes agree independently with clinical titles',
          np.array_equal(patients.patient_id, clinical['!Sample_title'].str.extract(r'\((\d+)\)')[0]))
    background = background.set_index('patient_id')
    tumor = tumor.set_index('GSM').loc[patients.index]
    check('All tumor OS times also equal GEO sample metadata',
          np.array_equal(patients.OS_time_days, tumor['survival time (days)'].astype(float)))
    check('All tumor OS events also equal GEO sample metadata',
          np.array_equal(patients.OS_event, tumor['survival status (0: alive or censored, 1: dead)'].astype(float)))
    linked = patients[patients.patient_id.isin(background.index)]
    bg = background.loc[linked.patient_id]
    check('Exactly 62 same-patient recurrence matches', len(linked) == 62)
    check('Every linked GSM is the same-patient source GSM', np.array_equal(linked.adjacent_GSM, bg.GSM))
    check('Every paired OS time agrees', np.array_equal(linked.OS_time_days, bg['survival time (days)'].astype(float)))
    check('Every paired OS event agrees', np.array_equal(linked.OS_event, bg['survival_status (0: alive_or_censored, 1:dead)'].astype(float)))
    check('Every recurrence time equals same-patient source metadata',
          np.array_equal(linked.TTR_time_days, bg['time to hcc recurrence (days)'].astype(float)))
    check('Every recurrence event equals same-patient source metadata',
          np.array_equal(linked.TTR_event, bg['recurrence status (0: no recurrence of censored, 1: recurrence)'].astype(float)))
    check('Unmatched tumors have no imputed recurrence time/event',
          patients.loc[~patients.patient_id.isin(background.index), ['TTR_time_days', 'TTR_event']].isna().all().all())

    annotation = pd.read_csv(I / 'GPL5474_annotation.tsv.gz', sep='\t', dtype=str)
    mapping = pd.read_csv(R / 'used_probe_mapping.csv', dtype=str)
    check('Mappings are official platform symbols', mapping.set_index('ID').gene.equals(
          annotation.set_index('ID').loc[mapping.ID, 'Symbol'].rename('gene')))
    check('No PGAM5 symbol, Entrez192111 or RefSeq NM138575 in supplied platform',
          not annotation.apply(lambda c: c.fillna('').str.contains(r'PGAM5|192111|NM_138575')).any(axis=1).any())
    mapped = raw.loc[mapping.ID].copy()
    mapped.index = mapping.gene.to_numpy()
    gene = mapped.groupby(level=0, sort=True).median()
    rebuilt = gene.div(gene.sum(axis=0), axis=1) * 10000
    saved = pd.read_csv(R / 'array_query_normalized_intensity.tsv.gz', sep='\t', index_col=0)
    check('All 33,280 query measurements reconstruct from linear source intensities',
          saved.shape == (416, 80) and np.allclose(saved, rebuilt.loc[saved.index, saved.columns], rtol=1e-12, atol=1e-12))
    ref = pd.read_csv(I / 'reference_850x17.tsv', sep='\t', index_col=0)
    scales = pd.read_csv(I / 'reference_row_scales.csv', index_col=0).row_scale
    used_ref = pd.read_csv(R / 'used_reference.tsv', sep='\t', index_col=0)
    used_scale = pd.read_csv(R / 'used_reference_row_scales.csv', index_col=0).row_scale
    check('Frozen 17 reference components unchanged', used_ref.columns.equals(ref.columns) and len(used_ref.columns) == 17)
    check('Measured reference rows are unmodified', np.allclose(used_ref, ref.loc[saved.index], atol=0, rtol=1e-14))
    check('Original feature weights unchanged', np.allclose(used_scale, scales.loc[saved.index], atol=0, rtol=1e-14))
    coverage = pd.read_csv(R / 'reference_gene_coverage.csv').set_index('gene').measured
    check('Feature coverage is exact intersection, no zero imputation', set(coverage[coverage].index) == set(ref.index).intersection(gene.index))
    A = ref.loc[saved.index].to_numpy() / scales.loc[saved.index].to_numpy()[:, None]
    check('Weighted reference full rank 17', np.linalg.matrix_rank(A) == 17)
    all_coeff = pd.read_csv(R / 'all_component_coefficients.csv', index_col=0)
    check('All 80 nonnegative coefficient sums equal one', (all_coeff.to_numpy() >= 0).all() and np.allclose(all_coeff.sum(axis=1), 1))
    check('Target coefficients match exported full reference solution', np.allclose(
          patients.target_coefficient, all_coeff.loc[patients.index, 'TAM_PGAM5_detected'], atol=1e-14))
    selection = np.unique(np.linspace(0, 79, 24, dtype=int))
    numerical = []
    for index in selection:
        sample = saved.columns[index]
        y = rebuilt.loc[saved.index, sample].to_numpy() / scales.loc[saved.index].to_numpy()
        fit = lsq_linear(A, y, bounds=(0, np.inf), method='bvls', tol=1e-12, max_iter=10000)
        estimate = fit.x / fit.x.sum()
        error = float(np.max(np.abs(estimate - all_coeff.loc[sample].to_numpy())))
        check(f'Independent BVLS 17-component solution: {sample}', fit.success and error < 1e-8)
        numerical.append({'GSM': sample, 'max_component_difference': error, 'independent_solver_success': fit.success})
    pd.DataFrame(numerical).to_csv(R / 'independent_bounded_solver_checks.csv', index=False, lineterminator='\n')
    cutoff = float(patients.target_coefficient.median())
    check('Outcome-blind median group assignments across all 80 tumors', np.array_equal(
          patients.group, np.where(patients.target_coefficient.gt(cutoff), 'High', 'Low')))
    check('Exactly 40 High and 40 Low tumors', patients.group.value_counts().to_dict() == {'High': 40, 'Low': 40})
    stats = pd.read_csv(R / 'survival_statistics.csv').set_index('endpoint')
    curves = pd.read_csv(R / 'KM_curve_values.csv')
    risk = pd.read_csv(R / 'KM_numbers_at_risk.csv')
    for endpoint, count in [('OS', 80), ('TTR', 62)]:
        data = pd.read_csv(R / f'{endpoint}_analysis_patients.csv', dtype={'patient_id': str}).set_index('GSM')
        t, e = patients[f'{endpoint}_time_days'], patients[f'{endpoint}_event']
        valid = t.gt(0) & np.isfinite(t) & e.isin([0, 1])
        if endpoint == 'TTR':
            valid &= ~(patients.OS_event.eq(1) & t.gt(patients.OS_time_days))
        check(f'{endpoint} cohort selected solely by valid source follow-up', len(data) == count and set(data.index) == set(patients.index[valid]))
        check(f'{endpoint} times/events unchanged', np.array_equal(data.time_days, patients.loc[data.index, f'{endpoint}_time_days']) and
              np.array_equal(data.event, patients.loc[data.index, f'{endpoint}_event']))
        check(f'{endpoint} grouping remains the all-80 fixed median', np.array_equal(data.group, patients.loc[data.index].group))
        high, low = data[data.group.eq('High')], data[data.group.eq('Low')]
        oe, variance = 0., 0.
        for time in sorted(data.loc[data.event.eq(1), 'time_days'].unique()):
            nh, nl = high.time_days.ge(time).sum(), low.time_days.ge(time).sum()
            dh = (high.time_days.eq(time) & high.event.eq(1)).sum()
            dl = (low.time_days.eq(time) & low.event.eq(1)).sum()
            n, d = int(nh + nl), int(dh + dl)
            oe += dh - d * nh / n
            if n > 1:
                variance += nh * nl * d * (n - d) / (n * n * (n - 1))
        check(f'{endpoint} logrank independently reconstructs from patient risk sets',
              np.isclose(chi2.sf(oe**2 / variance, 1), stats.loc[endpoint, 'logrank_p'], atol=1e-12))
        check(f'{endpoint} events and group sizes correctly reported', stats.loc[endpoint, 'events'] == data.event.sum() and
              stats.loc[endpoint, 'high_patients'] == len(high) and stats.loc[endpoint, 'low_patients'] == len(low))
        for group, q in [('High', high), ('Low', low)]:
            c = curves[curves.endpoint.eq(endpoint) & curves.group.eq(group)]
            product = 1.
            previous = -np.inf
            for row in c.itertuples():
                for time in sorted(q.loc[q.time_days.gt(previous) & q.time_days.le(row.time_days), 'time_days'].unique()):
                    product *= 1 - int((q.time_days.eq(time) & q.event.eq(1)).sum()) / int(q.time_days.ge(time).sum())
                assert np.isclose(product, row.survival, atol=1e-12)
                previous = row.time_days
            check(f'{endpoint} {group} KM values independently reconstruct', True)
            table = risk[risk.endpoint.eq(endpoint) & risk.group.eq(group)]
            check(f'{endpoint} {group} numbers at risk match source times', all(
                  row.at_risk == int(q.time_days.ge(row.time_days).sum()) for row in table.itertuples()))
    ps = stats.logrank_p.to_numpy()
    order = np.argsort(ps)
    expected_q = np.minimum.accumulate((ps[order] * 2 / np.arange(1, 3))[::-1])[::-1].clip(max=1)
    check('BH correction across the two endpoint logrank tests', np.allclose(stats.logrank_BH_q.to_numpy()[order], expected_q))
    check('KM PNG and PDF both exist', all((R / f'KM_OS_TTR.{ext}').stat().st_size > 10000 for ext in ['png', 'pdf']))
    summary = {'checks_passed': len(checks), 'checks_failed': 0, 'independent_solvers': len(numerical),
               'independent_solver_max_difference': max(x['max_component_difference'] for x in numerical),
               'rebuilt_query_values': int(saved.size), 'checks': checks,
               'limit': 'These checks verify computation and source joins, not biological identity or cell-abundance validity.'}
    (R / 'audit_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({key: value for key, value in summary.items() if key != 'checks'}, indent=2))


if __name__ == '__main__':
    main()
