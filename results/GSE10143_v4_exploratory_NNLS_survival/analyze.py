"""Frozen v4 NNLS projection and patient-linked survival analysis; run beside inputs/."""
from pathlib import Path
import csv, gzip, io, json, re, platform, importlib.metadata as md
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, proportional_hazard_test
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R = Path(__file__).resolve().parent
I = R / 'inputs'
TARGET = 'TAM_PGAM5_detected'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'pdf.fonttype': 42})


def write(df, name):
    df.to_csv(R / name, index=False, lineterminator='\n')


def load_series(accession):
    text = gzip.open(I / f'{accession}_series_matrix.txt.gz', 'rt').read()
    matrix = pd.read_csv(io.StringIO(text.split('!series_matrix_table_begin')[1]
                         .split('!series_matrix_table_end')[0].strip()), sep='\t', index_col=0)
    rows = [next(csv.reader([line], delimiter='\t')) for line in text.splitlines()
            if line.startswith('!Sample')]
    field = {row[0]: row[1:] for row in rows if row[0] != '!Sample_characteristics_ch1'}
    meta = pd.DataFrame({'GSM': field['!Sample_geo_accession'],
                         'sample_title': field['!Sample_title'],
                         'tissue': field['!Sample_source_name_ch1'],
                         'platform': field['!Sample_platform_id']})
    for row in rows:
        if row[0] == '!Sample_characteristics_ch1':
            # Descriptions contain colons (e.g. 0: alive); the value follows the last colon.
            key = next(value for value in row[1:] if ':' in value).rsplit(':', 1)[0].strip()
            meta[key] = [value.rsplit(':', 1)[1].strip() if ':' in value else '' for value in row[1:]]
    meta['patient_id'] = meta.sample_title.str.extract(r'\((?:bkg_|hcc_)(\d+)\)')[0]
    assert meta.patient_id.notna().all() and meta.patient_id.is_unique
    assert meta.GSM.is_unique and matrix.index.is_unique and matrix.columns.is_unique
    assert matrix.columns.tolist() == meta.GSM.tolist()
    assert meta.platform.eq('GPL5474').all() and np.isfinite(matrix.to_numpy()).all()
    return matrix, meta


def project(raw):
    annotation = pd.read_csv(I / 'GPL5474_annotation.tsv.gz', sep='\t', dtype=str)
    annotation['gene'] = annotation.Symbol.fillna('').str.strip()
    usable = annotation[annotation.ID.isin(raw.index) & annotation.gene.ne('') &
                        annotation.gene.ne('---') &
                        ~annotation.gene.str.contains(r'///|;|\|', regex=True)].copy()
    assert usable.ID.is_unique and raw.min().min() > 0
    gene = raw.loc[usable.ID].copy()
    gene.index = usable.gene.to_numpy()
    gene = gene.groupby(level=0, sort=True).median()
    query = gene.div(gene.sum(axis=0), axis=1) * 1e4
    reference = pd.read_csv(I / 'reference_850x17.tsv', sep='\t', index_col=0)
    scales = pd.read_csv(I / 'reference_row_scales.csv', index_col=0).row_scale
    assert reference.shape == (850, 17)
    genes = [g for g in reference.index if g in query.index]
    A = reference.loc[genes].to_numpy() / scales.loc[genes].to_numpy()[:, None]
    assert np.linalg.matrix_rank(A) == 17, 'Target coefficients not uniquely identifiable in measured reference.'
    coefficients, diagnostics = [], []
    for sample in raw.columns:
        y = query.loc[genes, sample].to_numpy() / scales.loc[genes].to_numpy()
        beta, residual = nnls(A, y, maxiter=10000)
        assert beta.sum() > 0
        gradient = A.T @ (A @ beta - y)
        scale = max(1., np.linalg.norm(A) * np.linalg.norm(y))
        assert gradient.min() / scale > -1e-9
        assert np.max(np.abs(gradient * beta)) / max(1., np.linalg.norm(y)**2) < 1e-9
        coefficients.append(beta / beta.sum())
        diagnostics.append({'GSM': sample, 'relative_fit_residual': residual / np.linalg.norm(y),
                            'raw_beta_sum': beta.sum(), 'KKT_min_gradient_scaled': gradient.min() / scale,
                            'normalization_denominator': gene[sample].sum(),
                            'positive_components': int((beta > 1e-12).sum())})
    coefficients = pd.DataFrame(coefficients, index=raw.columns, columns=reference.columns)
    write(coefficients.rename_axis('GSM').reset_index(), 'all_component_coefficients.csv')
    write(pd.DataFrame(diagnostics), 'fit_diagnostics.csv')
    write(usable[['ID', 'gene', 'GeneID', 'REFSEQ']], 'used_probe_mapping.csv')
    write(pd.DataFrame({'gene': reference.index, 'measured': reference.index.isin(genes)}),
          'reference_gene_coverage.csv')
    reference.loc[genes].to_csv(R / 'used_reference.tsv', sep='\t', lineterminator='\n')
    scales.loc[genes].to_csv(R / 'used_reference_row_scales.csv', lineterminator='\n')
    query.loc[genes].to_csv(R / 'array_query_normalized_intensity.tsv.gz', sep='\t', lineterminator='\n')
    summary = {'tumor_samples': len(raw.columns), 'original_reference_genes': 850,
               'measured_reference_genes': len(genes), 'coverage_fraction': len(genes) / 850,
               'components': 17, 'PGAM5_measured': 'PGAM5' in gene.index,
               'PGAM5_identity_in_platform': bool(annotation.apply(
                   lambda col: col.fillna('').str.contains(r'PGAM5|192111|NM_138575')).any(axis=1).any()),
               'mapped_probes': len(usable), 'mapped_genes': len(gene),
               'raw_intensity_min': float(raw.min().min()), 'raw_intensity_max': float(raw.max().max()),
               'weighted_reference_rank': int(np.linalg.matrix_rank(A)),
               'weighted_reference_condition': float(np.linalg.cond(A)),
               'median_relative_fit_residual': float(np.median([x['relative_fit_residual'] for x in diagnostics]))}
    return coefficients, summary


def join_clinical(tumor, adjacent):
    clinical = pd.read_excel(io.BytesIO(gzip.decompress((I / 'clinical.xls.gz').read_bytes())))
    clinical = clinical.rename(columns={'!Sample_geo_accession': 'GSM', '!Sample_title': 'clinical_sample_title',
        'Survival time (days)': 'OS_time_days',
        'Survival status (0: alive or censored, 1: dead)': 'OS_event'})
    assert clinical.GSM.is_unique and set(clinical.GSM) == set(tumor.GSM)
    patients = tumor.merge(clinical, on='GSM', how='left', validate='one_to_one')
    assert np.array_equal(patients.OS_time_days.to_numpy(),
                          pd.to_numeric(patients['survival time (days)']).to_numpy())
    assert np.array_equal(patients.OS_event.to_numpy(),
                          pd.to_numeric(patients['survival status (0: alive or censored, 1: dead)']).to_numpy())
    linked = adjacent[['GSM', 'sample_title', 'patient_id', 'survival time (days)',
                       'survival_status (0: alive_or_censored, 1:dead)',
                       'time to hcc recurrence (days)',
                       'recurrence status (0: no recurrence of censored, 1: recurrence)']].copy()
    linked.columns = ['adjacent_GSM', 'adjacent_sample_title', 'patient_id', 'adjacent_OS_time_days',
                      'adjacent_OS_event', 'TTR_time_days', 'TTR_event']
    patients = patients.merge(linked, on='patient_id', how='left', validate='one_to_one')
    for col in ['adjacent_OS_time_days', 'adjacent_OS_event', 'TTR_time_days', 'TTR_event']:
        patients[col] = pd.to_numeric(patients[col], errors='coerce')
    matched = patients.adjacent_GSM.notna()
    assert patients.loc[matched, 'OS_time_days'].eq(patients.loc[matched, 'adjacent_OS_time_days']).all()
    assert patients.loc[matched, 'OS_event'].eq(patients.loc[matched, 'adjacent_OS_event']).all()
    patients['recurrence_censor_exceeds_OS_censor'] = (
        patients.TTR_time_days.gt(patients.OS_time_days) & patients.TTR_event.eq(0) & patients.OS_event.eq(0))
    write(patients[['GSM', 'sample_title', 'patient_id', 'OS_time_days', 'OS_event',
                    'adjacent_GSM', 'adjacent_sample_title', 'adjacent_OS_time_days', 'adjacent_OS_event',
                    'TTR_time_days', 'TTR_event', 'recurrence_censor_exceeds_OS_censor']], 'clinical_link_audit.csv')
    return patients


def endpoint_frame(patients, endpoint):
    time = pd.to_numeric(patients[f'{endpoint}_time_days'], errors='coerce')
    event = pd.to_numeric(patients[f'{endpoint}_event'], errors='coerce')
    reason = pd.Series('', index=patients.index)
    reason.loc[~np.isfinite(time) | ~event.isin([0, 1])] = 'No linked valid time/event'
    reason.loc[np.isfinite(time) & time.le(0)] = 'Nonpositive endpoint time'
    if endpoint == 'TTR':
        impossible = patients.OS_event.eq(1) & time.gt(patients.OS_time_days)
        reason.loc[impossible] = 'Recurrence follow-up exceeds known death time'
    good = reason.eq('')
    data = patients.loc[good].copy()
    data['time_days'], data['event'] = time.loc[good], event.loc[good].astype(int)
    data['high'] = data.group.eq('High').astype(int)
    excluded = patients.loc[~good, ['GSM', 'patient_id', 'group', f'{endpoint}_time_days', f'{endpoint}_event']].copy()
    excluded['endpoint'], excluded['exclusion_reason'] = endpoint, reason.loc[~good]
    return data, excluded


def compare(data, endpoint, cutoff):
    high, low = data[data.high.eq(1)], data[data.high.eq(0)]
    assert len(high) > 0 and len(low) > 0 and data.event.sum() > 0
    lr = logrank_test(high.time_days, low.time_days, event_observed_A=high.event, event_observed_B=low.event)
    risksets, oe, variance = [], 0., 0.
    for time in sorted(data.loc[data.event.eq(1), 'time_days'].unique()):
        nh, nl = int(high.time_days.ge(time).sum()), int(low.time_days.ge(time).sum())
        dh = int((high.time_days.eq(time) & high.event.eq(1)).sum())
        dl = int((low.time_days.eq(time) & low.event.eq(1)).sum())
        n, d = nh + nl, dh + dl
        expected = d * nh / n
        v = nh * nl * d * (n - d) / (n * n * (n - 1)) if n > 1 else 0.
        oe += dh - expected
        variance += v
        risksets.append({'endpoint': endpoint, 'time_days': time, 'risk_high': nh, 'risk_low': nl,
                         'events_high': dh, 'events_low': dl, 'expected_high': expected, 'variance': v})
    statistic = oe**2 / variance
    assert np.isclose(statistic, lr.test_statistic, atol=1e-10)
    assert np.isclose(chi2.sf(statistic, 1), lr.p_value, atol=1e-12)
    result = {'endpoint': endpoint, 'patients': len(data), 'events': int(data.event.sum()),
              'high_patients': len(high), 'low_patients': len(low),
              'high_events': int(high.event.sum()), 'low_events': int(low.event.sum()),
              'median_cutoff_all80': cutoff, 'logrank_p': lr.p_value,
              'logrank_chi_square': statistic, 'independent_logrank_p': float(chi2.sf(statistic, 1))}
    cox_rows = []
    for model, covariate in [('binary', 'high'), ('continuous_per_SD', 'coefficient_per_SD')]:
        design = data[['time_days', 'event', covariate]]
        cox = CoxPHFitter().fit(design, 'time_days', 'event')
        cr = cox.summary.loc[covariate]
        ph_p = proportional_hazard_test(cox, design, time_transform='rank').summary.loc[covariate, 'p']
        cox_rows.append({'endpoint': endpoint, 'model': model, 'patients': len(data),
                         'events': int(data.event.sum()), 'HR': cr['exp(coef)'],
                         'HR_lower95': cr['exp(coef) lower 95%'], 'HR_upper95': cr['exp(coef) upper 95%'],
                         'p': cr.p, 'PH_rank_p': ph_p})
        if model == 'binary':
            result.update({'HR_high_vs_low': cr['exp(coef)'], 'HR_lower95': cr['exp(coef) lower 95%'],
                           'HR_upper95': cr['exp(coef) upper 95%'], 'Cox_p': cr.p, 'PH_rank_p': ph_p})
    return result, cox_rows, risksets


def plot(endpoint_data, statistics, cutoff):
    colors = {'High': '#B64F4A', 'Low': '#2B78A5'}
    ticks = np.arange(0, 6001, 1000)
    fig = plt.figure(figsize=(13.5, 8.2))
    grid = fig.add_gridspec(2, 2, height_ratios=[4.8, 1], hspace=.42, wspace=.25)
    curves, risk = [], []
    for col, endpoint in enumerate(['OS', 'TTR']):
        ax, tab = fig.add_subplot(grid[0, col]), fig.add_subplot(grid[1, col])
        data = endpoint_data[endpoint]
        row = statistics[statistics.endpoint.eq(endpoint)].iloc[0]
        for group in ['High', 'Low']:
            q = data[data.group.eq(group)]
            fit = KaplanMeierFitter().fit(q.time_days, q.event, label=group)
            sf, ci = fit.survival_function_, fit.confidence_interval_
            ax.step(sf.index, sf.iloc[:, 0], where='post', color=colors[group], lw=2,
                    label=f'{group} coefficient (n={len(q)})')
            ax.fill_between(sf.index, ci.iloc[:, 0], ci.iloc[:, 1], step='post', alpha=.13, color=colors[group])
            censor = q.loc[q.event.eq(0), 'time_days']
            ax.plot(censor, fit.predict(censor), '+', ms=5, mew=.8, color=colors[group])
            survival = 1.
            for time in sorted(q.time_days.unique()):
                survival *= 1 - int((q.time_days.eq(time) & q.event.eq(1)).sum()) / int(q.time_days.ge(time).sum())
                assert np.isclose(fit.predict(time), survival, atol=1e-12)
            z = pd.concat([sf, ci], axis=1).reset_index()
            z.columns = ['time_days', 'survival', 'CI95_lower', 'CI95_upper']
            z['endpoint'], z['group'] = endpoint, group
            curves.append(z)
            statistics.loc[statistics.endpoint.eq(endpoint), f'median_{group}_days'] = fit.median_survival_time_
            for time in ticks:
                risk.append({'endpoint': endpoint, 'group': group, 'time_days': int(time),
                             'at_risk': int(q.time_days.ge(time).sum())})
        ax.set(title='Overall survival (OS)' if endpoint == 'OS' else 'Time to HCC recurrence (TTR)',
               xlabel='Time (days)', ylabel='Survival probability' if endpoint == 'OS' else 'Recurrence-free probability',
               ylim=(0, 1.02), xlim=(0, 6000), xticks=ticks)
        ax.text(.04, .06, f'Log-rank p = {row.logrank_p:.4g}\nBH q (OS / TTR) = {row.logrank_BH_q:.4g}\n'
                f'HR high / low = {row.HR_high_vs_low:.2f} ({row.HR_lower95:.2f}–{row.HR_upper95:.2f})',
                transform=ax.transAxes, fontsize=10, bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .85})
        ax.legend(loc='upper right', bbox_to_anchor=(1, .90), frameon=False, fontsize=10)
        ax.spines[['top', 'right']].set_visible(False)
        tab.axis('off')
        counts = [[str(int(((data.group == group) & data.time_days.ge(t)).sum())) for t in ticks] for group in ['High', 'Low']]
        table = tab.table(cellText=counts, rowLabels=['High', 'Low'], colLabels=[str(t) for t in ticks], loc='center', cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 1.3)
        tab.set_title('Number at risk (start of time point)', fontsize=10)
    fig.suptitle('GSE10143 HCC tumor tissues: exploratory v4 NNLS projection', fontsize=16, y=.98)
    fig.text(.5, .922, f'416 / 850 genes; 17 fixed components; all-80 median cutoff = {cutoff:.5f}', ha='center', fontsize=11)
    fig.text(.5, .020, 'OS: official tumor clinical data. TTR: linked same-patient recurrence metadata; not PFS or DFS.\n'
             'No PGAM5 probe in platform annotation. Coefficients are exploratory projections, not validated cell abundance.', ha='center', fontsize=10)
    fig.subplots_adjust(left=.075, right=.975, bottom=.18, top=.85)
    for ext in ['png', 'pdf']:
        fig.savefig(R / f'KM_OS_TTR.{ext}', dpi=220, facecolor='white')
    plt.close(fig)
    write(pd.concat(curves, ignore_index=True), 'KM_curve_values.csv')
    write(pd.DataFrame(risk), 'KM_numbers_at_risk.csv')


def main():
    raw, tumor = load_series('GSE10141')
    _, adjacent = load_series('GSE10140')
    assert raw.shape == (6144, 80) and tumor.tissue.eq('hepatocellular carcinoma').all()
    coefficients, summary = project(raw)
    patients = join_clinical(tumor, adjacent)
    patients['target_coefficient'] = patients.GSM.map(coefficients[TARGET])
    cutoff = float(patients.target_coefficient.median())
    patients['group'] = np.where(patients.target_coefficient.gt(cutoff), 'High', 'Low')
    patients['coefficient_per_SD'] = (patients.target_coefficient - patients.target_coefficient.mean()) / patients.target_coefficient.std(ddof=1)
    write(patients, 'patient_coefficients_clinical.csv')
    endpoint_data, stats, models, excluded, risksets = {}, [], [], [], []
    for endpoint in ['OS', 'TTR']:
        data, ex = endpoint_frame(patients, endpoint)
        result, cox, risk = compare(data, endpoint, cutoff)
        endpoint_data[endpoint] = data
        stats.append(result)
        models.extend(cox)
        excluded.append(ex)
        risksets.extend(risk)
        write(data, f'{endpoint}_analysis_patients.csv')
    stats = pd.DataFrame(stats)
    stats['logrank_BH_q'] = multipletests(stats.logrank_p, method='fdr_bh')[1]
    plot(endpoint_data, stats, cutoff)
    write(stats, 'survival_statistics.csv')
    write(pd.DataFrame(models), 'Cox_binary_and_continuous.csv')
    write(pd.concat(excluded, ignore_index=True), 'endpoint_excluded_patients.csv')
    write(pd.DataFrame(risksets), 'independent_logrank_risk_sets.csv')
    summary.update({'median_cutoff': cutoff, 'high_patients': int(patients.group.eq('High').sum()),
                    'low_patients': int(patients.group.eq('Low').sum()),
                    'nonzero_coefficients': int(patients.target_coefficient.gt(0).sum()),
                    'target_min': float(patients.target_coefficient.min()),
                    'target_max': float(patients.target_coefficient.max()),
                    'linked_recurrence_patients_before_QC': int(patients.adjacent_GSM.notna().sum()),
                    'recurrence_censor_exceeds_OS_censor': patients.loc[
                        patients.recurrence_censor_exceeds_OS_censor, 'patient_id'].tolist(),
                    'OS_valid': len(endpoint_data['OS']), 'TTR_valid': len(endpoint_data['TTR'])})
    (R / 'preparation_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (R / 'environment_versions.json').write_text(json.dumps({'python': platform.python_version(),
        **{n: md.version(n) for n in ['numpy', 'pandas', 'scipy', 'matplotlib', 'lifelines', 'statsmodels', 'xlrd']}}, indent=2) + '\n')
    print(json.dumps(summary, indent=2))
    print(stats.to_string(index=False))


if __name__ == '__main__':
    main()
