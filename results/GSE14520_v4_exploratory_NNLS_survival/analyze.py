from pathlib import Path
import csv, gzip, io, json, hashlib, itertools, platform, importlib.metadata as md
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import chi2, spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, proportional_hazard_test
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R = Path(__file__).resolve().parent
I = R / 'inputs'
TARGET = 'TAM_PGAM5_detected'
PLATFORMS = ['GPL3921', 'GPL571']
EP_FIELDS = {'OS': ('Survival months', 'Survival status'), 'DFS': ('Recurr months', 'Recurr status')}
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'pdf.fonttype': 42})


def write(df, name):
    df.to_csv(R / name, index=False, lineterminator='\n')


def load_platform(gpl):
    text = gzip.open(I / f'matrix_{gpl}.txt.gz', 'rt').read()
    matrix = pd.read_csv(io.StringIO(text.split('!series_matrix_table_begin')[1]
                         .split('!series_matrix_table_end')[0].strip()), sep='\t', index_col=0)
    matrix.index = matrix.index.astype(str)
    metadata = {}
    for line in text.splitlines():
        if line.startswith(('!Sample_geo_accession', '!Sample_title', '!Sample_platform_id')):
            row = next(csv.reader([line], delimiter='\t'))
            metadata[row[0]] = row[1:]
    assert metadata['!Sample_geo_accession'] == matrix.columns.tolist()
    assert set(metadata['!Sample_platform_id']) == {gpl}
    meta = pd.DataFrame({'Affy_GSM': matrix.columns, 'platform': gpl,
                         'sample_title': metadata['!Sample_title']})
    annotation = pd.read_csv(I / f'{gpl}_annotation.tsv.gz', sep='\t', dtype=str)
    annotation['ID'] = annotation.ID.str.strip()
    annotation['gene'] = annotation['Gene Symbol'].fillna('').str.strip()
    annotation['in_expression_matrix'] = annotation.ID.isin(matrix.index)
    annotation['unambiguous'] = (annotation.gene.ne('') & annotation.gene.ne('---') &
                                 ~annotation.gene.str.contains('///', regex=False))
    used = annotation[annotation.in_expression_matrix & annotation.unambiguous].copy()
    assert used.ID.is_unique and matrix.index.is_unique and matrix.columns.is_unique
    assert np.isfinite(matrix.to_numpy()).all() and matrix.min().min() > 0
    summarized = matrix.loc[used.ID].copy()
    summarized.index = used.gene.to_numpy()
    gene_log2 = summarized.groupby(level=0, sort=True).median()
    gene_linear = np.exp2(gene_log2)
    denominator = gene_linear.sum(axis=0)
    normalized = gene_linear.div(denominator, axis=1) * 1e4
    used['platform'] = gpl
    write(used[['platform','ID','gene','ENTREZ_GENE_ID','RefSeq Transcript ID']], f'{gpl}_used_probe_mapping.csv')
    stats = {'platform': gpl, 'samples': matrix.shape[1], 'probes': matrix.shape[0],
             'RMA_min': float(matrix.min().min()), 'RMA_max': float(matrix.max().max()),
             'mapped_probes': len(used), 'mapped_genes': len(normalized),
             'PGAM5_any_annotation': bool(annotation['Gene Symbol'].fillna('').str.contains(r'\bPGAM5\b', regex=True).any()),
             'PGAM5_measured': 'PGAM5' in normalized.index}
    return normalized, denominator, meta, stats


def solve(reference, scales, genes, query):
    ref = reference.loc[genes]
    weight = scales.loc[genes].to_numpy()
    A = ref.to_numpy() / weight[:,None]
    Y = query.loc[genes].to_numpy() / weight[:,None]
    coefficients, diagnostics = [], []
    for j, sample in enumerate(query.columns):
        y = Y[:,j]
        beta, residual = nnls(A, y, maxiter=10000)
        assert beta.sum() > 0 and (beta >= 0).all()
        fractions = beta / beta.sum()
        assert np.isclose(fractions.sum(), 1)
        coefficients.append(fractions)
        gradient = A.T @ (A @ beta - y)
        scale = max(1., np.linalg.norm(A) * np.linalg.norm(y))
        assert gradient.min() / scale > -1e-9
        assert np.max(np.abs(gradient * beta)) / max(1., np.linalg.norm(y)**2) < 1e-9
        diagnostics.append({'Affy_GSM': sample, 'fit_relative_residual': residual / np.linalg.norm(y),
                            'raw_coefficient_sum': beta.sum(), 'positive_components': int((beta > 1e-12).sum()),
                            'KKT_scaled_min_gradient': gradient.min()/scale})
    fractions = pd.DataFrame(coefficients, index=query.columns, columns=reference.columns)
    return fractions, pd.DataFrame(diagnostics), A


def endpoint_frame(scores, endpoint):
    time_field, event_field = EP_FIELDS[endpoint]
    times = pd.to_numeric(scores[time_field], errors='coerce')
    events = pd.to_numeric(scores[event_field], errors='coerce')
    valid = np.isfinite(times) & times.gt(0) & events.isin([0,1])
    data = scores.loc[valid].copy()
    data['time_months'] = times[valid].to_numpy()
    data['event'] = events[valid].astype(int).to_numpy()
    data['high'] = data.group.eq('High').astype(int)
    return data, scores.loc[~valid].copy()


def compare(scores, endpoint, label, cutoff):
    data, excluded = endpoint_frame(scores, endpoint)
    high, low = data[data.high.eq(1)], data[data.high.eq(0)]
    if len(high) == 0 or len(low) == 0:
        return {'analysis': label, 'endpoint': endpoint, 'status': 'Not estimable: one empty group'}, data, excluded, []
    lr = logrank_test(high.time_months, low.time_months, event_observed_A=high.event, event_observed_B=low.event)
    observed_minus_expected, variance, risksets = 0., 0., []
    for time in sorted(data.loc[data.event.eq(1), 'time_months'].unique()):
        nh, nl = int(high.time_months.ge(time).sum()), int(low.time_months.ge(time).sum())
        dh = int((high.time_months.eq(time) & high.event.eq(1)).sum())
        dl = int((low.time_months.eq(time) & low.event.eq(1)).sum())
        n,d = nh+nl, dh+dl
        expected = d*nh/n
        v = nh*nl*d*(n-d)/(n*n*(n-1)) if n > 1 else 0.
        observed_minus_expected += dh-expected
        variance += v
        risksets.append({'analysis': label, 'endpoint': endpoint, 'time_months': time,
                         'risk_high': nh, 'risk_low': nl, 'events_high': dh, 'events_low': dl,
                         'expected_high': expected, 'variance': v})
    statistic = observed_minus_expected**2/variance
    independent_p = chi2.sf(statistic, 1)
    assert np.isclose(statistic, lr.test_statistic, atol=1e-10)
    assert np.isclose(independent_p, lr.p_value, atol=1e-12)
    out = {'analysis': label, 'endpoint': endpoint, 'status': 'Computed', 'patients': len(data),
           'events': int(data.event.sum()), 'high_patients': len(high), 'low_patients': len(low),
           'high_events': int(high.event.sum()), 'low_events': int(low.event.sum()),
           'cutoff': cutoff, 'excluded': len(excluded), 'logrank_p': lr.p_value,
           'logrank_chi_square': lr.test_statistic, 'independent_logrank_p': independent_p}
    cox = CoxPHFitter().fit(data[['time_months','event','high']], 'time_months','event')
    cr = cox.summary.loc['high']
    out.update({'HR_high_vs_low': cr['exp(coef)'], 'HR_lower95': cr['exp(coef) lower 95%'],
                'HR_upper95': cr['exp(coef) upper 95%'], 'Cox_p': cr.p,
                'PH_p': proportional_hazard_test(cox, data[['time_months','event','high']], time_transform='rank').summary.loc['high','p']})
    return out, data, excluded, risksets


def fit_covariate_models(scores, endpoint):
    data,_ = endpoint_frame(scores, endpoint)
    data['score_per_SD'] = (data.target_coefficient - scores.target_coefficient.mean()) / scores.target_coefficient.std(ddof=1)
    data['age_per_10y'] = pd.to_numeric(data.Age, errors='coerce') / 10.
    data['male'] = data.Gender.map({'M':1., 'F':0.})
    def stage_group(value):
        if pd.isna(value): return np.nan
        s = str(value).strip()
        if s.startswith(('III','IV')): return 'III-IV'
        if s.startswith('II'): return 'II'
        if s.startswith('I'): return 'I'
        return np.nan
    data['stage_group'] = data['TNM staging'].map(stage_group)
    data['platform_GPL571'] = data.platform.eq('GPL571').astype(float)
    rows = []
    for model, covariate, adjusted in [('continuous_unadjusted','score_per_SD',False),
                                      ('binary_adjusted','high',True), ('continuous_adjusted','score_per_SD',True)]:
        columns = ['time_months','event',covariate]
        if adjusted:
            complete = data.dropna(subset=['age_per_10y','male','stage_group'])
            design = complete[columns+['age_per_10y','male','platform_GPL571']].copy()
            design['TNM_II'] = complete.stage_group.eq('II').astype(float)
            design['TNM_III_IV'] = complete.stage_group.eq('III-IV').astype(float)
            # Remove a structurally constant adjustment term; record it rather than silently fit.
            removed = [c for c in design if c not in ['time_months','event',covariate] and design[c].nunique() <= 1]
            design = design.drop(columns=removed)
        else:
            design = data[columns].copy(); removed=[]
        write(pd.concat([data.loc[design.index,['Affy_GSM','ID']], design],axis=1), f'{endpoint}_{model}_patients.csv')
        cox = CoxPHFitter().fit(design, 'time_months','event')
        summary = cox.summary.copy().reset_index()
        summary.insert(0,'model',model);summary.insert(0,'endpoint',endpoint)
        summary['patients'] = len(design);summary['events'] = int(design.event.sum())
        summary['removed_constant_adjustment_terms'] = ';'.join(removed)
        ph = proportional_hazard_test(cox, design, time_transform='rank').summary
        summary['PH_rank_p'] = summary.covariate.map(ph.p)
        rows.append(summary)
    return pd.concat(rows,ignore_index=True)


def plot_km(endpoint_data, stats, ticks):
    fig = plt.figure(figsize=(13.5,7.8))
    grid = fig.add_gridspec(2,2,height_ratios=[4.8,1.0],hspace=.43,wspace=.23)
    curves,risk = [],[]
    colors = {'High':'#B64F4A','Low':'#2B78A5'}
    for col, endpoint in enumerate(['OS','DFS']):
        ax,tab = fig.add_subplot(grid[0,col]),fig.add_subplot(grid[1,col])
        data = endpoint_data[endpoint]
        row = stats[stats.endpoint.eq(endpoint)].iloc[0]
        for group in ['High','Low']:
            q = data[data.group.eq(group)]
            fit = KaplanMeierFitter().fit(q.time_months,q.event,label=group)
            sf,ci = fit.survival_function_,fit.confidence_interval_
            ax.step(sf.index,sf.iloc[:,0],where='post',lw=2,color=colors[group],label=f'{group} coefficient (n={len(q)})')
            ax.fill_between(sf.index,ci.iloc[:,0],ci.iloc[:,1],step='post',alpha=.13,color=colors[group])
            censored = q.loc[q.event.eq(0),'time_months']
            ax.plot(censored,fit.predict(censored),'+',ms=5,mew=.8,color=colors[group])
            independent_survival = 1.
            for time in sorted(q.time_months.unique()):
                n = q.time_months.ge(time).sum()
                d = (q.time_months.eq(time)&q.event.eq(1)).sum()
                independent_survival *= 1-d/n
                assert np.isclose(fit.predict(time),independent_survival,atol=1e-12)
            z = pd.concat([sf,ci],axis=1).reset_index()
            z.columns=['time_months','survival','CI95_lower','CI95_upper']
            z['endpoint'],z['group']=endpoint,group;curves.append(z)
            for time in ticks:
                risk.append({'endpoint':endpoint,'group':group,'time_months':time,'at_risk':int(q.time_months.ge(time).sum())})
            stats.loc[stats.endpoint.eq(endpoint),f'median_{group}_months'] = fit.median_survival_time_
        ax.set(title='Overall survival (OS)' if endpoint=='OS' else 'Disease-free survival (DFS)',
               xlabel='Time (months)',ylabel='Survival probability',ylim=(0,1.02),xlim=(0,ticks[-1]),xticks=ticks)
        ax.text(.04,.06,f'Log-rank p = {row.logrank_p:.4g}\nBH q (OS/DFS) = {row.logrank_BH_q:.4g}\n'
                f'HR high/low = {row.HR_high_vs_low:.2f} ({row.HR_lower95:.2f}–{row.HR_upper95:.2f})',
                transform=ax.transAxes,fontsize=10,va='bottom',bbox={'facecolor':'white','edgecolor':'none','alpha':.85})
        ax.legend(loc='upper right',bbox_to_anchor=(1, .87),frameon=False,fontsize=10)
        ax.spines[['top','right']].set_visible(False)
        tab.axis('off')
        counts=[[str(int(((data.group==group)&data.time_months.ge(time)).sum())) for time in ticks] for group in ['High','Low']]
        table=tab.table(cellText=counts,rowLabels=['High','Low'],colLabels=[str(int(x)) for x in ticks],loc='center',cellLoc='center')
        table.auto_set_font_size(False);table.set_fontsize(9);table.scale(1,1.3)
        tab.set_title('Number at risk (start of time point)',fontsize=10)
    fig.suptitle('GSE14520: exploratory v4 NNLS reference-projected coefficient',fontsize=17,y=.98)
    fig.text(.5,.925,'850-gene / 17-component reference; 569 common measured genes; global median = 0: High > 0, Low = 0',ha='center',fontsize=11)
    fig.text(.5,.02,'No PGAM5 probe identified in supplied annotations. Coefficients are exploratory projections, not validated cell abundance.',ha='center',fontsize=10)
    fig.subplots_adjust(top=.85,bottom=.11,left=.07,right=.98)
    for ext in ['png','pdf']:fig.savefig(R/f'KM_OS_DFS.{ext}',dpi=220,facecolor='white')
    plt.close(fig)
    write(pd.concat(curves,ignore_index=True),'KM_curve_values.csv')
    write(pd.DataFrame(risk),'KM_numbers_at_risk.csv')


def main():
    reference=pd.read_csv(I/'reference_850x17.tsv',sep='\t',index_col=0)
    scales=pd.read_csv(I/'reference_row_scales.csv',index_col=0).row_scale
    assert reference.shape==(850,17) and reference.index.is_unique and reference.columns.is_unique
    assert (scales.loc[reference.index]>0).all()
    clinical=pd.read_csv(I/'GSE14520_Extra_Supplement.txt.gz',sep='\t',na_values=['.',''])
    tumor=clinical[clinical['Tissue Type'].eq('Tumor')].copy()
    assert len(tumor)==247 and tumor.Affy_GSM.is_unique and tumor.ID.is_unique
    platforms={p:load_platform(p) for p in PLATFORMS}
    query=pd.concat([platforms[p][0] for p in PLATFORMS],axis=1,join='inner')
    meta=pd.concat([platforms[p][2] for p in PLATFORMS],ignore_index=True)
    meta['linear_gene_denominator']=meta.Affy_GSM.map(pd.concat([platforms[p][1] for p in PLATFORMS]))
    assert meta.Affy_GSM.is_unique and len(meta)==488
    tumor=tumor.merge(meta,on='Affy_GSM',how='left',validate='one_to_one')
    assert tumor.platform.notna().all() and tumor.sample_title.str.startswith('Liver Tumor Tissue').all()
    query=query[tumor.Affy_GSM]
    genes=[g for g in reference.index if g in query.index]
    coverage=pd.DataFrame({'gene':reference.index})
    for p in PLATFORMS:coverage[p+'_measured']=coverage.gene.isin(platforms[p][0].index)
    coverage['used_primary']=coverage.gene.isin(genes)
    write(coverage,'reference_gene_coverage.csv')
    reference.loc[genes].to_csv(R/'used_reference.tsv',sep='\t')
    scales.loc[genes].rename('row_scale').to_csv(R/'used_reference_row_scales.csv')
    query.loc[genes].to_csv(R/'primary_array_query_normalized_intensity.tsv.gz',sep='\t')
    fractions,diagnostics,A=solve(reference,scales,genes,query)
    write(fractions.rename_axis('Affy_GSM').reset_index(),'all_component_coefficients.csv')
    write(diagnostics,'fit_diagnostics.csv')
    scores=tumor.merge(fractions[TARGET].rename('target_coefficient'),left_on='Affy_GSM',right_index=True,validate='one_to_one')
    cutoff=float(scores.target_coefficient.median())
    scores['group']=np.where(scores.target_coefficient>cutoff,'High','Low')
    write(scores,'patient_coefficients_clinical.csv')
    summaries=[platforms[p][3] for p in PLATFORMS]
    for s in summaries:
        s['clinical_tumor_patients']=int(scores.platform.eq(s['platform']).sum())
        s['measurable_reference_genes']=int(coverage[s['platform']+'_measured'].sum())
    preparation={'original_reference_genes':850,'components':17,'common_measured_reference_genes':len(genes),
                 'coverage_fraction':len(genes)/850,'PGAM5_measured':False,'patients':len(scores),
                 'global_median':cutoff,'high_patients':int(scores.group.eq('High').sum()),
                 'low_patients':int(scores.group.eq('Low').sum()),'nonzero_coefficients':int(scores.target_coefficient.gt(1e-8).sum()),
                 'min_coefficient':float(scores.target_coefficient.min()),'max_coefficient':float(scores.target_coefficient.max()),
                 'matrix_rank':int(np.linalg.matrix_rank(A)),'condition_number':float(np.linalg.cond(A)),
                 'platforms':summaries,'fit_relative_residual_median':float(diagnostics.fit_relative_residual.median())}
    (R/'preparation_summary.json').write_text(json.dumps(preparation,indent=2)+'\n')
    primary,endpoint_data,excluded,risksets,secondary=[],{},[],[],[]
    for endpoint in EP_FIELDS:
        row,data,ex,risk=compare(scores,endpoint,'primary_common_gene_projection',cutoff)
        primary.append(row);endpoint_data[endpoint]=data
        write(data,f'{endpoint}_analysis_patients.csv')
        ex['endpoint']=endpoint;ex['reason']='missing/nonfinite/nonpositive time or nonbinary event';excluded.append(ex)
        risksets+=risk
    stats=pd.DataFrame(primary)
    stats['logrank_BH_q']=multipletests(stats.logrank_p,method='fdr_bh')[1]
    max_time=max(data.time_months.max() for data in endpoint_data.values())
    ticks=np.arange(0,np.ceil(max_time/12)*12+1,12)
    plot_km(endpoint_data,stats,ticks)
    write(stats,'survival_statistics.csv')
    write(pd.concat(excluded,ignore_index=True),'endpoint_excluded_patients.csv')
    cox=pd.concat([fit_covariate_models(scores,ep) for ep in EP_FIELDS],ignore_index=True)
    write(cox,'Cox_continuous_and_adjusted.csv')
    for p in PLATFORMS:
        subset=scores[scores.platform.eq(p)].copy()
        local_cut=float(subset.target_coefficient.median())
        subset['group']=np.where(subset.target_coefficient>local_cut,'High','Low')
        for endpoint in EP_FIELDS:
            row,_,_,risk=compare(subset,endpoint,f'{p}_primary_projection_platform_median',local_cut)
            secondary.append(row);risksets+=risk
        local_genes=[g for g in reference.index if g in platforms[p][0].index]
        fq,diag,_=solve(reference,scales,local_genes,platforms[p][0][subset.Affy_GSM])
        native=subset.drop(columns='target_coefficient').merge(fq[TARGET].rename('target_coefficient'),left_on='Affy_GSM',right_index=True,validate='one_to_one')
        native_cut=float(native.target_coefficient.median())
        native['group']=np.where(native.target_coefficient>native_cut,'High','Low')
        write(native,f'{p}_all_measurable_reference_genes_patients.csv')
        write(fq.rename_axis('Affy_GSM').reset_index(),f'{p}_all_measurable_component_coefficients.csv')
        for endpoint in EP_FIELDS:
            row,_,_,risk=compare(native,endpoint,f'{p}_all_measurable_genes_platform_median',native_cut)
            row['reference_genes']=len(local_genes)
            secondary.append(row);risksets+=risk
    sec=pd.DataFrame(secondary)
    sec['BH_q_8_secondary_tests']=multipletests(sec.logrank_p,method='fdr_bh')[1]
    write(sec,'platform_sensitivity_survival_statistics.csv')
    write(pd.DataFrame(risksets),'independent_logrank_risk_sets.csv')
    platform_values=[scores.loc[scores.platform.eq(p),'target_coefficient'] for p in PLATFORMS]
    u=mannwhitneyu(*platform_values,alternative='two-sided')
    (R/'platform_coefficient_comparison.json').write_text(json.dumps({'test':'Mann-Whitney, descriptive platform/cohort comparison; cannot separate platform and cohort',
        'p':float(u.pvalue),'medians':{p:float(scores.loc[scores.platform.eq(p),'target_coefficient'].median()) for p in PLATFORMS}},indent=2)+'\n')
    (R/'environment_versions.json').write_text(json.dumps({'python':platform.python_version(),**{n:md.version(n) for n in ['numpy','pandas','scipy','lifelines','matplotlib','statsmodels']}},indent=2)+'\n')
    print(json.dumps(preparation,indent=2),flush=True)
    print(stats.to_string(index=False),flush=True)
    print(sec.to_string(index=False),flush=True)


if __name__=='__main__':main()
