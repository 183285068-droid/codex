"""Exploratory coefficient-survival association; no abundance or causal claim."""
from pathlib import Path
import json
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import chi2
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, proportional_hazard_test
from statsmodels.stats.multitest import multipletests
from statsmodels.duration.hazard_regression import PHReg

ROOT = Path(__file__).resolve().parent
R = ROOT/'survival_exploratory'


def stage(value):
    if pd.isna(value):
        return None
    text = re.sub(r'^STAGE\s*', '', str(value).strip().upper())
    match = re.fullmatch(r'(IV|III|II|I)[ABC]?', text)
    if not match:
        return None
    return {'I':'I','II':'II','III':'III-IV','IV':'III-IV'}[match.group(1)]


def main():
    protocol = json.loads((R/'protocol.json').read_text())
    coeff = pd.read_csv(ROOT/'EXPLORATORY_TCGA_coefficients_NOT_CELL_FRACTIONS.csv')
    coeff['coefficient'] = coeff.exploratory_PGAM5_RNA_macrophage_coefficient
    assert np.isclose(coeff.coefficient.median(), protocol['global_median'], atol=1e-15)
    coeff['group'] = np.where(coeff.coefficient.gt(protocol['global_median']), 'High','Low')
    coeff['high'] = coeff.group.eq('High').astype(int)
    coeff['coefficient_per_IQR'] = coeff.coefficient/protocol['global_IQR']
    surv = pd.read_csv(R/'TCGA_CDR_mirror_survival.csv.gz')
    surv['patient'] = surv['sample'].str[:12]
    surv = surv[surv.patient.isin(coeff.patient)].copy()
    endpoints = ['OS','OS.time','DSS','DSS.time','DFI','DFI.time','PFI','PFI.time']
    assert surv.groupby('patient')[endpoints].nunique(dropna=False).le(1).all().all()
    surv.to_csv(R/'survival_alias_consistency_audit.csv',index=False)
    surv = surv[surv['sample'].str[13:15].eq('01')].rename(columns={'sample':'survival_sample'})
    assert surv.patient.is_unique
    data = coeff.merge(surv,on='patient',how='left',validate='one_to_one',indicator='survival_match')
    clinical = pd.read_csv(R/'TCGA_LIHC_clinical_primary.csv').rename(columns={'sample':'clinical_sample'})
    data = data.merge(clinical,on='patient',how='left',validate='one_to_one')
    data['age_years'] = pd.to_numeric(data.age_at_initial_pathologic_diagnosis,errors='coerce')
    data.loc[~data.age_years.between(18,120),'age_years'] = np.nan
    data['age_per_10_years'] = data.age_years/10
    data['male'] = data.gender.str.upper().map({'MALE':1.,'FEMALE':0.})
    data['stage_group'] = data.ajcc_pathologic_tumor_stage.map(stage)
    data['stage_II'] = data.stage_group.eq('II').astype(float)
    data['stage_III_IV'] = data.stage_group.eq('III-IV').astype(float)
    data['covariates_complete'] = data.age_years.notna() & data.male.notna() & data.stage_group.notna()
    data.to_csv(R/'patient_coefficients_survival_join.csv',index=False)
    statistics, exclusions, curves, risksets, risk, coxrows, phrows, checks = [],[],[],[],[],[],[],[]
    fits = {}
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
    for ep in ['OS','PFI']:
        time = pd.to_numeric(data[ep+'.time'],errors='coerce')
        event = pd.to_numeric(data[ep],errors='coerce')
        valid = np.isfinite(time) & time.gt(0) & event.isin([0,1])
        t = data.loc[valid].copy()
        t['time_days'] = time[valid].to_numpy()
        t['event'] = event[valid].astype(int).to_numpy()
        t.to_csv(R/f'{ep}_analysis_patients.csv',index=False)
        for _, row in data.loc[~valid].iterrows():
            reasons = []
            if row.survival_match != 'both': reasons.append('no matched survival record')
            if not np.isfinite(row[ep+'.time']) or row[ep+'.time']<=0: reasons.append('missing/nonpositive/nonfinite time')
            if row[ep] not in [0,1]: reasons.append('missing/nonbinary event')
            exclusions.append({'endpoint':ep,'patient':row.patient,'group':row.group,
                'time_days':row[ep+'.time'],'event':row[ep],'reason':'; '.join(reasons)})
        high, low = t[t.high.eq(1)], t[t.high.eq(0)]
        lr = logrank_test(high.time_days,low.time_days,event_observed_A=high.event,event_observed_B=low.event)
        oe, variance = 0.,0.
        for day in sorted(t.loc[t.event.eq(1),'time_days'].unique()):
            yh, yl = int(high.time_days.ge(day).sum()), int(low.time_days.ge(day).sum())
            dh = int((high.time_days.eq(day)&high.event.eq(1)).sum())
            dl = int((low.time_days.eq(day)&low.event.eq(1)).sum())
            n, d = yh+yl, dh+dl
            expected = d*yh/n
            var = yh*yl*d*(n-d)/(n*n*(n-1)) if n>1 else 0.
            oe += dh-expected
            variance += var
            risksets.append({'endpoint':ep,'day':day,'risk_high':yh,'risk_low':yl,'events_high':dh,
                             'events_low':dl,'expected_high':expected,'variance':var})
        stat = oe**2/variance
        assert np.isclose(stat,lr.test_statistic,atol=1e-10)
        assert np.isclose(chi2.sf(stat,1),lr.p_value,atol=1e-12)
        checks.append(ep+' independently reconstructed log-rank')
        out = {'endpoint':ep,'patients':len(t),'events':int(t.event.sum()),
               'high_patients':len(high),'high_events':int(high.event.sum()),
               'low_patients':len(low),'low_events':int(low.event.sum()),
               'excluded_patients':int((~valid).sum()),'median_cutpoint':protocol['global_median'],
               'logrank_chi_square':float(lr.test_statistic),'logrank_p':float(lr.p_value)}
        endpoint_fits = {}
        for label, d in [('High',high),('Low',low)]:
            km = KaplanMeierFitter().fit(d.time_days,event_observed=d.event,label=label)
            endpoint_fits[label] = km
            product, greenwood = 1.,0.
            for day in sorted(d.time_days.unique()):
                n = int(d.time_days.ge(day).sum())
                ne = int((d.time_days.eq(day)&d.event.eq(1)).sum())
                product *= 1-ne/n
                assert np.isclose(float(km.predict(day)),product,atol=1e-12)
                if ne and n>ne: greenwood += ne/(n*(n-ne))
                if 0 < product < 1:
                    loglog = np.log(-np.log(product))
                    se = np.sqrt(greenwood)/abs(np.log(product))
                    ci = [np.exp(-np.exp(loglog+1.959963984540054*se)),
                          np.exp(-np.exp(loglog-1.959963984540054*se))]
                    assert np.allclose(ci,km.confidence_interval_.loc[day].to_numpy(),atol=1e-11)
            checks.append(ep+' '+label+' all KM product-limit and log-log CI steps')
            cur = pd.concat([km.survival_function_,km.confidence_interval_],axis=1).reset_index()
            cur.columns = ['time_days','survival','CI95_lower','CI95_upper']
            cur['endpoint'],cur['group'] = ep,label
            curves.append(cur)
            out['median_'+label+'_days'] = float(km.median_survival_time_) if np.isfinite(km.median_survival_time_) else np.nan
        fits[ep] = (t,endpoint_fits)
        covariates = ['age_per_10_years','male','stage_II','stage_III_IV']
        for name, term, adjusted in [('high_low_unadjusted','high',False),
                ('continuous_unadjusted','coefficient_per_IQR',False),
                ('high_low_adjusted','high',True),('continuous_adjusted','coefficient_per_IQR',True)]:
            use = t[t.covariates_complete].copy() if adjusted else t.copy()
            variables = [term]+covariates if adjusted else [term]
            if adjusted:
                use.to_csv(R/f'{ep}_adjusted_analysis_patients.csv',index=False)
            modeldata = use[['time_days','event']+variables]
            cox = CoxPHFitter().fit(modeldata,duration_col='time_days',event_col='event',
                                   fit_options={'precision':1e-10,'max_steps':500})
            sm = PHReg(use.time_days.to_numpy(),use[variables].to_numpy(),status=use.event.to_numpy(),ties='efron').fit()
            assert np.allclose(cox.params_.to_numpy(),sm.params,rtol=1e-4,atol=1e-5)
            assert np.allclose(cox.variance_matrix_.to_numpy(),sm.cov_params(),rtol=1e-4,atol=1e-5)
            checks.append(ep+' '+name+' Efron Cox coefficients/covariance vs statsmodels PHReg')
            ph = proportional_hazard_test(cox,modeldata,time_transform='rank').summary
            for variable, row in cox.summary.iterrows():
                coxrows.append({'endpoint':ep,'model':name,'term':variable,'target_term':variable==term,
                    'patients':len(use),'events':int(use.event.sum()),'log_HR':row['coef'],
                    'HR':row['exp(coef)'],'CI95_lower':row['exp(coef) lower 95%'],
                    'CI95_upper':row['exp(coef) upper 95%'],'Wald_p':row['p'],
                    'PH_rank_time_p':ph.loc[variable,'p'],'secondary_BH_q_6_target_tests':np.nan})
                phrows.append({'endpoint':ep,'model':name,'term':variable,
                               'test_statistic':ph.loc[variable,'test_statistic'],'p':ph.loc[variable,'p']})
            if name=='high_low_unadjusted':
                h = cox.summary.loc[term]
                out.update({'HR_high_vs_low':h['exp(coef)'],'HR_CI95_lower':h['exp(coef) lower 95%'],
                            'HR_CI95_upper':h['exp(coef) upper 95%'],'Cox_high_low_p':h['p'],
                            'PH_high_low_p':ph.loc[term,'p']})
        for _, row in t.loc[~t.covariates_complete].iterrows():
            missing = [field for field in ['age_years','male','stage_group'] if pd.isna(row[field])]
            exclusions.append({'endpoint':ep+'_adjusted','patient':row.patient,'group':row.group,
                'time_days':row.time_days,'event':row.event,'reason':'incomplete covariates: '+', '.join(missing)})
        # Post-diagnostic sensitivity: stage terms showed PH departures for PFI.
        # Run both endpoints, retain the original analyses, and label this post hoc.
        use = t[t.covariates_complete].copy()
        variables = ['coefficient_per_IQR','age_per_10_years','male']
        modeldata = use[['time_days','event','stage_group']+variables]
        stratified = CoxPHFitter().fit(modeldata,duration_col='time_days',event_col='event',
                                      strata=['stage_group'],fit_options={'precision':1e-10,'max_steps':500})
        sm = PHReg(use.time_days.to_numpy(),use[variables].to_numpy(),status=use.event.to_numpy(),
                   strata=use.stage_group.to_numpy(),ties='efron').fit()
        assert np.allclose(stratified.params_.to_numpy(),sm.params,rtol=1e-4,atol=1e-5)
        assert np.allclose(stratified.variance_matrix_.to_numpy(),sm.cov_params(),rtol=1e-4,atol=1e-5)
        checks.append(ep+' posthoc stage-stratified Cox vs independent PHReg')
        ph = proportional_hazard_test(stratified,modeldata,time_transform='rank').summary
        for variable, row in stratified.summary.iterrows():
            coxrows.append({'endpoint':ep,'model':'posthoc_continuous_stage_stratified','term':variable,
                'target_term':variable=='coefficient_per_IQR','patients':len(use),'events':int(use.event.sum()),
                'log_HR':row['coef'],'HR':row['exp(coef)'],'CI95_lower':row['exp(coef) lower 95%'],
                'CI95_upper':row['exp(coef) upper 95%'],'Wald_p':row['p'],
                'PH_rank_time_p':ph.loc[variable,'p'],'secondary_BH_q_6_target_tests':np.nan})
            phrows.append({'endpoint':ep,'model':'posthoc_continuous_stage_stratified','term':variable,
                           'test_statistic':ph.loc[variable,'test_statistic'],'p':ph.loc[variable,'p']})
        statistics.append(out)
    stats = pd.DataFrame(statistics)
    stats['logrank_BH_q_2_endpoints'] = multipletests(stats.logrank_p,method='fdr_bh')[1]
    stats.to_csv(R/'survival_statistics.csv',index=False)
    cox = pd.DataFrame(coxrows)
    secondary = cox.target_term & cox.model.isin(['continuous_unadjusted','high_low_adjusted','continuous_adjusted'])
    assert secondary.sum()==6
    cox.loc[secondary,'secondary_BH_q_6_target_tests'] = multipletests(cox.loc[secondary,'Wald_p'],method='fdr_bh')[1]
    posthoc = cox.target_term & cox.model.eq('posthoc_continuous_stage_stratified')
    assert posthoc.sum()==2
    cox['posthoc_stage_stratified_BH_q_2'] = np.nan
    cox.loc[posthoc,'posthoc_stage_stratified_BH_q_2'] = multipletests(cox.loc[posthoc,'Wald_p'],method='fdr_bh')[1]
    cox.to_csv(R/'Cox_all_models.csv',index=False)
    pd.DataFrame(phrows).to_csv(R/'PH_diagnostics.csv',index=False)
    pd.DataFrame(exclusions).to_csv(R/'excluded_patients.csv',index=False)
    pd.concat(curves).to_csv(R/'KM_curve_values.csv',index=False)
    pd.DataFrame(risksets).to_csv(R/'independent_logrank_risk_sets.csv',index=False)
    for ep, (t,models) in fits.items():
        result = stats[stats.endpoint.eq(ep)].iloc[0]
        fig = plt.figure(figsize=(9.2,7.6))
        gs = fig.add_gridspec(2,1,height_ratios=[4.5,1.2],hspace=.32)
        ax,tab = fig.add_subplot(gs[0]),fig.add_subplot(gs[1])
        colors = {'High':'#c44940','Low':'#2673af'}
        max_month = t.time_days.max()/30.4375
        ticks = np.arange(0,np.ceil(max_month/24)*24+1,24)
        for label, km in models.items():
            sf,ci = km.survival_function_,km.confidence_interval_
            xx = sf.index.to_numpy()/30.4375
            groupdata = t[t.group.eq(label)]
            ax.step(xx,sf.iloc[:,0],where='post',lw=2,color=colors[label],
                    label=f'{label} coefficient (n={len(groupdata)})')
            ax.fill_between(xx,ci.iloc[:,0],ci.iloc[:,1],step='post',alpha=.15,color=colors[label])
            censored = groupdata.loc[groupdata.event.eq(0),'time_days'].to_numpy()
            ax.plot(censored/30.4375,km.predict(censored).to_numpy(),'+',ms=5,mew=.8,color=colors[label])
            for month in ticks:
                risk.append({'endpoint':ep,'group':label,'month':float(month),
                             'at_risk':int(groupdata.time_days.ge(month*30.4375).sum())})
        ax.set(title='TCGA-LIHC: '+('Overall survival (OS)' if ep=='OS' else 'Progression-free interval (PFI)')+
                     '\nExploratory PGAM5 RNA-detected macrophage reference coefficient',
               xlabel='Time (months)',ylabel='Survival probability',ylim=(0,1.02),xlim=(0,ticks[-1]),xticks=ticks)
        ax.text(.98,.98,f'Log-rank p = {result.logrank_p:.4g}\nBH q (2 endpoints) = {result.logrank_BH_q_2_endpoints:.4g}\n'
                f'HR high/low = {result.HR_high_vs_low:.2f}\n95% CI {result.HR_CI95_lower:.2f}-{result.HR_CI95_upper:.2f}',
                transform=ax.transAxes,ha='right',va='top',fontsize=10,
                bbox={'facecolor':'white','alpha':.8,'edgecolor':'none'})
        ax.legend(loc='lower left',frameon=False)
        ax.spines[['top','right']].set_visible(False)
        tab.axis('off')
        numbers = [[str(int(t[t.group.eq(label)].time_days.ge(m*30.4375).sum())) for m in ticks] for label in ['High','Low']]
        table = tab.table(cellText=numbers,rowLabels=['High','Low'],colLabels=[str(int(m)) for m in ticks],
                          loc='center',cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1,1.4)
        tab.set_title('Number at risk (start of time point)',fontsize=10,pad=4)
        footer = 'High > global median; Low <= median. Unvalidated reference: coefficients are not cell infiltration.\n'
        footer += 'PFI is the available TCGA-CDR progression endpoint; it is not strictly defined PFS.'
        if ep=='PFI':
            footer += f'\nHigh/low PH diagnostic p = {result.PH_high_low_p:.4g}; the fitted HR is not constant over time.'
        fig.text(.5,.025,footer,ha='center',fontsize=8.8)
        fig.subplots_adjust(top=.89,bottom=.12,left=.12,right=.96)
        fig.savefig(R/f'KM_{ep}.png',dpi=190)
        fig.savefig(R/f'KM_{ep}.pdf')
        plt.close(fig)
    pd.DataFrame(risk).to_csv(R/'KM_numbers_at_risk.csv',index=False)
    (R/'implementation_checks.json').write_text(json.dumps({'status':'PASS','checks':checks,
        'n_checks':len(checks),'grouping_fixed_before_endpoint_analysis':True,
        'prognostic_specificity_or_cell_abundance_validated':False},indent=2)+'\n')
    print(stats.to_string(index=False))
    print(cox[cox.target_term].to_string(index=False))


if __name__ == '__main__':
    main()
