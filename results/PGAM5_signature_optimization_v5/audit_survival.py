"""Independent audit of exported survival input, risk sets and test families."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm
from statsmodels.duration.hazard_regression import PHReg

ROOT = Path(__file__).resolve().parent
R = ROOT/'survival_exploratory'


def main():
    checks = []
    def check(name, condition):
        assert condition, name
        checks.append(name)
    def bh(p):
        p = np.asarray(p)
        order = np.argsort(p)
        q = np.empty(len(p))
        q[order] = np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
        return q
    protocol = json.loads((R/'protocol.json').read_text())
    original = pd.read_csv(ROOT/'EXPLORATORY_TCGA_coefficients_NOT_CELL_FRACTIONS.csv').set_index('patient')
    joined = pd.read_csv(R/'patient_coefficients_survival_join.csv').set_index('patient')
    original = original.loc[joined.index]
    check('371 unique patients',len(joined)==371 and joined.index.is_unique)
    check('unaltered exploratory target coefficients',np.allclose(joined.coefficient,
          original.exploratory_PGAM5_RNA_macrophage_coefficient,atol=1e-15))
    median = original.exploratory_PGAM5_RNA_macrophage_coefficient.median()
    check('same global outcome-blind median',abs(median-protocol['global_median'])<1e-15)
    check('exact group membership',np.array_equal(joined.group.eq('High'),joined.coefficient.gt(median)))
    check('all coefficient validation flags remain false',not joined.usable_as_cellular_infiltration.any())
    source = pd.read_csv(R/'TCGA_CDR_mirror_survival.csv.gz')
    source['patient'] = source['sample'].str[:12]
    source = source[source.patient.isin(joined.index)]
    fields = ['OS','OS.time','PFI','PFI.time','DFI','DFI.time','DSS','DSS.time']
    check('survival aliases agree in all eight fields',source.groupby('patient')[fields].nunique(dropna=False).le(1).all().all())
    source = source[source['sample'].str[13:15].eq('01')].set_index('patient').reindex(joined.index)
    check('joined survival values exactly reproduce source',np.allclose(joined[fields],source[fields],equal_nan=True))
    stats = pd.read_csv(R/'survival_statistics.csv').set_index('endpoint')
    cov = pd.read_csv(R/'Cox_all_models.csv')
    risks = pd.read_csv(R/'KM_numbers_at_risk.csv')
    curves = pd.read_csv(R/'KM_curve_values.csv')
    for ep in ['OS','PFI']:
        valid = np.isfinite(joined[ep+'.time']) & joined[ep+'.time'].gt(0) & joined[ep].isin([0,1])
        data = pd.read_csv(R/f'{ep}_analysis_patients.csv').set_index('patient')
        check(ep+' endpoint exact eligible patient set',set(data.index)==set(joined.index[valid]))
        check(ep+' global groups preserved',data.group.eq(joined.loc[data.index,'group']).all())
        check(ep+' sample and event counts',len(data)==stats.loc[ep,'patients'] and data.event.sum()==stats.loc[ep,'events'])
        oe,v = 0.,0.
        for day in np.unique(data.loc[data.event.eq(1),'time_days']):
            at = data[data.time_days.ge(day)]
            deaths = data[data.time_days.eq(day)&data.event.eq(1)]
            nh,nl = int(at.high.sum()),int((1-at.high).sum())
            n,d = nh+nl,len(deaths)
            oe += deaths.high.sum()-d*nh/n
            v += nh*nl*d*(n-d)/(n*n*(n-1)) if n>1 else 0.
        check(ep+' exported log-rank p from raw patients',np.isclose(chi2.sf(oe**2/v,1),stats.loc[ep,'logrank_p'],atol=1e-12))
        for group in ['High','Low']:
            d = data[data.group.eq(group)]
            z = curves[curves.endpoint.eq(ep)&curves.group.eq(group)].sort_values('time_days')
            survival = 1.
            for day in np.unique(d.time_days):
                survival *= 1-int((d.time_days.eq(day)&d.event.eq(1)).sum())/int(d.time_days.ge(day).sum())
                check(ep+' '+group+' saved KM at '+str(day),abs(z.loc[z.time_days.eq(day),'survival'].iloc[0]-survival)<1e-12)
            q = risks[risks.endpoint.eq(ep)&risks.group.eq(group)]
            check(ep+' '+group+' all displayed risk counts',all(
                int(d.time_days.ge(row.month*30.4375).sum())==row.at_risk for row in q.itertuples()))
        for model, q in cov[cov.endpoint.eq(ep)].groupby('model'):
            target = 'high' if model.startswith('high_low') else 'coefficient_per_IQR'
            adjusted = model.endswith('adjusted') and model!='high_low_unadjusted' and model!='continuous_unadjusted'
            stratified = model=='posthoc_continuous_stage_stratified'
            use = data[data.covariates_complete] if adjusted or stratified else data
            variables = [target]+(['age_per_10_years','male'] if adjusted or stratified else [])
            if adjusted: variables += ['stage_II','stage_III_IV']
            sm = PHReg(use.time_days.to_numpy(),use[variables].to_numpy(),status=use.event.to_numpy(),
                       strata=use.stage_group.to_numpy() if stratified else None,ties='efron').fit()
            q = q.set_index('term').loc[variables]
            se = np.sqrt(np.diag(sm.cov_params()))
            p = 2*norm.sf(abs(sm.params/se))
            check(ep+' '+model+' all HR/CI/p reproduced independently',
                np.allclose(q.HR,np.exp(sm.params),rtol=1e-5,atol=1e-7) and
                np.allclose(q.CI95_lower,np.exp(sm.params-1.959963984540054*se),rtol=1e-5) and
                np.allclose(q.CI95_upper,np.exp(sm.params+1.959963984540054*se),rtol=1e-5) and
                np.allclose(q.Wald_p,p,rtol=1e-4,atol=1e-9))
    check('two primary logrank BH values',np.allclose(bh(stats.logrank_p),stats.logrank_BH_q_2_endpoints))
    secondary = cov[cov.target_term & cov.model.isin(['continuous_unadjusted','high_low_adjusted','continuous_adjusted'])]
    check('six secondary Cox BH values',len(secondary)==6 and np.allclose(bh(secondary.Wald_p),secondary.secondary_BH_q_6_target_tests))
    posthoc = cov[cov.target_term & cov.model.eq('posthoc_continuous_stage_stratified')]
    check('two posthoc stratified Cox BH values',len(posthoc)==2 and np.allclose(bh(posthoc.Wald_p),posthoc.posthoc_stage_stratified_BH_q_2))
    result = {'status':'PASS','checks_passed':len(checks),'checks':checks,
              'interpretation':'Reproduction of exploratory coefficient-survival analyses; no cellular abundance or specificity validation'}
    (R/'independent_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))


if __name__ == '__main__':
    main()
