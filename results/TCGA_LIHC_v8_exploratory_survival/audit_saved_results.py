"""Independent numerical audits using only the downloadable package inputs."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import chi2
from scipy.optimize import lsq_linear
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.stats.multitest import multipletests
from guard_solver import fit

R = Path(__file__).resolve().parent
UNITS = ['library_RNA_contribution','equalized_cell_fraction']

def main():
    checks = []
    def check(label,value):
        checks.append({'check':label,'passed':bool(value)})
        assert value,label
    counts = pd.read_csv(R/'TCGA_reference_union_raw_counts.csv.gz',index_col=0)
    cp = pd.read_csv(R/'TCGA_reference_union_CP10k.csv.gz',index_col=0)
    meta = pd.read_csv(R/'TCGA_primary_aliquot_selection.csv').set_index('patient')
    check('unique gene and patient identifiers',counts.index.is_unique and counts.columns.is_unique and meta.index.is_unique)
    check('one primary tumor per patient',len(meta)==371 and meta['sample'].str[13:15].eq('01').all())
    check('finite nonnegative integer raw counts',np.isfinite(counts.to_numpy()).all() and (counts.to_numpy()>=0).all() and (counts.to_numpy()==np.floor(counts.to_numpy())).all())
    expected = counts.to_numpy()/meta.loc[counts.columns,'full_library_total_counts'].to_numpy()*10000
    check('CP10k uses full-library denominator',np.allclose(expected,cp.loc[counts.index,counts.columns].to_numpy(),atol=1e-12,rtol=1e-12))
    frozen = json.loads((R/'frozen_grouping.json').read_text())
    cuts = {x['unit']:x['global_median'] for x in frozen['units']}
    stats = pd.read_csv(R/'survival_statistics.csv')
    curves = pd.read_csv(R/'KM_curve_values.csv')
    risk = pd.read_csv(R/'KM_numbers_at_risk.csv')
    exclusions = pd.read_csv(R/'excluded_patients.csv')
    risksets, kkt = [], []
    for unit in UNITS:
        original = np.load(R/'reference'/unit/f'{unit}_model.npz')
        used = np.load(R/unit/'used_reference_model.npz')
        keep = pd.Index(original['genes']).isin(counts.index)
        check(unit+' exact measured-row projection',np.array_equal(used['genes'],original['genes'][keep]) and np.array_equal(used['A'],original['A'][keep]) and np.array_equal(used['scale'],original['scale'][keep]))
        coeff = pd.read_csv(R/unit/'patient_coefficients.csv')
        c = coeff[['component_'+str(g) for g in used['groups']]].to_numpy()
        check(unit+' nonnegative sum-to-one constraints',np.min(c)>=-1e-10 and np.max(np.abs(c.sum(axis=1)-1))<1e-7)
        check(unit+' saved target and macrophage sums',np.allclose(c[:,used['target_mask']].sum(axis=1),coeff.target_coefficient,atol=1e-12) and np.allclose(c[:,used['macrophage_mask']].sum(axis=1),coeff.total_macrophage_coefficient,atol=1e-12))
        check(unit+' target-only PGAM5 RNA budget',np.max(coeff.target_PGAM5_scaled_contribution-coeff.observed_PGAM5_scaled)<=1e-7)
        check(unit+' fixed full-cohort median grouping',abs(float(coeff.target_coefficient.median())-cuts[unit])<1e-14 and np.array_equal(coeff.group,np.where(coeff.target_coefficient>cuts[unit],'High','Low')))
        pg = int(np.flatnonzero(used['genes']=='PGAM5')[0])
        cases = pd.read_csv(R/unit/'solver_case_manifest.csv')
        for _,case in cases.iterrows():
            z = np.load(R/unit/case['file'])
            v = cp.loc[used['genes'],case.patient].to_numpy()/used['scale']
            check(unit+' '+case.patient+' actual query input',np.allclose(z['b'],v,atol=1e-12,rtol=1e-12))
            refit,weights,residual = fit(z['A'],z['b'],True,pg,z['target_mask'])
            check(unit+' '+case.patient+' solver refit',np.allclose(refit,z['c'],atol=1e-8,rtol=0))
            A,b,c,w = z['A'],z['b'],z['c'],z['weights']
            H = A.T@(w[:,None]*A)/len(b)+.001*np.eye(A.shape[1])
            rhs = A.T@(w*b)/len(b)
            grad = H@c-rhs
            p = A[pg]*z['target_mask']
            budget_active = b[pg]-p@c < 1e-7
            active_lower = np.flatnonzero(c<1e-7)
            columns = [np.ones(len(c))]
            lower = [-np.inf]
            if budget_active: columns.append(p); lower.append(0.)
            for j in active_lower:
                col = np.zeros(len(c)); col[j] = -1
                columns.append(col); lower.append(0.)
            D = np.stack(columns,axis=1)
            multipliers = lsq_linear(D,-grad,bounds=(lower,np.full(len(lower),np.inf)),method='bvls',tol=1e-12).x
            error = float(np.max(np.abs(grad+D@multipliers)))
            check(unit+' '+case.patient+' final-weight convex QP KKT',error<1e-4)
            kkt.append({'unit':unit,'patient':case.patient,'stationarity_max_abs_error':error,
                        'target_budget_active':bool(budget_active),'checked_problem':'last fixed-weight ridge QP, not global IRLS optimality'})
        merged = pd.read_csv(R/unit/'patient_coefficients_survival_join.csv')
        check(unit+' patient join did not expand',len(merged)==371 and merged.patient.is_unique)
        check(unit+' join preserved grouping and coefficients',np.array_equal(merged.patient,coeff.patient) and np.array_equal(merged.group,coeff.group) and np.allclose(merged.target_coefficient,coeff.target_coefficient,atol=1e-12))
        for ep in ['OS','PFI']:
            stat = stats[(stats.unit==unit)&(stats.endpoint==ep)].iloc[0]
            time = pd.to_numeric(merged[ep+'.time'],errors='coerce')
            event = pd.to_numeric(merged[ep],errors='coerce')
            valid = np.isfinite(time)&time.gt(0)&event.isin([0,1])
            t = pd.read_csv(R/unit/f'{ep}_analysis_patients.csv')
            check(unit+' '+ep+' endpoint inclusion',set(t.patient)==set(merged.loc[valid,'patient']))
            ex = exclusions[(exclusions.unit==unit)&(exclusions.endpoint==ep)]
            check(unit+' '+ep+' exhaustive exclusions',set(ex.patient)==set(merged.loc[~valid,'patient']) and len(t)+len(ex)==371)
            check(unit+' '+ep+' counts and events',len(t)==stat.patients and int(t.event.sum())==stat.events and int(t.group.eq('High').sum())==stat.high_patients)
            OminusE,V = 0.,0.
            for day in np.unique(t.loc[t.event.eq(1),'time_days']):
                at = t.time_days.ge(day)
                events = t.time_days.eq(day)&t.event.eq(1)
                nh = int((at&t.group.eq('High')).sum()); nl = int((at&t.group.eq('Low')).sum())
                dh = int((events&t.group.eq('High')).sum()); dl = int((events&t.group.eq('Low')).sum())
                n,d = nh+nl,dh+dl
                expected_high = d*nh/n
                variance = nh*nl*d*(n-d)/(n*n*(n-1)) if n>1 else 0
                OminusE += dh-expected_high; V += variance
                risksets.append({'unit':unit,'endpoint':ep,'day':day,'risk_high':nh,'risk_low':nl,
                                 'events_high':dh,'events_low':dl,'expected_high':expected_high,'variance':variance})
            chisq = OminusE**2/V
            check(unit+' '+ep+' independent logrank risk-set statistic',np.isclose(chisq,stat.logrank_chi_square,rtol=1e-10,atol=1e-10))
            check(unit+' '+ep+' independent logrank p',np.isclose(chi2.sf(chisq,1),stat.logrank_p,rtol=1e-10,atol=1e-12))
            for group in ['High','Low']:
                g = t[t.group.eq(group)]
                cur = curves[(curves.unit==unit)&(curves.endpoint==ep)&(curves.group==group)]
                S,greenwood = 1.,0.
                for day in np.unique(g.time_days):
                    n = int(g.time_days.ge(day).sum())
                    d = int((g.time_days.eq(day)&g.event.eq(1)).sum())
                    S *= 1-d/n
                    if d and n>d: greenwood += d/(n*(n-d))
                    step = cur[cur.time_days.eq(day)].iloc[0]
                    check(unit+' '+ep+' '+group+' KM step '+str(day),abs(S-step.survival)<1e-11)
                    if 0<S<1:
                        se = np.sqrt(greenwood)/abs(np.log(S)); loglog = np.log(-np.log(S))
                        ci = [np.exp(-np.exp(loglog+1.959963984540054*se)),np.exp(-np.exp(loglog-1.959963984540054*se))]
                        check(unit+' '+ep+' '+group+' Greenwood CI '+str(day),np.allclose(ci,[step.CI95_lower,step.CI95_upper],atol=1e-10))
                rows = risk[(risk.unit==unit)&(risk.endpoint==ep)&(risk.group==group)]
                check(unit+' '+ep+' '+group+' displayed risk table',all(int(g.time_days.ge(row.month*30.4375).sum())==row.at_risk for _,row in rows.iterrows()))
            cox = PHReg(t.time_days.to_numpy(),t.high.to_numpy()[:,None],status=t.event.to_numpy(),ties='efron').fit()
            HR = float(np.exp(cox.params[0]));ci = np.exp(cox.conf_int()[0])
            check(unit+' '+ep+' independent Cox HR/CI',np.isclose(HR,stat.descriptive_unadjusted_HR_high_low,rtol=1e-5) and np.allclose(ci,[stat.HR_CI95_lower,stat.HR_CI95_upper],rtol=1e-5))
        selected = stats.unit.eq(unit)
        check(unit+' two-endpoint BH correction',np.allclose(multipletests(stats.loc[selected,'logrank_p'],method='fdr_bh')[1],stats.loc[selected,'BH_q_within_unit'],atol=1e-12))
    check('four-test BH correction retained',np.allclose(multipletests(stats.logrank_p,method='fdr_bh')[1],stats.BH_q_joint_four_tests,atol=1e-12))
    check('strict PFS was not invented',not json.loads((R/'protocol.json').read_text())['PFS_available'] and set(stats.endpoint)=={'OS','PFI'})
    pd.DataFrame(risksets).to_csv(R/'independent_logrank_risk_sets.csv',index=False)
    pd.DataFrame(kkt).to_csv(R/'independent_solver_KKT_checks.csv',index=False)
    pd.DataFrame(checks).to_csv(R/'independent_audit_checks.csv',index=False)
    result = {'status':'PASS','checks':len(checks),'all_passed':all(x['passed'] for x in checks),
              'solver_audit_cases':len(kkt),'maximum_KKT_stationarity_error':max(x['stationarity_max_abs_error'] for x in kkt),
              'new_biological_validation_or_cell_abundance_calibration':False}
    (R/'independent_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__ == '__main__': main()
