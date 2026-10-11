"""Independent clinical linkage, KM, tied logrank and Efron-Cox verification."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import chi2, norm

ROOT=Path(__file__).resolve().parent
checks=[]


def check(name,condition):
    checks.append({'check':name,'passed':bool(condition)})
    assert condition,name


def bh(p):
    p=np.array(p,dtype=float);order=np.argsort(p)
    adjusted=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    out=np.empty(len(p));out[order]=np.minimum(1,adjusted)
    return out


def univariate_cox(data):
    time=data.time_months.to_numpy();event=data.event.to_numpy();x=data.high.to_numpy()
    times=np.unique(time[event==1])
    def score_information(beta):
        score,information=0.,0.
        weights=np.exp(beta*x)
        for t in times:
            risk=time>=t;deaths=(time==t)&(event==1);d=deaths.sum()
            r0=weights[risk].sum();r1=(weights[risk]*x[risk]).sum()
            d0=weights[deaths].sum();d1=(weights[deaths]*x[deaths]).sum()
            score+=x[deaths].sum()
            for fraction in np.arange(d)/d:
                denominator=r0-fraction*d0
                first=(r1-fraction*d1)/denominator
                score-=first
                information+=first-first**2  # x is binary, so x**2 == x
        return score,information
    beta=brentq(lambda b:score_information(b)[0],-20,20,xtol=1e-12)
    information=score_information(beta)[1];se=1/np.sqrt(information)
    return np.exp(beta),np.exp(beta-1.959963984540054*se),np.exp(beta+1.959963984540054*se),2*norm.sf(abs(beta/se))


def main():
    patient_source=pd.read_csv(ROOT/'inputs/lihc_amc_prv_data_clinical_patient.txt',sep='\t',comment='#').set_index('PATIENT_ID')
    sample_source=pd.read_csv(ROOT/'inputs/lihc_amc_prv_data_clinical_sample.txt',sep='\t',comment='#').set_index('PATIENT_ID')
    for kind,sha in [('patient','629c2d2994250b55f069ddd4a7a5900f5ac0f3ab1d9da37b96686f2c8e8729c1'),
                     ('sample','5ecf933d8834175841f59efa2f8b9f251531877fd89a58a1da81b8042762ac49')]:
        check('Official clinical SHA256 '+kind,
              hashlib.sha256((ROOT/'inputs'/f'lihc_amc_prv_data_clinical_{kind}.txt').read_bytes()).hexdigest()==sha)
    original=pd.read_csv(ROOT/'sample_coefficients_and_metadata.csv').set_index('sample_id')
    for row in pd.read_csv(ROOT/'preclinical_frozen_coefficients.csv').itertuples():
        check('Preclinical frozen file: '+row.file,
              hashlib.sha256((ROOT/row.file).read_bytes()).hexdigest()==row.sha256)
    link=pd.read_csv(ROOT/'clinical_link_audit.csv').set_index('sample_id')
    stats=pd.read_csv(ROOT/'survival_statistics.csv')
    risk=pd.read_csv(ROOT/'KM_numbers_at_risk.csv')
    curves=pd.read_csv(ROOT/'KM_curve_values.csv')
    exclusions=pd.read_csv(ROOT/'endpoint_excluded_patients.csv')
    check('231 clinical patients and unique source IDs',len(patient_source)==231 and patient_source.index.is_unique and sample_source.index.is_unique)
    check('206 audited tumor expression sample IDs',len(original)==206 and original.index.is_unique)
    check('Strict IDs and age counts',link.exact_ID_link.sum()==198 and link.main_identity_eligible.sum()==196)
    check('Two specified age conflicts',set(link.index[link.exact_ID_link & ~link.age_agreement])=={'H080196T','H080204T'})
    check('Eight L/H aliases excluded from clinical matching',len(link[~link.exact_ID_link])==8 and link.index[~link.exact_ID_link].str.startswith('L').all())
    check('196 main records retain original median split',original.TAM_PGAM5_detected.median()==0)
    rebuilt=[]
    for ep,labels in [('OS',{'0:LIVING':0,'1:DECEASED':1}),('DFS',{'0:DiseaseFree':0,'1:Recurred':1})]:
        data=pd.read_csv(ROOT/f'{ep}_analysis_patients.csv')
        row=stats[stats.endpoint.eq(ep)].iloc[0]
        source=patient_source.loc[data.PATIENT_ID]
        check(ep+' exact full ID correspondence',np.array_equal(data.sample_id.str.removesuffix('T'),data.PATIENT_ID))
        check(ep+' clinical sample ID agrees',np.array_equal(sample_source.loc[data.PATIENT_ID].SAMPLE_ID,data.PATIENT_ID))
        check(ep+' independent age check',np.array_equal(original.loc[data.sample_id,'age'],source.AGE))
        check(ep+' strict earlier GEO full tumor-title match',data.earlier_sample_title.eq('Tumor ('+data.sample_id+')').all())
        check(ep+' exact unmodified source time and event',np.allclose(data.time_months,source[f'{ep}_MONTHS']) and np.array_equal(data.event,source[f'{ep}_STATUS'].map(labels)))
        check(ep+' groups and coefficients not retuned',np.array_equal(data.high,original.loc[data.sample_id].TAM_PGAM5_detected.gt(0).astype(int)) and np.allclose(data.TAM_PGAM5_detected,original.loc[data.sample_id].TAM_PGAM5_detected))
        check(ep+' all-206 coefficient scaling retained',np.allclose(data.coefficient_per_SD,
              (original.loc[data.sample_id].TAM_PGAM5_detected-original.TAM_PGAM5_detected.mean())/original.TAM_PGAM5_detected.std(ddof=1)))
        check(ep+' correct population and events',len(data)==196 and data.high.sum()==61 and data.event.sum()==(29 if ep=='OS' else 85))
        check(ep+' 10 exclusions are disjoint',len(exclusions[exclusions.endpoint.eq(ep)])==10 and not set(data.sample_id)&set(exclusions.loc[exclusions.endpoint.eq(ep),'sample_id']))
        high=data[data.high==1];low=data[data.high==0];oe,var=0.,0.
        for t in sorted(data.loc[data.event.eq(1),'time_months'].unique()):
            nh=(high.time_months>=t).sum();nl=(low.time_months>=t).sum()
            dh=((high.time_months==t)&(high.event==1)).sum();dl=((low.time_months==t)&(low.event==1)).sum()
            n,d=nh+nl,dh+dl
            oe+=dh-d*nh/n
            if n>1:var+=nh*nl*d*(n-d)/(n*n*(n-1))
        p=chi2.sf(oe*oe/var,1);rebuilt.append(p)
        check(ep+' independently rebuilt tied logrank',np.isclose(p,row.logrank_p,atol=1e-12) and np.isclose(oe*oe/var,row.logrank_chi_square,atol=1e-10))
        hr,lower,upper,cox_p=univariate_cox(data)
        check(ep+' independent Efron-Cox HR and CI',np.allclose([hr,lower,upper],[row.HR_high_vs_low,row.HR_lower95,row.HR_upper95],rtol=1e-5,atol=1e-7))
        check(ep+' independent Efron-Cox Wald p',np.isclose(cox_p,row.Cox_p,atol=1e-6))
        for group in ['High','Low']:
            q=data[data.group==group];saved=curves[(curves.endpoint==ep)&(curves.group==group)]
            survival=1.;values={0.:1.}
            for t in sorted(q.time_months.unique()):
                survival*=1-((q.time_months==t)&(q.event==1)).sum()/(q.time_months>=t).sum()
                values[t]=survival
            check(ep+' '+group+' full independent KM curve',all(np.isclose(row.survival,values[row.time_months],atol=1e-12) for row in saved.itertuples()))
            check(ep+' '+group+' CI within bounds',saved.CI95_lower.ge(0).all() and saved.CI95_upper.le(1).all() and saved.CI95_lower.le(saved.survival+1e-12).all() and saved.CI95_upper.ge(saved.survival-1e-12).all())
            r=risk[(risk.endpoint==ep)&(risk.group==group)]
            check(ep+' '+group+' complete plotted risk table',all(row.at_risk==(q.time_months>=row.time_months).sum() for row in r.itertuples()))
    check('Primary two-endpoint BH independently rebuilt',np.allclose(stats.logrank_BH_q,bh(rebuilt),atol=1e-12))
    cox=pd.read_csv(ROOT/'Cox_binary_continuous_age_adjusted.csv')
    target=cox.covariate.isin(['high','coefficient_per_SD'])
    check('Eight target Cox tests BH independently rebuilt',target.sum()==8 and np.allclose(cox.loc[target,'target_test_BH_q_across_8'],bh(cox.loc[target,'p']),atol=1e-12))
    sensitivity=pd.read_csv(ROOT/'survival_sensitivity.csv')
    full=sensitivity[sensitivity.analysis=='all198_exact_IDs_including_age_conflicts']
    check('Both 198-person ID sensitivities are present',len(full)==2 and full.patients.eq(198).all())
    drop=sensitivity[sensitivity.analysis=='exclude_living_followup_discordance']
    check('Living follow-up sensitivity excludes exactly one DFS patient',len(drop)==1 and drop.iloc[0].patients==195 and drop.iloc[0].events==84)
    check('Only H104572T has living follow-up discordance',set(link.index[link.DFS_later_than_OS_living_followup])=={'H104572T'})
    check('No DFS event beyond known death',not ((patient_source.OS_STATUS=='1:DECEASED') & (patient_source.DFS_MONTHS>patient_source.OS_MONTHS)).any())
    check('All erroneous source LIAD labels explicitly flagged',sample_source.ONCOTREE_CODE.eq('LIAD').all() and link.source_sample_diagnosis_conflict.sum()==198)
    check('Original reference and query still frozen',json.loads((ROOT/'audit_results.json').read_text())['passed'])
    check('Endpoints are OS / source DFS, never relabelled PFS',set(stats.endpoint)=={'OS','DFS'} and json.loads((ROOT/'survival_summary.json').read_text())['PFS_available'] is False)
    check('Both PNG and PDF KM exports exist',(ROOT/'KM_OS_DFS.png').stat().st_size>10000 and (ROOT/'KM_OS_DFS.pdf').stat().st_size>10000)
    result={'passed':all(x['passed'] for x in checks),'checks':checks,
            'note':'Source identity and numerical validation do not validate cellular abundance or retrospectively establish a true PFS endpoint.'}
    (ROOT/'survival_audit_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'{len(checks)} independent survival checks passed.')


if __name__=='__main__':
    main()
