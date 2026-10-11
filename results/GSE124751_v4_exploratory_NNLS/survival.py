"""Patient-linked OS and source-defined DFS with fixed prior v4 coefficients."""
from pathlib import Path
import json, importlib.metadata as md
import numpy as np
import pandas as pd
from scipy.stats import chi2
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, proportional_hazard_test
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
INPUT=ROOT/'inputs'
TARGET='TAM_PGAM5_detected'
COLORS={'High':'#B64F4A','Low':'#2B78A5'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42})


def write(df,name):
    df.to_csv(ROOT/name,index=False,lineterminator='\n')


def link():
    meta=pd.read_csv(ROOT/'sample_coefficients_and_metadata.csv')
    clinical=pd.read_csv(INPUT/'lihc_amc_prv_data_clinical_patient.txt',sep='\t',comment='#')
    sample=pd.read_csv(INPUT/'lihc_amc_prv_data_clinical_sample.txt',sep='\t',comment='#')
    assert clinical.PATIENT_ID.is_unique and sample.PATIENT_ID.is_unique and sample.SAMPLE_ID.is_unique
    assert set(clinical.PATIENT_ID)==set(sample.PATIENT_ID) and len(clinical)==231
    meta['PATIENT_ID']=meta.sample_id.str.removesuffix('T')
    assert meta.PATIENT_ID.is_unique and meta.sample_id.str.endswith('T').all()
    patients=meta.merge(clinical,on='PATIENT_ID',how='left',validate='one_to_one').merge(
        sample,on='PATIENT_ID',how='left',validate='one_to_one')
    earlier=pd.read_csv(ROOT/'earlier_GEO_exact_title_link_candidates.csv')
    patients=patients.merge(earlier[['sample_id','earlier_GSM','earlier_sample_title']],
                            on='sample_id',how='left',validate='one_to_one')
    patients['exact_ID_link']=patients.AGE.notna() & patients.SAMPLE_ID.eq(patients.PATIENT_ID) & patients.earlier_GSM.notna()
    patients['age_agreement']=patients.age.eq(patients.AGE)
    patients['link_exclusion_reason']=''
    patients.loc[~patients.exact_ID_link,'link_exclusion_reason']='No exact clinical and earlier GEO tumor-ID link'
    patients.loc[patients.exact_ID_link & ~patients.age_agreement,'link_exclusion_reason']='GEO age differs from clinical AGE'
    patients['main_identity_eligible']=patients.exact_ID_link & patients.age_agreement
    patients['source_sample_diagnosis_conflict']=patients.exact_ID_link & patients.CANCER_TYPE_DETAILED.ne('Hepatocellular Carcinoma')
    for endpoint,labels in [('OS',{'0:LIVING':0,'1:DECEASED':1}),
                            ('DFS',{'0:DiseaseFree':0,'1:Recurred':1})]:
        patients[f'{endpoint}_event']=patients[f'{endpoint}_STATUS'].map(labels)
        patients[f'{endpoint}_time_months']=pd.to_numeric(patients[f'{endpoint}_MONTHS'],errors='coerce')
    patients['DFS_later_than_OS_living_followup']=patients.OS_event.eq(0) & patients.DFS_time_months.gt(patients.OS_time_months)
    # Preserve the preclinical all-206 cutoff and scale, never recalibrate on follow-up.
    cutoff=float(meta[TARGET].median())
    assert cutoff==0 and np.array_equal(meta.group.eq('High'),meta[TARGET].gt(cutoff))
    patients['coefficient_per_SD']=(patients[TARGET]-meta[TARGET].mean())/meta[TARGET].std(ddof=1)
    patients['age_per_10_years']=patients.AGE/10.
    patients['high']=patients.group.eq('High').astype(int)
    write(patients,'patient_coefficients_clinical.csv')
    write(patients[['GSM','sample_id','GEO_sample_title','PATIENT_ID','SAMPLE_ID','earlier_GSM',
                    'earlier_sample_title','age','AGE','exact_ID_link','age_agreement','main_identity_eligible',
                    'link_exclusion_reason','source_sample_diagnosis_conflict',
                    'DFS_later_than_OS_living_followup']], 'clinical_link_audit.csv')
    return patients,cutoff


def endpoint_frame(patients,endpoint,strict_age=True,drop_living_discordance=False):
    reason=patients.link_exclusion_reason.copy() if strict_age else pd.Series('',index=patients.index)
    if not strict_age:
        reason.loc[~patients.exact_ID_link]='No exact clinical and earlier GEO tumor-ID link'
    time=patients[f'{endpoint}_time_months']
    event=patients[f'{endpoint}_event']
    reason.loc[reason.eq('') & (~np.isfinite(time) | time.le(0) | ~event.isin([0,1]))]='Invalid time or event'
    if endpoint=='DFS':
        reason.loc[reason.eq('') & patients.OS_event.eq(1) & time.gt(patients.OS_time_months)]='DFS later than known death'
        if drop_living_discordance:
            reason.loc[reason.eq('') & patients.DFS_later_than_OS_living_followup]='Sensitivity: DFS exceeds living OS follow-up'
    keep=reason.eq('')
    data=patients.loc[keep].copy()
    data['time_months'],data['event']=time.loc[keep],event.loc[keep].astype(int)
    excluded=patients.loc[~keep,['GSM','sample_id','PATIENT_ID','group',f'{endpoint}_time_months',f'{endpoint}_event']].copy()
    excluded['endpoint'],excluded['exclusion_reason']=endpoint,reason.loc[~keep]
    return data,excluded


def compare(data,endpoint,analysis,cutoff):
    high,low=data[data.high.eq(1)],data[data.high.eq(0)]
    lr=logrank_test(high.time_months,low.time_months,event_observed_A=high.event,event_observed_B=low.event)
    risksets,oe,var=[],0.,0.
    for time in sorted(data.loc[data.event.eq(1),'time_months'].unique()):
        nh,nl=int(high.time_months.ge(time).sum()),int(low.time_months.ge(time).sum())
        dh=int((high.time_months.eq(time)&high.event.eq(1)).sum())
        dl=int((low.time_months.eq(time)&low.event.eq(1)).sum())
        n,d=nh+nl,dh+dl
        expected=d*nh/n
        v=nh*nl*d*(n-d)/(n*n*(n-1)) if n>1 else 0.
        oe+=dh-expected;var+=v
        risksets.append({'analysis':analysis,'endpoint':endpoint,'time_months':time,'risk_high':nh,
                         'risk_low':nl,'events_high':dh,'events_low':dl,'expected_high':expected,'variance':v})
    statistic=oe**2/var
    assert np.isclose(statistic,lr.test_statistic,atol=1e-10)
    assert np.isclose(chi2.sf(statistic,1),lr.p_value,atol=1e-12)
    row={'analysis':analysis,'endpoint':endpoint,'patients':len(data),'events':int(data.event.sum()),
         'high_patients':len(high),'low_patients':len(low),'high_events':int(high.event.sum()),
         'low_events':int(low.event.sum()),'median_cutoff_all206':cutoff,
         'logrank_p':float(lr.p_value),'logrank_chi_square':float(statistic)}
    design=data[['time_months','event','high']]
    cox=CoxPHFitter().fit(design,'time_months','event')
    cr=cox.summary.loc['high']
    row.update({'HR_high_vs_low':cr['exp(coef)'],'HR_lower95':cr['exp(coef) lower 95%'],
                'HR_upper95':cr['exp(coef) upper 95%'],'Cox_p':cr.p,
                'PH_rank_p':proportional_hazard_test(cox,design,time_transform='rank').summary.loc['high','p']})
    return row,risksets


def models(data,endpoint):
    result=[]
    for model,covariate,adjusted in [('binary','high',False),('continuous_per_SD','coefficient_per_SD',False),
                                    ('binary_age_adjusted','high',True),('continuous_per_SD_age_adjusted','coefficient_per_SD',True)]:
        cols=['time_months','event',covariate]+(['age_per_10_years'] if adjusted else [])
        design=data[cols].copy()
        cox=CoxPHFitter().fit(design,'time_months','event')
        ph=proportional_hazard_test(cox,design,time_transform='rank').summary
        for name,cr in cox.summary.iterrows():
            result.append({'endpoint':endpoint,'model':model,'covariate':name,'patients':len(data),
                           'events':int(data.event.sum()),'HR':cr['exp(coef)'],
                           'HR_lower95':cr['exp(coef) lower 95%'],'HR_upper95':cr['exp(coef) upper 95%'],
                           'p':cr.p,'PH_rank_p':ph.loc[name,'p']})
    return result


def plot(frames,stats):
    ticks=np.arange(0,91,15)
    fig=plt.figure(figsize=(13.5,8.6))
    grid=fig.add_gridspec(2,2,height_ratios=[4.7,1],hspace=.43,wspace=.25)
    curves,risks=[] , []
    for col,endpoint in enumerate(['OS','DFS']):
        ax,tab=fig.add_subplot(grid[0,col]),fig.add_subplot(grid[1,col])
        data=frames[endpoint];row=stats[stats.endpoint.eq(endpoint)].iloc[0]
        for group in ['High','Low']:
            q=data[data.group.eq(group)]
            km=KaplanMeierFitter().fit(q.time_months,q.event,label=group)
            sf,ci=km.survival_function_,km.confidence_interval_
            ax.step(sf.index,sf.iloc[:,0],where='post',color=COLORS[group],lw=2,
                    label=f'{group}: '+('nonzero' if group=='High' else 'zero')+f' coefficient (n={len(q)})')
            ax.fill_between(sf.index,ci.iloc[:,0],ci.iloc[:,1],step='post',color=COLORS[group],alpha=.13)
            censor=q.loc[q.event.eq(0),'time_months']
            ax.plot(censor,km.predict(censor),'+',ms=5,mew=.8,color=COLORS[group])
            z=pd.concat([sf,ci],axis=1).reset_index()
            z.columns=['time_months','survival','CI95_lower','CI95_upper'];z['endpoint'],z['group']=endpoint,group
            curves.append(z)
            for t in ticks:
                risks.append({'endpoint':endpoint,'group':group,'time_months':int(t),'at_risk':int(q.time_months.ge(t).sum())})
            stats.loc[stats.endpoint.eq(endpoint),f'median_{group}_months']=km.median_survival_time_
        ax.set(title='Overall survival (OS)' if endpoint=='OS' else 'Disease-free survival (source DFS)',
               xlabel='Time (months)',ylabel='Survival probability' if endpoint=='OS' else 'Source DFS probability',
               xlim=(0,90),ylim=(0,1.02),xticks=ticks)
        ax.text(.04,.065,f'Log-rank p = {row.logrank_p:.4g}\nBH q (OS / DFS) = {row.logrank_BH_q:.4g}\n'
                f'HR high / low = {row.HR_high_vs_low:.2f} ({row.HR_lower95:.2f}–{row.HR_upper95:.2f})',
                transform=ax.transAxes,fontsize=10,bbox={'facecolor':'white','edgecolor':'none','alpha':.85})
        ax.legend(loc='lower left',bbox_to_anchor=(0,.24),frameon=False,fontsize=9)
        ax.spines[['top','right']].set_visible(False)
        tab.axis('off')
        counts=[[str(int(((data.group==group)&data.time_months.ge(t)).sum())) for t in ticks] for group in ['High','Low']]
        table=tab.table(cellText=counts,rowLabels=['High','Low'],colLabels=[str(t) for t in ticks],loc='center',cellLoc='center')
        table.auto_set_font_size(False);table.set_fontsize(9);table.scale(1,1.3)
        tab.set_title('Number at risk (start of time point)',fontsize=10)
    fig.suptitle('GSE124751 HCC tumor tissues: exploratory v4 NNLS survival analysis',fontsize=15,y=.975)
    fig.text(.5,.924,'196 patients with exact IDs and concordant age; High > 0 vs Low = 0; cutoff fixed on all 206 samples',ha='center',fontsize=10.5)
    fig.text(.5,.035,'763 / 850 genes; 17 fixed components. Coefficients are exploratory, not validated cell fractions.\n'
             'DFS uses the source recurrence labels; it is not PFS. Age-conflict and follow-up checks are reported separately.',ha='center',fontsize=10)
    fig.subplots_adjust(left=.075,right=.975,bottom=.185,top=.855)
    for ext in ['png','pdf']:fig.savefig(ROOT/f'KM_OS_DFS.{ext}',dpi=220,facecolor='white')
    plt.close(fig)
    write(pd.concat(curves,ignore_index=True),'KM_curve_values.csv')
    write(pd.DataFrame(risks),'KM_numbers_at_risk.csv')


def main():
    patients,cutoff=link()
    stats,frames,exclusions,risks,cox=[],{},[],[],[]
    for ep in ['OS','DFS']:
        data,excluded=endpoint_frame(patients,ep)
        assert len(data)==196 and int(data.high.sum())==61
        frames[ep]=data
        row,risk=compare(data,ep,'main_age_concordant',cutoff)
        stats.append(row);risks.extend(risk);cox.extend(models(data,ep));exclusions.append(excluded)
        write(data,f'{ep}_analysis_patients.csv')
    stats=pd.DataFrame(stats)
    stats['logrank_BH_q']=multipletests(stats.logrank_p,method='fdr_bh')[1]
    plot(frames,stats)
    write(stats,'survival_statistics.csv')
    cox=pd.DataFrame(cox)
    target=cox.covariate.isin(['high','coefficient_per_SD'])
    assert int(target.sum())==8
    cox.loc[target,'target_test_BH_q_across_8']=multipletests(cox.loc[target,'p'],method='fdr_bh')[1]
    write(cox,'Cox_binary_continuous_age_adjusted.csv')
    write(pd.concat(exclusions,ignore_index=True),'endpoint_excluded_patients.csv')
    sensitivity=[]
    for ep in ['OS','DFS']:
        data,_=endpoint_frame(patients,ep,strict_age=False)
        row,risk=compare(data,ep,'all198_exact_IDs_including_age_conflicts',cutoff)
        sensitivity.append(row);risks.extend(risk)
    data,_=endpoint_frame(patients,'DFS',drop_living_discordance=True)
    row,risk=compare(data,'DFS','exclude_living_followup_discordance',cutoff)
    sensitivity.append(row);risks.extend(risk)
    for variant in ['zero_preserving_log2','unshifted_log2']:
        score=pd.read_csv(ROOT/f'coefficients_{variant}.csv',index_col=0)[TARGET]
        alternative=patients.copy()
        alternative['high']=alternative.sample_id.map(score).gt(score.median()).astype(int)
        alternative['group']=np.where(alternative.high,'High','Low')
        for ep in ['OS','DFS']:
            data,_=endpoint_frame(alternative,ep)
            row,_=compare(data,ep,'expression_'+variant,float(score.median()))
            sensitivity.append(row)
    write(pd.DataFrame(sensitivity),'survival_sensitivity.csv')
    write(pd.DataFrame(risks),'independent_logrank_risk_sets.csv')
    summary={'source_clinical_patients':231,'exact_ID_linked':int(patients.exact_ID_link.sum()),
             'age_concordant':int(patients.main_identity_eligible.sum()),
             'age_conflicts':patients.loc[patients.exact_ID_link & ~patients.age_agreement,'sample_id'].tolist(),
             'no_exact_clinical_ID':patients.loc[~patients.exact_ID_link,'sample_id'].tolist(),
             'main_high':61,'main_low':135,'main_OS_patients':196,'main_DFS_patients':196,
             'DFS_later_than_living_OS_followup':patients.loc[patients.DFS_later_than_OS_living_followup,'sample_id'].tolist(),
             'source_LIAD_label_conflicts':int(patients.source_sample_diagnosis_conflict.sum()),
             'PFS_available':False,'survival_analysis_performed':True,
             'source_patient_sha256':'629c2d2994250b55f069ddd4a7a5900f5ac0f3ab1d9da37b96686f2c8e8729c1',
             'source_sample_sha256':'5ecf933d8834175841f59efa2f8b9f251531877fd89a58a1da81b8042762ac49',
             'packages':{n:md.version(n) for n in ['numpy','pandas','scipy','matplotlib','lifelines','statsmodels']}}
    (ROOT/'survival_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    projection_summary=json.loads((ROOT/'summary.json').read_text())
    projection_summary['survival_analysis_performed']=True
    projection_summary['clinical_reason']='User-uploaded clinical RTF decoded; exact SHA256 match to both official pinned clinical payloads.'
    (ROOT/'summary.json').write_text(json.dumps(projection_summary,indent=2)+'\n')
    print(stats.to_string(index=False))
    print(pd.DataFrame(sensitivity).to_string(index=False))


if __name__=='__main__':
    main()
