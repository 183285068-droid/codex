"""OS/source-RFS methods adapted from the prior frozen NNLS workflow."""
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
         'low_events':int(low.event.sum()),'median_cutoff_all115':cutoff,
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
                                    ('binary_age_sex_adjusted','high',True),('continuous_per_SD_age_sex_adjusted','coefficient_per_SD',True)]:
        cols=['time_months','event',covariate]+(['age_per_10_years','male'] if adjusted else [])
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
    for col,endpoint in enumerate(['OS','RFS']):
        ax,tab=fig.add_subplot(grid[0,col]),fig.add_subplot(grid[1,col])
        data=frames[endpoint];row=stats[stats.endpoint.eq(endpoint)].iloc[0]
        for group in ['High','Low']:
            q=data[data.group.eq(group)]
            km=KaplanMeierFitter().fit(q.time_months,q.event,label=group)
            sf,ci=km.survival_function_,km.confidence_interval_
            ax.step(sf.index,sf.iloc[:,0],where='post',color=COLORS[group],lw=2,
                    label=f'{group} coefficient (n={len(q)})')
            ax.fill_between(sf.index,ci.iloc[:,0],ci.iloc[:,1],step='post',color=COLORS[group],alpha=.13)
            censor=q.loc[q.event.eq(0),'time_months']
            ax.plot(censor,km.predict(censor),'+',ms=5,mew=.8,color=COLORS[group])
            z=pd.concat([sf,ci],axis=1).reset_index()
            z.columns=['time_months','survival','CI95_lower','CI95_upper'];z['endpoint'],z['group']=endpoint,group
            curves.append(z)
            for t in ticks:
                risks.append({'endpoint':endpoint,'group':group,'time_months':int(t),'at_risk':int(q.time_months.ge(t).sum())})
            stats.loc[stats.endpoint.eq(endpoint),f'median_{group}_months']=km.median_survival_time_
        ax.set(title='Overall survival (OS)' if endpoint=='OS' else 'Recurrence-free survival (source RFS)',
               xlabel='Time (months)',ylabel='Survival probability' if endpoint=='OS' else 'Recurrence-free probability',
               xlim=(0,96),ylim=(0,1.02),xticks=ticks)
        ax.text(.04,.065,f'Log-rank p = {row.logrank_p:.4g}\nBH q (OS / RFS) = {row.logrank_BH_q:.4g}\n'
                f'HR high / low = {row.HR_high_vs_low:.2f} ({row.HR_lower95:.2f}–{row.HR_upper95:.2f})',
                transform=ax.transAxes,fontsize=10,bbox={'facecolor':'white','edgecolor':'none','alpha':.85})
        ax.legend(loc='lower left',bbox_to_anchor=(0,.24),frameon=False,fontsize=9)
        ax.spines[['top','right']].set_visible(False)
        tab.axis('off')
        counts=[[str(int(((data.group==group)&data.time_months.ge(t)).sum())) for t in ticks] for group in ['High','Low']]
        table=tab.table(cellText=counts,rowLabels=['High','Low'],colLabels=[str(t) for t in ticks],loc='center',cellLoc='center')
        table.auto_set_font_size(False);table.set_fontsize(9);table.scale(1,1.3)
        tab.set_title('Number at risk (start of time point)',fontsize=10)
    fig.suptitle('GSE76427 HCC tumor tissues: exploratory v4 NNLS survival analysis',fontsize=15,y=.975)
    fig.text(.5,.924,f'115 HCC tumor expression samples; fixed global median cutoff = {CUTOFF:.5g}; OS n=115 / RFS n=107',ha='center',fontsize=10.5)
    fig.text(.5,.035,'755 / 850 genes; 17 fixed components. Coefficients are exploratory, not validated cell fractions.\n'
             'RFS uses GEO recurrence follow-up; seven missing records and one record beyond death are excluded.\n'
             'OS Cox proportional-hazards test p = 0.046; interpret the pooled OS HR cautiously.',ha='center',fontsize=9.5)
    fig.subplots_adjust(left=.075,right=.975,bottom=.185,top=.855)
    for ext in ['png','pdf']:fig.savefig(ROOT/f'KM_OS_RFS.{ext}',dpi=220,facecolor='white')
    plt.close(fig)
    write(pd.concat(curves,ignore_index=True),'KM_curve_values.csv')
    write(pd.DataFrame(risks),'KM_numbers_at_risk.csv')
