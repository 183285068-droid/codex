from pathlib import Path
import json,hashlib,platform,importlib.metadata as md
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter,CoxPHFitter
from lifelines.statistics import logrank_test,proportional_hazard_test
from scipy.stats import chi2,spearmanr
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parent;df=pd.read_csv(R/'patient_scores_survival_join.csv');prep=json.loads((R/'preparation_summary.json').read_text());results=[];exclusions=[];risk=[];audits=[];curve_data=[];fitdata={}
for ep in ['OS','PFI']:
 times=pd.to_numeric(df[ep+'.time'],errors='coerce');events=pd.to_numeric(df[ep],errors='coerce');valid=times.notna()&np.isfinite(times)&times.gt(0)&events.isin([0,1]);t=df.loc[valid,['patient','exploratory_target_fraction','group']].copy();t['time_days']=times[valid].to_numpy();t['event']=events[valid].astype(int).to_numpy();t['high']=t.group.eq('High').astype(int);t.to_csv(R/(ep+'_analysis_patients.csv'),index=False)
 for _,r in df.loc[~valid].iterrows():exclusions.append({'endpoint':ep,'patient':r.patient,'group':r.group,'time_days':r[ep+'.time'],'event':r[ep],'reason':'missing/invalid/nonpositive time or nonbinary/missing event'})
 hi=t[t.high.eq(1)];lo=t[t.high.eq(0)];lr=logrank_test(hi.time_days,lo.time_days,event_observed_A=hi.event,event_observed_B=lo.event)
 # Independent tied-event logrank computation from risk sets, without lifelines helper.
 O_E=0.;V=0.;steps=[]
 for time in sorted(t.loc[t.event.eq(1),'time_days'].unique()):
  yH=int(hi.time_days.ge(time).sum());yL=int(lo.time_days.ge(time).sum());n=yH+yL;dH=int((hi.time_days.eq(time)&hi.event.eq(1)).sum());dL=int((lo.time_days.eq(time)&lo.event.eq(1)).sum());d=dH+dL;ex=d*yH/n;var=yH*yL*d*(n-d)/(n*n*(n-1)) if n>1 else 0.;O_E+=dH-ex;V+=var;steps.append({'endpoint':ep,'time_days':time,'risk_high':yH,'risk_low':yL,'events_high':dH,'events_low':dL,'expected_high':ex,'variance':var})
 statistic=O_E*O_E/V;pmanual=chi2.sf(statistic,1);assert np.isclose(statistic,lr.test_statistic,atol=1e-10) and np.isclose(pmanual,lr.p_value,atol=1e-12);audits.extend(steps)
 cox=CoxPHFitter().fit(t[['time_days','event','high']],duration_col='time_days',event_col='event');cr=cox.summary.loc['high'];ph=proportional_hazard_test(cox,t[['time_days','event','high']],time_transform='rank').summary.loc['high','p']
 out={'endpoint':ep,'endpoint_label':'Overall survival (OS)' if ep=='OS' else 'Progression-free interval (PFI; not PFS)','patients':len(t),'events':int(t.event.sum()),'high_patients':len(hi),'high_events':int(hi.event.sum()),'low_patients':len(lo),'low_events':int(lo.event.sum()),'logrank_chi_square':float(lr.test_statistic),'logrank_p':float(lr.p_value),'independent_logrank_p':float(pmanual),'HR_high_vs_low':float(cr['exp(coef)']),'HR_95CI_lower':float(cr['exp(coef) lower 95%']),'HR_95CI_upper':float(cr['exp(coef) upper 95%']),'Cox_p':float(cr.p),'PH_test_p':float(ph),'estimate_cutpoint':prep['score_median'],'excluded_patients':int((~valid).sum())}
 models={}
 for label,q in [('High',hi),('Low',lo)]:
  fit=KaplanMeierFitter(alpha=.05).fit(q.time_days,event_observed=q.event,label=label+' estimate');models[label]=fit
  # Independent KM product-limit steps check event/censor ties are handled with full pre-event risk set.
  survival=1.
  for day in sorted(q.time_days.unique()):
   n=int(q.time_days.ge(day).sum());d=int((q.time_days.eq(day)&q.event.eq(1)).sum());survival*=1-d/n;assert np.isclose(float(fit.predict(day)),survival,atol=1e-12)
  sf=fit.survival_function_;ci=fit.confidence_interval_;z=pd.concat([sf,ci],axis=1).reset_index();z.columns=['time_days','survival','CI95_lower','CI95_upper'];z['endpoint']=ep;z['group']=label;curve_data.append(z)
  median=fit.median_survival_time_;out['median_'+label+'_days']=float(median) if np.isfinite(median) else None
 results.append(out);fitdata[ep]=(t,models)
pd.DataFrame(exclusions).to_csv(R/'endpoint_excluded_patients.csv',index=False);pd.DataFrame(audits).to_csv(R/'independent_logrank_risk_sets.csv',index=False);pd.concat(curve_data).to_csv(R/'KM_curve_values.csv',index=False);res=pd.DataFrame(results);res['logrank_BH_q_2_endpoints']=multipletests(res.logrank_p,method='fdr_bh')[1];res.to_csv(R/'survival_statistics.csv',index=False)
# Time scale in months for display; computations use original days. Risk counts are at start of each tick.
for ep,(t,models) in fitdata.items():
 rr=res[res.endpoint.eq(ep)].iloc[0];fig=plt.figure(figsize=(8.5,7));gs=fig.add_gridspec(2,1,height_ratios=[4.5,1.25],hspace=.25);ax=fig.add_subplot(gs[0]);tab=fig.add_subplot(gs[1]);colors={'High':'#c94b40','Low':'#2775b6'};maxmonths=float(t.time_days.max()/30.4375);ticks=np.arange(0,np.ceil(maxmonths/24)*24+1,24)
 for label,fit in models.items():
  sf=fit.survival_function_;ci=fit.confidence_interval_;xx=sf.index.to_numpy()/30.4375;ax.step(xx,sf.iloc[:,0],where='post',label=f'{label} estimate (n={int(t.group.eq(label).sum())})',color=colors[label],lw=2);ax.fill_between(xx,ci.iloc[:,0],ci.iloc[:,1],step='post',alpha=.15,color=colors[label]);q=t[t.group.eq(label)];c=q.loc[q.event.eq(0),'time_days'].to_numpy();ax.plot(c/30.4375,fit.predict(c).to_numpy(),'+',color=colors[label],ms=5,mew=.8)
  for m in ticks:risk.append({'endpoint':ep,'group':label,'month':float(m),'at_risk':int(q.time_days.ge(m*30.4375).sum())})
 ax.set(xlabel='Time (months)',ylabel='Survival probability',ylim=(0,1.02),xlim=(0,ticks[-1]),xticks=ticks,title='TCGA-LIHC: '+('OS' if ep=='OS' else 'PFI (not PFS)')+'\nGSE151530-only exploratory reference');ax.legend(loc='lower left',frameon=False);ax.text(.97,.97,f'Log-rank p = {rr.logrank_p:.4g}\nHR high/low = {rr.HR_high_vs_low:.2f} ({rr.HR_95CI_lower:.2f}–{rr.HR_95CI_upper:.2f})',transform=ax.transAxes,ha='right',va='top',fontsize=11);ax.spines[['top','right']].set_visible(False)
 tab.axis('off');rows=[[str(int(t[(t.group==label)&t.time_days.ge(m*30.4375)].shape[0])) for m in ticks] for label in ['High','Low']];table=tab.table(cellText=rows,rowLabels=['High','Low'],colLabels=[str(int(m)) for m in ticks],cellLoc='center',loc='center');table.auto_set_font_size(False);table.set_fontsize(10);table.scale(1,1.4);tab.set_title('Number at risk (start of time point)',fontsize=10,pad=3)
 fig.text(.5,.012,'Exploratory relative mixture estimate; High > 0, Low = 0. Abundance unvalidated.',ha='center',fontsize=9);fig.subplots_adjust(bottom=.08,top=.92,left=.12,right=.97);fig.savefig(R/('KM_'+ep+'.png'),dpi=220);fig.savefig(R/('KM_'+ep+'.pdf'));plt.close(fig)
pd.DataFrame(risk).to_csv(R/'KM_numbers_at_risk.csv',index=False)
score=pd.read_csv(R/'patient_scores.csv');rho,p=spearmanr(score.exploratory_target_fraction,score.PGAM5_log2_CPM_plus1);(R/'estimate_PGAM5_association.json').write_text(json.dumps({'patients':len(score),'Spearman_rho':float(rho),'p':float(p),'interpretation':'Bulk association only; does not establish PGAM5 specificity or macrophage state identity'},indent=2))
(R/'implementation_checks.json').write_text(json.dumps({'independent_logrank_tests':2,'logrank_agreement':True,'all_group_KM_steps_checked':True,'score_group_cutpoint':'single outcome-blind global median of NNLS target mixture coefficient','PFS_available':False,'PFI_available':True,'status':'Passed within scope; GSE151530-only reference failed internal diagnostic; abundance specificity unvalidated'},indent=2))
(R/'environment_versions.json').write_text(json.dumps({'python':platform.python_version(),**{n:md.version(n) for n in ['numpy','pandas','scipy','lifelines','matplotlib','pyreadr','statsmodels']}},indent=2));print(res.to_string(index=False),flush=True)
