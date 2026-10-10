"""Fixed-median KM/logrank analysis; PFI is never relabeled PFS."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, proportional_hazard_test
from statsmodels.stats.multitest import multipletests

R = Path(__file__).resolve().parent
UNITS = ['library_RNA_contribution','equalized_cell_fraction']
ENDPOINTS = ['OS','PFI']
COLORS = {'High':'#b24a3c','Low':'#2376ab'}
DAYS_PER_MONTH = 30.4375

def draw_panel(ax, tab, ep, t, models, stat, unit, ticks):
    for group in ['High','Low']:
        km = models[group]
        sf, ci = km.survival_function_, km.confidence_interval_
        x = sf.index.to_numpy()/DAYS_PER_MONTH
        subset = t[t.group.eq(group)]
        ax.step(x,sf.iloc[:,0],where='post',color=COLORS[group],lw=2,
                linestyle='-' if group == 'High' else '--',
                label=f'{group} (n={len(subset)}, events={int(subset.event.sum())})')
        ax.fill_between(x,ci.iloc[:,0],ci.iloc[:,1],step='post',color=COLORS[group],alpha=.12,lw=0)
        censor = subset.loc[subset.event.eq(0),'time_days'].to_numpy()
        ax.plot(censor/DAYS_PER_MONTH,km.predict(censor).to_numpy(),'+',color=COLORS[group],ms=5,mew=.8)
    ax.set(title='Overall survival (OS)' if ep == 'OS' else 'Progression-free interval (PFI)',
           xlabel='Time (months)',ylabel='Survival probability',ylim=(0,1.02),
           xlim=(0,ticks[-1]),xticks=ticks)
    ax.text(.98,.97,f'Logrank p = {stat.logrank_p:.4g}\nBH q (2 endpoints) = {stat.BH_q_within_unit:.4g}',
            transform=ax.transAxes,ha='right',va='top',fontsize=10,
            bbox={'facecolor':'white','alpha':.88,'edgecolor':'none','pad':5})
    ax.legend(loc='lower left',frameon=False,fontsize=10)
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y',color='#d7dce1',lw=.5,alpha=.6)
    tab.axis('off')
    labels = ['High','Low']
    numbers = [[str(int(t.loc[t.group.eq(label),'time_days'].ge(month*DAYS_PER_MONTH).sum()))
                for month in ticks] for label in labels]
    table = tab.table(cellText=numbers,rowLabels=labels,colLabels=[str(int(m)) for m in ticks],
                      loc='center',cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(9.5)
    table.scale(1,1.4)
    for (row,col),cell in table.get_celld().items():
        cell.set_edgecolor('#d7dce1'); cell.set_linewidth(.4)
        if col == -1 and row > 0: cell.get_text().set_color(COLORS[labels[row-1]])
    tab.set_title('Number at risk (start of month)',fontsize=9.5,pad=8)

def main():
    plt.rcParams.update({'font.size':11,'font.family':'DejaVu Sans','pdf.fonttype':42})
    frozen = json.loads((R/'frozen_grouping.json').read_text())
    cutpoints = {x['unit']:x['global_median'] for x in frozen['units']}
    surv = pd.read_csv(R/'TCGA_CDR_mirror_survival.csv.gz')
    surv['patient'] = surv['sample'].str[:12]
    patients = pd.read_csv(R/UNITS[0]/'patient_coefficients.csv').patient
    surv = surv[surv.patient.isin(patients)].copy()
    endpoint_columns = ['OS','OS.time','DSS','DSS.time','DFI','DFI.time','PFI','PFI.time']
    assert surv.groupby('patient')[endpoint_columns].nunique(dropna=False).le(1).all().all()
    surv.to_csv(R/'survival_alias_consistency_audit.csv',index=False)
    surv = surv[surv['sample'].str[13:15].eq('01')].rename(columns={'sample':'survival_sample'})
    assert surv.patient.is_unique
    stats, exclusions, curves, risk = [], [], [], []
    fits = {}
    for unit in UNITS:
        coeff = pd.read_csv(R/unit/'patient_coefficients.csv')
        assert coeff.patient.is_unique and len(coeff) == 371
        assert np.isclose(coeff.target_coefficient.median(),cutpoints[unit],atol=1e-14,rtol=0)
        assert (coeff.group == np.where(coeff.target_coefficient > cutpoints[unit],'High','Low')).all()
        merged = coeff.merge(surv,on='patient',how='left',validate='one_to_one',indicator='survival_match')
        assert len(merged) == 371
        merged.to_csv(R/unit/'patient_coefficients_survival_join.csv',index=False)
        for ep in ENDPOINTS:
            time = pd.to_numeric(merged[ep+'.time'],errors='coerce')
            event = pd.to_numeric(merged[ep],errors='coerce')
            valid = np.isfinite(time)&time.gt(0)&event.isin([0,1])
            t = merged.loc[valid].copy()
            t['time_days'] = time[valid].to_numpy()
            t['event'] = event[valid].to_numpy(dtype=int)
            t['high'] = t.group.eq('High').astype(int)
            t.to_csv(R/unit/f'{ep}_analysis_patients.csv',index=False)
            for idx,row in merged.loc[~valid].iterrows():
                why = []
                if row.survival_match != 'both': why.append('no matched survival record')
                if not np.isfinite(time[idx]) or time[idx] <= 0: why.append('missing/nonpositive/nonfinite time')
                if event[idx] not in [0,1]: why.append('missing/nonbinary event')
                exclusions.append({'unit':unit,'endpoint':ep,'patient':row.patient,'group':row.group,
                                   'time_days':time[idx],'event':event[idx],'reason':'; '.join(why)})
            high,low = t[t.high.eq(1)],t[t.high.eq(0)]
            lr = logrank_test(high.time_days,low.time_days,event_observed_A=high.event,event_observed_B=low.event)
            out = {'unit':unit,'analysis_role':'primary' if unit == UNITS[0] else 'sensitivity',
                   'endpoint':ep,'PFS_available':False,'patients':len(t),'events':int(t.event.sum()),
                   'high_patients':len(high),'high_events':int(high.event.sum()),
                   'low_patients':len(low),'low_events':int(low.event.sum()),
                   'excluded_patients':int((~valid).sum()),'global_median_cutpoint':cutpoints[unit],
                   'logrank_chi_square':float(lr.test_statistic),'logrank_p':float(lr.p_value)}
            endpoint_models = {}
            for label,subset in [('High',high),('Low',low)]:
                km = KaplanMeierFitter().fit(subset.time_days,event_observed=subset.event,label=label)
                endpoint_models[label] = km
                cur = pd.concat([km.survival_function_,km.confidence_interval_],axis=1).reset_index()
                cur.columns = ['time_days','survival','CI95_lower','CI95_upper']
                cur['unit'],cur['endpoint'],cur['group'] = unit,ep,label
                curves.append(cur)
                out[f'median_{label}_days'] = float(km.median_survival_time_) if np.isfinite(km.median_survival_time_) else np.nan
            modeldata = t[['time_days','event','high']]
            cox = CoxPHFitter().fit(modeldata,duration_col='time_days',event_col='event',
                                   fit_options={'precision':1e-10,'max_steps':500})
            h = cox.summary.loc['high']
            ph = proportional_hazard_test(cox,modeldata,time_transform='rank').summary.loc['high','p']
            out.update({'descriptive_unadjusted_HR_high_low':h['exp(coef)'],
                        'HR_CI95_lower':h['exp(coef) lower 95%'],'HR_CI95_upper':h['exp(coef) upper 95%'],
                        'descriptive_Cox_p':h['p'],'PH_rank_time_p':float(ph)})
            stats.append(out)
            fits[(unit,ep)] = (t,endpoint_models)
    summary = pd.DataFrame(stats)
    summary['BH_q_within_unit'] = np.nan
    for unit in UNITS:
        idx = summary.unit.eq(unit)
        summary.loc[idx,'BH_q_within_unit'] = multipletests(summary.loc[idx,'logrank_p'],method='fdr_bh')[1]
    summary['BH_q_joint_four_tests'] = multipletests(summary.logrank_p,method='fdr_bh')[1]
    summary.to_csv(R/'survival_statistics.csv',index=False)
    pd.DataFrame(exclusions).to_csv(R/'excluded_patients.csv',index=False)
    pd.concat(curves).to_csv(R/'KM_curve_values.csv',index=False)
    # Shared time scale across both endpoints and both references.
    max_days = max(t.time_days.max() for t,_ in fits.values())
    ticks = np.arange(0,np.ceil(max_days/DAYS_PER_MONTH/24)*24+1,24)
    for (unit,ep),(t,models) in fits.items():
        stat = summary[(summary.unit == unit)&(summary.endpoint == ep)].iloc[0]
        fig = plt.figure(figsize=(8.6,7.4))
        grid = fig.add_gridspec(2,1,height_ratios=[4.3,1.0],hspace=.38)
        ax,tab = fig.add_subplot(grid[0]),fig.add_subplot(grid[1])
        draw_panel(ax,tab,ep,t,models,stat,unit,ticks)
        kind = 'RNA contribution (primary)' if unit == UNITS[0] else 'Equalized cell RNA (sensitivity)'
        n_used = int(np.load(R/unit/'used_reference_model.npz')['genes'].size)
        fig.suptitle('TCGA-LIHC: PGAM5 RNA-detected macrophage coefficient\n'+kind+f'; {n_used}/395 genes matched',fontsize=12,y=.975)
        footer = 'High > global 371-patient median; Low <= median. Shading: 95% CI. +: censored.\n'
        footer += 'Coefficients are uncalibrated; they are not validated cell fractions. PFI is not PFS.'
        fig.text(.5,.026,footer,ha='center',fontsize=8.5)
        fig.subplots_adjust(top=.85,bottom=.13,left=.12,right=.97)
        for ext in ['png','pdf']: fig.savefig(R/unit/f'KM_{ep}.{ext}',dpi=200)
        plt.close(fig)
        for label in ['High','Low']:
            subset = t[t.group.eq(label)]
            for month in ticks:
                risk.append({'unit':unit,'endpoint':ep,'group':label,'month':float(month),
                             'at_risk':int(subset.time_days.ge(month*DAYS_PER_MONTH).sum())})
    unit = UNITS[0]
    fig = plt.figure(figsize=(13.2,7.3))
    grid = fig.add_gridspec(2,2,height_ratios=[4.3,1.0],hspace=.4,wspace=.23)
    for j,ep in enumerate(ENDPOINTS):
        t,models = fits[(unit,ep)]
        stat = summary[(summary.unit == unit)&(summary.endpoint == ep)].iloc[0]
        draw_panel(fig.add_subplot(grid[0,j]),fig.add_subplot(grid[1,j]),ep,t,models,stat,unit,ticks)
    fig.suptitle('TCGA-LIHC: PGAM5 RNA-detected macrophage coefficient\nExploratory v8 RNA contribution; 364/395 genes matched; 20 components',fontsize=13,y=.975)
    fig.text(.5,.027,'High > global median; Low <= median. Shading: 95% CI. +: censored.\n'
             'Uncalibrated coefficients, not validated cell fractions. Available progression endpoint: PFI, not PFS.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.845,bottom=.135,left=.07,right=.98)
    for ext in ['png','pdf']: fig.savefig(R/f'KM_OS_PFI_primary.{ext}',dpi=200)
    plt.close(fig)
    pd.DataFrame(risk).to_csv(R/'KM_numbers_at_risk.csv',index=False)
    print(summary.to_string(index=False))

if __name__ == '__main__': main()
