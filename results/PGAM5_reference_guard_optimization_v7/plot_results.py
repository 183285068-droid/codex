"""Standalone source-backed false-signal, dose and lineage figures."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R=Path(__file__).resolve().parent
C=json.loads((R/'protocol.json').read_text())['primary_cohorts']
U=['equalized_cell_fraction','library_RNA_contribution']
M=['v6_baseline','background_prototypes','target_only_budget','prototypes_target_only_budget']
LABELS=['v6 baseline','Background prototypes','Target-only PGAM5 bound','Prototypes + target-only bound']
COLORS=['#336b9b','#d58929','#478562','#9b4c8d']
plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})

def save(fig,name):
    fig.savefig(R/f'{name}.png',dpi=170,bbox_inches='tight',facecolor='white')
    fig.savefig(R/f'{name}.pdf',bbox_inches='tight',facecolor='white');plt.close(fig)

def bars(summary,value,name,ylabel,title,gate=None):
    fig,axes=plt.subplots(2,1,figsize=(12,8))
    for i,(unit,ax) in enumerate(zip(U,axes)):
        subset=summary[summary.unit.eq(unit)].set_index(['dataset','method'])
        for j,(m,l,color) in enumerate(zip(M,LABELS,COLORS)):
            y=[subset.loc[c,m][value] for c in C]
            b=ax.bar(np.arange(5)+(j-1.5)*.19,y,width=.18,color=color,label=l)
            ax.bar_label(b,fmt='%.2f',padding=2,fontsize=7)
        if gate is not None:ax.axhline(gate,color='#bc3338',ls='--',lw=1.2,label=f'Fixed gate <= {gate}')
        ax.set_xticks(range(5),C);ax.set_ylabel(ylabel)
        ax.set_title('Idealized equal-cell unit' if i==0 else 'Observed library RNA-contribution unit',fontsize=11)
        ax.set_ylim(0,max(subset.loc[(c,m)][value] for c in C for m in M)*1.3)
        if i==0:ax.legend(frameon=False,ncol=3,fontsize=8,loc='upper left')
    fig.suptitle(title,fontsize=14,y=.99)
    fig.subplots_adjust(hspace=.38,bottom=.17)
    fig.text(.07,.025,'Same actual held-cohort mixtures; coefficients are not calibrated TCGA cell fractions.\n'
        'Target-only bound is a documented retrospective repair. Original total-PGAM5 constraints were rejected.\n'
        'Lower false signal must be assessed together with dose recovery and total-macrophage error.',fontsize=9)
    save(fig,name)

def main():
    summary=pd.read_csv(R/'all_validation_summary.csv')
    bars(summary,'all_zero_P95_percent','false_signal_comparison','Zero-target coefficient P95 x100',
         'False target signal after reference and RNA-bound changes',.5)
    bars(summary,'standard_macrophage_MAE_pp','macrophage_tradeoff','Total-macrophage MAE (percentage points)',
         'Total-macrophage recovery checks the cost of target-specific changes')
    d=pd.read_csv(R/'all_predictions.csv.gz');d=d[d.scenario.eq('standard')]
    fig,axes=plt.subplots(2,5,figsize=(14,7),sharex=True)
    for i,unit in enumerate(U):
        for j,c in enumerate(C):
            ax=axes[i,j];x=d[d.unit.eq(unit)&d.dataset.eq(c)]
            truth=x[x.method.eq('v6_baseline')].groupby('nominal_target_cell_fraction').truth.median()*100
            ax.plot(truth.index*100,truth.values,color='black',ls='--',lw=1.5,label='Median actual truth')
            for m,l,color in zip(M,LABELS,COLORS):
                y=x[x.method.eq(m)].groupby('nominal_target_cell_fraction').predicted.median()*100
                ax.plot(y.index*100,y.values,color=color,marker='o',ms=3,lw=1.2,label=l)
            ax.set_title(c,fontsize=10);ax.set_xticks([0,1,2,5]);ax.set_ylim(bottom=0)
            if j==0:ax.set_ylabel(('Equal-cell' if i==0 else 'Library RNA')+'\nMedian coefficient x100')
            if i==1:ax.set_xlabel('Nominal target cell %',fontsize=9)
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.96),ncol=3,frameon=False,fontsize=9)
    fig.suptitle('Dose recovery exposes zero-output and background-confounding failure',fontsize=14,y=1.02)
    fig.subplots_adjust(top=.81,bottom=.19,hspace=.4,wspace=.34)
    fig.text(.06,.025,'Each line aggregates donor-labelled and technical mixture repetitions; lines are descriptive, not confidence intervals.\n'
        'RNA-contribution truth differs from nominal cell proportion. GSE202642 labels are sample proxies.\n'
        'GSE189903 has only two eligible independent donors. Exact per-dose values are in all_dose_recovery.csv.',fontsize=9)
    save(fig,'dose_recovery_comparison')

if __name__=='__main__':main()
