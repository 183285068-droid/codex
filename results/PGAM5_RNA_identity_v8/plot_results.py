"""Plots for pooled internal validation; protein/cross-cohort proof not required."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R=Path(__file__).resolve().parent


def save(fig,name):
    fig.savefig(R/f'{name}.png',dpi=180,bbox_inches='tight')
    fig.savefig(R/f'{name}.pdf',bbox_inches='tight')
    plt.close(fig)


def main():
    summary=pd.read_csv(R/'pooled_internal_validation_summary.csv')
    doses=pd.read_csv(R/'pooled_internal_dose_recovery.csv')
    colors={'continuous_rank':'#2A9D8F','log2FC_1':'#355C7D','log2FC_0p5':'#C06C84'}
    labels={'continuous_rank':'No FDR/log2FC cutoff','log2FC_1':'log2FC >= 1; FDR < 0.05','log2FC_0p5':'log2FC >= 0.5; FDR < 0.05'}
    units=['equalized_cell_fraction','library_RNA_contribution']
    fig,axes=plt.subplots(1,2,figsize=(11,4.3))
    for ax,unit in zip(axes,units):
        d=summary[summary.unit.eq(unit)]
        for i,row in enumerate(d.itertuples()):
            ax.bar(i,row.all_zero_P95_percent,color=colors[row.threshold_version],width=.55)
            ax.text(i,row.all_zero_P95_percent,f'{row.all_zero_P95_percent:.2f}%',ha='center',va='bottom')
        ax.axhline(.5,color='#A32020',ls='--',label='Project acceptance limit: 0.5%')
        ax.set_xticks(range(len(d)),[labels[x] for x in d.threshold_version],rotation=12,ha='right',fontsize=9)
        ax.set_title('Idealized equal-cell mixtures' if unit==units[0] else 'Observed-library RNA mixtures')
        ax.set_ylabel('Target-absent estimate P95 (coefficient x100)')
        ax.set_ylim(0,max(.8,d.all_zero_P95_percent.max()*1.3));ax.legend(fontsize=8)
    fig.suptitle('Pooled internal patient/sample holdout: target-absent false signal')
    fig.text(.5,-.03,'Existing computational mixtures; no validated TCGA cell-fraction interpretation.',ha='center',fontsize=9)
    fig.tight_layout();save(fig,'pooled_false_signal')
    fig,axes=plt.subplots(1,2,figsize=(11,4.3))
    for ax,unit in zip(axes,units):
        for level in colors:
            d=doses[doses.unit.eq(unit)&doses.threshold_version.eq(level)].sort_values('nominal_cell_fraction')
            x=d.median_truth.to_numpy()*100;y=d.median_prediction.to_numpy()*100
            lo=d.prediction_q25.to_numpy()*100;hi=d.prediction_q75.to_numpy()*100
            ax.errorbar(x,y,yerr=np.array([np.maximum(y-lo,0),np.maximum(hi-y,0)]),marker='o',capsize=4,color=colors[level],label=labels[level])
        upper=max(5,float(doses[doses.unit.eq(unit)].median_truth.max()*100))*1.05
        ax.plot([0,upper],[0,upper],ls='--',color='gray',label='Identity')
        ax.set_title('Idealized equal-cell mixtures' if unit==units[0] else 'Observed-library RNA mixtures')
        ax.set_xlabel('Median actual target fraction (%)' if unit==units[0] else 'Median actual target RNA contribution (%)')
        ax.set_ylabel('Median estimated coefficient x100')
        ax.legend(fontsize=8)
    fig.suptitle('Pooled internal dose recovery: median and interquartile range')
    fig.text(.5,-.03,'Technical mixtures are repeated draws; they are not independent patients or physical experiments.',ha='center',fontsize=9)
    fig.tight_layout();save(fig,'pooled_dose_recovery')
    a=pd.read_csv(R/'annotation_coverage.csv');a=a[a.scope.eq('primary_all_QC_HCC')]
    fig,ax=plt.subplots(figsize=(8,4))
    ax.bar(a.dataset,a.PGAM5_RNA_detected_macrophages,color='#9AA9B1',label='PGAM5 RNA-detected macrophages')
    ax.bar(a.dataset,a.marker_corroborated_PGAM5_RNA_detected_macrophages,color='#355C7D',label='With descriptive lineage-marker corroboration')
    for row in a.itertuples():
        ax.text(row.dataset,row.PGAM5_RNA_detected_macrophages,f'{row.marker_corroborated_PGAM5_RNA_detected_macrophages}/{row.PGAM5_RNA_detected_macrophages}',ha='center',va='bottom')
    ax.set_ylabel('Observed cells');ax.set_ylim(0,a.PGAM5_RNA_detected_macrophages.max()*1.3)
    ax.set_title('RNA detection labels and macrophage-lineage evidence')
    ax.legend(fontsize=8);fig.tight_layout();save(fig,'RNA_lineage_evidence')


if __name__=='__main__':main()
