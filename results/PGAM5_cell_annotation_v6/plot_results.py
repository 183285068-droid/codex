"""Evidence figures from saved counts, selected held-cohort states and recovery."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R = Path(__file__).resolve().parent
COHORTS = json.loads((R/'protocol.json').read_text())['primary_cohorts']
plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})


def save(fig,name):
    fig.savefig(R/f'{name}.png',dpi=180,bbox_inches='tight',facecolor='white')
    fig.savefig(R/f'{name}.pdf',bbox_inches='tight',facecolor='white')
    plt.close(fig)


def main():
    d = pd.read_csv(R/'raw_count_threshold_sensitivity.csv')
    values = np.array([[int(d[d.dataset.eq(c)&d.raw_count_cutoff.eq(t)].cells_above_cutoff.iloc[0])
                        for t in [1,2,3,5]] for c in COHORTS])
    fig,ax = plt.subplots(figsize=(9,5.5))
    im = ax.imshow(np.log1p(values),aspect='auto',cmap='Blues',vmin=0)
    for i in range(5):
        for j in range(4):
            ax.text(j,i,str(values[i,j]),ha='center',va='center',color='white' if values[i,j]>100 else 'black',fontsize=12)
    ax.set_xticks(range(4),['>=1','>=2','>=3','>=5'])
    ax.set_yticks(range(5),COHORTS)
    ax.set(xlabel='Observed PGAM5 raw RNA count threshold',title='Increasing the raw-count threshold sharply reduces available cells')
    fig.colorbar(im,ax=ax,label='log(1 + cell count)',fraction=.035,pad=.03)
    fig.subplots_adjust(bottom=.24)
    fig.text(.08,.025,'All 30,533 primary macrophages retained. Counts >=3 give only 21 cells across five cohorts.\n'
             'A raw-count cutoff denotes RNA capture; it does not establish a stable biological subtype.',fontsize=9)
    save(fig,'raw_count_thresholds')

    donors = pd.read_csv(R/'program_PGAM5_association_by_donor.csv')
    assoc = pd.read_csv(R/'program_PGAM5_association_by_cohort.csv')
    fig,ax = plt.subplots(figsize=(10,5.7))
    for j,c in enumerate(COHORTS):
        target = json.loads((R/c/'definition.json').read_text())['target_program']
        a = assoc[assoc.fold.eq(c)&assoc.role.eq('heldout')&assoc.program.eq(target)].iloc[0]
        subset = donors[donors.fold.eq(c)&donors.role.eq('heldout')&donors.program.eq(target)&donors.eligible]
        y = np.log2(subset.PGAM5_CP10k_ratio.to_numpy())
        if len(y):
            ax.scatter(j+np.linspace(-.10,.10,len(y)),y,s=30,color='#236ca4',alpha=.8)
            ax.hlines(np.log2(np.median(subset.PGAM5_CP10k_ratio)),j-.18,j+.18,color='black',lw=2)
        ax.text(j,2.6,f'n={len(y)}\nSupport failed',ha='center',fontsize=9)
    ax.axhline(0,color='#777777',linestyle='--',label='No normalized-expression enrichment')
    ax.axhline(np.log2(1.5),color='#bc3338',linestyle=':',label='Required median ratio >=1.5')
    ax.set_xticks(range(5),COHORTS)
    ax.set(ylabel='Donor log2 PGAM5 CP10k ratio\n(candidate / other macrophages)',
           title='Selected training programs fail PGAM5 association in held cohorts',ylim=(-2,3.3))
    ax.legend(loc='lower left',frameon=False,fontsize=9)
    fig.subplots_adjust(bottom=.22)
    fig.text(.08,.025,'Dots: evaluable donor/sample labels; black lines: medians. GSE202642 uses sample proxies.\n'
             'Each program and threshold is selected using other cohorts only. No patient-paired DE or depth regression.\n'
             'Enrichment, coverage and donor consistency must pass together; no held-cohort candidate passes.',fontsize=9)
    save(fig,'heldout_program_PGAM5_association')

    summary = pd.read_csv(R/'reference_validation_summary.csv')
    primary = summary[summary.estimand.eq('all_PGAM5_RNA')].set_index(['dataset','unit'])
    fig,ax = plt.subplots(figsize=(10,5.6))
    x = np.arange(5)
    max_y = 0
    for i,(unit,label,color) in enumerate([('equalized_cell_fraction','Idealized equal-cell normalization','#236ca4'),
                                         ('library_RNA_contribution','Observed library RNA contribution','#d47c32')]):
        y = [primary.loc[c,unit].all_zero_P95_percent for c in COHORTS]
        max_y = max(max_y,max(y))
        bars = ax.bar(x+(i-.5)*.36,y,width=.36,color=color,label=label)
        ax.bar_label(bars,fmt='%.2f',padding=3,fontsize=9)
    ax.axhline(.5,color='#bc3338',linestyle='--',label='Fixed gate: <=0.5%')
    ax.set_xticks(x,COHORTS)
    ax.set(ylabel='Zero-target estimated coefficient P95 x 100',ylim=(0,max_y*1.3),
           title='False RNA-positive macrophage signal in independent target-absent mixtures')
    ax.legend(frameon=False,loc='upper left',fontsize=9)
    fig.subplots_adjust(bottom=.24)
    fig.text(.08,.025,'Target is observed PGAM5 raw-count>0 macrophages; inferred RNA-zero state members remain separate.\n'
             'Sum the two RNA-detected state components and any retained cycling counterparts.\n'
             'All zero-target stress scenarios included. Coefficients are not calibrated TCGA cellular abundance.',fontsize=9)
    save(fig,'RNA_positive_reference_false_signal')


if __name__=='__main__':
    main()
