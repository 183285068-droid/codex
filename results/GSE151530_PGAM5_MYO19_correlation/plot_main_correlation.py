"""Plot the primary all-macrophage comparison from saved correlation results."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def render(root):
    root=Path(root)
    d=pd.read_csv(root/'macrophage_expression.csv.gz')
    stats=pd.read_csv(root/'correlation_statistics.csv').set_index('analysis')
    sp=stats.loc['Spearman_all_cells_log1p_CP10k']
    pe=stats.loc['Pearson_all_cells_log1p_CP10k']
    both=d.PGAM5_detected & d.MYO19_detected
    fig,ax=plt.subplots(figsize=(6.7,5.8))
    fig.subplots_adjust(left=.14,right=.97,bottom=.21,top=.86)
    ax.scatter(d.PGAM5_log1p_CP10k,d.MYO19_log1p_CP10k,s=17,alpha=.3,color='#538aa2',edgecolors='none')
    ax.scatter(d.loc[both,'PGAM5_log1p_CP10k'],d.loc[both,'MYO19_log1p_CP10k'],s=32,alpha=.85,color='#a36837',edgecolors='none',label='Both genes detected (n=22)')
    label=(f'Spearman rho={sp.correlation:.3f}\n'
           f'Asymptotic p={sp.asymptotic_p:.3g}\n'
           f'Permutation p={sp.permutation_p:.3g} (200,000 shuffles)\n\n'
           f'Pearson r={pe.correlation:.3f}\n'
           f'p={pe.asymptotic_p:.3g}; secondary BH q={pe.secondary_BH_q:.3g}')
    ax.text(.98,.97,label,transform=ax.transAxes,ha='right',va='top',fontsize=9,
            bbox={'facecolor':'white','alpha':.9,'edgecolor':'none'})
    ax.legend(loc='center right',bbox_to_anchor=(1.,.46),fontsize=8,frameon=False)
    ax.set(xlabel='PGAM5 log1p(CP10k)',ylabel='MYO19 log1p(CP10k)',
           title=f'All HCC macrophages, including zeros (n={len(d):,})')
    ax.spines[['top','right']].set_visible(False)
    ax.grid(alpha=.12)
    fig.suptitle('GSE151530: PGAM5 and MYO19 RNA expression',fontsize=12,y=.98)
    fig.text(.14,.06,'Detected: PGAM5 161; MYO19 189; both 22; neither 3,165.\n'
             'Cell-level p-values do not establish independent patient-level association.',
             fontsize=8,ha='left',va='bottom')
    fig.savefig(root/'PGAM5_MYO19_all_macrophages.png',dpi=220)
    fig.savefig(root/'PGAM5_MYO19_all_macrophages.pdf')
    plt.close(fig)


if __name__=='__main__':
    render(Path(__file__).resolve().parent)
