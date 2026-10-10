"""Saved-result plots; every displayed cell has both genes detected."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def render(root):
    root=Path(root)
    d=pd.read_csv(root/'codetected_macrophage_expression.csv.gz')
    summary=pd.read_csv(root/'correlation_summary.csv')
    assert d[['PGAM5_raw_count','MYO19_raw_count']].gt(0).all().all()
    for scale,suffix,stem in [('log1p','_log1p_CP10k','codetected_normalized_correlations'),('raw','_raw_count','codetected_raw_count_correlations')]:
        fig,axs=plt.subplots(3,3,figsize=(13.5,12),constrained_layout=True)
        for ax,(_,row) in zip(axs.flat,summary.iterrows()):
            t=d[d.cohort.eq(row.cohort)]
            x=t['PGAM5'+suffix];y=t['MYO19'+suffix]
            if scale=='raw':
                multiplicity=t.groupby(['PGAM5_raw_count','MYO19_raw_count']).size().reset_index(name='cells')
                ax.scatter(multiplicity.PGAM5_raw_count,multiplicity.MYO19_raw_count,s=25+12*multiplicity.cells,color='#548f75',alpha=.8,edgecolors='none')
                for z in multiplicity.itertuples():
                    if z.cells>1:ax.annotate(f'n={z.cells}',(z.PGAM5_raw_count,z.MYO19_raw_count),xytext=(5,5),textcoords='offset points',fontsize=8)
            else:ax.scatter(x,y,s=24,color='#a36837',alpha=.75,edgecolors='none')
            prefix='Spearman_'+scale
            if row[prefix+'_status']=='estimable':
                detail=f"rho={row[prefix+'_r']:.3f}; permutation p={row[prefix+'_p']:.3g}\nBH q={row[prefix+'_q']:.3g}"
                if scale=='log1p' and row.Spearman_raw_status=='not_estimable_constant_vector':detail+='; raw count vector constant'
            else:detail='Not estimable: '+('only 1 co-detected cell' if row.both_detected<3 else 'constant raw-count vector')
            name=row.cohort.replace('_',' / ')
            if row.cohort=='GSE140228_Droplet':name+=' (Tumor incl. CC)'
            ax.set_title(f'{name}\nBoth detected n={int(row.both_detected)}\n{detail}',fontsize=9.2)
            unit='log1p(CP10k)' if scale=='log1p' else f'{row.count_type}'
            ax.set_xlabel(f'PGAM5 {unit}',fontsize=9);ax.set_ylabel(f'MYO19 {unit}',fontsize=9)
            ax.tick_params(labelsize=8);ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.12)
        fig.suptitle('PGAM5 and MYO19: co-detected macrophages only\n'+('Normalized RNA expression' if scale=='log1p' else 'Original raw counts; bubble size reflects overlapping cells'),fontsize=13)
        fig.savefig(root/(stem+'.png'),dpi=190);fig.savefig(root/(stem+'.pdf'));plt.close(fig)
    # Contrast scales directly: a large normalized coefficient does not imply
    # raw counts covary, especially when both genes have only one count.
    fig,ax=plt.subplots(figsize=(8.8,5.3))
    yy=np.arange(len(summary))
    ax.scatter(summary.Spearman_log1p_r,yy-.10,color='#a36837',s=45,label='log1p(CP10k)')
    ax.scatter(summary.Spearman_raw_r,yy+.10,color='#548f75',s=45,label='Raw counts')
    for i,row in summary.iterrows():
        if pd.isna(row.Spearman_log1p_r):ax.text(-.95,i,'n=1; not estimable',va='center',fontsize=8)
        elif pd.isna(row.Spearman_raw_r):ax.text(-.95,i+.12,'Raw vector constant',va='center',fontsize=8,color='#39735c')
    ax.axvline(0,color='#999999',linewidth=.8)
    ax.set(yticks=yy,yticklabels=[f'{r.cohort} (n={int(r.both_detected)})' for r in summary.itertuples()],xlabel='Spearman rho',xlim=(-1.02,1.02))
    ax.invert_yaxis();ax.spines[['top','right']].set_visible(False);ax.grid(axis='x',alpha=.12);ax.legend(fontsize=9,loc='lower left',bbox_to_anchor=(0,1.01),ncol=2,frameon=False);fig.suptitle('Co-detected macrophages: compare expression scales',fontsize=12,y=.98)
    fig.tight_layout(rect=(0,0,1,.9));fig.savefig(root/'correlation_scale_comparison.png',dpi=190);fig.savefig(root/'correlation_scale_comparison.pdf');plt.close(fig)


if __name__=='__main__':render(Path(__file__).resolve().parent)
