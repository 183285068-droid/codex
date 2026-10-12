from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/gse245908_mpl')
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
R=Path(__file__).resolve().parent;F=R/'figures';F.mkdir(exist_ok=True)
TYPES=['All_CD8','CD8_exhaustion','CD8_activation','CD8_effector','CD8_naive_memory','CD4_other','Treg','T_unresolved'];LABELS=['All CD8 candidates','CD8 exhaustion','CD8 activation','CD8 effector','CD8 naive/memory','CD4 other','Treg candidates','T unresolved'];GROUPS=['PGAM5_high','Other_detected','Undetected'];COLORS=['#bd5124','#276c9b','#929292']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,name):
 fig.savefig(F/(name+'.png'),dpi=180,bbox_inches='tight');fig.savefig(F/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
def heat(ax,x,rows,cols,title,vmax=None):
 im=ax.imshow(np.ma.masked_invalid(x),aspect='auto',vmin=0,vmax=vmax,cmap='Blues');ax.set_xticks(range(len(cols)),cols,rotation=35,ha='right');ax.set_yticks(range(len(rows)),rows);ax.set_title(title)
 for i in range(x.shape[0]):
  for j in range(x.shape[1]):ax.text(j,i,'NA' if np.isnan(x[i,j]) else f'{x[i,j]:.1f}',ha='center',va='center',fontsize=8,color='white' if not np.isnan(x[i,j]) and vmax and x[i,j]>.55*vmax else 'black')
 return im
summary=pd.read_csv(R/'subtype_neighbor_summary.csv');stats=pd.read_csv(R/'spatial_subtype_tests.csv');lr=pd.read_csv(R/'LR_by_subtype.csv')
fig,axs=plt.subplots(1,2,figsize=(14,6));width=.24
for ax,s in zip(axs,['CHC20','CHC23']):
 for j,g in enumerate(GROUPS):
  d=summary[(summary['sample']==s)&(summary.source_group==g)].set_index('subtype').loc[TYPES];ax.barh(np.arange(8)+(j-1)*width,d.mean_neighbor_fraction*100,height=width,color=COLORS[j],label=f'{g} (n={int(d.n_sources.iloc[0])})')
 ax.set_yticks(range(8),LABELS);ax.invert_yaxis();ax.set_xlim(0,20);ax.set_xlabel('Mean proportion of retained neighbors (%)');ax.set_title(s);ax.legend(fontsize=8,loc='lower right')
fig.suptitle('T subtype candidate neighborhoods | source-weighted averages\nCD8 states overlap; barcodes are mixed tissue spots');fig.tight_layout();save(fig,'subtype_neighbor_coverage')
fig,axs=plt.subplots(1,2,figsize=(14,7))
for ax,s in zip(axs,['CHC20','CHC23']):
 d=stats[(stats['sample']==s)&(stats.kind=='subtype_coverage')].set_index('feature').loc[TYPES];ax.axvline(0,color='#555',lw=1);ax.scatter(d.delta_high_minus_other*100,range(8),color=COLORS[0]);ax.set_yticks(range(8),LABELS);ax.invert_yaxis();ax.set_xlim(-16,16);ax.set_xlabel('High − other detected (percentage points)');ax.set_title(s)
 for i,row in enumerate(d.itertuples()):ax.text(.99,i,f'p={row.spatial_conditional_p:.3f}; q={row.BH_q:.3f}',transform=ax.get_yaxis_transform(),ha='right',va='bottom',fontsize=8)
fig.suptitle('Unadjusted neighborhood differences with spatial-conditional p/q\nSingle BH family of 24 prespecified comparisons');fig.tight_layout();save(fig,'subtype_spatial_conditional_comparison')
fig,axs=plt.subplots(2,4,figsize=(17,9))
for row,s in enumerate(['CHC20','CHC23']):
 o=pd.read_csv(R/f'{s}_spot_subtypes.csv');
 for ax,k,title in zip(axs[row],['All_CD8','CD8_exhaustion','CD8_naive_memory','Treg'],['CD8 candidates','CD8 exhaustion','CD8 naive/memory','Treg candidates']):
  ax.scatter(o.pixel_col,o.pixel_row,s=2,color='#ddd');z=o['gate_'+k];ax.scatter(o.loc[z,'pixel_col'],o.loc[z,'pixel_row'],s=13,color='#276c9b',label=f'{title}: {z.sum()}');h=o.source_group.eq('PGAM5_high');ax.scatter(o.loc[h,'pixel_col'],o.loc[h,'pixel_row'],s=30,facecolors='none',edgecolors=COLORS[0],linewidths=1.1,label=f'High sources: {h.sum()}');ax.invert_yaxis();ax.set_aspect('equal');ax.set_title(s+' | '+title);ax.set_xlabel('Image pixel column');ax.set_ylabel('Image pixel row');ax.legend(fontsize=7,loc='upper right')
fig.suptitle('True tissue coordinates | overlapping RNA proxy labels');fig.tight_layout();save(fig,'spatial_T_subtypes')
pairs=['SPP1--CD44','PVR--TIGIT','NECTIN2--TIGIT','TGFB1--TGFBR1_TGFBR2','LGALS9--HAVCR2','CXCL9--CXCR3','CXCL10--CXCR3','CD274--PDCD1'];fig,axs=plt.subplots(1,2,figsize=(17,6))
for ax,s in zip(axs,['CHC20','CHC23']):
 d=lr[(lr['sample']==s)&(lr.source_group=='PGAM5_high')].copy();d['pair']=d.ligand+'--'+d.receptor;mat=d.pivot(index='subtype',columns='pair',values='codetection_pct').reindex(index=TYPES,columns=pairs).to_numpy();rowlabels=[f'{label} (edges={int(d[d.subtype==k].n_edges.iloc[0])})' for label,k in zip(LABELS,TYPES)];im=heat(ax,mat,rowlabels,pairs,s+' | High source → subtype neighbors',100)
fig.suptitle('Selected candidate LR pairs | NA means no qualifying edges\nOverlapping states; percentages are not communication probabilities');fig.subplots_adjust(top=.84,bottom=.28,wspace=.72,right=.9);cax=fig.add_axes([.94,.3,.015,.5]);fig.colorbar(im,cax=cax,label='RNA codetected edges (%)');save(fig,'LR_subtype_codetection')
# Show candidate edges supported by SPP1-CD44, without ranking by apparent significance.
fig,axs=plt.subplots(1,2,figsize=(12,6))
for ax,s in zip(axs,['CHC20','CHC23']):
 o=pd.read_csv(R/f'{s}_spot_subtypes.csv',index_col=0);e=pd.read_csv(R/f'{s}_subtype_LR_edges.csv.gz');e=e[e.source_group.eq('PGAM5_high')&e.target_All_CD8&(e['SPP1--CD44']>0)];xy=o[['pixel_col','pixel_row']];segments=[np.array([xy.loc[x].to_numpy(),xy.loc[y].to_numpy()]) for x,y in zip(e.source_barcode,e.target_barcode)];ax.scatter(o.pixel_col,o.pixel_row,s=2,c='#ddd');ax.add_collection(LineCollection(segments,colors='#276c9b',linewidths=1.2));h=o.source_group.eq('PGAM5_high');ax.scatter(o.loc[h,'pixel_col'],o.loc[h,'pixel_row'],s=18,c=COLORS[0]);ax.invert_yaxis();ax.set_aspect('equal');ax.set_title(f'{s}: SPP1–CD44, {len(e)} codetected CD8 edges');ax.set_xlabel('Image pixel column');ax.set_ylabel('Image pixel row')
fig.suptitle('Candidate edges at genuine tissue coordinates | no causal inference');fig.tight_layout();save(fig,'SPP1_CD44_CD8_spatial_edges')
