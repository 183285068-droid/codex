from pathlib import Path
import numpy as np,pandas as pd,json
from scTenifold.core._networks import make_networks
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent
x=pd.read_csv(R/'inputs/network_raw_counts.csv.gz',index_col=0);m=pd.read_csv(R/'inputs/target_cell_metadata.csv.gz',index_col=0);cp=x.div(m.total_counts,axis='columns')*1e6;j=x.index.get_loc('PGAM5')
# Identical bootstrap sequence, no sparsification: inspect what q=.95 would discard.
with threadpool_limits(limits=4):nets=make_networks(cp,n_nets=30,n_samp_cells=171,n_comp=3,q=0,random_state=151530,replace=True)
rows=[]
for i,n in enumerate(nets):
 w=n.toarray();active=np.any(w!=0,axis=0)|np.any(w!=0,axis=1);threshold=float(np.quantile(np.abs(w[np.ix_(active,active)]),.95));pg=np.abs(w[:,j]);rows.append(dict(network=i+1,max_unfiltered_PGAM5_outgoing_abs_weight=float(pg.max()),global_q95_cutoff=threshold,PGAM5_edges_before_sparsification=int((pg>0).sum()),PGAM5_edges_meeting_q95=int((pg>=threshold).sum()),maximum_to_cutoff_ratio=float(pg.max()/threshold)))
del nets
with threadpool_limits(limits=4):filtered=make_networks(cp,n_nets=30,n_samp_cells=171,n_comp=3,q=.95,random_state=151530,replace=True)
for row,n in zip(rows,filtered):
 actual=int(n.tocsr()[:,j].count_nonzero());assert actual==row['PGAM5_edges_meeting_q95'];row['actual_retained_edges_verified']=actual
del filtered
t=pd.DataFrame(rows);t.to_csv(R/'PGAM5_network_sparsification_diagnostic.csv',index=False)
fig,axs=plt.subplots(1,2,figsize=(12,4));axs[0].plot(t.network,t.max_unfiltered_PGAM5_outgoing_abs_weight,marker='o',markersize=3,label='Maximum PGAM5 outgoing |weight|');axs[0].plot(t.network,t.global_q95_cutoff,label='Global 95% sparsification cutoff');axs[0].set_xlabel('Bootstrap network (seed 151530)');axs[0].set_ylabel('Scaled absolute PC-regression weight');axs[0].legend(fontsize=8);axs[0].set_title('PGAM5 outgoing weights relative to cutoff');axs[1].bar(t.network,t.PGAM5_edges_meeting_q95,color='#277ea7');axs[1].set_xlabel('Bootstrap network');axs[1].set_ylabel('PGAM5 edges retained at q=.95');axs[1].set_ylim(0,max(1,float(t.PGAM5_edges_meeting_q95.max())*1.2));axs[1].set_title('Pre-tensor intervention feasibility');fig.tight_layout();fig.savefig(R/'06_PGAM5_edge_filtering_diagnostic.png',dpi=180);fig.savefig(R/'06_PGAM5_edge_filtering_diagnostic.pdf');plt.close(fig)
print(t[['PGAM5_edges_meeting_q95','maximum_to_cutoff_ratio']].describe().to_string())
