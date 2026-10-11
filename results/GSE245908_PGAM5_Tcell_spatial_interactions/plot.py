from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE245908_Tcell_mpl')
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Patch
from matplotlib.collections import PatchCollection,LineCollection
from matplotlib.colors import Normalize
import numpy as np,pandas as pd,json,gzip
from PIL import Image
from io import BytesIO
from analyze import ROOT,BASE,SAMPLES,PAIRS
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'figure.facecolor':'white'})
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
DATA={s:pd.read_csv(ROOT/f'{s}_spot_states.csv',index_col=0) for s in SAMPLES}
COL={'PGAM5_high':'#B2182B','Other_detected':'#EF8A62','Undetected':'#7B8DA1','Outside_mac_gate':'#DADDE1'}
T_COL={'Activation_supported':'#2166AC','Exhaustion_supported':'#762A83','Both_programs':'#E6AB02','T_other':'#80B6B1','Outside_T_gate':'#DADDE1'}
def save(fig,name):
 fig.savefig(OUT/f'{name}.png',dpi=210,bbox_inches='tight');fig.savefig(OUT/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
def setup(ax,s):
 gsm={'CHC20':'GSM7850822','CHC23':'GSM7850823'}[s];p=BASE/'inputs'/s;sf=json.loads(gzip.decompress((p/f'{gsm}_{s}_scalefactors_json.json.gz').read_bytes()));image=np.asarray(Image.open(BytesIO(gzip.decompress((p/f'{gsm}_{s}_tissue_hires_image.png.gz').read_bytes()))));xy=DATA[s][['pixel_col','pixel_row']].to_numpy()*sf['tissue_hires_scalef'];ax.imshow(image,alpha=.5);ax.set_xlim(xy[:,0].min()-30,xy[:,0].max()+30);ax.set_ylim(xy[:,1].max()+30,xy[:,1].min()-30);ax.set_aspect('equal');ax.axis('off');return xy,sf['spot_diameter_fullres']*sf['tissue_hires_scalef']/2

def circles(ax,s,colors):
 xy,r=setup(ax,s);collection=PatchCollection([Circle(p,r) for p in xy],facecolor=colors,linewidth=0,alpha=.9,rasterized=True);ax.add_collection(collection);return xy

def states():
 fig,axs=plt.subplots(2,4,figsize=(20,10))
 for row,s in enumerate(SAMPLES):
  d=DATA[s];circles(axs[row,0],s,d.source_group.map(COL));axs[row,0].set_title(f'{s}: PGAM5 source RNA groups\nHigh n={(d.source_group=="PGAM5_high").sum()}')
  circles(axs[row,1],s,d.T_state.map(T_COL));axs[row,1].set_title(f'{s}: marker-supported T states\nActivation={d.Activation_early_supported_high.sum()}, exhaustion={d.Exhaustion_core_supported_high.sum()}',fontsize=11)
  for col,k in [(2,'Activation_early'),(3,'Exhaustion_core')]:
   mx=max(DATA[x].loc[DATA[x].T_enriched,k].max() for x in SAMPLES);norm=Normalize(0,mx);cm=plt.get_cmap('viridis');color=cm(norm(d[k].to_numpy()));color[~d.T_enriched]=[.84,.86,.88,.85];circles(axs[row,col],s,color);axs[row,col].set_title(f'{s}: {k.replace("_"," ")} RNA');fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cm),ax=axs[row,col],fraction=.032,pad=.02)
 legends=[Patch(fc=COL[k],label=k.replace('_',' ')) for k in ['PGAM5_high','Other_detected','Undetected']]+[Patch(fc=T_COL[k],label=k.replace('_',' ')) for k in ['Activation_supported','Exhaustion_supported','Both_programs','T_other']]
 fig.legend(handles=legends,loc='lower center',ncol=4,bbox_to_anchor=(.5,.02),frameon=False,fontsize=10);fig.suptitle('GSE245908: spatial RNA groups and T activation/exhaustion programs',fontsize=17)
 fig.text(.025,.0,'Mixed spots; RNA marker support does not establish exhausted/activated single-cell identity. Grey = outside the displayed lineage gate.',fontsize=11);fig.subplots_adjust(top=.9,bottom=.13,wspace=.12,hspace=.19);save(fig,'spatial_PGAM5_sources_T_states')

def neighbor():
 tests=pd.read_csv(ROOT/'spatial_conditional_tests.csv');keys=['T_neighbor_fraction','Activation_early','Effector_cytotoxic','Exhaustion_core'];groups=['PGAM5_high','Other_detected','Undetected'];fig,axs=plt.subplots(2,4,figsize=(18,9));rng=np.random.default_rng(245908)
 for row,s in enumerate(SAMPLES):
  d=DATA[s]
  for col,k in enumerate(keys):
   key=k if k=='T_neighbor_fraction' else 'T_neighbor_'+k;ax=axs[row,col];vals=[d.loc[d.source_group==g,key].dropna().to_numpy() for g in groups]
   bp=ax.boxplot(vals,positions=range(3),showfliers=False,widths=.5,patch_artist=True,medianprops={'color':'#333333','linewidth':1.2})
   for box,g in zip(bp['boxes'],groups):box.set_facecolor(COL[g]);box.set_alpha(.5)
   for i,(v,g) in enumerate(zip(vals,groups)):ax.scatter(rng.uniform(i-.17,i+.17,len(v)),v,s=8,alpha=.42,c=COL[g],linewidths=0,rasterized=True)
   r=tests[(tests['sample']==s)&(tests.feature==k)].iloc[0];ax.set_title(f'{s}: {k.replace("_"," ")}\np={r.spatial_conditional_p:.3f}, q={r.BH_q:.3f}',fontsize=11);ax.set_xticks(range(3),labels=[f'High\nn={len(vals[0])}',f'Other detected\nn={len(vals[1])}',f'Undetected\nn={len(vals[2])}'],fontsize=10);ax.set_ylabel('Adjacent T fraction' if k=='T_neighbor_fraction' else 'Mean neighbor T RNA program');ax.spines[['top','right']].set_visible(False)
 for col,k in enumerate(keys):
  key=k if k=='T_neighbor_fraction' else 'T_neighbor_'+k;maximum=max(d.loc[d.Mac_enriched,key].max() for d in DATA.values())
  for row in range(2):axs[row,col].set_ylim(0,maximum*1.07)
 fig.suptitle('PGAM5-high vs other detected sources: adjacent T signals',fontsize=17);fig.text(.02,.01,'p: spatial-block/depth conditional randomization; q: BH across all 50 primary tests. Undetected group is descriptive.\nProgram comparisons require at least one T-enriched neighbor; plotted observations are source spots, not individual T cells.',fontsize=11);fig.tight_layout(rect=[0,.065,1,.95]);save(fig,'neighbor_T_program_comparison')

def LRheat():
 detail=pd.read_csv(ROOT/'LR_by_receiver_state.csv');cols=[(s,g) for s in SAMPLES for g in ['PGAM5_high','Other_detected','Undetected']];pairs=[l+'--'+r for l,r,_ in PAIRS];matrix=np.zeros((len(pairs),len(cols)));counts=np.zeros_like(matrix);den=np.zeros_like(matrix)
 for i,(l,r,_) in enumerate(PAIRS):
  for j,(s,g) in enumerate(cols):
   z=detail[(detail['sample']==s)&(detail.source_group==g)&(detail.receiver_T_state=='All_T')&(detail.ligand==l)&(detail.receptor_complex==r)].iloc[0];matrix[i,j]=z.codetection_pct;counts[i,j]=z.n_codetected_edges;den[i,j]=z.n_edges
 fig,ax=plt.subplots(figsize=(13,11));im=ax.imshow(matrix,cmap='YlOrRd',vmin=0,vmax=100,aspect='auto');ax.set_xticks(range(len(cols)),labels=[s+'\n'+g.replace('_',' ') for s,g in cols],fontsize=10);ax.set_yticks(range(len(pairs)),labels=[p.replace('--',' → ').replace('_',' + ') for p in pairs],fontsize=10)
 for i in range(len(pairs)):
  for j in range(len(cols)):ax.text(j,i,f'{int(counts[i,j])}/{int(den[i,j])}\n{matrix[i,j]:.1f}%',ha='center',va='center',fontsize=8,color='white' if matrix[i,j]>65 else '#222222')
 fig.colorbar(im,ax=ax,pad=.03,fraction=.025,label='RNA co-detected adjacent edges (%)');ax.set_title('Candidate ligand–receptor RNA support: myeloid source → T-enriched neighbor',fontsize=14,pad=15);fig.text(.015,.02,'Ligand and all receptor subunits must be detected. Zero-support edges retained. Edges may share spots.\nRNA co-detection does not verify cellular origin, protein interaction, communication direction or immunosuppression.',fontsize=10);fig.tight_layout(rect=[0,.065,1,1]);save(fig,'LR_edge_codetection_heatmap')
 # Specific activation/exhaustion groups: descriptive denominators remain explicit.
 cols=[(s,state) for s in SAMPLES for state in ['Activation_supported','Exhaustion_supported','Both_programs','T_other']];m=np.full((len(PAIRS),len(cols)),np.nan);text={}
 for i,(l,r,_) in enumerate(PAIRS):
  for j,(s,st) in enumerate(cols):
   z=detail[(detail['sample']==s)&(detail.source_group=='PGAM5_high')&(detail.receiver_T_state==st)&(detail.ligand==l)&(detail.receptor_complex==r)].iloc[0];m[i,j]=z.codetection_pct;text[i,j]=f'{z.n_codetected_edges}/{z.n_edges}' if z.n_edges else 'No edges'
 fig,ax=plt.subplots(figsize=(15,11));cm=plt.get_cmap('YlOrRd').copy();cm.set_bad('#DADDE1');im=ax.imshow(m,cmap=cm,vmin=0,vmax=100,aspect='auto');ax.set_xticks(range(len(cols)),labels=[s+'\n'+st.replace('_supported','').replace('_',' ') for s,st in cols],fontsize=10);ax.set_yticks(range(len(PAIRS)),labels=[l+' → '+r.replace('_',' + ') for l,r,_ in PAIRS],fontsize=10)
 for i in range(len(PAIRS)):
  for j in range(len(cols)):ax.text(j,i,text[i,j],ha='center',va='center',fontsize=8,color='white' if m[i,j]>65 else '#222222')
 ax.set_title('PGAM5-high sources: RNA edge support by provisional T state',fontsize=15,pad=15);fig.colorbar(im,ax=ax,fraction=.025,pad=.03,label='Co-detection (%)');fig.text(.015,.02,'Numbers = supported edges / eligible edges. Grey = no eligible edge, not measured absence of signaling.\nState categories are mixed-spot RNA proxies; receptor selection can overlap state markers. Descriptive only.',fontsize=10);fig.tight_layout(rect=[0,.06,1,1]);save(fig,'LR_by_activation_exhaustion_state')

def edge_examples():
 fig,axs=plt.subplots(2,2,figsize=(12,11))
 for row,s in enumerate(SAMPLES):
  d=DATA[s];edges=pd.read_csv(ROOT/f'{s}_directed_LR_edges.csv.gz')
  for col,key in enumerate(['SPP1--CD44','PVR--TIGIT']):
   ax=axs[row,col];xy=circles(ax,s,['#DADDE1']*len(d));ix={b:i for i,b in enumerate(d.index)};e=edges[(edges.source_group=='PGAM5_high')&(edges[key]>0)]
   seg=[np.array([xy[ix[r.source_barcode]],xy[ix[r.target_barcode]]]) for r in e.itertuples()];ax.add_collection(LineCollection(seg,colors='#D95F02',linewidths=1.7,alpha=.9))
   high=d.source_group=='PGAM5_high';ax.scatter(xy[high,0],xy[high,1],s=15,c=COL['PGAM5_high'],marker='D',edgecolors='white',linewidths=.2);ax.scatter(xy[d.T_enriched,0],xy[d.T_enriched,1],s=9,facecolors='none',edgecolors='#2166AC',linewidths=.6);ax.set_title(f'{s}: {key.replace("--"," → ")}\n{len(e)} RNA-supported high-source edges',fontsize=13)
 fig.suptitle('Spatial positions of candidate RNA-supported edges',fontsize=17);fig.legend(handles=[Patch(fc=COL['PGAM5_high'],label='PGAM5-high source RNA point'),Patch(fc='#2166AC',label='T-enriched RNA point'),Patch(fc='#D95F02',label='Adjacent co-detection')],loc='lower center',ncol=3,frameon=False,bbox_to_anchor=(.5,.035));fig.text(.025,.01,'Source/target is the analytic orientation of a mixed-spot edge, not a demonstrated molecular signaling direction.',fontsize=10);fig.subplots_adjust(top=.9,bottom=.1,wspace=.12,hspace=.17);save(fig,'candidate_LR_spatial_edges')
if __name__=='__main__':states();neighbor();LRheat();edge_examples();print('Figures complete')
