import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE281759_function_mpl')
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze import ROOT,PROGRAMS,COMPS,SAMPLES,PANELS
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
LABELS={'MHCII_presentation':'MHC-II antigen presentation','Inflammatory_myeloid':'Inflammatory myeloid RNA','CCL_CXCL8_chemokines':'CCL / CXCL8 chemokines','CXCR3_ligands':'CXCR3 ligands','Checkpoint_ligands':'Checkpoint ligands','Regulatory_mediators':'IL10 / TGFB1 / LGALS9','Costimulation_ligands':'Costimulation ligands','Scavenger_LPL_state':'Scavenger / LPL state','C1Q_complement':'C1Q complement (gate overlap)','Uptake_clearance_proxy':'Uptake proxy (gate overlap)','SPP1_MMP12_remodeling':'SPP1 / MMP12 remodeling','Hypoxia_response':'Hypoxia response','Cycling':'Cycling'}
def save(fig,name):
 fig.savefig(ROOT/'figures'/f'{name}.png',dpi=180,bbox_inches='tight');fig.savefig(ROOT/'figures'/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
def heat(mat,title,name,limit=1):
 fig,ax=plt.subplots(figsize=(max(8,len(mat.columns)*.6+4),max(5,len(mat)*.43+1)))
 im=ax.imshow(mat.to_numpy(float),cmap='RdBu_r',vmin=-limit,vmax=limit,aspect='auto');ax.set_xticks(range(len(mat.columns)),mat.columns,rotation=45 if len(mat.columns)>5 else 0,ha='right' if len(mat.columns)>5 else 'center');ax.set_yticks(range(len(mat)),[LABELS.get(k,k) for k in mat.index])
 for i in range(len(mat)):
  for j in range(len(mat.columns)):
   v=mat.iloc[i,j];label=(f'{v:+.3f}' if abs(v)<.01 else f'{v:+.2f}') if np.isfinite(v) else 'NA';ax.text(j,i,label,ha='center',va='center',fontsize=8,color='white' if abs(v)>limit*.65 else '#252525')
 ax.set_title(title,pad=15);fig.colorbar(im,ax=ax,shrink=.7,label='Adjusted coefficient / within-section outcome SD');fig.tight_layout();save(fig,name)
def run():
 s=pd.read_csv(ROOT/'program_statistics.csv');p=pd.read_csv(ROOT/'patient_equal_section_effects.csv');g=pd.read_csv(ROOT/'marker_gene_descriptive_contrasts.csv')
 primary=s[s.comparison.eq(COMPS[0])&s.status.eq('tested')];pp=p[p.comparison.eq(COMPS[0])];mat=pp.pivot(index='program',columns='patient',values='median_standardized_beta').reindex(PROGRAMS)
 heat(mat,'PGAM5-high vs all other macrophage-RNA spots\nEqual-section median adjusted effects; 4 independent patients','patient_adjusted_function_effects',.8)
 section=primary.pivot(index='program',columns='sample',values='standardized_beta').reindex(index=PROGRAMS,columns=[i['sample'] for i in SAMPLES]);heat(section,'Section-level adjusted RNA-program contrasts\nPositive = higher in PGAM5-high; not a measured functional activity','section_adjusted_function_effects',1)
 fig,axes=plt.subplots(1,2,figsize=(14,7),sharey=True);patients=['A009','B013','C004','C008'];colors=['#285A8E','#B5822C','#A8536A','#607C44']
 for c,comp in enumerate(COMPS[:2]):
  z=p[p.comparison.eq(comp)]
  for j,patient in enumerate(patients):
   v=z[z.patient.eq(patient)].set_index('program').reindex(PROGRAMS).median_standardized_beta;axes[c].scatter(v,np.arange(len(PROGRAMS))+(j-1.5)*.13,s=32,color=colors[j],marker=['o','s','^','D'][j],label=patient)
  axes[c].axvline(0,color='#555555',lw=.8);axes[c].set_title('All other macrophage-RNA' if c==0 else 'Other PGAM5-detected macrophage-RNA\nC004: not testable');axes[c].set_xlabel('Adjusted coefficient / outcome SD');axes[c].set_xlim(-1.2,1.2);axes[c].spines[['top','right']].set_visible(False)
 axes[0].set_yticks(np.arange(len(PROGRAMS)),[LABELS[k] for k in PROGRAMS]);axes[0].invert_yaxis();axes[1].legend(loc='lower right');fig.suptitle('Comparison choice and patient heterogeneity | Points are patient effects, not confidence intervals');fig.tight_layout();save(fig,'patient_comparison_sensitivity')
 # Unadjusted and adjusted section effects summarized identically, avoiding pooled-spot inference.
 means=pd.read_csv(ROOT/'program_group_means.csv');combined=[]
 for row in primary.itertuples():
  z=pd.read_csv(ROOT/f'{row.sample}_spot_program_scores.csv.gz',index_col=0);use=z.Mac_RNA_enriched.astype(bool);sd=z.loc[use,row.program].std()
  combined.append(dict(patient=row.patient,program=row.program,unadjusted=row.unadjusted_delta/sd,adjusted=row.standardized_beta))
 delta=pd.DataFrame(combined).groupby(['patient','program'])[['unadjusted','adjusted']].median().reset_index();columns=[]
 for patient in patients:
  for kind in ['unadjusted','adjusted']:columns.append(delta[delta.patient.eq(patient)].set_index('program')[kind].rename(patient+' '+kind))
 mm=pd.concat(columns,axis=1).reindex(PROGRAMS);heat(mm,'Unadjusted vs adjusted contrasts (equal-section patient medians)\nAdjustment includes depth, MT, spatial tiles and mixing RNA proxies','unadjusted_adjusted_comparison',1)
 for comp,name in [(COMPS[2],'nonT_function_sensitivity'),(COMPS[3],'MT5_function_sensitivity')]:
  z=p[p.comparison.eq(comp)];mm=z.pivot(index='program',columns='patient',values='median_standardized_beta').reindex(index=PROGRAMS,columns=patients);heat(mm,'Fixed PGAM5/source labels | '+('Exclude sources satisfying T-RNA gate' if comp.startswith('nonT') else 'Restrict to author QC MT<=5% subset')+'\nNA = insufficient testable observations',''+name,1)
 genes=[]
 for k in PROGRAMS:
  for gene in PANELS[k]:
   if gene not in genes:genes.append(gene)
 gg=g[g.comparison.eq(COMPS[0])].groupby(['gene','patient']).standardized_beta.median().unstack().reindex(genes);heat(gg,'Selected marker-level adjusted effects | Descriptive only\n179-marker input; no genome-wide differential-expression claims','marker_gene_effects',1)
 qc=pd.read_csv(ROOT/'group_quality_mixing.csv');variables=['total_counts','pct_mt','Macrophage_core','T','Hepatocyte','Plasma'];fig,axes=plt.subplots(2,3,figsize=(13,7))
 colors2=['#285A8E','#B5822C','#777777'];source=['PGAM5_high','Other_detected','Undetected']
 for ax,var in zip(axes.flat,variables):
  for j,group in enumerate(source):
   zz=qc[qc.variable.eq(var)&qc.source_group.eq(group)].groupby('patient')['mean'].median().reindex(patients);ax.bar(np.arange(4)+(j-1)*.24,zz,width=.23,color=colors2[j],label=group)
  ax.set_xticks(range(4),patients);ax.set_title(var);ax.spines[['top','right']].set_visible(False)
 axes[0,0].legend(fontsize=7);fig.suptitle('Quality and mixing proxies | Median of section means; not cell fractions');fig.tight_layout();save(fig,'quality_and_mixing_proxies')
if __name__=='__main__':run()
