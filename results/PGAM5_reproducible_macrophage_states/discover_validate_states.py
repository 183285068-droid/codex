from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
from scipy import sparse
from scipy.spatial.distance import cdist
from scipy.stats import fisher_exact
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parent;V=Path('/workspace/scratch/PGAM5_myeloid_reference_v2');names=['GSE151530','GSE149614','GSE189903','GSE242889','GSE202642'];train_names=names[:2];validation=['GSE189903','GSE242889'];supplementary='GSE202642'
protocol={'goal':'Re-establish an RNA-expression-defined reproducible macrophage state associated with PGAM5, not assume prior count>0 group is a subtype','discovery':train_names,'primary_external':['GSE189903','GSE242889'],'supplementary_external':'GSE202642 HCC tumors only; sample identifiers are proxies, independent patients unverified','history':'All cohorts explored in earlier DE. GSE189903 also used in prior model checks. GSE242889/GSE202642 not used for prior state performance. This is retrospective, not prospectively blinded. No external expression/PGAM5 outcome used in state feature or state selection','macrophage_identity':'Use existing harmonized Macrophage labels for151/149/189.242/202 old marker-inferred candidates reclustered PGAM5-free, accept C1QA/B/C fractions>=.4,CD68>=.3,CSF1R>=.2,TYROBP>=.5, exclude dual CD1C>=.3 andFCER1A>=.2. No ambient/doublet correction claimed','normalization':'Raw symbol counts duplicate-summed, CP10k using original full-cell library count, log1p. No depth matching/covariates, no paired donor DE','cycling_state':'Fixed separate Cycling_macrophage: >=2 detected genes among MKI67/TOP2A/UBE2C/CENPF/CDK1/BIRC5. Defined without PGAM5. Never inferred as PGAM5-specific without association checks','noncycling_features':'PGAM5, finite recorded cell-cycle list,MT-/RPL/RPS excluded. Macrophage detection>=10%; exclude strong nonmacrophage background (max supported other-cell mean>4x macrophage donor-equal mean) and explicit hepatocyte/Ig/T/stromal marker list. Retain shared nuclear metabolic genes: within-macrophage state discovery does not demand all features be macrophage-specific. Nonmacrophage donor means only groups>=20 cells, fine-groups>=50 cells and>=2 donors. Rank training noncycling log-expression variance/mean; max1000 genes','clustering':'Patients with>=30 noncycling macrophages, random cap200 cells/patient seed0. Feature z-score training means,SD floor.1,clip[-10,10]; PCA30 randomized seed0;20NN;Leiden.6 igraph2iterations seed0. No batch regression/integration. Seeds1/2 sensitivity on same graph, not biological independent validation','projection':'Fixed nearest PCA centroid. Reject if squared distance>1.5*training cluster95th percentile or nearest/second margin<.02; labelUnknown_state. Cycling rule applied first. Freeze before external association results','donor_association':'For each state vs all other retained macrophages, donor evaluable if state>=10, rest>=20,totalPGAM5 detected>=3. CP10k mean log2 ratio uses+.05, rate ratio usesJeffreys+.5/+1. Zero PGAM5 donors retained in state presence, excluded only from informative association calculation','cohort_support':'State>=30cells,>=3 evaluable donor units,>=5PGAM5 detected state cells; pooled raw detection rate ratio>=1.5; median donor normalizedPGAM5 ratio>=1.5;>=70% evaluable donors normalizedPGAM5 ratio>1; classifier assigned coverage>=70%','selection':'Only states satisfying BOTH discovery cohorts selected. Ifnone, no PGAM5-related state named. All state-by-cohort results exported regardless. Primary validation requires same selected state passes BOTH189903 and242889;202642 supplementary only','technical_support':'For noncycling state best-match graph-seed Jaccard>=.6 atbothseeds; cycling fixed rule has no random clustering seed. Protein/function/true patient independence and depth confounding unresolved. RNA gatepass does not establish bulk-deconvolution usability','statistics':'Fisher and BH across state/cohort rows are descriptive cell-level checks (cells not independent biological replicates); primary evidence is donor consistency and cross-cohort reproducibility. Prior strictDEtables unchanged. Marker lists are descriptive training mean contrasts, not new FDR-qualified DE.'}
(R/'protocol.json').write_text(json.dumps(protocol,indent=2));loaded={}
for n in names:
 p=V/(n+'_harmonized_counts_QC.h5ad') if n in names[:3] else R/(n+'_marker_consistent_macrophages.h5ad');a=ad.read_h5ad(p);a.obs['cell_id']=a.obs_names.astype(str);a.obs.index.name='barcode_index';a.obs['dataset']=n
 if n in names[:3]:a.obs['donor_kind']='patient_or_author_patient_proxy'
 loaded[n]=a
common=sorted(set.intersection(*(set(a.var_names) for a in loaded.values())));assert 'PGAM5' in common and len(common)>12000;genes=np.array(common);pg=common.index('PGAM5');cycle_genes=['MKI67','TOP2A','UBE2C','CENPF','CDK1','BIRC5'];cy_ix=[common.index(g) for g in cycle_genes];D={};fullprofiles=[]
for n,a in loaded.items():
 x=a[:,common].X.astype(np.float32).tocsr();norm=x.multiply((1e4/a.obs.total_counts.to_numpy(dtype=np.float32))[:,None]).tocsr();o=a.obs.copy()
 if n in train_names:
  for g,t in o.groupby('fine_group',observed=True):
   if g.startswith('TAM_') or g=='Cycling_TAM_PGAM5_undetected':continue
   means=[]
   for d,z in t.groupby('donor_id',observed=True):
    if len(z)>=20:means.append(np.asarray(norm[o.index.isin(z.index)].mean(0)).ravel())
   if len(t)>=50 and len(means)>=2:fullprofiles.append(np.mean(means,axis=0))
 if n in names[:3]:sel=o.harmonized_group.eq('Macrophage').to_numpy();o=o.loc[sel].copy();x=x[sel];norm=norm[sel]
 o['PGAM5_raw']=x[:,pg].toarray().ravel();o['PGAM5_CP10k']=norm[:,pg].toarray().ravel();o['cycling_flag']=np.asarray((x[:,cy_ix]>0).sum(1)).ravel()>=2;D[n]={'raw':x,'norm':norm,'obs':o};print(n,'macrophages',len(o),'PGAM5detected',int(o.PGAM5_raw.gt(0).sum()),flush=True)
del loaded
obs=pd.concat([D[n]['obs'] for n in train_names]);obs.index=obs.dataset.astype(str)+':'+obs.cell_id.astype(str);raw=sparse.vstack([D[n]['raw'] for n in train_names],format='csr');norm=sparse.vstack([D[n]['norm'] for n in train_names],format='csr');log=norm.copy();log.data=np.log1p(log.data)
macro_means=[]
for d,t in obs.groupby('donor_id',observed=True):
 if len(t)>=30:macro_means.append(np.asarray(norm[obs.index.isin(t.index)].mean(0)).ravel())
mp=np.mean(macro_means,axis=0);maxother=np.max(np.vstack(fullprofiles),axis=0);freq=np.asarray((raw>0).mean(0)).ravel();cycles=set((R/'cell_cycle_genes.txt').read_text().split());background=set('ALB APOA1 APOA2 APOC3 TTR RBP4 APOH AHSG FGB FGA FGG ORM1 ORM2 PRAP1 TF GSTA1 CES1 CD3D CD3E TRAC CD79A MS4A1 COL1A1 COL1A2 DCN PECAM1 VWF'.split());(R/'excluded_background_genes.txt').write_text('\n'.join(sorted(background))+'\n');eligible=np.array([g!='PGAM5' and g not in cycles and g not in background and not g.startswith(('MT-','RPL','RPS','IGH','IGK','IGL','HBA','HBB')) for g in genes])&(mp*4>=maxother)&(freq>=.1)
rng=np.random.default_rng(0);chosen=[]
for d,t in obs.loc[~obs.cycling_flag].groupby('donor_id',observed=True):
 if len(t)>=30:
  ii=np.flatnonzero(obs.index.isin(t.index));chosen.extend(rng.choice(ii,min(200,len(ii)),replace=False).tolist())
chosen=np.array(sorted(chosen));mu=np.asarray(log[chosen].mean(0)).ravel();var=np.maximum(0,np.asarray(log[chosen].multiply(log[chosen]).mean(0)).ravel()-mu*mu);rank=var/(mu+.05);ix=np.flatnonzero(eligible&(mu>.01));ix=ix[np.argsort(-rank[ix],kind='stable')[:1000]];assert len(ix)>=50
pd.DataFrame({'gene':genes,'macro_donor_equal_CP10k':mp,'max_nonmacro_donor_equal_CP10k':maxother,'macro_detection':freq,'eligible':eligible,'chosen_feature':np.isin(np.arange(len(genes)),ix)}).to_csv(R/'training_feature_evidence.csv',index=False)
xx=log[chosen][:,ix].toarray();mean=xx.mean(0);sd=np.maximum(xx.std(0),.1);scaled=np.clip((xx-mean)/sd,-10,10);pca=PCA(n_components=min(30,len(ix)-1),svd_solver='randomized',random_state=0);pc=pca.fit_transform(scaled);b=ad.AnnData(scaled,obs=obs.iloc[chosen].copy());b.obsm['X_pca']=pc;sc.pp.neighbors(b,n_neighbors=20,use_rep='X_pca',random_state=0);sc.tl.leiden(b,resolution=.6,flavor='igraph',directed=False,n_iterations=2,random_state=0,key_added='state');lab=b.obs.state.astype(str).to_numpy();labels=sorted(set(lab),key=int);cent=np.vstack([pc[lab==k].mean(0) for k in labels]);radius=np.array([np.quantile(np.sum((pc[lab==k]-cent[j])**2,axis=1),.95) for j,k in enumerate(labels)]);seedrows=[];aris=[]
for seed in [1,2]:
 sc.tl.leiden(b,resolution=.6,flavor='igraph',directed=False,n_iterations=2,random_state=seed,key_added='seed'+str(seed));alt=b.obs['seed'+str(seed)].astype(str).to_numpy();aris.append({'seed':seed,'ARI_vs_seed0':float(adjusted_rand_score(lab,alt))})
 for k in labels:
  jac=[((lab==k)&(alt==z)).sum()/((lab==k)|(alt==z)).sum() for z in set(alt)];seedrows.append({'state':'State_'+k,'seed':seed,'best_match_Jaccard':float(max(jac))})
pd.DataFrame(seedrows).to_csv(R/'graph_state_seed_stability.csv',index=False);pd.DataFrame(aris).to_csv(R/'graph_global_seed_stability.csv',index=False)
model={'features':genes[ix].tolist(),'feature_mean':mean.tolist(),'feature_SD':sd.tolist(),'PCA_mean':pca.mean_.tolist(),'PCA_components':pca.components_.tolist(),'centroid_states':['State_'+k for k in labels],'PCA_centroids':cent.tolist(),'squared_distance_q95':radius.tolist(),'squared_distance_multiplier':1.5,'minimum_margin':.02,'cycling_genes':cycle_genes,'cycling_detected_gene_min':2,'status':'LOCKED_BEFORE_EXTERNAL_PGAM5_ASSOCIATION'};(R/'locked_state_classifier.json').write_text(json.dumps(model,indent=2));pd.DataFrame(pca.components_.T,index=genes[ix],columns=['PC'+str(j+1) for j in range(pc.shape[1])]).to_csv(R/'PCA_gene_loadings.tsv',sep='\t')
for n in names:
 a=D[n];x=a['norm'][:,ix].toarray();q=pca.transform(np.clip((np.log1p(x)-mean)/sd,-10,10));dist=cdist(q,cent,'sqeuclidean');order=np.argsort(dist,axis=1);best=order[:,0];d1=dist[np.arange(len(q)),best];d2=dist[np.arange(len(q)),order[:,1]];margin=(d2-d1)/np.maximum(d2,1e-12);state=np.array(['State_'+labels[j] for j in best],dtype=object);state[(d1>1.5*radius[best])|(margin<.02)]='Unknown_state';state[a['obs'].cycling_flag.to_numpy()]='Cycling_macrophage';a['obs']['state']=state;a['obs']['mapping_distance_squared']=d1;a['obs']['mapping_margin']=margin;a['obs'].to_csv(R/(n+'_macrophage_state_metadata.csv.gz'));print(n,'statecounts',a['obs'].state.value_counts().to_dict(),flush=True)
# Freeze candidates based on discovery, then separately show validation; no state identity/config changes.
state_names=['State_'+k for k in labels]+['Cycling_macrophage'];cohort=[];donorrows=[];presence=[]
for n in names:
 o=D[n]['obs'];coverage=float(o.state.ne('Unknown_state').mean())
 for state in state_names:
  inside=o.state.eq(state);s=o[inside];rest=o[~inside];a=int(s.PGAM5_raw.gt(0).sum());bb=len(s)-a;c=int(rest.PGAM5_raw.gt(0).sum());dd=len(rest)-c;pval=fisher_exact([[a,bb],[c,dd]],alternative='two-sided').pvalue if len(s) and len(rest) else np.nan;dr=[]
  for donor,t in o.groupby('donor_id',observed=True):
   tt=t[t.state.eq(state)];rr=t[t.state.ne(state)];p1=int(tt.PGAM5_raw.gt(0).sum());p0=int(rr.PGAM5_raw.gt(0).sum());valid=len(tt)>=10 and len(rr)>=20 and p1+p0>=3;ratio=(tt.PGAM5_CP10k.mean()+.05)/(rr.PGAM5_CP10k.mean()+.05) if len(tt) and len(rr) else np.nan;rate=((p1+.5)/(len(tt)+1))/((p0+.5)/(len(rr)+1)) if len(tt) and len(rr) else np.nan
   row={'dataset':n,'donor':donor,'donor_kind':t.donor_kind.iloc[0],'state':state,'state_cells':len(tt),'other_macrophages':len(rr),'state_PGAM5_detected':p1,'other_PGAM5_detected':p0,'state_mean_PGAM5_CP10k':tt.PGAM5_CP10k.mean(),'other_mean_PGAM5_CP10k':rr.PGAM5_CP10k.mean(),'PGAM5_CP10k_ratio':ratio,'PGAM5_detection_rate_ratio_Jeffreys':rate,'association_evaluable':valid};donorrows.append(row);presence.append({'dataset':n,'donor':donor,'state':state,'cells':len(tt)})
   if valid:dr.append(row)
  median=float(np.median([z['PGAM5_CP10k_ratio'] for z in dr])) if dr else np.nan;frac=float(np.mean([z['PGAM5_CP10k_ratio']>1 for z in dr])) if dr else 0;rate=(a/len(s))/(c/len(rest)) if len(s) and len(rest) and c>0 else np.nan
  cohort.append({'dataset':n,'role':'discovery' if n in train_names else 'supplementary_sample_proxies' if n==supplementary else 'external_validation','state':state,'state_cells':len(s),'other_cells':len(rest),'state_PGAM5_detected':a,'other_PGAM5_detected':c,'state_PGAM5_detected_fraction':a/len(s) if len(s) else np.nan,'other_PGAM5_detected_fraction':c/len(rest) if len(rest) else np.nan,'pooled_detection_rate_ratio':rate,'state_mean_PGAM5_CP10k':s.PGAM5_CP10k.mean(),'other_mean_PGAM5_CP10k':rest.PGAM5_CP10k.mean(),'evaluable_donor_units':len(dr),'median_donor_PGAM5_CP10k_ratio':median,'fraction_evaluable_donors_enriched':frac,'state_mapping_coverage':coverage,'Fisher_cell_level_p':pval,'cohort_support_passed':bool(len(s)>=30 and len(dr)>=3 and a>=5 and pd.notna(rate) and rate>=1.5 and pd.notna(median) and median>=1.5 and frac>=.7 and coverage>=.7)})
co=pd.DataFrame(cohort);co['Fisher_cell_level_BH_q']=np.nan;m=co.Fisher_cell_level_p.notna();co.loc[m,'Fisher_cell_level_BH_q']=multipletests(co.loc[m,'Fisher_cell_level_p'],method='fdr_bh')[1];co.to_csv(R/'state_PGAM5_association_by_cohort.csv',index=False);pd.DataFrame(donorrows).to_csv(R/'state_PGAM5_association_by_donor.csv',index=False);pd.DataFrame(presence).to_csv(R/'state_presence_by_donor.csv',index=False)
selected=[];accepted=[]
for state in state_names:
 st=co[(co.state==state)&co.dataset.isin(train_names)];tech=True if state=='Cycling_macrophage' else min(z['best_match_Jaccard'] for z in seedrows if z['state']==state)>=.6
 if len(st)==2 and st.cohort_support_passed.all() and tech:
  selected.append(state);ext=co[(co.state==state)&co.dataset.isin(validation)]
  if len(ext)==2 and ext.cohort_support_passed.all():accepted.append(state)
# Training-only descriptive marker profiles, not statistical DE.
profiles=[]
for n in train_names:
 o=D[n]['obs'];norm=D[n]['norm']
 for state in state_names:
  mask=o.state.eq(state).to_numpy()
  if not mask.sum() or mask.all():continue
  me=np.asarray(norm[mask][:,ix].mean(0)).ravel();other=np.asarray(norm[~mask][:,ix].mean(0)).ravel();det=np.asarray((D[n]['raw'][mask][:,ix]>0).mean(0)).ravel();contrast=np.log2((me+.05)/(other+.05));top=np.argsort(-contrast*np.minimum(1,me))[:30]
  profiles.extend({'dataset':n,'state':state,'gene':genes[ix[j]],'state_mean_CP10k':float(me[j]),'other_mean_CP10k':float(other[j]),'state_detection_fraction':float(det[j]),'descriptive_log2_contrast':float(contrast[j])} for j in top)
pd.DataFrame(profiles).to_csv(R/'training_state_descriptive_markers.csv',index=False);b.obs.to_csv(R/'balanced_discovery_graph_cell_labels.csv.gz')
result={'discovery_macrophages':sum(len(D[n]['obs']) for n in train_names),'graph_training_cells':len(chosen),'classifier_features':len(ix),'noncycling_states':len(labels),'discovery_selected_PGAM5_related_states':selected,'cross_cohort_RNA_support_states':accepted,'PGAM5_related_group_established_at_RNA_level':bool(accepted),'supplementary202642_not_patient_independence_proof':True,'limits':'RetrospectiveRNA inference. No protein/function gold standard, no ambient/doublet correction, no depth controls. Neither RNA association nor reproducible state implies PGAM5-specific bulk deconvolution signature. Prior DE and TCGA analyses unchanged.'};(R/'result.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
