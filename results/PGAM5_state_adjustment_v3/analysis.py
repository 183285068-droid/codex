from pathlib import Path
import importlib.util,json,hashlib
import numpy as np,pandas as pd
from scipy.optimize import nnls,lsq_linear
from scipy.stats import rankdata,pearsonr,spearmanr
R=Path(__file__).resolve().parent;OLD=Path('/workspace/scratch/PGAM5_myeloid_reference_v2')
spec=importlib.util.spec_from_file_location('v2',OLD/'train_joint_reference.py');v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
TARGET=v.TARGET;NEG=v.NEG
# Fixed thresholds; panels/folds specified before this run's external performance.
protocol={'target':'PGAM5-related continuous macrophage state alongside retained RNA count>0 abundance target','training':'GSE151530+GSE149614','retest':'GSE189903 previously inspected; retrospective reuse, not new independent validation','patient_pairing':'No paired DE. Donor summaries and leave-donor-out prediction validation only. No depth covariates/matching.','candidate_filter':'Exclude PGAM5/MT-/RPL/RPS. Eligible training donors >=5 detected and >=20 undetected macrophages. Detection>=10% in target, median per-donor log2((meanCP10k+.05)/(negative+.05))>=1, >=70% donors log2 contrast>=.5, each cohort median>=.5. This is reference feature screening, not new FDR-qualified DE.','specific_filter':'Training target donor-equal mean >=2x each competing fine-group donor-equal mean with pseudocount .05; includes other macrophage states, mono/DC/APC/tumor.','panels':['stable','stable_no_cycle','specific','specific_no_cycle','prior_consensus','prior_consensus_no_cycle'],'cycle':'Fixed finite S/G2M list + documented extension (cell_cycle_provenance.json), not complete GO ontology. Excludes features only, never removes cells.','models':['direct_NNLS','hierarchical_NNLS'],'state_score':'Mean standardized log1p(CP10k) of training-selected genes; training macrophage cell means/SD. PGAM5 excluded. No score when <5 genes. Positive/zero classification AUC and PGAM5 expression correlation are descriptive and not biological gold standard.','gates':'Same MAE<=1pp,r>=.7,zero95q<=.5%,external>=3; independent pureDC>=2 remains necessary. State recommendation requires >=5 genes, testable >=3 external donors, >=70% external donors AUC>=.7 and no pure nonmacrophage pool median score exceeding external macrophage median. Exploratory retrospective protocol; cannot establish clinical utility.'}
(R/'protocol.json').write_text(json.dumps(protocol,indent=2))
cycle=set((R/'cell_cycle_genes.txt').read_text().split());gi={g:i for i,g in enumerate(v.genes)};cy=np.array([g in cycle for g in v.genes]);base=v.feature_eligible.copy();base[v.pg]=False
mac_types=[g for g in v.types if g in [TARGET,NEG,'Cycling_'+NEG]]
def evidence(ds):
 ix=[v.donors.index(d) for d in ds];n=v.num[ix];s=v.sums[ix];de=v.det[ix];p=v.typeindex[TARGET];nn=[v.typeindex[g] for g in [NEG,'Cycling_'+NEG] if g in v.typeindex];pc=n[:,p];nc=n[:,nn].sum(1);eligible=(pc>=5)&(nc>=20);px=s[eligible,p]/pc[eligible,None];nx=s[eligible][:,nn].sum(1)/nc[eligible,None];lf=np.log2((px+.05)/(nx+.05));coh=np.array([d.split(':')[0] for d in ds])[eligible];med=np.median(lf,axis=0);fraction=(lf>=.5).mean(0);freq=de[:,p].sum(0)/pc.sum();pm=(s[pc>0,p]/pc[pc>0,None]).mean(0)
 refs=[];labels=[]
 for j,g in enumerate(v.types):
  if g==TARGET or n[:,j].sum()<50 or (n[:,j]>0).sum()<2:continue
  z=n[:,j]>0;refs.append((s[z,j]/n[z,j,None]).mean(0));labels.append(g)
 comps=np.column_stack(refs);maxc=comps.max(1);contrast=np.log2((pm+.05)/(maxc+.05));cm=[np.median(lf[coh==c],axis=0) for c in sorted(set(coh))];valid=base&(freq>=.1)&(med>=1)&(fraction>=.7)&(np.min(cm,axis=0)>=.5)
 ev=pd.DataFrame({'gene':v.genes,'median_donor_log2_contrast':med,'donor_fraction_log2_contrast_ge_0_5':fraction,'target_detection_fraction':freq,'target_mean_CP10k':pm,'max_competitor_mean_CP10k':maxc,'max_competitor':np.array(labels)[comps.argmax(1)],'specific_log2_contrast':contrast,'cell_cycle':cy,'stable_candidate':valid,'specific_candidate':valid&(contrast>=1)})
 return ev,lf,np.array(ds)[eligible]
full,lf,ed=evidence(v.donors);full.to_csv(R/'training_gene_stability_specificity.csv',index=False);pd.DataFrame(lf,index=ed,columns=v.genes).rename_axis('donor').to_csv(R/'per_donor_gene_contrasts.csv.gz',compression='gzip')
prior_genes='LMNB1 TK1 BIRC5 CDK1 CENPH CENPK CLSPN FEN1 GINS2 HCFC1 MAD2L1 MKI67 NUSAP1 PRC1 RFC4 RRM2 TOP2A TPX2 TYMS UBE2T ZWINT'.split()
protocol['prior_consensus']='Exploratory comparator: fixed previous 21-gene >=3/5-cohort consensus. Not required to pass new stability/specificity screens; previous analyses included external GSE189903, therefore this panel has prior external information and is NOT held-out feature discovery. Added after training-only screening showed insufficient candidates, before this run external performance. No validated state/abundance claim.'
(R/'protocol.json').write_text(json.dumps(protocol,indent=2))
panels={};scores=[];cell_scores=[];pred=[];coverage=[];coefchecks=[];audit=[];selected=[]
def panel(ev,name):
 if name.startswith('prior_consensus'):
  zz=np.array([gi[g] for g in prior_genes if g in gi],dtype=int)
  return zz[~cy[zz]] if name.endswith('no_cycle') else zz
 mask=ev.stable_candidate.to_numpy().copy() if name.startswith('stable') else ev.specific_candidate.to_numpy().copy()
 if name.endswith('no_cycle'):mask &= ~ev.cell_cycle.to_numpy()
 idx=np.flatnonzero(mask);return idx[np.argsort(-ev.median_donor_log2_contrast.to_numpy()[idx],kind='stable')[:100]]
def collapsed_model(ds):
 ix=[v.donors.index(d) for d in ds];n=v.num[ix];s=v.sums[ix];labels=[('Macrophage' if g in mac_types else v.parent(g)) for g in v.types];active=[];refs=[];freq=[]
 for g in sorted(set(labels)):
  j=[k for k,z in enumerate(labels) if z==g];nc=n[:,j].sum(1);su=s[:,j].sum(1);de=v.det[ix][:,j].sum(1)
  if nc.sum()<50 or (nc>0).sum()<2:continue
  active.append(g);refs.append((su[nc>0]/nc[nc>0,None]).mean(0));freq.append(de.sum(0)/nc.sum())
 M=np.column_stack(refs);fr=np.column_stack(freq);features=set()
 for j,g in enumerate(active):
  contrast=np.log2((M[:,j]+.05)/(np.delete(M,j,1).max(1)+.05));mask=v.feature_eligible&(fr[:,j]>=.1)&(M[:,j]>=.05);z=np.flatnonzero(mask);features.update(z[np.argsort(-(contrast*np.minimum(1,M[:,j]))[z],kind='stable')[:100]].tolist())
 ix=np.array(sorted(features));scale=np.maximum(np.sqrt(np.mean(M[ix]**2,axis=1)),1e-6)
 return active,M,ix,scale

def state_eval(ds,name,d,ix,phase):
 if len(ix)<5:return
 train=v.D['joint_discovery'];o=train['obs'];mask=o.donor.isin(ds)&o.harmonized_group.eq('Macrophage');x=train['norm'][mask.to_numpy()][:,ix].toarray();x=np.log1p(x);mu=x.mean(0);sd=np.maximum(x.std(0),.1)
 data=v.D[name];o=data['obs'];sel=o.donor.eq(d).to_numpy();x=np.log1p(data['norm'][sel][:,ix].toarray());sc=((x-mu)/sd).mean(1);t=o.loc[sel].copy();t['state_score']=sc
 cell_scores.append(t[['cell_id','donor','harmonized_group','PGAM5_counts','state_score']].assign(dataset=name,phase=phase,panel=current_panel))
 for g,a in t.groupby('harmonized_group',observed=True):scores.append({'phase':phase,'dataset':name,'donor':d,'panel':current_panel,'group':g,'cells':len(a),'median_score':a.state_score.median(),'mean_score':a.state_score.mean(),'p95_score':a.state_score.quantile(.95)})
 a=t[t.harmonized_group.eq('Macrophage')];y=a.fine_group.eq(TARGET).to_numpy();z=a.state_score.to_numpy();np_=y.sum();nn=(~y).sum();auc=(rankdata(z)[y].sum()-np_*(np_+1)/2)/(np_*nn) if np_ and nn else np.nan
 pg=np.log1p(data['norm'][sel][:,v.pg].toarray().ravel());rho=spearmanr(z,pg[t.harmonized_group.eq('Macrophage').to_numpy()]).statistic if len(a)>1 else np.nan
 coverage.append({'phase':phase,'dataset':name,'donor':d,'panel':current_panel,'selected_genes':len(ix),'macrophages':len(a),'PGAM5_detected':int(np_),'AUC_RNA_detection':auc,'Spearman_score_PGAM5':rho})

for phase,name,dd in [('internal_leave_donor_out','joint_discovery',v.internal),('external_retest','GSE189903',v.external)]:
 for d in dd:
  ds=[k for k in v.donors if k!=d] if name=='joint_discovery' else v.donors
  ev,_,_=evidence(ds);m=v.model(ds,'donor_equal_50');groups,M,i1,sc1=collapsed_model(ds);mix=v.mixtures(name,d,replicates=5)
  for current_panel in protocol['panels']:
   ix=panel(ev,current_panel);selected.extend({'phase':phase,'donor':d,'panel':current_panel,'gene':v.genes[k]} for k in ix)
   if phase=='external_retest':panels[current_panel]=v.genes[ix].tolist()
   audit.append({'phase':phase,'donor':d,'panel':current_panel,'genes':len(ix),'status':'evaluated' if len(ix)>=5 else 'UNTESTABLE: fewer than 5 training-supported genes'})
   if len(ix)<5:continue
   state_eval(ds,name,d,ix,phase)
   # Direct model adds broad lineage features; target-only signature must be tested alone inside macrophages in stage2.
   ii=np.union1d(m['ix'],ix);MM=m['ref'][ii];ss=np.maximum(np.sqrt(np.mean(MM**2,axis=1)),1e-6)
   # Fair cycle sensitivity removes cycle genes from BOTH broad and state features for no_cycle panels.
   if current_panel.endswith('no_cycle'):
    kk=~cy[ii];ii=ii[kk];MM=MM[kk];ss=ss[kk];jj=i1[~cy[i1]];s1=sc1[~cy[i1]]
   else:jj=i1;s1=sc1
   pos=m['ref'][:,m['groups'].index(TARGET)];di=[v.donors.index(k) for k in ds];ni=[v.typeindex[g] for g in [NEG,'Cycling_'+NEG] if g in v.typeindex];nc=v.num[di][:,ni].sum(1);su=v.sums[di][:,ni].sum(1);neg=(su[nc>0]/nc[nc>0,None]).mean(0);stage=np.column_stack([pos[ix],neg[ix]]);s2=np.maximum(np.sqrt(np.mean(stage**2,axis=1)),1e-6)
   for z in mix:
    y=z['y'];b,r=nnls(MM/ss[:,None],y[ii]/ss,maxiter=1000);direct=b[m['groups'].index(TARGET)]/max(b.sum(),1e-12)
    b1,r1=nnls(M[jj]/s1[:,None],y[jj]/s1,maxiter=1000);mi=groups.index('Macrophage');totalmac=b1[mi]/max(b1.sum(),1e-12);other=[j for j in range(len(groups)) if j!=mi];residual=y[ix]-M[ix][:,other]@b1[other];b2,r2=nnls(stage/s2[:,None],residual/s2,maxiter=1000);conditional=b2[0]/max(b2.sum(),1e-12);hier=totalmac*conditional
    info={k:w for k,w in z.items() if k!='y'}|{'phase':phase,'dataset':name,'donor':d,'panel':current_panel,'state_genes':len(ix)}
    pred.append(info|{'model':'direct_NNLS','predicted':direct,'residual':r});pred.append(info|{'model':'hierarchical_NNLS','predicted':hier,'estimated_total_macrophage':totalmac,'estimated_conditional_state_fraction':conditional,'residual':r2})
    if phase=='external_retest' and z['scenario']=='standard' and z['replicate']==0 and z['unit']=='equal_RNA_cell_fraction':
     alt=lsq_linear(stage/s2[:,None],residual/s2,bounds=(0,np.inf),method='bvls',tol=1e-10);assert alt.success;assert np.allclose(alt.x,b2,atol=1e-6);assert np.isclose(np.linalg.norm(stage/s2[:,None]@alt.x-residual/s2),r2,atol=1e-6);coefchecks.append({'donor':d,'panel':current_panel,'truth':z['truth'],'NNLS_residual':r2,'bounded_residual':float(np.linalg.norm(stage/s2[:,None]@alt.x-residual/s2))})
  print('Completed',phase,d,flush=True)
# All external donors can inform macrophage state discrimination even if ineligible for mixture dose-response.
for current_panel in protocol['panels']:
 ix=panel(full,current_panel)
 for d in sorted(v.D['GSE189903']['obs'].donor.unique()):
  if d not in v.external:state_eval(v.donors,'GSE189903',d,ix,'external_state_low_positive_support')
pd.DataFrame(audit).to_csv(R/'panel_fold_testability.csv',index=False);pd.DataFrame(selected).to_csv(R/'fold_selected_genes.csv',index=False);pd.DataFrame(pred).to_csv(R/'predictions.csv',index=False);pd.DataFrame(scores).to_csv(R/'state_score_by_donor_lineage.csv',index=False);pd.DataFrame(coverage).to_csv(R/'state_discrimination.csv',index=False);pd.DataFrame(coefchecks).to_csv(R/'bounded_solver_checks.csv',index=False);(R/'locked_state_panels.json').write_text(json.dumps(panels,indent=2))
summary=[]
if pred:
 t=pd.DataFrame(pred)
 for keys,a in t.groupby(['phase','panel','model','unit','scenario']):
  zero=a[a.truth.eq(0)];r=pearsonr(a.truth,a.predicted).statistic if a.truth.nunique()>1 and a.predicted.nunique()>1 else None
  summary.append(dict(zip(['phase','panel','model','unit','scenario'],keys))|{'donors':a.donor.nunique(),'mixtures':len(a),'MAE_pp':abs(a.predicted-a.truth).mean()*100,'Pearson_r':r,'zero_p95_percent':zero.predicted.quantile(.95)*100 if len(zero) else None})
pd.DataFrame(summary).to_csv(R/'validation_summary.csv',index=False)
if cell_scores:pd.concat(cell_scores,ignore_index=True).to_csv(R/'heldout_cell_state_scores.csv.gz',index=False,compression='gzip')
params=[]
o=v.D['joint_discovery']['obs'];sel=o.harmonized_group.eq('Macrophage').to_numpy()
for name in protocol['panels']:
 ix=panel(full,name)
 if len(ix)<5:continue
 x=np.log1p(v.D['joint_discovery']['norm'][sel][:,ix].toarray());mu=x.mean(0);sd=np.maximum(x.std(0),.1)
 params.extend({'panel':name,'gene':v.genes[k],'training_mean_log1p_CP10k':float(a),'training_SD_floor_0_1':float(b),'mean_weight':1/len(ix)} for k,a,b in zip(ix,mu,sd))
pd.DataFrame(params).to_csv(R/'locked_state_score_parameters.csv',index=False)
# Synthetic matched two-state reference controls verify numerical recovery, not biological validity.
exact=[];rng=np.random.default_rng(20261009)
for j in range(20):
 w=rng.dirichlet([1,1]);b,res=nnls(stage/s2[:,None],(stage@w)/s2);err=float(abs(b-w).max());assert err<1e-6;exact.append({'case':j,'max_coefficient_error':err,'residual':res})
pd.DataFrame(exact).to_csv(R/'exact_two_state_recovery.csv',index=False)

result={'training_stable_genes':int(full.stable_candidate.sum()),'training_specific_genes':int(full.specific_candidate.sum()),'training_stable_noncycle_genes':int((full.stable_candidate&~full.cell_cycle).sum()),'training_specific_noncycle_genes':int((full.specific_candidate&~full.cell_cycle).sum()),'internal_donors':len(v.internal),'external_mixture_donors':len(v.external),'external_DC_challenge_donors':0,'models_tested':len(pred)>0,'TCGA_abundance':'NOT_RUN: external coverage gate remains unmet','limits':'Retrospective external reuse. Scores proxy PGAM5 RNA detection, not protein validation. No depth adjustment; depth-associated signal remains possible. No paired DE.'}
(R/'result.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
