from pathlib import Path
import importlib.util,json,warnings
import numpy as np,pandas as pd
from sklearn.decomposition import NMF
from scipy.stats import spearmanr,rankdata
from sklearn.exceptions import ConvergenceWarning
R=Path(__file__).resolve().parent;OLD=Path('/workspace/scratch/PGAM5_myeloid_reference_v2');spec=importlib.util.spec_from_file_location('v2',OLD/'train_joint_reference.py');v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
cy=set((R/'cell_cycle_genes.txt').read_text().split());eligible=np.array([g not in cy and g!='PGAM5' and not g.startswith(('MT-','RPL','RPS')) for g in v.genes]);outputs=[];genes_out=[];lineages=[];factor_rows=[];cellout=[];convergence=[]
protocol={'purpose':'Exploratory noncycling macrophage metaprograms potentially related to continuous PGAM5 expression; not abundance','fixed_parameters':'CP10k log1p; training-only top500 variance/mean genes excluding PGAM5, fixed cell-cycle list, MT/RPL/RPS; NMF rank4 MU Frobenius,seed0,max_iter300,tol1e-4; choose factor with greatest training Spearman correlation to log1p PGAM5 CP10k; top30 nonnegative loading genes; score=mean training-standardized log1p CP10k (SD floor .1)','validation':'Two reciprocal dataset-held-out fits plus joint fit retested on previously inspected GSE189903. Fixed all parameters; no tuning against test outcomes. Program fitted on all training macrophages without PGAM5 labels; factor selection uses training PGAM5. No depth covariate/matching, no paired DE.','limitation':'NMF rank4 exploratory, finite cycle exclusion not complete; prior dataset reuse not new blinded validation. Scores not validated protein/state identity or cell abundance.','support_rule':'At least70% eligible test patients (>=5 detected) score-PGAM5 rho>=.1, at least3 external eligible donors, and external nonmacrophage group median never exceeds macrophage median (>=20 cells/group). Otherwise not recommended for bulk PGAM5-specific state inference.'}
(R/'program_discovery_protocol.json').write_text(json.dumps(protocol,indent=2))
for train,test,phase in [('GSE151530','GSE149614','cohort_holdout_151_to_149'),('GSE149614','GSE151530','cohort_holdout_149_to_151'),('joint_discovery','GSE189903','joint_external_retest')]:
 data=v.D[train];o=data['obs'];mac=o.harmonized_group.eq('Macrophage').to_numpy();X=data['norm'][mac].copy();X.data=np.log1p(X.data);mu=np.asarray(X.mean(0)).ravel();var=np.asarray(X.multiply(X).mean(0)).ravel()-mu*mu;rank=var/(mu+.05);ix=np.flatnonzero(eligible&(mu>.01));ix=ix[np.argsort(-rank[ix],kind='stable')[:500]];x=X[:,ix].toarray().astype(np.float32);model=NMF(n_components=4,init='nndsvda',solver='mu',random_state=0,max_iter=300,tol=1e-4)
 with warnings.catch_warnings(record=True) as ww:
  warnings.simplefilter('always');W=model.fit_transform(x)
 conv=not any(issubclass(w.category,ConvergenceWarning) for w in ww);convergence.append({'phase':phase,'iterations':model.n_iter_,'converged':conv,'reconstruction_error':float(model.reconstruction_err_)})
 pg=np.log1p(data['norm'][mac][:,v.pg].toarray().ravel());rhos=[spearmanr(W[:,j],pg).statistic for j in range(4)];selected=int(np.nanargmax(rhos));top=np.argsort(-model.components_[selected],kind='stable')[:30];gi=ix[top];gx=np.log1p(data['norm'][mac][:,gi].toarray());m=gx.mean(0);sd=np.maximum(gx.std(0),.1)
 for j in range(4):
  tt=np.argsort(-model.components_[j],kind='stable')[:30]
  factor_rows.append({'phase':phase,'factor':j,'training_Spearman_PGAM5':rhos[j],'selected':j==selected,'top30_genes':';'.join(v.genes[ix[tt]])})
 for g,a,b,w in zip(v.genes[gi],m,sd,model.components_[selected,top]):genes_out.append({'phase':phase,'gene':g,'training_mean_log1p_CP10k':float(a),'training_SD_floor_0_1':float(b),'NMF_loading':float(w),'score_weight':1/30})
 # Export full learned basis so factor construction can be audited/reused.
 pd.DataFrame(model.components_.T,index=v.genes[ix],columns=['factor_'+str(j) for j in range(4)]).to_csv(R/(phase+'_NMF_basis.tsv'),sep='\t')
 for name,role in [(train,'training'),(test,'test')]:
  dd=v.D[name];ob=dd['obs'];xx=np.log1p(dd['norm'][:,gi].toarray());score=((xx-m)/sd).mean(1);pgl=np.log1p(dd['norm'][:,v.pg].toarray().ravel());t=ob.copy();t['score']=score;t['PGAM5_log_CP10k']=pgl
  cellout.append(t[['cell_id','donor','harmonized_group','PGAM5_counts','score']].assign(phase=phase,role=role,dataset=name))
  for d,a in t.groupby('donor',observed=True):
   b=a[a.harmonized_group.eq('Macrophage')];y=b.fine_group.eq(v.TARGET).to_numpy();z=b.score.to_numpy();n=y.sum();nn=(~y).sum();auc=(rankdata(z)[y].sum()-n*(n+1)/2)/(n*nn) if n and nn else np.nan;rho=spearmanr(z,b.PGAM5_log_CP10k).statistic if len(b)>1 and b.PGAM5_log_CP10k.nunique()>1 else np.nan
   outputs.append({'phase':phase,'role':role,'dataset':name,'donor':d,'macrophages':len(b),'PGAM5_detected':int(n),'eligible_positive_ge5':n>=5,'AUC_RNA_detection':auc,'Spearman_score_PGAM5':rho,'selected_training_factor_Spearman':rhos[selected],'NMF_converged':conv})
   for g,b in a.groupby('harmonized_group',observed=True):lineages.append({'phase':phase,'role':role,'dataset':name,'donor':d,'lineage':g,'cells':len(b),'median_score':b.score.median(),'p95_score':b.score.quantile(.95)})
 print('Completed NMF',phase,'training factor rhos',rhos,'top genes',v.genes[gi].tolist(),flush=True)
pd.DataFrame(outputs).to_csv(R/'program_validation_by_patient.csv',index=False);pd.DataFrame(genes_out).to_csv(R/'program_score_parameters.csv',index=False);pd.DataFrame(lineages).to_csv(R/'program_lineage_specificity.csv',index=False);pd.DataFrame(factor_rows).to_csv(R/'all_programs_training_association.csv',index=False);pd.concat(cellout,ignore_index=True).to_csv(R/'program_cell_scores.csv.gz',index=False,compression='gzip');pd.DataFrame(convergence).to_csv(R/'NMF_convergence.csv',index=False)
