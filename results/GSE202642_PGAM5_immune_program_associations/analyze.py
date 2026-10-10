from pathlib import Path
import json,numpy as np,pandas as pd
from scipy.stats import rankdata,spearmanr
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parent; B=9999

def rho_batch(x,y):
 a=rankdata(x,axis=-1);b=rankdata(y,axis=-1);a-=a.mean(axis=-1,keepdims=True);b-=b.mean(axis=-1,keepdims=True)
 return (a*b).sum(axis=-1)/np.sqrt((a*a).sum(axis=-1)*(b*b).sum(axis=-1))

def main():
 a=pd.read_csv(R/'cell_metadata.csv.gz',index_col=0);genes=pd.read_csv(R/'program_genes.csv').gene.tolist();raw=np.load(R/'program_raw_counts.npz')['counts'];z=np.log1p(raw/a.total_counts.to_numpy()[:,None]*10000);x=z[:,genes.index('PGAM5')];panels=json.loads((R/'program_definitions.json').read_text());scores=pd.DataFrame({p:z[:,[genes.index(g) for g in gs]].mean(axis=1) for p,gs in panels.items()},index=a.index);scores.to_csv(R/'cell_program_scores.csv.gz');rows=[];patients=[]
 for c,idx in a.groupby('primary_cluster').indices.items():
  xx=x[idx];n=len(idx);pos=int((xx>0).sum());rng=np.random.default_rng(202642+int(c));nulls={p:[] for p in panels}
  if pos:
   if pos==1:
    perms=np.zeros((n,n));value=xx[xx>0][0];perms[np.arange(n),np.arange(n)]=value;btotal=n
   else:perms=None;btotal=B
   for start in range(0,btotal,128):
    xp=perms[start:start+128] if perms is not None else np.array([rng.permutation(xx) for _ in range(min(128,btotal-start))])
    for p,gs in panels.items():
     other=z[idx][:,[genes.index(g) for g in gs if g!='PGAM5']].sum(axis=1);yp=(xp+other)/len(gs);nulls[p].extend(rho_batch(xp,yp).tolist())
  for p,gs in panels.items():
   yy=scores[p].to_numpy()[idx];r={'cluster':int(c),'functional_name':a.iloc[idx[0]].functional_name_EN,'program':p,'cells':n,'PGAM5_detected':pos,'mixed_RNA_flag':bool(a.iloc[idx[0]].mixed_RNA_flag),'low_information':pos<5,'rho':np.nan,'null_median':np.nan,'delta_rho':np.nan,'null_2.5pct':np.nan,'null_97.5pct':np.nan,'p_calibrated':np.nan,'permutations':0}
   if pos:
    obs=float(spearmanr(xx,yy).statistic);nu=np.array(nulls[p]);med=np.median(nu);pval=min(1,2*min((int((nu>=obs-1e-12).sum())+(pos!=1))/(len(nu)+(pos!=1)),(int((nu<=obs+1e-12).sum())+(pos!=1))/(len(nu)+(pos!=1))));r.update(rho=obs,null_median=med,delta_rho=obs-med,p_calibrated=pval,permutations=len(nu));r['null_2.5pct'],r['null_97.5pct']=np.quantile(nu,[.025,.975])
   rows.append(r)
   for patient,pidx in a.iloc[idx].groupby('patient_id').indices.items():
    ii=idx[pidx];px=x[ii];py=scores[p].to_numpy()[ii];pp=int((px>0).sum());eligible=len(ii)>=20 and pp>=3 and (px==0).sum()>=3
    # Full scores retain PGAM5; subtract only its known covariance contribution for a descriptive direction diagnostic.
    excess=float(np.mean((px-px.mean())*(py-py.mean()))-np.var(px)/len(gs)) if eligible else np.nan
    patients.append(dict(cluster=int(c),program=p,patient_id=patient,cells=len(ii),PGAM5_detected=pp,eligible=eligible,excess_covariance=excess))
  print('Finished C'+str(c),flush=True)
 out=pd.DataFrame(rows);valid=out.p_calibrated.notna();out['q_BH']=np.nan;out.loc[valid,'q_BH']=multipletests(out.loc[valid,'p_calibrated'],method='fdr_bh')[1];out.to_csv(R/'subgroup_program_associations.csv',index=False);pd.DataFrame(patients).to_csv(R/'patient_direction_diagnostics.csv',index=False)
 # Self-only program sanity: correlation is identically one under observed and permuted scores.
 t=np.array([0.,0.,1.,2.]);assert np.allclose(rho_batch(np.array([t,t[::-1]]),np.array([t,t[::-1]])/4),1)
 manifest=[dict(program=p,gene=g,role='shared_PGAM5_anchor' if g=='PGAM5' else 'immune_RNA_member',detected_cells=int((raw[:,genes.index(g)]>0).sum())) for p,gs in panels.items() for g in gs];pd.DataFrame(manifest).to_csv(R/'program_gene_manifest.csv',index=False)
 print(out.loc[out.q_BH<.05,['cluster','program','rho','delta_rho','q_BH']].to_string(index=False),flush=True)
if __name__=='__main__':main()
