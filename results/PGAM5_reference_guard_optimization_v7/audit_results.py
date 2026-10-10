"""Independent truth/summary/prototype and final convex KKT certificates."""
from pathlib import Path
import argparse
import json
import hashlib
import numpy as np
import pandas as pd
from scipy.optimize import linprog

R=Path(__file__).resolve().parent

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--saved-only',action='store_true');args=parser.parse_args()
    p=json.loads((R/'protocol.json').read_text());S=Path(p['source_root']);V=Path(p['v6_root'])
    d=pd.read_csv(R/'all_predictions.csv.gz');summary=pd.read_csv(R/'all_validation_summary.csv')
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    check('19920 rows,3320 queries,six fixed/documented methods',len(d)==19920 and d.mixture_id.nunique()==3320 and d.method.nunique()==6 and not d.duplicated(['mixture_id','method']).any())
    check('All interpretation flags false',not d.usable_as_cellular_abundance.any())
    truth=d[d.method.eq('v6_baseline')].set_index('mixture_id')
    for method,x in d.groupby('method'):
        check(method+' same actual truth and donor identities',np.array_equal(x.set_index('mixture_id').sort_index().truth,truth.sort_index().truth) and
              np.array_equal(x.set_index('mixture_id').sort_index().donor,truth.sort_index().donor))
        if method in ['PGAM5_budget','prototypes_plus_budget']:
            check(method+' total PGAM5 budget satisfied',(x.predicted_PGAM5_scaled<=x.observed_PGAM5_scaled+1e-7).all())
        if 'target_only' in method:
            check(method+' target-only budget satisfied',(x.target_PGAM5_scaled<=x.observed_PGAM5_scaled+1e-7).all())
        if 'budget' in method:
            check(method+' zero PGAM5 gives zero RNA-target coefficient',x.loc[x.observed_PGAM5_scaled.eq(0),'predicted'].abs().max()<1e-8)
    for row in summary.itertuples():
        x=d[d.dataset.eq(row.dataset)&d.unit.eq(row.unit)&d.method.eq(row.method)]
        s=x[x.scenario.eq('standard')];z=x[x.truth.eq(0)]
        a,b=s.truth.to_numpy(),s.predicted.to_numpy()
        corr=np.dot(a-a.mean(),b-b.mean())/(np.linalg.norm(a-a.mean())*np.linalg.norm(b-b.mean())) if np.std(a)>0 and np.std(b)>0 else np.nan
        check(row.dataset+row.unit+row.method+' independent summary',
              np.isclose(abs(a-b).mean()*100,row.standard_MAE_pp) and np.isclose(corr,row.standard_r,equal_nan=True) and
              np.isclose(np.quantile(z.predicted,.95)*100,row.all_zero_P95_percent) and
              np.isclose(abs(s.truth_macrophage-s.predicted_macrophage).mean()*100,row.standard_macrophage_MAE_pp))
    case_index=pd.concat([pd.read_csv(R/'fit_cases.csv'),pd.read_csv(R/'target_only_fit_cases.csv')]).set_index('case')
    for f in sorted(R.glob('fit_case_*.npz')):
        m=np.load(f)
        if 'A' in m:
            A=m['A']
        else:
            path=R/str(m['model_relative_path'])
            assert hashlib.sha256(path.read_bytes()).hexdigest()==str(m['model_SHA256'])
            A=np.load(path)['A'][:,:int(m['columns'])]
        b,c,w=m['b'],m['c'],m['w'];pg=int(m['pg_index']);budget=bool(m['budget'])
        vector=m['budget_vector'] if 'budget_vector' in m else A[pg]
        grad=A.T@(w*(A@c-b))/len(b)+.001*c
        free=vector<=1e-12 if budget and b[pg]<=1e-12 else np.ones(len(c),bool)
        active=c>1e-6
        # Multipliers certify nonnegative simplex + optional linear RNA budget,
        # independent of the SLSQP implementation used to obtain the coefficients.
        use_budget=budget and b[pg]>1e-12 and abs(vector@c-b[pg])<1e-6
        X=np.column_stack([np.ones(len(c)),-vector])
        inequalities=np.vstack([X[free],-X[active]])
        rhs=np.r_[grad[free]+2e-5,-grad[active]+2e-5]
        solved=linprog([0,0],A_ub=inequalities,b_ub=rhs,bounds=[(None,None),(0,None) if use_budget else (0,0)],method='highs')
        check(f.name+' fixed-weight QP certificate',solved.success and abs(c.sum()-1)<1e-7 and c.min()>=-1e-9 and
              (not budget or vector@c<=b[pg]+1e-7) and (c[~free].max(initial=0)<1e-9))
        key=f.stem[len('fit_case_'):];row=case_index.loc[key]
        actual=d[d.mixture_id.eq(row.mixture_id)&d.method.eq(row.method)].iloc[0]
        check(f.name+' actual prediction attachment',np.isclose(c[m['target_mask']].sum(),actual.predicted) and
              np.isclose(c[m['macrophage_mask']].sum(),actual.predicted_macrophage))
    for fold in p['primary_cohorts']+['ALL_TRAINING']:
        for unit in ['equalized_cell_fraction','library_RNA_contribution']:
            m=np.load(R/fold/f'{unit}_model.npz')
            membership=pd.read_csv(R/fold/f'{unit}_prototype_membership.csv')
            check(fold+unit+' prototype support and no RNA-target members',
                  membership.donor_labels.ge(2).all() and membership.cells.ge(50).all() and
                  not membership.original_fine_group.str.contains('PGAM5_detected').any() and
                  (fold=='ALL_TRAINING' or not membership.donors.str.contains(fold+':',regex=False).any()))
            check(fold+unit+' nonnegative reference and exact PGAM5 feature',
                  (m['A']>=0).all() and np.sum(m['genes']=='PGAM5')==1)
            if not args.saved_only:
                original=pd.read_csv(V/fold/f'EXPLORATORY_reference_{unit}.tsv',sep='\t',index_col=0)
                check(fold+unit+' unchanged baseline feature profiles',np.allclose(m['A'][:,:int(m['baseline_columns'])]*m['scale'][:,None],original.to_numpy()))
                allrows=[];allvalues=[]
                common=pd.Index(pd.read_csv(S/'common_genes.csv').gene);ix=common.get_indexer(m['genes'])
                training=membership.iloc[0].training_cohorts.split(',')
                for cohort in training:
                    table=pd.read_csv(S/f'{cohort}_group_stats_index.csv');stats=np.load(S/f'{cohort}_group_statistics.npz')['values']
                    for j,row in table.iterrows():
                        if row.cells<20 or 'PGAM5_detected' in row.group:continue
                        denom=row.cells if unit=='equalized_cell_fraction' else row.total_counts/10000
                        allrows.append(row.to_dict());allvalues.append(stats[j,0 if unit=='equalized_cell_fraction' else 1,ix].astype(float)/denom)
                table=pd.DataFrame(allrows);values=np.asarray(allvalues)
                for row in membership.itertuples():
                    group=table[table.group.eq(row.original_fine_group)]
                    subset=group[group.donor.isin(row.donors.split(','))]
                    weights=np.array([1/(group.dataset.nunique()*sum(group.dataset.eq(c))) for c in subset.dataset])
                    mean=np.average(values[subset.index],axis=0,weights=weights)
                    k=np.flatnonzero(m['groups']==row.prototype)[0]
                    check(fold+unit+row.prototype+' source centroid recomputed',np.allclose(mean,m['A'][:,k]*m['scale'],rtol=1e-7,atol=1e-7))
    if not args.saved_only:
        old=pd.read_csv(V/'heldout_mixture_predictions.csv.gz')
        old=old[old.estimand.eq('all_PGAM5_RNA')].set_index('mixture_id').sort_index()
        base=truth.sort_index()
        check('All3320 baseline target coefficients reproduced',np.max(abs(base.predicted.to_numpy()-old.predicted.to_numpy()))<1e-6)
        check('All3320 RNA truths unchanged',np.array_equal(base.truth.to_numpy(),old.truth.to_numpy()))
    check('Biological/TCGA validity remains false',not json.loads((R/'validation_result.json').read_text())['usable_as_TCGA_cellular_abundance'])
    out={'status':'PASS','source_audit':not args.saved_only,'checks_passed':len(checks),'checks':checks,
         'interpretation':'Saved calculation/constraint/centroid/last convex subproblem verification, not biological or cellular-abundance validation'}
    (R/('saved_audit.json' if args.saved_only else 'source_audit.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='checks'},indent=2))

if __name__=='__main__':main()
