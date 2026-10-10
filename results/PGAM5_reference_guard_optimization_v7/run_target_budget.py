"""Documented post-diagnostic repair: constrain only target PGAM5 contribution."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import json
import hashlib
import numpy as np
import pandas as pd
from guard_solver import fit
from optimize_reference import summarize,is_target,is_mac

R=Path(__file__).resolve().parent
P=json.loads((R/'protocol.json').read_text());V=Path(P['v6_root']);S=Path(P['source_root'])

def worker(task):
    held,unit=task
    m=np.load(R/held/f'{unit}_model.npz');A=m['A'];groups=m['groups']
    base=int(m['baseline_columns']);pg=int(np.flatnonzero(m['genes']=='PGAM5')[0])
    common=pd.Index(pd.read_csv(S/'common_genes.csv').gene)
    x=np.load(S/'all_mixture_expression.npz')['values'];ix=common.get_indexer(m['genes'])
    truth=pd.read_csv(V/'reconstructed_mixture_truths.csv')
    truth=truth[truth.dataset.eq(held)&truth.unit.eq(unit)]
    targets=np.array([is_target(g) for g in groups]);mac=np.array([is_mac(g) for g in groups])
    rows,cases=[],[]
    for j,row in enumerate(truth.itertuples()):
        b=x[int(row.mixture_id),ix].astype(float)/m['scale']
        for name,size in [('target_only_budget',base),('prototypes_target_only_budget',len(groups))]:
            mask=targets[:size]
            c,w,res=fit(A[:,:size],b,True,pg,mask)
            rows.append({'dataset':held,'unit':unit,'method':name,'mixture_id':row.mixture_id,'donor':row.donor,
                'scenario':row.scenario,'replicate':row.replicate,'nominal_target_cell_fraction':row.nominal_target_cell_fraction,
                'truth':row.truth_all_PGAM5_RNA,'predicted':float(c[mask].sum()),
                'truth_macrophage':row.macrophage_cell_fraction if unit=='equalized_cell_fraction' else row.macrophage_RNA_contribution,
                'predicted_macrophage':float(c[mac[:size]].sum()),'relative_residual':res,
                'observed_PGAM5_scaled':b[pg],'predicted_PGAM5_scaled':float(A[pg,:size]@c),
                'target_PGAM5_scaled':float((A[pg,:size]*mask)@c),'components':size,
                'usable_as_cellular_abundance':False})
            if row.scenario=='standard' and row.replicate==0 and row.nominal_target_cell_fraction in [0,.05] and len(cases)<4:
                key=f'targetOnly_{held}_{unit}_{len(cases)}'
                model_path=R/held/f'{unit}_model.npz'
                np.savez_compressed(R/f'fit_case_{key}.npz',b=b,c=c,w=w,budget=True,
                    pg_index=pg,target_mask=mask,macrophage_mask=mac[:size],budget_vector=A[pg,:size]*mask,
                    model_relative_path=np.array(str(model_path.relative_to(R))),columns=size,
                    model_SHA256=np.array(hashlib.sha256(model_path.read_bytes()).hexdigest()))
                cases.append({'case':key,'dataset':held,'unit':unit,'method':name,'mixture_id':row.mixture_id})
        if (j+1)%150==0:print('Target-only',held,unit,j+1,'of',len(truth),flush=True)
    print('Completed target-only',held,unit,flush=True)
    return pd.DataFrame(rows),pd.DataFrame(cases)

def main():
    provenance={'amendment_SHA256':hashlib.sha256((R/'diagnostic_amendment.json').read_bytes()).hexdigest(),
                'initial_prediction_SHA256':hashlib.sha256((R/'heldout_predictions.csv.gz').read_bytes()).hexdigest(),
                'timing':'Captured before new target-only runs'}
    (R/'amendment_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    tasks=[(c,u) for c in P['primary_cohorts'] for u in ['equalized_cell_fraction','library_RNA_contribution']]
    results=[]
    with ProcessPoolExecutor(max_workers=3) as pool:
        for f in as_completed([pool.submit(worker,t) for t in tasks]):results.append(f.result())
    d=pd.concat([x[0] for x in results],ignore_index=True).sort_values(['dataset','unit','method','mixture_id'])
    d.to_csv(R/'target_only_predictions.csv.gz',index=False)
    pd.concat([x[1] for x in results],ignore_index=True).to_csv(R/'target_only_fit_cases.csv',index=False)
    combined=pd.concat([pd.read_csv(R/'heldout_predictions.csv.gz'),d],ignore_index=True)
    combined.to_csv(R/'all_predictions.csv.gz',index=False)
    summary=summarize(combined);summary.to_csv(R/'all_validation_summary.csv',index=False)
    validity=json.loads((R/'validation_result.json').read_text())
    validity['diagnostic_amendment_methods']=['target_only_budget','prototypes_target_only_budget']
    validity['amendment_is_fresh_independent_validation']=False
    validity['all_component_budget_recommended']=False
    validity['usable_as_TCGA_cellular_abundance']=False
    (R/'validation_result.json').write_text(json.dumps(validity,indent=2)+'\n')
    print(summary[summary.method.str.contains('target_only')].to_string(index=False),flush=True)

if __name__=='__main__':main()
