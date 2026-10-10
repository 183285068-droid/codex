"""Fixed factorial optimization, held-cohort recovery and nuisance diagnostics."""
from pathlib import Path
import json
import hashlib
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from scipy.optimize import nnls
from guard_solver import fit

R = Path(__file__).resolve().parent
P = json.loads((R/'protocol.json').read_text())
V = Path(P['v6_root']);S = Path(P['source_root'])
COHORTS = P['primary_cohorts']

def prototype_profiles(training,unit,genes,scale):
    common = pd.Index(pd.read_csv(S/'common_genes.csv').gene)
    rows,values = [],[]
    for cohort in training:
        table = pd.read_csv(S/f'{cohort}_group_stats_index.csv')
        stats = np.load(S/f'{cohort}_group_statistics.npz')['values']
        for j,row in table.iterrows():
            if row.cells<20 or 'PGAM5_detected' in row.group: continue
            denominator = row.cells if unit=='equalized_cell_fraction' else row.total_counts/10000
            values.append(stats[j,0 if unit=='equalized_cell_fraction' else 1,common.get_indexer(genes)].astype(float)/denominator)
            rows.append(row.to_dict())
    table = pd.DataFrame(rows);values = np.asarray(values)
    columns,names,roles = [],[],[]
    for group,d in table.groupby('group',sort=True):
        if len(d)<6: continue
        weights = np.array([1/(d.dataset.nunique()*sum(d.dataset.eq(c))) for c in d.dataset])
        x = values[d.index]
        model = KMeans(n_clusters=min(3,len(d)//3),random_state=P['seed'],n_init=10)
        labels = model.fit_predict(np.log1p(x/scale),sample_weight=weights)
        for k in range(model.n_clusters):
            flag = labels==k
            cluster = d.iloc[np.flatnonzero(flag)]
            if cluster.donor.nunique()<2 or cluster.cells.sum()<50: continue
            name = 'NuisancePrototype::'+group+'::'+str(k)
            columns.append(np.average(x[flag],axis=0,weights=weights[flag]))
            names.append(name)
            roles.append({'prototype':name,'original_fine_group':group,'unit':unit,'donor_labels':cluster.donor.nunique(),
                          'cohorts':cluster.dataset.nunique(),'cells':int(cluster.cells.sum()),
                          'donors':','.join(cluster.donor),'training_cohorts':','.join(training)})
    return np.column_stack(columns),names,roles

def is_target(name):
    return not name.startswith('NuisancePrototype::') and 'PGAM5_detected_' in name

def is_mac(name):
    group = name.split('::')[1] if name.startswith('NuisancePrototype::') else name
    return group.startswith(('TAM_','Cycling_TAM_'))

def build(held,unit):
    folder = R/held;folder.mkdir(exist_ok=True)
    ref = pd.read_csv(V/held/f'EXPLORATORY_reference_{unit}.tsv',sep='\t',index_col=0)
    scales = pd.read_csv(V/held/f'{unit}_row_scales.csv').set_index('gene').loc[ref.index,'scale']
    definition = json.loads((V/held/'definition.json').read_text())
    proto,names,roles = prototype_profiles(definition['training_cohorts'],unit,ref.index,scales.to_numpy())
    augmented = np.column_stack([ref.to_numpy(),proto])
    groups = ref.columns.tolist()+names
    pd.DataFrame(roles).to_csv(folder/f'{unit}_prototype_membership.csv',index=False)
    pd.DataFrame(augmented,index=ref.index,columns=groups).to_csv(folder/f'EXPLORATORY_augmented_reference_{unit}.tsv',sep='\t')
    scales.rename('scale').to_csv(folder/f'{unit}_row_scales.csv')
    np.savez_compressed(folder/f'{unit}_model.npz',A=augmented/scales.to_numpy()[:,None],
        scale=scales.to_numpy(),genes=ref.index.to_numpy(dtype=str),groups=np.array(groups,dtype=str),baseline_columns=len(ref.columns))
    return {'fold':held,'unit':unit,'baseline_components':len(ref.columns),'added_prototypes':len(names),
            'genes':len(ref),'training_cohorts':','.join(definition['training_cohorts'])}

def run_hold(task):
    held,unit = task
    data = np.load(R/held/f'{unit}_model.npz')
    A,groups = data['A'],data['groups']
    base = int(data['baseline_columns'])
    common = pd.Index(pd.read_csv(S/'common_genes.csv').gene)
    measured = np.load(S/'all_mixture_expression.npz')['values']
    index = common.get_indexer(data['genes'])
    pg = int(np.flatnonzero(data['genes']=='PGAM5')[0])
    truth = pd.read_csv(V/'reconstructed_mixture_truths.csv')
    truth = truth[truth.dataset.eq(held)&truth.unit.eq(unit)]
    previous = pd.read_csv(V/'heldout_mixture_predictions.csv.gz')
    previous = previous[previous.estimand.eq('all_PGAM5_RNA')&previous.dataset.eq(held)&previous.unit.eq(unit)].set_index('mixture_id')
    rows,cases,geometry = [],[],[]
    targets = np.array([is_target(g) for g in groups]);mac = np.array([is_mac(g) for g in groups])
    for size,label in [(base,'baseline'),(len(groups),'augmented')]:
        B = A[:,:size];mask = targets[:size]
        for k in np.flatnonzero(mask):
            rival = B[:,~mask]
            c,res = nnls(rival,B[:,k],maxiter=6000)
            geometry.append({'dataset':held,'unit':unit,'reference':label,'target_component':groups[k],
                'relative_cone_residual':float(res/max(np.linalg.norm(B[:,k]),1e-12)),
                'top_competitor':str(groups[:size][~mask][np.argmax(c)]),
                'meaning':'Geometry diagnostic only; cone coefficients do not have to sum to1'})
    for j,row in enumerate(truth.itertuples()):
        b = measured[int(row.mixture_id),index].astype(float)/data['scale']
        for name,size,budget in [('v6_baseline',base,False),('PGAM5_budget',base,True),
                                 ('background_prototypes',len(groups),False),('prototypes_plus_budget',len(groups),True)]:
            c,w,residual = fit(A[:,:size],b,budget,pg)
            target = float(c[targets[:size]].sum())
            if name=='v6_baseline':
                assert abs(target-previous.loc[row.mixture_id,'predicted'])<1e-6
            truthmac = row.macrophage_cell_fraction if unit=='equalized_cell_fraction' else row.macrophage_RNA_contribution
            rows.append({'dataset':held,'unit':unit,'method':name,'mixture_id':row.mixture_id,
                'donor':row.donor,'scenario':row.scenario,'replicate':row.replicate,
                'nominal_target_cell_fraction':row.nominal_target_cell_fraction,'truth':row.truth_all_PGAM5_RNA,
                'predicted':target,'truth_macrophage':truthmac,'predicted_macrophage':float(c[mac[:size]].sum()),
                'relative_residual':residual,'observed_PGAM5_scaled':b[pg],'predicted_PGAM5_scaled':float(A[pg,:size]@c),
                'components':size,'usable_as_cellular_abundance':False})
            if row.scenario=='standard' and row.replicate==0 and row.nominal_target_cell_fraction in [0,.05] and len(cases)<8:
                key=f'{held}_{unit}_{len(cases)}'
                model_path=R/held/f'{unit}_model.npz'
                np.savez_compressed(R/f'fit_case_{key}.npz',b=b,c=c,w=w,
                    budget=budget,pg_index=pg,target_mask=targets[:size],macrophage_mask=mac[:size],
                    model_relative_path=np.array(str(model_path.relative_to(R))),columns=size,
                    model_SHA256=np.array(hashlib.sha256(model_path.read_bytes()).hexdigest()))
                cases.append({'case':key,'dataset':held,'unit':unit,'method':name,'mixture_id':row.mixture_id})
        if (j+1)%100==0: print('Progress',held,unit,j+1,'of',len(truth),flush=True)
    print('Completed',held,unit,len(truth),'mixture rows x4',flush=True)
    return pd.DataFrame(rows),pd.DataFrame(cases),pd.DataFrame(geometry)

def correlation(x,y):
    return float(np.corrcoef(x,y)[0,1]) if np.std(x)>0 and np.std(y)>0 else np.nan

def summarize(d):
    rows=[]
    for (cohort,unit,method),x in d.groupby(['dataset','unit','method']):
        s=x[x.scenario.eq('standard')];z=x[x.truth.eq(0)]
        rows.append({'dataset':cohort,'unit':unit,'method':method,'standard_donor_labels':s.donor.nunique(),
            'standard_MAE_pp':float(abs(s.truth-s.predicted).mean()*100),'standard_r':correlation(s.truth,s.predicted),
            'all_zero_P95_percent':float(z.predicted.quantile(.95)*100),'all_zero_max_percent':float(z.predicted.max()*100),
            'standard_macrophage_MAE_pp':float(abs(s.truth_macrophage-s.predicted_macrophage).mean()*100),
            'median_relative_residual':float(x.relative_residual.median())})
    return pd.DataFrame(rows)

def main():
    provenance={'protocol_SHA256':hashlib.sha256((R/'protocol.json').read_bytes()).hexdigest(),
                'inputs':{},'v6_cell_annotation_and_state_validation_unchanged':True}
    for cohort in COHORTS:
        for file in [S/f'{cohort}_group_statistics.npz',S/f'{cohort}_group_stats_index.csv']:
            provenance['inputs'][str(file)]={'bytes':file.stat().st_size,'SHA256':hashlib.sha256(file.read_bytes()).hexdigest()}
    for file in [S/'all_mixture_expression.npz',V/'heldout_mixture_predictions.csv.gz',V/'reconstructed_mixture_truths.csv',V/'validation_result.json']:
        provenance['inputs'][str(file)]={'bytes':file.stat().st_size,'SHA256':hashlib.sha256(file.read_bytes()).hexdigest()}
    (R/'source_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    models=[]
    for held in COHORTS+['ALL_TRAINING']:
        for unit in ['equalized_cell_fraction','library_RNA_contribution']:
            models.append(build(held,unit));print('Built',models[-1],flush=True)
    pd.DataFrame(models).to_csv(R/'model_dimensions.csv',index=False)
    tasks=[(c,u) for c in COHORTS for u in ['equalized_cell_fraction','library_RNA_contribution']]
    results=[]
    with ProcessPoolExecutor(max_workers=3) as pool:
        futures=[pool.submit(run_hold,t) for t in tasks]
        for future in as_completed(futures): results.append(future.result())
    d=pd.concat([x[0] for x in results],ignore_index=True).sort_values(['dataset','unit','method','mixture_id'])
    d.to_csv(R/'heldout_predictions.csv.gz',index=False)
    pd.concat([x[1] for x in results],ignore_index=True).to_csv(R/'fit_cases.csv',index=False)
    pd.concat([x[2] for x in results],ignore_index=True).to_csv(R/'reference_geometry.csv',index=False)
    summary=summarize(d);summary.to_csv(R/'validation_summary.csv',index=False)
    scenarios=[]
    doses=[]
    for key,x in d.groupby(['dataset','unit','method','scenario']):
        z=x[x.truth.eq(0)]
        scenarios.append(dict(zip(['dataset','unit','method','scenario'],key))|{'mixtures':len(x),
            'MAE_pp':float(abs(x.truth-x.predicted).mean()*100),'zero_P95_percent':float(z.predicted.quantile(.95)*100) if len(z) else np.nan})
    for key,x in d[d.scenario.eq('standard')].groupby(['dataset','unit','method','nominal_target_cell_fraction']):
        doses.append(dict(zip(['dataset','unit','method','nominal_target_cell_fraction'],key))|{'mixtures':len(x),
            'mean_truth_percent':float(x.truth.mean()*100),'median_prediction_percent':float(x.predicted.median()*100),
            'fraction_prediction_over_0p5_percent':float(x.predicted.gt(.005).mean()),
            'MAE_pp':float(abs(x.truth-x.predicted).mean()*100)})
    pd.DataFrame(scenarios).to_csv(R/'scenario_diagnostics.csv',index=False)
    pd.DataFrame(doses).to_csv(R/'dose_recovery.csv',index=False)
    gates=[]
    for row in summary.itertuples():
        for name,val,passed in [('MAE_pp<=1',row.standard_MAE_pp,row.standard_MAE_pp<=1),
            ('r>=0.7',row.standard_r,row.standard_r>=.7),('zero_P95_percent<=0.5',row.all_zero_P95_percent,row.all_zero_P95_percent<=.5),
            ('independent_donors>=3',row.standard_donor_labels,row.standard_donor_labels>=3 and row.dataset!='GSE202642')]:
            gates.append({'dataset':row.dataset,'unit':row.unit,'method':row.method,'criterion':name,'value':val,'passed':passed})
    pd.DataFrame(gates).to_csv(R/'validation_gates.csv',index=False)
    result={'v6_biological_state_evidence_passed':False,'independent_protein_validation_available':False,
            'real_known_mixture_validation_available':False,'bulk_platform_cell_RNA_calibrated':False,
            'usable_as_TCGA_cellular_abundance':False,'fixed_methods_reported':P['comparison'],
            'interpretation':'Computational guard/reference comparison only; do not convert coefficients into validated PGAM5-positive cell abundance'}
    (R/'validation_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(summary.to_string(index=False),flush=True)

if __name__=='__main__':
    main()
