"""Inspect error reductions together with dose, lineage and residual costs."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import nnls

R=Path(__file__).resolve().parent

def main():
    d=pd.read_csv(R/'all_predictions.csv.gz')
    rows=[]
    for key,x in d.groupby(['dataset','unit','method','scenario']):
        zero=x[x.truth_macrophage.eq(0)]
        rows.append(dict(zip(['dataset','unit','method','scenario'],key))|{
            'mixtures':len(x),'macrophage_MAE_pp':float(abs(x.truth_macrophage-x.predicted_macrophage).mean()*100),
            'median_predicted_macrophage_percent':float(x.predicted_macrophage.median()*100),
            'zero_macrophage_P95_percent':float(zero.predicted_macrophage.quantile(.95)*100) if len(zero) else np.nan,
            'median_relative_residual':float(x.relative_residual.median())})
    pd.DataFrame(rows).to_csv(R/'lineage_residual_tradeoffs.csv',index=False)
    geometry=[]
    for fold in json.loads((R/'protocol.json').read_text())['primary_cohorts']+['ALL_TRAINING']:
        for unit in ['equalized_cell_fraction','library_RNA_contribution']:
            m=np.load(R/fold/f'{unit}_model.npz')
            targets=np.array([not g.startswith('NuisancePrototype::') and 'PGAM5_detected_' in g for g in m['groups']])
            for size,label in [(int(m['baseline_columns']),'baseline'),(len(targets),'augmented')]:
                for pg_mode in ['included','excluded']:
                    genes=np.ones(len(m['genes']),bool) if pg_mode=='included' else m['genes']!='PGAM5'
                    A=m['A'][genes,:size];flag=targets[:size]
                    for j in np.flatnonzero(flag):
                        c,res=nnls(A[:,~flag],A[:,j],maxiter=6000)
                        geometry.append({'fold':fold,'unit':unit,'reference':label,'PGAM5_feature':pg_mode,
                            'target_component':m['groups'][j],'relative_cone_residual':float(res/max(np.linalg.norm(A[:,j]),1e-12)),
                            'meaning':'Reference-column geometry; not an independent classification or simplex abundance recovery'})
    pd.DataFrame(geometry).to_csv(R/'geometry_PGAM5_sensitivity.csv',index=False)
    summary=pd.read_csv(R/'all_validation_summary.csv')
    base=summary[summary.method.eq('v6_baseline')].set_index(['dataset','unit'])
    changes=[]
    for row in summary.itertuples():
        before=base.loc[row.dataset,row.unit]
        changes.append({'dataset':row.dataset,'unit':row.unit,'method':row.method,
            'zero_P95_change_pp':row.all_zero_P95_percent-before.all_zero_P95_percent,
            'standard_MAE_change_pp':row.standard_MAE_pp-before.standard_MAE_pp,
            'standard_r_change':row.standard_r-before.standard_r,
            'standard_macrophage_MAE_change_pp':row.standard_macrophage_MAE_pp-before.standard_macrophage_MAE_pp})
    pd.DataFrame(changes).to_csv(R/'changes_vs_v6.csv',index=False)
    doses=[]
    for key,x in d[d.scenario.eq('standard')].groupby(['dataset','unit','method','nominal_target_cell_fraction']):
        doses.append(dict(zip(['dataset','unit','method','nominal_target_cell_fraction'],key))|{'mixtures':len(x),
            'mean_truth_percent':float(x.truth.mean()*100),'median_prediction_percent':float(x.predicted.median()*100),
            'fraction_prediction_over_0p5_percent':float(x.predicted.gt(.005).mean()),
            'MAE_pp':float(abs(x.truth-x.predicted).mean()*100)})
    pd.DataFrame(doses).to_csv(R/'all_dose_recovery.csv',index=False)
    gates=[]
    for row in summary.itertuples():
        for name,val,passed in [('MAE_pp<=1',row.standard_MAE_pp,row.standard_MAE_pp<=1),
            ('r>=0.7',row.standard_r,row.standard_r>=.7),('zero_P95_percent<=0.5',row.all_zero_P95_percent,row.all_zero_P95_percent<=.5),
            ('independent_donors>=3',row.standard_donor_labels,row.standard_donor_labels>=3 and row.dataset!='GSE202642')]:
            gates.append({'dataset':row.dataset,'unit':row.unit,'method':row.method,'criterion':name,'value':val,'passed':passed})
    pd.DataFrame(gates).to_csv(R/'all_validation_gates.csv',index=False)
    print(summary.groupby('method')[['standard_MAE_pp','standard_r','all_zero_P95_percent','standard_macrophage_MAE_pp']].mean().to_string())

if __name__=='__main__':main()
