"""Independent RNA labels, BVLS projections, mixture truth, summaries and KKT audit."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import lsq_linear

R = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--saved-only',action='store_true')
    args = parser.parse_args()
    p = json.loads((R/'protocol.json').read_text())
    source = Path(p['source_root'])
    checks = []
    def check(name,value):
        assert value,name
        checks.append(name)
    assoc = pd.read_csv(R/'program_PGAM5_association_by_cohort.csv')
    selected = pd.read_csv(R/'candidate_selection_by_fold.csv')
    diagnostics = pd.read_csv(R/'NMF_fit_diagnostics.csv')
    check('All NMF fits converged',diagnostics.iterations.lt(2000).all()&diagnostics.warnings.eq('[]').all())
    cycle = set((R/'cell_cycle_genes.txt').read_text().split())
    for held in p['primary_cohorts']+['ALL_TRAINING']:
        folder = R/held
        definition = json.loads((folder/'definition.json').read_text())
        features = pd.read_csv(folder/'features.csv')
        h = pd.read_csv(folder/'program_loadings.tsv',sep='\t',index_col=0)
        check(held+' exact feature rows',np.array_equal(features.gene,h.index))
        check(held+' PGAM5 and listed cycle genes excluded', 'PGAM5' not in set(features.gene) and not (cycle&set(features.gene)))
        check(held+' nonnegative normalized programs',(h.to_numpy()>=0).all() and np.allclose(np.linalg.norm(h.to_numpy(),axis=0),1))
        check(held+' held cohort excluded from training',held=='ALL_TRAINING' or held not in definition['training_cohorts'])
        a = assoc[assoc.fold.eq(held)&assoc.role.eq('training')]
        support = a.groupby('program').cohort_support_passed.sum()
        ranked = selected[selected.fold.eq(held)].sort_values(['supported_training_cohorts','median_cohort_donor_ratio','program'],ascending=[False,False,True])
        check(held+' candidate training-only ranked selection',int(ranked.iloc[0].program)==definition['target_program'])
        check(held+' training support recomputed',int(support.loc[definition['target_program']])==int(ranked.iloc[0].supported_training_cohorts))
        check(held+' eligibility exact training support rule',definition['training_selection_eligible']==(support.loc[definition['target_program']]>=np.ceil(.6*len(definition['training_cohorts']))))
        if not args.saved_only:
            cohort = held if held!='ALL_TRAINING' else p['primary_cohorts'][0]
            exported = pd.read_csv(folder/'heldout_macrophage_annotations.csv.gz') if held!='ALL_TRAINING' else pd.read_csv(R/f'{cohort}_macrophage_annotations.csv.gz')
            original = pd.read_csv(source/f'{cohort}_metadata.csv.gz')
            raw = sparse.load_npz(source/f'{cohort}_raw_common.npz')
            genes = pd.Index(pd.read_csv(source/'common_genes.csv').gene)
            rows = pd.Index(original.cell_id).get_indexer(exported.cell_id)
            sample = np.random.default_rng(20261010).choice(len(rows),min(25,len(rows)),replace=False)
            # Match the actual float32 full-library preprocessing; solve with BVLS,
            # independently of the discovery/project NNLS implementation.
            counts = raw[rows[sample]][:,genes.get_indexer(features.gene)].astype(np.float32)
            total = original.iloc[rows[sample]].total_counts.to_numpy(dtype=np.float32)
            x = counts.multiply((10000/total)[:,None]).toarray()
            x = np.minimum(np.log1p(x)/features.SD.to_numpy(),10)
            for j,row in enumerate(x):
                fit = lsq_linear(h.to_numpy(),row,bounds=(0,np.inf),method='bvls',tol=1e-12,max_iter=1000)
                w = fit.x
                order = np.argsort(w)
                winner = order[-1]
                margin = (w[winner]-w[order[-2]])/max(w[winner],1e-12)
                target = definition['target_program']
                member = winner==target and w[target]>=definition['activation_thresholds'][target] and margin>=.05
                observed = exported.iloc[sample[j]]
                check(held+' independent program projection '+str(j),winner==observed.dominant_program and
                      np.isclose(w[target],observed.candidate_program_score,rtol=1e-5,atol=1e-6) and member==observed.candidate_state_member)
    totalmac,totalpositive = 0,0
    sensitivity = pd.read_csv(R/'raw_count_threshold_sensitivity.csv')
    for cohort in p['primary_cohorts']:
        annotation = pd.read_csv(R/f'{cohort}_all_cell_annotation.csv.gz')
        mac = annotation.cell_type.eq('Macrophage')
        positive = annotation.PGAM5_counts.gt(0)
        check(cohort+' exact RNA-status labels',np.array_equal(annotation.PGAM5_RNA_status.eq('RNA_detected'),positive))
        check(cohort+' exact PGAM5 macrophage labels',np.array_equal(annotation.PGAM5_macrophage_annotation.eq('Macrophage_PGAM5_RNA_detected'),mac&positive))
        check(cohort+' cell IDs unique and biological validation flag false',annotation.cell_id.is_unique and not annotation.usable_for_validated_TCGA_PGAM5_cell_abundance.any())
        totalmac += int(mac.sum())
        totalpositive += int((mac&positive).sum())
        for cutoff in [1,2,3,5]:
            row = sensitivity[sensitivity.dataset.eq(cohort)&sensitivity.raw_count_cutoff.eq(cutoff)].iloc[0]
            check(cohort+' threshold '+str(cutoff)+' independently counted',int((mac&annotation.PGAM5_counts.ge(cutoff)).sum())==row.cells_above_cutoff)
        if not args.saved_only:
            original = pd.read_csv(source/f'{cohort}_metadata.csv.gz').set_index('cell_id')
            raw = sparse.load_npz(source/f'{cohort}_raw_common.npz')
            genes = pd.Index(pd.read_csv(source/'common_genes.csv').gene)
            raw_pg = raw[:,genes.get_loc('PGAM5')].toarray().ravel()
            check(cohort+' all raw PGAM5 source counts',np.array_equal(raw_pg,original.PGAM5_counts))
            aligned = annotation.set_index('cell_id').loc[original.index]
            check(cohort+' all annotation counts and lineage from source',np.array_equal(raw_pg,aligned.PGAM5_counts) and original.harmonized_group.eq(aligned.cell_type).all())
    check('primary30533 macrophages and1413 RNA detections',totalmac==30533 and totalpositive==1413)
    truth = pd.read_csv(R/'reconstructed_mixture_truths.csv').set_index('mixture_id')
    for path in sorted(R.glob('membership_audit_*.npz')):
        a = np.load(path)
        pg = a['macrophage'] & (a['PGAM5_raw']>0)
        intersection = pg&a['candidate_state']
        library = a['full_library_counts']
        m = int(a['equalized_mixture_id'])
        for offset,flag in [(0,None),(1,None)]:
            row = truth.loc[m+offset]
            weight = np.ones(len(library)) if offset==0 else library
            check(path.name+' unit '+str(offset)+' exact target memberships',
                  np.isclose(np.sum(weight[pg])/weight.sum(),row.truth_all_PGAM5_RNA) and
                  np.isclose(np.sum(weight[intersection])/weight.sum(),row.truth_PGAM5_RNA_in_candidate_state))
    pred = pd.read_csv(R/'heldout_mixture_predictions.csv.gz')
    summary = pd.read_csv(R/'reference_validation_summary.csv')
    check('3320 heldout mixtures,3 distinct estimands',len(pred)==9960 and not pred.duplicated(['mixture_id','estimand']).any())
    for row in summary.itertuples():
        d = pred[pred.dataset.eq(row.dataset)&pred.unit.eq(row.unit)&pred.estimand.eq(row.estimand)]
        s = d[d.scenario.eq('standard')]
        z = d[d.truth.eq(0)]
        x,y = s.truth.to_numpy(),s.predicted.to_numpy()
        r = np.dot(x-x.mean(),y-y.mean())/(np.linalg.norm(x-x.mean())*np.linalg.norm(y-y.mean())) if np.std(x)>0 and np.std(y)>0 else np.nan
        check(row.dataset+row.unit+row.estimand+' independent metrics',
              np.isclose(np.mean(abs(x-y))*100,row.standard_MAE_pp) and
              np.isclose(r,row.standard_r,equal_nan=True) and
              np.isclose(np.quantile(z.predicted,.95)*100,row.all_zero_P95_percent,equal_nan=True))
    for path in sorted(R.glob('fit_audit_case_*.npz')):
        a = np.load(path)
        A,b,c,w = a['A'],a['b'],a['coefficients'],a['weights']
        grad = A.T@(w*(A@c-b))/len(b)+.001*c
        active = c>1e-6
        multiplier = grad[active].mean()
        violation = max(np.max(abs(grad[active]-multiplier)),
                        np.max(np.maximum(multiplier-grad[~active],0)) if (~active).any() else 0.)
        check(path.name+' fixed-weight convex KKT certificate',abs(c.sum()-1)<1e-7 and c.min()>=-1e-9 and violation<2e-5)
    gates = pd.read_csv(R/'validation_gates.csv')
    result = json.loads((R/'validation_result.json').read_text())
    check('failed validation retained',int((~gates.passed).sum())==result['failed_gates'] and not result['usable_as_TCGA_cellular_abundance'])
    out = {'status':'PASS','source_data_audited':not args.saved_only,'checks_passed':len(checks),'checks':checks,
           'meaning':'Implementation/source/split/last convex solver checks; not proof of stable PGAM5-positive subtype or cellular abundance'}
    (R/('saved_audit.json' if args.saved_only else 'source_audit.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='checks'},indent=2))


if __name__=='__main__':
    main()
