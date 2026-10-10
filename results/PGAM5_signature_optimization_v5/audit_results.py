"""Independent source/count, rank-test, BH, split, summary and KKT certificates."""
from pathlib import Path
import argparse
import json
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata, norm, pearsonr
from scipy.optimize import lsq_linear

R = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--saved-results-only', action='store_true')
    args = p.parse_args()
    checks = []
    def check(name, condition):
        assert condition, name
        checks.append(name)
    def close(a, b, rtol=1e-7, atol=1e-10):
        return np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True)
    config = json.loads((R/'config.json').read_text())
    primary = config['primary_cohorts']
    genes = pd.read_csv(R/'common_genes.csv').gene.to_numpy()
    de = pd.read_csv(R/'primary_all_gene_DE.csv.gz')
    votes = []
    for cohort in primary:
        d = de[de.cohort.eq(cohort)].copy()
        check(cohort+' full gene universe and finite p/FDR', np.array_equal(d.gene, genes)
              and np.isfinite(d.Mann_Whitney_two_sided_p).all() and np.isfinite(d.BH_FDR).all())
        pv = d.Mann_Whitney_two_sided_p.to_numpy()
        order = np.argsort(pv)
        q = np.empty(len(pv))
        q[order] = np.minimum(1, np.minimum.accumulate((pv[order]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1])
        check(cohort+' all BH FDR independently reconstructed', close(q, d.BH_FDR, rtol=1e-10, atol=1e-14))
        fc = np.log2((d.mean_CP10k_positive.to_numpy()+.001)/(d.mean_CP10k_undetected.to_numpy()+.001))
        passes = (q<.05)&(fc>=1)&d.positive_detection.ge(.1).to_numpy()
        check(cohort+' exact linear mean log2FC and strict criteria', close(fc,d.log2FC_linear_CP10k)
              and np.array_equal(passes,d.strict_upregulated))
        zero = d.mean_CP10k_positive.eq(0)&d.mean_CP10k_undetected.eq(0)
        check(cohort+' noninformative constant-zero genes have p=1', d.loc[zero,'Mann_Whitney_two_sided_p'].eq(1).all())
        votes.append(passes)
        if not args.saved_results_only:
            obs = pd.read_csv(R/f'{cohort}_metadata.csv.gz')
            raw = sparse.load_npz(R/f'{cohort}_raw_common.npz')
            mac = obs.harmonized_group.eq('Macrophage').to_numpy()
            pos = mac & obs.PGAM5_counts.gt(0).to_numpy()
            neg = mac & ~pos
            chosen_genes = list(dict.fromkeys(['PGAM5','DHFR','MKI67','NUSAP1','TOP2A','CENPK','TYMS']
                           + d.sort_values('BH_FDR').gene.head(5).tolist()+d.loc[zero,'gene'].head(1).tolist()))
            for gene in chosen_genes:
                k = int(np.flatnonzero(genes==gene)[0])
                expression = raw[:,k].toarray().ravel()/obs.total_counts.to_numpy()*10000
                a,b = np.log1p(expression[pos]),np.log1p(expression[neg])
                n1,n0 = len(a),len(b)
                combined = np.r_[a,b]
                u1 = rankdata(combined)[:n1].sum()-n1*(n1+1)/2
                u2 = n1*n0-u1
                counts = np.unique(combined,return_counts=True)[1].astype(float)
                variance = n1*n0/12*((n1+n0+1)-np.sum(counts**3-counts)/((n1+n0)*(n1+n0-1)))
                pvalue = min(1.,2*norm.sf((max(u1,u2)-n1*n0/2-.5)/np.sqrt(variance))) if variance>0 else 1.
                expected = d.set_index('gene').loc[gene]
                check(cohort+' '+gene+' manual tie-corrected U tail and raw means',
                      close(pvalue,expected.Mann_Whitney_two_sided_p,rtol=1e-7,atol=1e-300)
                      and close([expression[pos].mean(),expression[neg].mean()],
                                [expected.mean_CP10k_positive,expected.mean_CP10k_undetected]))
            provenance = json.loads((R/'source_provenance.json').read_text())
            record = next(x for x in provenance['primary'] if x['cohort']==cohort)
            source = next(Path(path) for path in record['inputs'] if path.endswith('.h5ad'))
            with h5py.File(source,'r') as f:
                def decode(node):
                    if isinstance(node,h5py.Group):
                        if 'values' in node:
                            values = decode(node['values']).astype(object)
                            values[node['mask'][:]] = 'None'
                            return values
                        cats=decode(node['categories']);codes=node['codes'][:]
                        return np.array([cats[x] if x>=0 else 'None' for x in codes],dtype=object)
                    return node.asstr()[:] if h5py.check_string_dtype(node.dtype) is not None else node[:]
                source_genes = decode(f['var/gene'])
                pg_columns = np.flatnonzero(source_genes=='PGAM5')
                pointer = f['X/indptr'][:]
                recovered = np.zeros(len(pointer)-1,dtype=float)
                for start in range(0,len(f['X/indices']),1000000):
                    stop = min(start+1000000,len(f['X/indices']))
                    indices = f['X/indices'][start:stop]
                    at = np.flatnonzero(np.isin(indices,pg_columns))
                    if len(at):
                        rows = np.searchsorted(pointer,start+at,side='right')-1
                        np.add.at(recovered,rows,f['X/data'][start:stop][at])
                source_ids = decode(f['obs/cell_id']) if 'cell_id' in f['obs'] else decode(f['obs'][f['obs'].attrs['_index']])
                mapping = pd.Index(source_ids).get_indexer(obs.cell_id)
                check(cohort+' all selected cell IDs match original HDF5', (mapping>=0).all() and obs.cell_id.is_unique)
                check(cohort+' all source PGAM5 raw counts and label memberships match',
                      np.array_equal(recovered[mapping],obs.PGAM5_counts)
                      and np.array_equal(raw[:,int(np.flatnonzero(genes=='PGAM5')[0])].toarray().ravel(),recovered[mapping]))
                rng=np.random.default_rng(1010)
                audit_rows=np.r_[rng.choice(np.flatnonzero(pos),10,replace=False),
                                 rng.choice(np.flatnonzero(neg),10,replace=False),rng.choice(len(obs),10,replace=False)]
                for row in audit_rows:
                    begin,end=pointer[mapping[row]:mapping[row]+2]
                    columns=f['X/indices'][begin:end];values=f['X/data'][begin:end]
                    for gene in ['DHFR','MKI67','NUSAP1','TOP2A','CENPK','TYMS']:
                        locations=np.flatnonzero(source_genes==gene)
                        expected=values[np.isin(columns,locations)].sum()
                        k=int(np.flatnonzero(genes==gene)[0])
                        check(cohort+' original '+gene+' count cell '+str(row),raw[row,k]==expected)
    evidence=pd.read_csv(R/'all_gene_candidate_evidence.csv.gz')
    support=np.column_stack(votes).sum(1)
    check('All cross-cohort recurrence votes independently recomputed',close(support,evidence.support_cohorts))
    stable=pd.read_csv(R/'patient_consistent_state_candidates.csv')
    check('PGAM5 included as model anchor and excluded from independent marker evidence',
          'PGAM5' not in set(stable.gene) and 'PGAM5' in set(pd.read_csv(R/'EXPLORATORY_signature_gene_roles.csv').gene))
    check('Final independent candidates are the three supported exploratory genes',set(stable.gene)=={'MKI67','NUSAP1','DHFR'})
    inner=pd.read_csv(R/'inner_configuration_scores.csv')
    check('Every inner split excludes both held cohorts from training',all(row.outer_held not in row.training_inner.split(',')
          and row.held_inner not in row.training_inner.split(',') for _,row in inner.iterrows()))
    out=pd.read_csv(R/'nested_outer_predictions.csv.gz')
    check('Every outer split excludes its complete query cohort',all(row.dataset not in row.training_cohorts.split(',')
          and len(row.training_cohorts.split(','))==4 for _,row in out.iterrows()))
    check('Exactly3320 predictions with unique mixture IDs and valid nonnegative coefficients',
          len(out)==3320 and out.mixture_id.is_unique and out.predicted.between(-1e-8,1+1e-8).all())
    summary=pd.read_csv(R/'nested_validation_summary.csv')
    for _,row in summary.iterrows():
        d=out[out.dataset.eq(row.dataset)&out.unit.eq(row.unit)]
        standard=d[d.scenario.eq('standard')];zero=d[d.truth.eq(0)]
        actual=[abs(standard.predicted-standard.truth).mean()*100,
                pearsonr(standard.truth,standard.predicted).statistic,zero.predicted.quantile(.95)*100]
        check(row.dataset+' '+row.unit+' MAE/r/zero95 independently recomputed',
              close(actual,[row.standard_MAE_pp,row.standard_Pearson_r,row.all_zero_target_p95_percent]))
    for prefix in ['solver_audit_case','TCGA_solver_audit_case']:
        manifest=pd.read_csv(R/f'{prefix}_manifest.csv')
        for _,row in manifest.iterrows():
            data=np.load(R/f'{prefix}_{int(row.case)}.npz')
            A,b,c,w=data['A'],data['b'],data['coefficients'],data['weights']
            solver=str(data['solver'])
            check(prefix+' '+str(row.case)+' coefficients feasible',c.min()>=-1e-8 and close(c.sum(),1,atol=1e-7))
            if solver=='NNLS':
                alt=lsq_linear(A,b,bounds=(0,np.inf),method='bvls',tol=1e-10,max_iter=2000)
                check(prefix+' '+str(row.case)+' independent BVLS coefficients',alt.success and close(alt.x/alt.x.sum(),c,rtol=1e-5,atol=1e-5))
            else:
                grad=A.T@(w*(A@c-b))/len(b)+.001*c
                active=c>1e-6
                multiplier=grad[active].mean()
                violation=max(np.max(np.abs(grad[active]-multiplier)),
                              max(0,float(multiplier-grad[~active].min())) if (~active).any() else 0)
                check(prefix+' '+str(row.case)+' independent convex final-step KKT certificate',violation<2e-5)
    tcga=pd.read_csv(R/'EXPLORATORY_TCGA_coefficients_NOT_CELL_FRACTIONS.csv')
    components=pd.read_csv(R/'EXPLORATORY_TCGA_all_component_coefficients.csv').set_index('patient')
    check('Exactly371 unique TCGA primary patients, all outputs flagged unvalidated',
          len(tcga)==371 and tcga.patient.is_unique and (~tcga.usable_as_cellular_infiltration).all())
    check('TCGA coefficients sum to one and target column matches published table',
          close(components.sum(axis=1),1,atol=1e-7) and close(components.loc[tcga.patient,'TAM_PGAM5_detected'],tcga.exploratory_PGAM5_RNA_macrophage_coefficient))
    gates=pd.read_csv(R/'validation_gate_checks.csv')
    result=json.loads((R/'validation_result.json').read_text())
    check('Overall validation status agrees with every reported gate',not gates.passed.all()
          and result['failed_gates']==int((~gates.passed).sum()) and not result['valid_for_TCGA_cellular_abundance'])
    output={'status':'PASS','checks_passed':len(checks),'source_data_audited':not args.saved_results_only,
            'checks':checks,'meaning':'Implementation/split/unit checks, not biological validation; KKT certifies only the fixed-weight convex last IRLS step'}
    filename='audit_saved_results.json' if args.saved_results_only else 'audit.json'
    (R/filename).write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k!='checks'},indent=2))


if __name__=='__main__':
    main()
