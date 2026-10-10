"""Full competitors, RNA-status-preserving state splits and original RNA truth."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import pearsonr
from solver_core import fit

R = Path(__file__).resolve().parent
P = json.loads((R/'protocol.json').read_text())
S = Path(P['source_root'])
COHORTS = P['primary_cohorts']
POS = 'TAM_PGAM5_detected_'
NEG = 'TAM_PGAM5_undetected_'


def mixture_truths(cohort,obs,fold_annotations):
    table = pd.read_csv(S/'mixture_index.csv')
    original = table[table.dataset.eq(cohort)]
    candidate = obs.cell_id.map(fold_annotations.set_index('cell_id').candidate_state_member).eq(True).to_numpy(bool)
    total = obs.total_counts.to_numpy(float)
    rng = np.random.default_rng(20261010+1000*COHORTS.index(cohort))
    records,checks,position,audited = [],[],0,0
    for donor,d in obs.groupby('donor',sort=True):
        ix = d.index.to_numpy()
        mac,pg = d.harmonized_group.eq('Macrophage').to_numpy(),d.PGAM5_counts.gt(0).to_numpy()
        pools = {'positive':ix[mac&pg],'negative':ix[mac&~pg],'other':ix[~mac],
                 'PGAM5_nonmac':ix[~mac&pg],'cycling_nonmac':ix[~mac&d.cycling_flag.to_numpy()]}
        eligible = len(pools['positive'])>=5 and len(pools['negative'])>=20 and len(pools['other'])>=200
        cases = []
        if eligible:
            cases += [('standard',f,.25,'other',rep) for f in [0,.005,.01,.02,.05] for rep in range(3)]
            cases += [('macrophage_fraction_'+str(m),f,m,'other',rep)
                      for m in [.1,.5] for f in [0,.02,.05] for rep in range(2)]
        if len(pools['negative'])>=20:
            for group in ['PGAM5_nonmac','cycling_nonmac']:
                if len(pools[group])>=20:
                    cases += [('zero_'+group,0,.25,group,rep) for rep in range(3)]
        for group,part in d.groupby('fine_group',sort=True):
            j = part.index.to_numpy()
            if group=='TAM_PGAM5_detected' or len(j)<20: continue
            pools['pure_'+group] = j
            cases += [('pure_'+group,0,0,'pure_'+group,rep) for rep in range(2)]
        for scenario,f,m,other,rep in cases:
            if scenario.startswith('pure_'):
                chosen = rng.choice(pools[other],2000,replace=True)
            else:
                npg,nneg = int(round(f*2000)),int(round((m-f)*2000))
                pos = rng.choice(pools['positive'],npg,replace=True) if npg else np.array([],int)
                neg = rng.choice(pools['negative'],nneg,replace=True) if nneg else np.array([],int)
                chosen = np.r_[pos,neg,rng.choice(pools[other],2000-npg-nneg,replace=True)]
            ismac = obs.loc[chosen,'harmonized_group'].eq('Macrophage').to_numpy()
            ispg = ismac & obs.loc[chosen,'PGAM5_counts'].gt(0).to_numpy()
            isstate = ismac & candidate[chosen]
            isintersection = ispg & isstate
            assert np.isclose(ispg.mean(),f)
            if scenario=='standard' and rep==0 and f in [0,.05] and audited<2:
                np.savez_compressed(R/f'membership_audit_{cohort}_{audited}.npz',
                    chosen_indices=chosen,full_library_counts=total[chosen],
                    PGAM5_raw=obs.loc[chosen,'PGAM5_counts'].to_numpy(),macrophage=ismac,
                    candidate_state=isstate, equalized_mixture_id=int(original.iloc[position].mixture_id),
                    dataset=np.array(cohort))
                audited += 1
            for unit in ['equalized_cell_fraction','library_RNA_contribution']:
                row = original.iloc[position]
                position += 1
                assert row.donor==donor and row.scenario==scenario and row.replicate==rep and row.unit==unit
                truth = lambda flag:float(flag.mean()) if unit=='equalized_cell_fraction' else float(total[chosen][flag].sum()/total[chosen].sum())
                assert np.isclose(truth(ispg),row.truth,atol=1e-12)
                records.append(row.to_dict()|{'truth_all_PGAM5_RNA':truth(ispg),
                    'truth_PGAM5_RNA_in_candidate_state':truth(isintersection),'truth_candidate_state':truth(isstate)})
            checks.append({'dataset':cohort,'donor':donor,'scenario':scenario,'replicate':rep,
                'sampled_cells':2000,'all_PGAM5_RNA_cells':int(ispg.sum()),
                'PGAM5_RNA_in_candidate_cells':int(isintersection.sum()),'candidate_state_cells':int(isstate.sum()),
                'library_total':float(total[chosen].sum())})
    assert position==len(original)
    return pd.DataFrame(records),pd.DataFrame(checks)


def new_group_labels(obs,annotation):
    labels = obs.fine_group.to_numpy(dtype=object).copy()
    mac = obs.harmonized_group.eq('Macrophage').to_numpy()
    subset = obs.loc[mac].copy()
    state = subset.cell_id.map(annotation.set_index('cell_id').candidate_state_member).to_numpy()
    assert pd.notna(state).all()
    pg = subset.PGAM5_counts.gt(0).to_numpy()
    base = np.where(pg,POS,NEG)
    suffix = np.where(state,'candidate_state','other_state')
    names = np.char.add(base.astype(str),suffix)
    cycling = subset.cycling_flag.to_numpy(bool)
    names = np.where(cycling,np.char.add('Cycling_',names),names)
    labels[mac] = names
    return labels


def build(training,annotations,full,genes,unit,program_genes):
    entries = []
    for cohort in training:
        obs,raw,norm = full[cohort]
        labels = new_group_labels(obs,annotations[cohort])
        work = obs[['donor','total_counts']].copy()
        work['group'] = labels
        for (group,donor),d in work.groupby(['group','donor'],sort=True):
            ix = d.index.to_numpy()
            expression = np.asarray(norm[ix].sum(0)).ravel() if unit=='equalized_cell_fraction' else np.asarray(raw[ix].sum(0)).ravel().astype(float)
            denominator = len(ix) if unit=='equalized_cell_fraction' else d.total_counts.sum()/10000
            detection = np.asarray((raw[ix]>0).sum(0)).ravel()
            entries.append({'cohort':cohort,'donor':donor,'group':group,'cells':len(ix),
                            'expression_sum':expression,'denominator':denominator,'detection_sum':detection})
    info = pd.DataFrame([{k:v for k,v in e.items() if not k.endswith('_sum')} for e in entries])
    support = info.groupby('group').agg(cells=('cells','sum'),donors=('donor','nunique'))
    active = set(support[(support.cells>=50)&(support.donors>=2)].index)
    def mapped(group):
        if group in active: return group
        parent = group[8:] if group.startswith('Cycling_') else group
        if parent in active: return parent
        if parent.startswith(POS): return POS+'other_state'
        if parent.startswith(NEG): return NEG+'other_state'
        return 'Unknown_other'
    for e in entries: e['mapped_group'] = mapped(e['group'])
    groups = sorted(set(e['mapped_group'] for e in entries))
    profiles,frequencies,support_rows = [],[],[]
    for group in groups:
        cohort_means,cohort_freq = [],[]
        used = [e for e in entries if e['mapped_group']==group]
        for cohort in training:
            subset = [e for e in used if e['cohort']==cohort]
            if not subset: continue
            donor_profiles,donor_freq = [],[]
            for donor in sorted(set(e['donor'] for e in subset)):
                d = [e for e in subset if e['donor']==donor]
                donor_profiles.append(np.stack([e['expression_sum'] for e in d]).sum(0)/sum(e['denominator'] for e in d))
                donor_freq.append(np.stack([e['detection_sum'] for e in d]).sum(0)/sum(e['cells'] for e in d))
            cohort_means.append(np.stack(donor_profiles).mean(0))
            cohort_freq.append(np.stack(donor_freq).mean(0))
        profiles.append(np.stack(cohort_means).mean(0))
        frequencies.append(np.stack(cohort_freq).mean(0))
        support_rows.append({'group':group,'cells':sum(e['cells'] for e in used),
            'donor_labels':len(set(e['donor'] for e in used)),'cohorts':len(cohort_means)})
    reference = np.column_stack(profiles)
    frequency = np.column_stack(frequencies)
    base = np.array([g!='PGAM5' and not g.startswith(('MT-','RPL','RPS')) for g in genes])
    choices,markers = set(),[]
    for j,group in enumerate(groups):
        other = np.delete(reference,j,axis=1).max(1)
        contrast = np.log2((reference[:,j]+.05)/(other+.05))
        eligible = base & (frequency[:,j]>=.1)&(reference[:,j]>=.05)
        ix = np.flatnonzero(eligible)
        ix = ix[np.argsort(-(contrast[ix]*np.minimum(1,reference[ix,j])),kind='stable')[:25]]
        choices.update(ix)
        markers.extend({'gene':genes[k],'group':group,'log2_specificity':contrast[k],
                        'detection_fraction':frequency[k,j],'role':'descriptive reference discrimination'} for k in ix)
    choices.update(np.flatnonzero(np.isin(genes,program_genes)).tolist())
    choices.add(int(np.flatnonzero(genes=='PGAM5')[0]))
    features = np.asarray(sorted(choices))
    scale = np.maximum(np.sqrt(np.mean(reference[features]**2,axis=1)),.05)
    return {'ref':reference[features],'groups':groups,'features':features,'scale':scale,
            'A':reference[features]/scale[:,None],'support':support_rows,'markers':markers}


def summarize(table,keys):
    rows = []
    for key,d in table.groupby(keys,sort=True):
        key = key if isinstance(key,tuple) else (key,)
        standard = d[d.scenario.eq('standard')]
        zero = d[d.truth.eq(0)]
        r = pearsonr(standard.truth,standard.predicted).statistic if len(standard) and standard.truth.nunique()>1 and standard.predicted.nunique()>1 else np.nan
        rows.append(dict(zip(keys,key))|{'standard_donors':standard.donor.nunique(),
            'mixtures':len(d),'standard_MAE_pp':float(abs(standard.truth-standard.predicted).mean()*100),
            'standard_r':r,'all_zero_P95_percent':float(zero.predicted.quantile(.95)*100) if len(zero) else np.nan,
            'all_zero_max_percent':float(zero.predicted.max()*100) if len(zero) else np.nan})
    return pd.DataFrame(rows)


def main():
    genes = pd.read_csv(S/'common_genes.csv').gene.to_numpy()
    expression = np.load(S/'all_mixture_expression.npz')['values']
    full = {}
    for cohort in COHORTS:
        obs = pd.read_csv(S/f'{cohort}_metadata.csv.gz')
        raw = sparse.load_npz(S/f'{cohort}_raw_common.npz').astype(np.float32)
        norm = raw.multiply((10000/obs.total_counts.to_numpy(dtype=np.float32))[:,None]).tocsr()
        full[cohort] = obs,raw,norm
    rows,truth_rows,member_checks,cases = [],[],[],[]
    top_genes = pd.read_csv(R/'program_top_gene_roles.csv')
    for held in COHORTS+['ALL_TRAINING']:
        folder = R/held
        definition = json.loads((folder/'definition.json').read_text())
        training = definition['training_cohorts']
        train = pd.read_csv(folder/'training_macrophage_annotations.csv.gz')
        annotations = {c:train[train.dataset.eq(c)] for c in training}
        if held!='ALL_TRAINING':
            annotation = pd.read_csv(folder/'heldout_macrophage_annotations.csv.gz')
            truths,check = mixture_truths(held,full[held][0],annotation)
            truth_rows.append(truths)
            member_checks.append(check)
        program = top_genes[top_genes.fold.eq(held)&top_genes.selected_target].gene.tolist()
        for unit in ['equalized_cell_fraction','library_RNA_contribution']:
            model = build(training,annotations,full,genes,unit,program)
            pd.DataFrame(model['ref'],index=genes[model['features']],columns=model['groups']).to_csv(folder/f'EXPLORATORY_reference_{unit}.tsv',sep='\t')
            pd.DataFrame({'gene':genes[model['features']],'scale':model['scale']}).to_csv(folder/f'{unit}_row_scales.csv',index=False)
            pd.DataFrame(model['support']).to_csv(folder/f'{unit}_group_support.csv',index=False)
            pd.DataFrame(model['markers']).to_csv(folder/f'{unit}_marker_roles.csv',index=False)
            if held=='ALL_TRAINING': continue
            pg = np.array(['PGAM5_detected_' in g for g in model['groups']])
            intersection = np.array(['PGAM5_detected_candidate_state' in g for g in model['groups']])
            state = np.array([g.endswith('candidate_state') for g in model['groups']])
            assert pg.any()
            for _,row in truths[truths.unit.eq(unit)].iterrows():
                b = expression[int(row.mixture_id),model['features']].astype(float)/model['scale']
                coefficients,weights,residual = fit(model['A'],b)
                for estimand,mask,truth in [('all_PGAM5_RNA',pg,row.truth_all_PGAM5_RNA),
                        ('PGAM5_RNA_in_candidate_state',intersection,row.truth_PGAM5_RNA_in_candidate_state),
                        ('candidate_state_algorithm_label',state,row.truth_candidate_state)]:
                    rows.append(row.to_dict()|{'estimand':estimand,'truth':truth,'predicted':float(coefficients[mask].sum()),
                        'training_cohorts':','.join(training),'features':len(b),'components':len(coefficients),'residual':residual})
                if row.scenario=='standard' and row.replicate==0 and row.nominal_target_cell_fraction in [0,.05] and sum(c['dataset']==held and c['unit']==unit for c in cases)<4:
                    case = len(cases)
                    np.savez_compressed(R/f'fit_audit_case_{case}.npz',A=model['A'],b=b,coefficients=coefficients,
                                        weights=weights,PGAM5_RNA_mask=pg)
                    cases.append({'case':case,'dataset':held,'unit':unit,'mixture_id':int(row.mixture_id)})
            print('Reference validation',held,unit,'features',len(model['features']),'groups',len(model['groups']),flush=True)
    pred = pd.DataFrame(rows)
    pred.to_csv(R/'heldout_mixture_predictions.csv.gz',index=False)
    pd.concat(truth_rows,ignore_index=True).to_csv(R/'reconstructed_mixture_truths.csv',index=False)
    pd.concat(member_checks,ignore_index=True).to_csv(R/'mixture_membership_checks.csv',index=False)
    pd.DataFrame(cases).to_csv(R/'fit_audit_cases.csv',index=False)
    summary = summarize(pred,['dataset','unit','estimand'])
    summary.to_csv(R/'reference_validation_summary.csv',index=False)
    summarize(pred,['dataset','unit','estimand','scenario']).to_csv(R/'reference_validation_by_scenario.csv',index=False)
    gates = []
    for row in summary[summary.estimand.eq('all_PGAM5_RNA')].itertuples():
        for name,value,passed in [('MAE_pp<=1',row.standard_MAE_pp,row.standard_MAE_pp<=1),
                                 ('r>=0.7',row.standard_r,row.standard_r>=.7),
                                 ('zero_P95_percent<=0.5',row.all_zero_P95_percent,row.all_zero_P95_percent<=.5),
                                 ('independent_donors>=3',row.standard_donors,row.standard_donors>=3 and row.dataset!='GSE202642')]:
            gates.append({'dataset':row.dataset,'unit':row.unit,'criterion':name,'value':value,'passed':passed})
    assoc = pd.read_csv(R/'program_PGAM5_association_by_cohort.csv')
    stability = pd.read_csv(R/'bootstrap_stability.csv')
    selection = pd.read_csv(R/'candidate_selection_by_fold.csv')
    for held in COHORTS:
        definition = json.loads((R/held/'definition.json').read_text())
        a = assoc[assoc.fold.eq(held)&assoc.role.eq('heldout')&assoc.program.eq(definition['target_program'])].iloc[0]
        st = stability[stability.fold.eq(held)]
        for name,value,passed in [('training_selection_supported',definition['training_selection_eligible'],definition['training_selection_eligible']),
                                 ('held_cohort_PGAM5_support',a.cohort_support_passed,bool(a.cohort_support_passed)),
                                 ('bootstrap_cosine>=0.8',st.aligned_cosine.min(),st.aligned_cosine.min()>=.8),
                                 ('bootstrap_Jaccard>=0.6',st.target_cell_Jaccard.min(),st.target_cell_Jaccard.min()>=.6)]:
            gates.append({'dataset':held,'unit':'state','criterion':name,'value':value,'passed':passed})
    gate = pd.DataFrame(gates)
    gate.to_csv(R/'validation_gates.csv',index=False)
    result = {'computational_gates_all_passed':bool(gate.passed.all()),
        'failed_gates':int((~gate.passed).sum()),'total_gates':len(gate),
        'biological_PGAM5_positive_gold_standard_available':False,
        'bulk_platform_and_cell_RNA_calibrated':False,'usable_as_TCGA_cellular_abundance':False,
        'interpretation':'Exact RNA-detection labels and exploratory algorithmic programs/reference. Do not name inferred RNA-zero cells PGAM5-positive or claim validated TCGA cell abundance.'}
    (R/'validation_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    main()
