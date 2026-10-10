"""Training-only NMF programs and exact PGAM5 RNA status; all failures retained."""
from pathlib import Path
import json
import hashlib
import warnings
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import nnls, linear_sum_assignment
from scipy.stats import fisher_exact
from sklearn.decomposition import NMF
from statsmodels.stats.multitest import multipletests

R = Path(__file__).resolve().parent
P = json.loads((R/'protocol.json').read_text())
S = Path(P['source_root'])
COHORTS = P['primary_cohorts']


def load():
    genes = pd.read_csv(S/'common_genes.csv').gene.to_numpy()
    data = {}
    for cohort in COHORTS:
        obs = pd.read_csv(S/f'{cohort}_metadata.csv.gz')
        macro = obs.harmonized_group.eq('Macrophage').to_numpy()
        raw = sparse.load_npz(S/f'{cohort}_raw_common.npz')[macro].tocsr()
        obs = obs[macro].reset_index(drop=True)
        assert np.array_equal(raw[:,np.flatnonzero(genes=='PGAM5')[0]].toarray().ravel(),obs.PGAM5_counts)
        norm = raw.astype(np.float32).multiply((10000/obs.total_counts.to_numpy(dtype=np.float32))[:,None]).tocsr()
        log = norm.copy()
        log.data = np.log1p(log.data)
        data[cohort] = {'obs':obs,'raw':raw,'norm':norm,'log':log}
    return genes,data


def features(training,genes,data,seed):
    obs = pd.concat([data[c]['obs'] for c in training],ignore_index=True)
    log = sparse.vstack([data[c]['log'] for c in training],format='csr')
    norm = sparse.vstack([data[c]['norm'] for c in training],format='csr')
    rng = np.random.default_rng(seed)
    selected = []
    for donor,d in obs.groupby('donor',sort=True):
        if len(d)>=30:
            selected.extend(rng.choice(d.index,min(200,len(d)),replace=False).tolist())
    selected = np.asarray(sorted(selected))
    mean = np.asarray(log[selected].mean(0)).ravel()
    variance = np.maximum(np.asarray(log[selected].power(2).mean(0)).ravel()-mean**2,0)
    freq = np.asarray((log[selected]>0).mean(0)).ravel()
    macro_mean = np.asarray(norm[selected].mean(0)).ravel()
    others = []
    for cohort in training:
        stats = np.load(S/f'{cohort}_group_statistics.npz')['values']
        index = pd.read_csv(S/f'{cohort}_group_stats_index.csv')
        for group,d in index.groupby('group',sort=True):
            if group.startswith('TAM_') or group.startswith('Cycling_TAM_'):
                continue
            supported = d[d.cells.ge(20)]
            if d.cells.sum()>=50 and supported.donor.nunique()>=2:
                profile = (stats[supported.index,0].astype(float)/supported.cells.to_numpy()[:,None]).mean(0)
                others.append(profile)
    max_other = np.stack(others).max(0)
    cycle = set((S/'cell_cycle_genes.txt').read_text().split())
    background = set('ALB APOA1 APOA2 APOC3 TTR RBP4 APOH AHSG FGB FGA FGG ORM1 ORM2 PRAP1 TF GSTA1 CES1 CD3D CD3E TRAC CD79A MS4A1 COL1A1 COL1A2 DCN PECAM1 VWF'.split())
    eligible = np.array([g!='PGAM5' and g not in cycle and g not in background and
                         not g.startswith(('MT-','RPL','RPS','IGH','IGK','IGL','HBA','HBB')) for g in genes])
    eligible &= (freq>=.1)&(macro_mean*4>=max_other)&(variance>0)
    rank = variance/np.maximum(mean,.01)
    ix = np.flatnonzero(eligible)
    ix = ix[np.argsort(-rank[ix],kind='stable')[:1000]]
    sd = np.maximum(np.sqrt(variance[ix]),.1)
    dense = np.minimum(log[:,ix].toarray()/sd,10).astype(np.float64)
    return obs,selected,ix,sd,dense


def fit_nmf(x,rank,seed):
    model = NMF(n_components=rank,init='nndsvda',solver='cd',beta_loss='frobenius',
                max_iter=2000,tol=.0001,random_state=seed)
    with warnings.catch_warnings(record=True) as caught:
        model.fit_transform(x)
    assert model.n_iter_<2000, 'NMF did not converge within the expanded numerical limit'
    h = model.components_.copy()
    h /= np.maximum(np.linalg.norm(h,axis=1),1e-12)[:,None]
    return h,{'iterations':int(model.n_iter_),'reconstruction_error':float(model.reconstruction_err_),
              'warnings':[str(w.message) for w in caught]}


def project(x,h):
    a = np.ascontiguousarray(h.T)
    return np.stack([nnls(a,row,maxiter=3000)[0] for row in x])


def state_masks(w,threshold):
    order = np.argsort(w,axis=1)
    winner = order[:,-1]
    first = w[np.arange(len(w)),winner]
    second = w[np.arange(len(w)),order[:,-2]]
    margin = (first-second)/np.maximum(first,1e-12)
    masks = (winner[:,None]==np.arange(w.shape[1])) & (w>=threshold) & (margin[:,None]>=.05)
    return masks,winner,margin


def association(obs,masks,fold,role):
    rows,donors = [],[]
    pg = obs.PGAM5_counts.to_numpy()>0
    normpg = obs.PGAM5_counts.to_numpy()/obs.total_counts.to_numpy()*10000
    for cohort,d in obs.groupby('dataset',sort=True):
        ci = d.index.to_numpy()
        for program in range(masks.shape[1]):
            inside = masks[:,program]
            subset = inside[ci]
            a,c = int(pg[ci][subset].sum()),int(pg[ci][~subset].sum())
            n1,n0 = int(subset.sum()),int((~subset).sum())
            rate = (a/n1)/(c/n0) if n1 and n0 and c else np.nan
            p = fisher_exact([[a,n1-a],[c,n0-c]]).pvalue if n1 and n0 else 1.
            local = []
            for donor,part in d.groupby('donor',sort=True):
                ix = part.index.to_numpy()
                flag = inside[ix]
                target,rest = int(flag.sum()),int((~flag).sum())
                positive = int(pg[ix].sum())
                ratio = (normpg[ix][flag].mean()+.05)/(normpg[ix][~flag].mean()+.05) if target and rest else np.nan
                eligible = target>=10 and rest>=20 and positive>=3
                row = {'fold':fold,'role':role,'dataset':cohort,'program':program,'donor':donor,
                    'state_cells':target,'other_macrophages':rest,'state_PGAM5_detected':int(pg[ix][flag].sum()),
                    'PGAM5_CP10k_ratio':ratio,'eligible':eligible}
                donors.append(row)
                if eligible: local.append(row)
            median = float(np.median([x['PGAM5_CP10k_ratio'] for x in local])) if local else np.nan
            direction = float(np.mean([x['PGAM5_CP10k_ratio']>1 for x in local])) if local else np.nan
            passed = bool(n1>=30 and a>=5 and len(local)>=3 and pd.notna(rate) and rate>=1.5 and
                          pd.notna(median) and median>=1.5 and direction>=.7)
            rows.append({'fold':fold,'role':role,'dataset':cohort,'program':program,
                'state_cells':n1,'other_macrophages':n0,'state_PGAM5_detected':a,'other_PGAM5_detected':c,
                'detection_rate_ratio':rate,'evaluable_donors':len(local),
                'median_donor_PGAM5_CP10k_ratio':median,'fraction_donors_positive':direction,
                'cycling_fraction':float(obs.loc[ci[subset],'cycling_flag'].mean()) if n1 else np.nan,
                'cell_level_Fisher_p_descriptive':p,'cohort_support_passed':passed})
    table = pd.DataFrame(rows)
    for cohort,ix in table.groupby('dataset').groups.items():
        table.loc[ix,'cell_level_BH_q_descriptive'] = multipletests(table.loc[ix,'cell_level_Fisher_p_descriptive'],method='fdr_bh')[1]
    return table,pd.DataFrame(donors)


def select_program(association_table,rank,ntraining):
    records = []
    for j in range(rank):
        d = association_table[association_table.program.eq(j)]
        median = d.median_donor_PGAM5_CP10k_ratio.median()
        records.append({'program':j,'supported_training_cohorts':int(d.cohort_support_passed.sum()),
                        'median_cohort_donor_ratio':float(median) if pd.notna(median) else 0.,
                        'selection_eligible':bool(d.cohort_support_passed.sum()>=np.ceil(.6*ntraining))})
    ranked = pd.DataFrame(records).sort_values(['supported_training_cohorts','median_cohort_donor_ratio','program'],ascending=[False,False,True])
    return int(ranked.iloc[0].program),ranked


def main():
    genes,data = load()
    pd.DataFrame({'gene':genes}).to_csv(R/'common_genes.csv',index=False)
    threshold_rows = []
    for cohort,part in data.items():
        o = part['obs']
        for cutoff in [1,2,3,5]:
            mask = o.PGAM5_counts.ge(cutoff)
            donor_count = o[mask].groupby('donor').size()
            threshold_rows.append({'dataset':cohort,'raw_count_cutoff':cutoff,'macrophages':len(o),
                'cells_above_cutoff':int(mask.sum()),'donors_with_any':len(donor_count),
                'donors_with_ge5':int(donor_count.ge(5).sum()),
                'median_full_library_above':float(o.loc[mask,'total_counts'].median()),
                'median_full_library_below':float(o.loc[~mask,'total_counts'].median()),
                'nonzero_RNA_cells_below_cutoff':int((~mask & o.PGAM5_counts.gt(0)).sum()),
                'interpretation':'RNA count threshold, not a validated biological subtype'})
    pd.DataFrame(threshold_rows).to_csv(R/'raw_count_threshold_sensitivity.csv',index=False)
    all_assoc,all_donor,selection_rows,stability_rows,fit_rows,program_genes = [],[],[],[],[],[]
    for fi,held in enumerate(COHORTS+['ALL_TRAINING']):
        training = [c for c in COHORTS if c!=held]
        folder = R/held
        folder.mkdir(exist_ok=True)
        obs,selected,ix,sd,x = features(training,genes,data,P['seed']+fi)
        h,info = fit_nmf(x[selected],8,P['seed']+fi)
        w = project(x,h)
        thresholds = np.quantile(w[selected],.7,axis=0)
        masks,winner,margin = state_masks(w,thresholds)
        assoc,donor = association(obs,masks,held,'training')
        target,ranking = select_program(assoc,8,len(training))
        ranking['fold'] = held
        selection_rows.append(ranking)
        all_assoc.append(assoc)
        all_donor.append(donor)
        fit_rows.append({'fold':held,'kind':'primary','rank':8,**info})
        print('NMF',held,'training cells',len(obs),'balanced',len(selected),'features',len(ix),
              'target',target,'supported cohorts',int(ranking.iloc[0].supported_training_cohorts),flush=True)
        bootstrap_confidence = np.zeros(len(obs))
        for repeat in [1,2]:
            rng = np.random.default_rng(P['seed']+100*fi+repeat)
            sampled = []
            for donor,d in obs.iloc[selected].groupby('donor',sort=True):
                sampled.extend(rng.choice(d.index,len(d),replace=True).tolist())
            hb,binfo = fit_nmf(x[np.asarray(sampled)],8,P['seed']+100*fi+repeat)
            similarity = h@hb.T
            a,b = linear_sum_assignment(-similarity)
            hb = hb[b[np.argsort(a)]]
            wb = project(x,hb)
            bt = np.quantile(wb[selected],.7,axis=0)
            mb,_,_ = state_masks(wb,bt)
            bootstrap_confidence += (mb[:,target]==masks[:,target])/2
            union = (mb[:,target]|masks[:,target]).sum()
            jaccard = (mb[:,target]&masks[:,target]).sum()/union if union else 0.
            stability_rows.append({'fold':held,'repeat':repeat,'target_program':target,
                'aligned_cosine':float(h[target]@hb[target]),'target_cell_Jaccard':float(jaccard)})
            fit_rows.append({'fold':held,'kind':'bootstrap'+str(repeat),'rank':8,**binfo})
            print('Bootstrap',held,repeat,'cos',h[target]@hb[target],'Jaccard',jaccard,flush=True)
        feature_table = pd.DataFrame({'gene':genes[ix],'SD':sd})
        feature_table.to_csv(folder/'features.csv',index=False)
        pd.DataFrame(h.T,index=genes[ix],columns=['Program_'+str(j) for j in range(8)]).to_csv(folder/'program_loadings.tsv',sep='\t')
        for j in range(8):
            top = np.argsort(-h[j],kind='stable')[:30]
            program_genes.extend({'fold':held,'program':j,'rank':rank+1,'gene':genes[ix[k]],'loading':h[j,k],
                                  'selected_target':j==target} for rank,k in enumerate(top))
        definition = {'training_cohorts':training,'held_cohort':held,'target_program':target,
            'rank':8,'activation_thresholds':thresholds.tolist(),'minimum_margin':.05,
            'training_selection_eligible':bool(ranking.iloc[0].selection_eligible),
            'definition':'Dominant program, fixed activation threshold and margin; RNA positivity is a separate raw-count field',
            'status':'CANDIDATE_NOT_BIOLOGICALLY_VALIDATED'}
        (folder/'definition.json').write_text(json.dumps(definition,indent=2)+'\n')
        obs['dominant_program'] = winner
        obs['mapping_margin'] = margin
        obs['candidate_program_score'] = w[:,target]
        obs['candidate_state_member'] = masks[:,target]
        obs['bootstrap_membership_agreement'] = bootstrap_confidence
        obs.to_csv(folder/'training_macrophage_annotations.csv.gz',index=False)
        if held!='ALL_TRAINING':
            o = data[held]['obs'].copy()
            xx = np.minimum(data[held]['log'][:,ix].toarray()/sd,10)
            ww = project(xx,h)
            mm,win,mar = state_masks(ww,thresholds)
            a,d = association(o,mm,held,'heldout')
            all_assoc.append(a)
            all_donor.append(d)
            o['dominant_program'],o['mapping_margin'] = win,mar
            o['candidate_program_score'],o['candidate_state_member'] = ww[:,target],mm[:,target]
            o.to_csv(folder/'heldout_macrophage_annotations.csv.gz',index=False)
        else:
            for cohort,o in obs.groupby('dataset',sort=True):
                o.to_csv(R/f'{cohort}_macrophage_annotations.csv.gz',index=False)
            for rank in [6,10]:
                hs,sinfo = fit_nmf(x[selected],rank,P['seed']+rank)
                ws = project(x,hs)
                ms,_,_ = state_masks(ws,np.quantile(ws[selected],.7,axis=0))
                a,d = association(obs,ms,'rank'+str(rank),'rank_sensitivity')
                all_assoc.append(a)
                all_donor.append(d)
                fit_rows.append({'fold':'ALL_TRAINING','kind':'rank_sensitivity','rank':rank,**sinfo})
                pd.DataFrame(hs.T,index=genes[ix]).to_csv(R/f'rank{rank}_program_loadings.tsv',sep='\t')
    pd.concat(all_assoc,ignore_index=True).to_csv(R/'program_PGAM5_association_by_cohort.csv',index=False)
    pd.concat(all_donor,ignore_index=True).to_csv(R/'program_PGAM5_association_by_donor.csv',index=False)
    pd.concat(selection_rows,ignore_index=True).to_csv(R/'candidate_selection_by_fold.csv',index=False)
    pd.DataFrame(stability_rows).to_csv(R/'bootstrap_stability.csv',index=False)
    pd.DataFrame(fit_rows).to_csv(R/'NMF_fit_diagnostics.csv',index=False)
    pd.DataFrame(program_genes).to_csv(R/'program_top_gene_roles.csv',index=False)
    hashes = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [S/'source_provenance.json',S/'common_genes.csv',S/'cell_cycle_genes.txt']}
    (R/'source_provenance.json').write_text(json.dumps({'base_source_provenance':json.loads((S/'source_provenance.json').read_text()),
        'base_artifact_hashes':hashes,'protocol_SHA256':hashlib.sha256((R/'protocol.json').read_bytes()).hexdigest()},indent=2)+'\n')
    print('Discovery completed',flush=True)


if __name__ == '__main__':
    main()
