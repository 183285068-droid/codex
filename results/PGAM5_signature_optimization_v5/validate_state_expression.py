"""Test independent state genes on held-cohort macrophages and competitors."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata, spearmanr
from reference_core import ReferenceBuilder, COHORTS

R = Path(__file__).resolve().parent


def evaluate(cohort, obs, raw, genes, state_genes, phase):
    shared = [g for g in state_genes if g in genes]
    score = np.full(len(obs), np.nan)
    if shared:
        ix = genes.get_indexer(shared)
        x = raw[:, ix].toarray().astype(float)/obs.total_counts.to_numpy()[:, None]*10000
        score = np.log1p(x).mean(1)
    obs = obs.copy()
    obs['score'] = score
    donor_results, lineage_results = [], []
    for donor, d in obs.groupby('donor', sort=True):
        macro = d[d.harmonized_group.eq('Macrophage')]
        target = macro.PGAM5_counts.gt(0).to_numpy()
        n1, n0 = int(target.sum()), int((~target).sum())
        ranks = rankdata(macro.score.to_numpy())
        auc = (ranks[target].sum()-n1*(n1+1)/2)/(n1*n0) if n1 and n0 and shared else np.nan
        pg = np.log1p(macro.PGAM5_counts.to_numpy()/macro.total_counts.to_numpy()*10000)
        rho = spearmanr(macro.score, pg).statistic if n1 and macro.score.nunique()>1 else np.nan
        donor_results.append({'dataset': cohort, 'phase': phase, 'donor': donor,
            'state_genes': ','.join(state_genes), 'genes_measurable': len(shared), 'genes_selected': len(state_genes),
            'n_macrophages': len(macro), 'PGAM5_detected': n1, 'undetected': n0,
            'RNA_detection_AUC': auc, 'Spearman_score_PGAM5': rho,
            'eligible_support': bool(n1>=5 and n0>=20 and len(state_genes)>=5),
            'descriptive_support_n1_ge5': bool(n1>=5 and n0>=20)})
        group_column = 'fine_group' if 'fine_group' in d else 'harmonized_group'
        for group, part in d.groupby(group_column, sort=True):
            lineage_results.append({'dataset': cohort, 'phase': phase, 'donor': donor, 'group': group,
                'cells': len(part), 'mean_state_score': part.score.mean(),
                'p90_state_score': part.score.quantile(.9), 'state_genes': ','.join(state_genes)})
    saved = obs[obs.harmonized_group.eq('Macrophage')][['cell_id','donor','PGAM5_counts','total_counts','score']].copy()
    saved['dataset'] = cohort
    saved['phase'] = phase
    return donor_results, lineage_results, saved


def main():
    builder = ReferenceBuilder()
    genes = pd.Index(builder.genes)
    results, lineages, cells = [], [], []
    for held in COHORTS:
        train = tuple(c for c in COHORTS if c != held)
        _, ix, _ = builder.state_evidence(train)
        obs = pd.read_csv(R/f'{held}_metadata.csv.gz')
        raw = sparse.load_npz(R/f'{held}_raw_common.npz')
        d, l, saved = evaluate(held, obs, raw, genes, builder.genes[ix].tolist(), 'primary_cohort_holdout')
        results += d
        lineages += l
        if saved is not None:
            cells.append(saved)
    _, selected, _ = builder.state_evidence(tuple(COHORTS))
    for held in ['GSE125449','GSE140228_Droplet_HCC','GSE140228_Smartseq2_HCC','GSE146115']:
        obs = pd.read_csv(R/f'{held}_metadata.csv.gz')
        raw = sparse.load_npz(R/f'{held}_raw_measured.npz')
        measured = pd.Index(pd.read_csv(R/f'{held}_measured_genes.csv').gene)
        d, l, saved = evaluate(held, obs, raw, measured, builder.genes[selected].tolist(), 'supplementary_technology_or_overlap')
        results += d
        lineages += l
        if saved is not None:
            cells.append(saved)
    table = pd.DataFrame(results)
    table.to_csv(R/'state_discrimination_by_donor.csv', index=False)
    pd.DataFrame(lineages).to_csv(R/'state_score_competitor_profiles.csv', index=False)
    pd.concat(cells, ignore_index=True).to_csv(R/'heldout_macrophage_state_scores.csv.gz', index=False)
    summary = []
    for (phase, cohort), d in table.groupby(['phase','dataset']):
        descriptive = d[d.descriptive_support_n1_ge5]
        summary.append({'dataset': cohort, 'phase': phase, 'donor_labels': len(d),
            'descriptive_testable_donors': len(descriptive), 'selected_genes': int(d.genes_selected.max()),
            'support_eligible_donors': int(d.eligible_support.sum()),
            'median_AUC_descriptive': descriptive.RNA_detection_AUC.median(),
            'fraction_AUC_ge0.7_descriptive': descriptive.RNA_detection_AUC.ge(.7).mean(),
            'median_rho_descriptive': descriptive.Spearman_score_PGAM5.median(),
            'status': 'SUPPORT_UNTESTABLE_TOO_FEW_GENES' if d.genes_selected.max()<5 else 'EVALUATED'})
    pd.DataFrame(summary).to_csv(R/'state_discrimination_summary.csv', index=False)
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == '__main__':
    main()
