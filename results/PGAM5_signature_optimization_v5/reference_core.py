"""Cohort-balanced references, training-only state selection, and fixed solvers."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import nnls, minimize

R = Path(__file__).resolve().parent
CONFIG = json.loads((R/'config.json').read_text())
COHORTS = CONFIG['primary_cohorts']
T = 'TAM_PGAM5_detected'
N = 'TAM_PGAM5_undetected'
RIDGE = .001


class ReferenceBuilder:
    def __init__(self, root=R):
        self.root = Path(root)
        self.genes = pd.read_csv(self.root/'common_genes.csv').gene.to_numpy()
        self.cycle = set((self.root/'cell_cycle_genes.txt').read_text().split())
        self.base = np.array([not (g.startswith(('MT-', 'RPL', 'RPS')) or g == 'PGAM5')
                              for g in self.genes])
        self.pg = int(np.flatnonzero(self.genes == 'PGAM5')[0])
        self.statistics, self.group_index, self.de, self.donor_means, self.donor_index = {}, {}, {}, {}, {}
        self.cache = {}
        for cohort in COHORTS:
            self.statistics[cohort] = np.load(self.root/f'{cohort}_group_statistics.npz')['values']
            self.group_index[cohort] = pd.read_csv(self.root/f'{cohort}_group_stats_index.csv')
            self.de[cohort] = pd.read_csv(self.root/f'{cohort}_tie_corrected_DE.csv.gz')
            assert np.array_equal(self.de[cohort].gene, self.genes)
            self.donor_means[cohort] = np.load(self.root/f'{cohort}_donor_means.npz')['values']
            self.donor_index[cohort] = pd.read_csv(self.root/f'{cohort}_donor_contrasts_index.csv')

    def state_evidence(self, training):
        votes = np.column_stack([self.de[c].strict_upregulated.to_numpy() for c in training])
        support = votes.sum(1)
        minimum = max(2, int(np.ceil(.6 * len(training))))
        contrasts, donors, cohorts = [], [], []
        for cohort in training:
            ind, values = self.donor_index[cohort], self.donor_means[cohort]
            for donor, sub in ind.groupby('donor', sort=True):
                pos = sub[sub.status.eq('detected')].index[0]
                neg = sub[sub.status.eq('undetected')].index[0]
                if ind.loc[pos, 'cells'] >= 5 and ind.loc[neg, 'cells'] >= 20:
                    contrasts.append(np.log2((values[pos].astype(float)+.05)/(values[neg].astype(float)+.05)))
                    donors.append(donor)
                    cohorts.append(cohort)
        assert len(contrasts)
        contrasts = np.stack(contrasts)
        fraction = (contrasts > 0).mean(0)
        median = np.median(contrasts, axis=0)
        eligible = self.base & (support >= minimum) & (fraction >= .7)
        result = pd.DataFrame({'gene': self.genes, 'support_cohorts': support,
            'minimum_required_cohorts': minimum, 'eligible_donors': len(donors),
            'donor_positive_effect_fraction': fraction, 'median_donor_log2_contrast': median,
            'cell_cycle_list_member': [g in self.cycle for g in self.genes],
            'stable_candidate': eligible})
        for i, cohort in enumerate(training):
            result[cohort+'_strict_up'] = votes[:, i]
            result[cohort+'_log2FC'] = self.de[cohort].log2FC_linear_CP10k.to_numpy()
            result[cohort+'_FDR'] = self.de[cohort].BH_FDR.to_numpy()
        ranks = support.astype(float) * np.maximum(median, 0) * fraction
        choices = np.flatnonzero(eligible)
        choices = choices[np.argsort(-ranks[choices], kind='stable')[:50]]
        return result, choices, pd.DataFrame(contrasts, index=donors, columns=self.genes)

    def model(self, training, top, unit, include_pgam5=True, remove_cycle=False):
        training = tuple(sorted(training))
        key = (training, top, unit, include_pgam5, remove_cycle)
        if key in self.cache:
            return self.cache[key]
        index = pd.concat([self.group_index[c] for c in training], ignore_index=True)
        statistics = np.concatenate([self.statistics[c] for c in training])
        support = index.groupby('group').agg(cells=('cells', 'sum'), donors=('donor', 'nunique'))
        active = set(support[(support.cells >= 50) & (support.donors >= 2)].index)
        assert T in active and N in active
        def mapped(group):
            if group in active:
                return group
            if group.startswith('Cycling_') and group[8:] in active:
                return group[8:]
            return 'Unknown_other'
        index['mapped_group'] = index.group.map(mapped)
        groups = sorted(index.mapped_group.unique())
        refs, frequencies, variability, support_rows = [], [], [], []
        for group in groups:
            cohort_profiles, cohort_frequency, profiles = [], [], []
            donor_count = 0
            for cohort in training:
                sub = index[index.dataset.eq(cohort) & index.mapped_group.eq(group)]
                donor_profiles = []
                for donor, part in sub.groupby('donor', sort=True):
                    ix = part.index.to_numpy()
                    stats = statistics[ix].astype(float).sum(0)
                    denominator = part.cells.sum() if unit == 'equalized_cell_fraction' else part.total_counts.sum()
                    profile = stats[0] / denominator if unit == 'equalized_cell_fraction' else stats[1] / denominator * 10000
                    donor_profiles.append(profile)
                    profiles.append(profile)
                    donor_count += 1
                if donor_profiles:
                    cohort_profiles.append(np.stack(donor_profiles).mean(0))
                    cohort_frequency.append(statistics[sub.index, 2].astype(float).sum(0) / sub.cells.sum())
            reference = np.stack(cohort_profiles).mean(0)
            refs.append(reference)
            frequencies.append(np.stack(cohort_frequency).mean(0))
            variability.append(np.stack(profiles).var(0))
            support_rows.append({'group': group, 'training_cohorts': len(cohort_profiles),
                                 'donor_labels': donor_count, 'cells': int(index[index.mapped_group.eq(group)].cells.sum())})
        ref = np.column_stack(refs)
        freq = np.column_stack(frequencies)
        variance = np.column_stack(variability)
        selected, marker_rows = set(), []
        for j, group in enumerate(groups):
            other = np.delete(ref, j, axis=1).max(1)
            contrast = np.log2((ref[:, j]+.05)/(other+.05))
            rank = contrast * np.minimum(1, ref[:, j])
            eligible = self.base & (freq[:, j] >= .1) & (ref[:, j] >= .05)
            choices = np.flatnonzero(eligible)
            choices = choices[np.argsort(-rank[choices], kind='stable')[:top]]
            selected.update(choices)
            marker_rows.extend({'gene': self.genes[k], 'group': group, 'role': 'reference_discrimination',
                'mean_expression': float(ref[k, j]), 'maximum_competitor_expression': float(other[k]),
                'log2_specificity': float(contrast[k]), 'detection_fraction': float(freq[k, j]),
                'cell_cycle_list_member': self.genes[k] in self.cycle} for k in choices)
        evidence, state_choices, _ = self.state_evidence(training)
        selected.update(state_choices)
        if include_pgam5:
            selected.add(self.pg)
        if remove_cycle:
            selected = {k for k in selected if self.genes[k] not in self.cycle}
        ix = np.array(sorted(selected), dtype=int)
        scale = np.maximum(np.sqrt(np.mean(ref[ix]**2, axis=1)), .05)
        model = {'training': training, 'unit': unit, 'top': top, 'groups': groups, 'ref': ref,
            'features': ix, 'scale': scale, 'matrix': ref[ix]/scale[:, None],
            'state_features': state_choices, 'evidence': evidence, 'markers': marker_rows,
            'support': support_rows, 'variance': variance[ix], 'include_pgam5': include_pgam5,
            'remove_cycle': remove_cycle}
        self.cache[key] = model
        return model


def simplex_fit(A, b, weights=None, initial=None):
    n = len(b)
    if weights is None:
        weights = np.ones(n)
    gram = A.T @ (weights[:, None]*A) / n + RIDGE*np.eye(A.shape[1])
    rhs = A.T @ (weights*b) / n
    if initial is None:
        initial = nnls(A, b, maxiter=3000)[0]
        initial = initial / initial.sum() if initial.sum() else np.ones(A.shape[1])/A.shape[1]
    objective = lambda c: .5*c@gram@c - rhs@c
    jacobian = lambda c: gram@c-rhs
    result = minimize(objective, initial, jac=jacobian, method='SLSQP',
        bounds=[(0, 1)]*A.shape[1], constraints={'type': 'eq', 'fun': lambda c: c.sum()-1,
        'jac': lambda c: np.ones(len(c))}, options={'ftol': 1e-11, 'maxiter': 500})
    assert result.success, result.message
    assert result.x.min() >= -1e-8 and np.isclose(result.x.sum(), 1, atol=1e-7)
    return np.maximum(result.x, 0), weights


def fit(model, y, solver):
    A = model['matrix']
    b = y[model['features']]/model['scale']
    weights = np.ones(len(b))
    if solver == 'NNLS':
        coefficients, _ = nnls(A, b, maxiter=3000)
        assert coefficients.sum() > 0
        coefficients = coefficients/coefficients.sum()
    elif solver == 'FCLS':
        coefficients, weights = simplex_fit(A, b)
    elif solver == 'robust_FCLS':
        coefficients, weights = simplex_fit(A, b)
        for _ in range(3):
            weights = 1/(1+(A@coefficients-b)**2)
            coefficients, _ = simplex_fit(A, b, weights, coefficients)
    else:
        raise ValueError(solver)
    return coefficients, float(np.linalg.norm(A@coefficients-b)/max(np.linalg.norm(b), 1e-12)), weights
