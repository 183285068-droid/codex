"""Frozen exploratory v4 projection. Does not fabricate missing follow-up."""
from pathlib import Path
import gzip, json, platform, importlib.metadata as md
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / 'inputs'
TARGET = 'TAM_PGAM5_detected'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'pdf.fonttype': 42, 'savefig.facecolor': 'white'})


def metadata(path):
    result, row = [], None
    for line in gzip.open(path, 'rt').read().splitlines():
        if line.startswith('^SAMPLE = '):
            if row is not None:
                result.append(row)
            row = {'GSM': line.split(' = ', 1)[1]}
        elif row is not None and line.startswith('!Sample_'):
            key, value = line.split(' = ', 1)
            row.setdefault(key, []).append(value)
    if row is not None:
        result.append(row)
    return result


def load():
    log = pd.read_csv(INPUT / 'rna.processed.hcc_I.log2.txt.gz', sep='\t', index_col=0)
    assert log.shape == (20501, 206) and log.index.is_unique and log.columns.is_unique
    assert np.isfinite(log.to_numpy()).all() and log.min().min() >= 0
    meta = []
    for row in metadata(INPUT / 'sample_metadata.soft.gz'):
        chars = dict(x.split(': ', 1) for x in row['!Sample_characteristics_ch1'])
        assert row['!Sample_platform_id'] == ['GPL11154']
        meta.append({'GSM': row['GSM'], 'GEO_sample_title': row['!Sample_title'][0],
                     'age': int(chars['age']), 'tissue': chars['tissue'],
                     'subtype': row['!Sample_description'][0]})
    meta = pd.DataFrame(meta)
    # Eight HCC-IS matrix IDs have L prefixes whereas GEO titles have H.
    # Audit the complete ordered list and preserve both names explicitly.
    canonical = [s if not s.startswith('L') else 'H' + s[1:] for s in log.columns]
    assert canonical == meta.GEO_sample_title.tolist()
    meta['sample_id'] = log.columns
    meta['GEO_metadata_link'] = np.where(meta.sample_id == meta.GEO_sample_title,
                                         'exact_title', 'L/H_prefix_alias_order_and_suffix')
    assert (meta.GEO_metadata_link != 'exact_title').sum() == 8
    assert meta.loc[meta.GEO_metadata_link != 'exact_title', 'subtype'].eq('HCC-IS').all()
    assert meta.GSM.is_unique and meta.sample_id.is_unique and len(meta) == 206
    assert meta.subtype.value_counts().to_dict() == {'C-HCC': 195, 'HCC-IS': 11}
    reference = pd.read_csv(INPUT / 'reference_850x17.tsv', sep='\t', index_col=0)
    scales = pd.read_csv(INPUT / 'reference_row_scales.csv', index_col=0).row_scale
    assert reference.shape == (850, 17) and (scales > 0).all()
    genes = reference.index[reference.index.isin(log.index)]
    assert len(genes) == 763 and 'PGAM5' in genes
    return log, meta, reference, scales, genes


def transform(log, variant):
    # GEO documents log2 upper-quartile normalized values, but not an offset.
    # Primary is a predeclared assumption, checked against two alternatives.
    if variant == 'log2_plus1_primary':
        linear = np.expm1(log * np.log(2))
    elif variant == 'zero_preserving_log2':
        linear = np.where(log > 0, np.exp2(log), 0.)
    elif variant == 'unshifted_log2':
        linear = np.exp2(log)
    else:
        raise ValueError(variant)
    linear = pd.DataFrame(linear, index=log.index, columns=log.columns)
    denominator = linear.sum(axis=0)
    return linear / denominator * 1e4, denominator


def solve(query, denominator, reference, scales, genes, variant):
    A = reference.loc[genes].to_numpy() / scales.loc[genes].to_numpy()[:, None]
    assert np.linalg.matrix_rank(A) == 17
    coefs, fits = [], []
    for i, sample in enumerate(query.columns):
        y = query.loc[genes, sample].to_numpy() / scales.loc[genes].to_numpy()
        beta, residual = nnls(A, y, maxiter=10000)
        assert beta.sum() > 0
        gradient = A.T @ (A @ beta - y)
        norm = max(1., np.linalg.norm(A) * np.linalg.norm(y))
        assert gradient.min() / norm > -1e-9
        assert np.max(np.abs(gradient * beta)) / max(1., np.linalg.norm(y)**2) < 1e-9
        coefs.append(beta / beta.sum())
        fits.append({'sample_id': sample, 'variant': variant,
                     'relative_fit_residual': residual / np.linalg.norm(y),
                     'raw_beta_sum': beta.sum(), 'normalization_denominator': denominator.iloc[i],
                     'KKT_min_gradient_scaled': gradient.min() / norm})
    coefs = pd.DataFrame(coefs, index=query.columns, columns=reference.columns)
    return coefs, pd.DataFrame(fits), A


def main():
    log, meta, ref, scales, genes = load()
    primary = None
    sensitivities, summaries = [], []
    for variant in ['log2_plus1_primary', 'zero_preserving_log2', 'unshifted_log2']:
        query, denominator = transform(log, variant)
        coefs, fits, A = solve(query, denominator, ref, scales, genes, variant)
        cutoff = float(coefs[TARGET].median())
        high = coefs[TARGET] > cutoff
        coefs.rename_axis('sample_id').to_csv(ROOT / f'coefficients_{variant}.csv')
        fits.to_csv(ROOT / f'fit_diagnostics_{variant}.csv', index=False)
        if primary is None:
            primary = coefs[TARGET].copy()
            primary_high = high.copy()
            main_cutoff = cutoff
            meta[TARGET] = primary.to_numpy()
            meta['group'] = np.where(high, 'High', 'Low')
            meta['PGAM5_log2_upper_quartile_value'] = log.loc['PGAM5'].to_numpy()
            query.loc[genes].to_csv(ROOT / 'used_query_primary.tsv.gz', sep='\t')
            summaries = {'samples': 206, 'subtypes': meta.subtype.value_counts().to_dict(),
                         'expression_genes': len(log), 'reference_genes': 850,
                         'GEO_exact_title_matches': 198, 'GEO_prefix_aliases': 8,
                         'covered_reference_genes': len(genes), 'PGAM5_measured': True,
                         'reference_components': 17, 'weighted_reference_rank': int(np.linalg.matrix_rank(A)),
                         'weighted_reference_condition': float(np.linalg.cond(A)),
                         'primary_transform_assumption': 'log2(x+1); invert as 2^value - 1',
                         'transform_offset_documented_by_GEO': False,
                         'primary_cutoff': cutoff,
                         'primary_nonzero': int((primary > 0).sum()),
                         'primary_zero': int((primary == 0).sum()),
                         'primary_high': int(high.sum()), 'primary_low': int((~high).sum()),
                         'primary_median_relative_fit_residual': float(fits.relative_fit_residual.median()),
                         'survival_analysis_performed': False,
                         'clinical_reason': 'GEO lacks time/event; linked cBioPortal LFS payload blocked by runtime network policy.',
                         'cellular_abundance_validated': False,
                         'source_datahub_commit': '18763f8eee86d4886d7210068f745f1e509498ee',
                         'python': platform.python_version(),
                         'packages': {x: md.version(x) for x in ['numpy','pandas','scipy','matplotlib']}}
        sensitivities.append({'variant': variant, 'median_cutoff': cutoff,
                              'nonzero': int((coefs[TARGET] > 0).sum()),
                              'high': int(high.sum()), 'low': int((~high).sum()),
                              'Spearman_r_vs_primary': float(spearmanr(primary, coefs[TARGET]).statistic),
                              'max_absolute_difference_vs_primary': float(np.abs(primary - coefs[TARGET]).max()),
                              'group_changes_vs_primary': int((primary_high != high).sum()),
                              'median_relative_fit_residual': float(fits.relative_fit_residual.median())})
    meta.to_csv(ROOT / 'sample_coefficients_and_metadata.csv', index=False)
    pd.DataFrame(sensitivities).to_csv(ROOT / 'transformation_sensitivity.csv', index=False)
    pd.DataFrame({'gene': ref.index, 'measured': ref.index.isin(genes)}).to_csv(
        ROOT / 'reference_gene_coverage.csv', index=False)
    ref.loc[genes].to_csv(ROOT / 'used_reference.tsv', sep='\t')
    scales.loc[genes].to_csv(ROOT / 'used_reference_row_scales.csv')
    # Same-patient earlier GEO records are evidence of a linkage opportunity,
    # never substitutes for unobserved survival time/event.
    earlier = []
    for row in metadata(INPUT / 'GSE54504_sample_metadata.soft.gz'):
        title = row['!Sample_title'][0]
        earlier.append({'earlier_GSM': row['GSM'], 'earlier_sample_title': title,
                        'sample_id': title.removeprefix('Tumor (').removesuffix(')')})
    earlier = pd.DataFrame(earlier)
    assert earlier.sample_id.is_unique
    link = meta[['sample_id','GSM','subtype','age']].merge(earlier, on='sample_id', how='left', validate='one_to_one')
    link.to_csv(ROOT / 'earlier_GEO_exact_title_link_candidates.csv', index=False)
    summaries['earlier_GEO_title_matches'] = int(link.earlier_GSM.notna().sum())
    summaries['earlier_GEO_title_matches_by_subtype'] = link.dropna(subset=['earlier_GSM']).subtype.value_counts().to_dict()
    (ROOT / 'summary.json').write_text(json.dumps(summaries, indent=2) + '\n')
    plot(meta, main_cutoff)
    print(json.dumps(summaries, indent=2))
    print(pd.DataFrame(sensitivities).to_string(index=False))


def plot(meta, cutoff):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={'width_ratios': [1, 1.1]})
    ax = axes[0]
    ax.hist(meta[TARGET], bins=24, color='#3688A9', edgecolor='white')
    ax.axvline(cutoff, color='#B5483B', ls='--', label=f'Global median = {cutoff:.4g}')
    ax.set(xlabel='Exploratory PGAM5 macrophage NNLS coefficient', ylabel='Tumor samples',
           title='206 HCC tumor expression samples')
    ax.legend(frameon=False)
    ax = axes[1]
    primary = meta[TARGET].to_numpy()
    for variant, color, label in [('zero_preserving_log2','#278A72','Zero-preserving 2^x'),
                                  ('unshifted_log2','#BF7735','Unshifted 2^x')]:
        alt = pd.read_csv(ROOT / f'coefficients_{variant}.csv')[TARGET]
        ax.scatter(primary, alt, s=19, alpha=.7, color=color, label=label)
    bound = max(primary.max(), alt.max()) * 1.07
    ax.plot([0,bound], [0,bound], color='grey', lw=1, ls='--')
    ax.set(xlabel='Primary coefficient: 2^x - 1', ylabel='Alternative coefficient',
           title='Predeclared transformation sensitivity')
    ax.legend(frameon=False)
    for ax in axes:
        ax.spines[['top','right']].set_visible(False)
    fig.suptitle('GSE124751: fixed 850-gene / 17-component reference', fontsize=14)
    fig.text(.06,.025,'Exploratory coefficients; not validated cell fractions. Survival follow-up not yet retrieved.',
             fontsize=10, color='#555555')
    fig.tight_layout(rect=[0,.09,1,.92])
    fig.savefig(ROOT / 'NNLS_distribution_and_sensitivity.png', dpi=220)
    fig.savefig(ROOT / 'NNLS_distribution_and_sensitivity.pdf')
    plt.close(fig)


if __name__ == '__main__':
    main()
