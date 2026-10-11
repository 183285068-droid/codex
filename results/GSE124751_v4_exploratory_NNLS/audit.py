"""Independent reconstruction, BVLS checks and source integrity audit."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear

ROOT = Path(__file__).resolve().parent
checks = []


def check(name, condition):
    checks.append({'check': name, 'passed': bool(condition)})
    assert condition, name


def main():
    inputs = ROOT / 'inputs'
    manifest = pd.read_csv(ROOT / 'input_manifest.csv')
    for row in manifest.itertuples():
        check('Input SHA256: ' + row.file, hashlib.sha256((ROOT / row.file).read_bytes()).hexdigest() == row.sha256)
    ref = pd.read_csv(inputs / 'reference_850x17.tsv', sep='\t', index_col=0)
    scales = pd.read_csv(inputs / 'reference_row_scales.csv', index_col=0).row_scale
    log = pd.read_csv(inputs / 'rna.processed.hcc_I.log2.txt.gz', sep='\t', index_col=0)
    genes = ref.index[ref.index.isin(log.index)]
    check('763 measured genes; PGAM5 included', len(genes) == 763 and 'PGAM5' in genes)
    check('No duplicate/missing expression values', log.index.is_unique and log.columns.is_unique and np.isfinite(log.to_numpy()).all())
    linear = (2 ** log) - 1
    query = linear.div(linear.sum(axis=0), axis=1) * 10000
    saved = pd.read_csv(ROOT / 'used_query_primary.tsv.gz', sep='\t', index_col=0)
    check('Rebuild all 157178 used query values', np.allclose(query.loc[genes], saved, rtol=1e-12, atol=1e-12))
    check('Query order preserved', saved.index.equals(genes) and saved.columns.equals(log.columns))
    A = ref.loc[genes].to_numpy() / scales.loc[genes].to_numpy()[:, None]
    check('All 17 components retained and rank 17', A.shape == (763,17) and np.linalg.matrix_rank(A) == 17)
    coef = pd.read_csv(ROOT / 'coefficients_log2_plus1_primary.csv', index_col=0)
    check('Coefficient rows sum to one and are nonnegative', np.allclose(coef.sum(axis=1), 1.) and coef.min().min() >= 0)
    indices = np.unique(np.r_[np.linspace(0, 205, 24, dtype=int),
                             np.flatnonzero(coef.TAM_PGAM5_detected.to_numpy() > 0)[:4]])
    errors = []
    for i in indices:
        sample = log.columns[i]
        y = query.loc[genes, sample].to_numpy() / scales.loc[genes].to_numpy()
        independent = lsq_linear(A, y, bounds=(0,np.inf), method='bvls', tol=1e-12, max_iter=10000)
        check('Independent BVLS convergence: ' + sample, independent.success)
        scaled = independent.x / independent.x.sum()
        error = np.abs(scaled - coef.loc[sample]).max()
        errors.append(float(error))
        check('Independent BVLS coefficient agreement: ' + sample, error < 1e-8)
    meta = pd.read_csv(ROOT / 'sample_coefficients_and_metadata.csv')
    check('206 expression-to-metadata records', len(meta)==206 and meta.sample_id.is_unique and meta.GSM.is_unique)
    check('195 C-HCC and 11 HCC-IS', meta.subtype.value_counts().to_dict() == {'C-HCC':195,'HCC-IS':11})
    check('Exact 198 plus flagged 8 prefix aliases', meta.GEO_metadata_link.value_counts().to_dict() == {'exact_title':198,'L/H_prefix_alias_order_and_suffix':8})
    check('Flagged aliases preserve suffix and column order',
          meta.sample_id.tolist()==log.columns.tolist() and
          [x if not x.startswith('L') else 'H'+x[1:] for x in meta.sample_id] == meta.GEO_sample_title.tolist())
    cutoff = coef.TAM_PGAM5_detected.median()
    check('No optimized cutpoint: global median equals zero', cutoff == 0)
    check('68 High and 138 Low', meta.group.value_counts().to_dict() == {'Low':138,'High':68})
    check('Group assignments reconstruct', np.array_equal(meta.group.eq('High'), coef.TAM_PGAM5_detected.gt(cutoff)))
    check('PGAM5 observed values independently reconstructed', np.allclose(meta.PGAM5_log2_upper_quartile_value, log.loc['PGAM5']))
    for variant in ['zero_preserving_log2','unshifted_log2']:
        alt = pd.read_csv(ROOT / f'coefficients_{variant}.csv',index_col=0).TAM_PGAM5_detected
        check('Transform sensitivity group consistency: ' + variant,
              np.array_equal(alt.gt(alt.median()), meta.group.eq('High')))
    link = pd.read_csv(ROOT / 'earlier_GEO_exact_title_link_candidates.csv')
    matched = link.earlier_GSM.notna()
    check('198 exact earlier GEO tumor links', matched.sum()==198 and link.loc[matched,'earlier_GSM'].is_unique)
    check('Links do not use L/H alias for clinical matching',
          link.loc[matched,'earlier_sample_title'].str.removeprefix('Tumor (').str.removesuffix(')').equals(link.loc[matched,'sample_id']))
    summary = json.loads((ROOT/'summary.json').read_text())
    check('No fabricated survival statistics', summary['survival_analysis_performed'] is False and not list(ROOT.glob('KM*')))
    output = {'checks':checks,'passed':all(x['passed'] for x in checks),
              'independent_BVLS_samples':len(indices), 'BVLS_max_coefficient_error':max(errors),
              'reconstructed_query_values':int(saved.size),
              'survival_analysis_audited':False,
              'note':'Numerical checks do not establish cell-fraction validity or clinical endpoint availability.'}
    (ROOT/'audit_results.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k!='checks'},indent=2))


if __name__ == '__main__':
    main()
