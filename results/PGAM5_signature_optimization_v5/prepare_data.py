"""Prepare full competitors and macrophage DE with tie-corrected statistics."""
from pathlib import Path
import json
import hashlib
import gc
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

R = Path(__file__).resolve().parent
S = Path('/workspace/scratch')
ANNOTATIONS = Path('/workspace/codex/results/PGAM5_RNA_annotation_deconvolution_v4')
CONFIG = json.loads((R/'config.json').read_text())
SOURCES = {n: S/'PGAM5_myeloid_reference_v2'/f'{n}_harmonized_counts_QC.h5ad'
           for n in CONFIG['primary_cohorts'][:3]}
SOURCES.update(GSE242889=S/'GSE242889/all_tumor_counts_QC.h5ad',
               GSE202642=S/'GSE202642/all_library_counts_QC.h5ad')
T = 'TAM_PGAM5_detected'
N = 'TAM_PGAM5_undetected'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def symbols(a):
    labels = a.var.gene.astype(str).to_numpy()
    invalid = np.isin(labels, ['nan', '', 'None'])
    labels[invalid] = a.var_names.astype(str).to_numpy()[invalid]
    return labels


def remap(x, labels, genes):
    index = genes.get_indexer(labels)
    good = index >= 0
    mapping = sparse.csr_matrix((np.ones(good.sum(), dtype=np.float32),
              (np.flatnonzero(good), index[good])), shape=(len(labels), len(genes)))
    raw = (x @ mapping).tocsr()
    raw.sum_duplicates()
    assert np.isfinite(raw.data).all() and (raw.data >= 0).all()
    assert (raw.data == np.floor(raw.data)).all()
    return raw.astype(np.int32)


def de(raw, total, obs, genes, cohort):
    norm = raw.astype(float).multiply(1e4/total[:, None]).tocsr()
    mac = obs.harmonized_group.eq('Macrophage').to_numpy()
    positive = obs.PGAM5_counts.gt(0).to_numpy() & mac
    negative = ~positive & mac
    n1, n0 = positive.sum(), negative.sum()
    means1 = np.asarray(norm[positive].mean(0)).ravel()
    means0 = np.asarray(norm[negative].mean(0)).ravel()
    freq1 = np.asarray((raw[positive] > 0).mean(0)).ravel()
    freq0 = np.asarray((raw[negative] > 0).mean(0)).ravel()
    p = np.ones(len(genes))
    for start in range(0, len(genes), 256):
        stop = min(start+256, len(genes))
        xp = np.log1p(norm[positive, start:stop].toarray())
        xn = np.log1p(norm[negative, start:stop].toarray())
        p[start:stop] = mannwhitneyu(xp, xn, axis=0, alternative='two-sided',
                                   method='asymptotic', use_continuity=True).pvalue
    # A constant-zero gene carries no comparison information; some scipy
    # versions return NaN when the tie-corrected rank variance is exactly zero.
    p[(means1 == 0) & (means0 == 0)] = 1.0
    assert np.isfinite(p).all() and (p >= 0).all() and (p <= 1).all()
    q = multipletests(p, method='fdr_bh')[1]
    fc = np.log2((means1+.001)/(means0+.001))
    passes = (q < .05) & (fc >= 1) & (freq1 >= .1)
    result = pd.DataFrame({'cohort': cohort, 'gene': genes, 'positive_cells': n1,
        'undetected_cells': n0, 'mean_CP10k_positive': means1, 'mean_CP10k_undetected': means0,
        'positive_detection': freq1, 'undetected_detection': freq0,
        'log2FC_linear_CP10k': fc, 'Mann_Whitney_two_sided_p': p, 'BH_FDR': q,
        'strict_upregulated': passes})
    result.to_csv(R/f'{cohort}_tie_corrected_DE.csv.gz', index=False)
    return norm, result


def main():
    common = None
    for name, path in SOURCES.items():
        a = ad.read_h5ad(path, backed='r')
        known = set(symbols(a))
        common = known if common is None else common & known
        a.file.close()
    genes = pd.Index(sorted(common), name='gene')
    assert 'PGAM5' in genes and len(genes) > 15000
    pd.Series(genes).to_csv(R/'common_genes.csv', index=False, header=['gene'])
    summary, provenance, de_tables = [], [], []
    for name, path in SOURCES.items():
        a = ad.read_h5ad(path)
        annotation_path = ANNOTATIONS/f'{name}_all_cell_annotation.csv.gz'
        annotation = pd.read_csv(annotation_path)
        ids = a.obs.cell_id.astype(str) if 'cell_id' in a.obs else a.obs_names.astype(str)
        assert pd.Index(ids).is_unique and annotation.cell_id.is_unique
        index = pd.Index(ids).get_indexer(annotation.cell_id)
        assert (index >= 0).all()
        obs = annotation[['cell_id', 'donor', 'harmonized_group', 'fine_group',
                          'PGAM5_counts', 'total_counts', 'cycling_flag']].copy()
        obs['dataset'] = name
        raw = remap(a.X[index], symbols(a), genes)
        total = obs.total_counts.to_numpy(float)
        pg = raw[:, genes.get_loc('PGAM5')].toarray().ravel()
        assert np.array_equal(pg, obs.PGAM5_counts) and (total > 0).all()
        assert np.allclose(a.obs.total_counts.to_numpy()[index], total)
        norm, result = de(raw, total, obs, genes, name)
        de_tables.append(result)
        obs.to_csv(R/f'{name}_metadata.csv.gz', index=False)
        sparse.save_npz(R/f'{name}_raw_common.npz', raw)
        donor_stats, donor_rows, group_stats, group_rows = [], [], [], []
        for donor, d in obs.groupby('donor', sort=True):
            for group, v in d.groupby('fine_group', sort=True):
                ix = v.index.to_numpy()
                n = len(ix)
                linear_sum = np.asarray(norm[ix].sum(0)).ravel()
                raw_sum = np.asarray(raw[ix].sum(0)).ravel()
                detection = np.asarray((raw[ix] > 0).sum(0)).ravel()
                sq_sum = np.asarray(norm[ix].power(2).sum(0)).ravel()
                group_stats.append(np.stack([linear_sum, raw_sum, detection, sq_sum]).astype(np.float32))
                group_rows.append({'dataset': name, 'donor': donor, 'group': group,
                                   'cells': n, 'total_counts': float(total[ix].sum())})
            mac = d.harmonized_group.eq('Macrophage').to_numpy()
            pos = mac & d.PGAM5_counts.gt(0).to_numpy()
            neg = mac & ~pos
            for status, mask in [('detected', pos), ('undetected', neg)]:
                ix = d.index.to_numpy()[mask]
                mean = np.asarray(norm[ix].mean(0)).ravel() if len(ix) else np.zeros(len(genes))
                donor_stats.append(mean.astype(np.float32))
                donor_rows.append({'dataset': name, 'donor': donor, 'status': status, 'cells': len(ix)})
        pd.DataFrame(group_rows).to_csv(R/f'{name}_group_stats_index.csv', index=False)
        np.savez_compressed(R/f'{name}_group_statistics.npz', values=np.stack(group_stats))
        pd.DataFrame(donor_rows).to_csv(R/f'{name}_donor_contrasts_index.csv', index=False)
        np.savez_compressed(R/f'{name}_donor_means.npz', values=np.stack(donor_stats))
        mac = obs.harmonized_group.eq('Macrophage')
        summary.append({'dataset': name, 'role': 'primary', 'all_QC_cells': len(obs),
            'macrophages': int(mac.sum()), 'PGAM5_detected': int((mac & obs.PGAM5_counts.gt(0)).sum()),
            'donor_labels': obs.donor.nunique(), 'strict_upregulated_excluding_PGAM5':
                int((result.strict_upregulated & result.gene.ne('PGAM5')).sum())})
        provenance.append({'cohort': name, 'inputs': {str(p): {'bytes': p.stat().st_size, 'sha256': digest(p)}
                                                   for p in [path, annotation_path]}, 'genes': len(genes)})
        print('Prepared', summary[-1], flush=True)
        del a, raw, norm, group_stats, donor_stats
        gc.collect()
    combined = pd.concat(de_tables, ignore_index=True)
    combined.to_csv(R/'primary_all_gene_DE.csv.gz', index=False)
    pd.DataFrame(summary).to_csv(R/'dataset_summary.csv', index=False)
    (R/'source_provenance.json').write_text(json.dumps({'primary': provenance,
        'config_sha256': digest(R/'config.json'), 'scope': 'Retrospective optimization, not new untouched cohorts'}, indent=2)+'\n')


if __name__ == '__main__':
    main()
