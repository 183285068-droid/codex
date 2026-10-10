"""HCC-only supplementary gene evidence, retaining platform and overlap limits."""
from pathlib import Path
import json
import hashlib
import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from scipy import sparse
from prepare_data import remap, symbols, de, digest

R = Path(__file__).resolve().parent
S = Path('/workspace/scratch')
SPECS = [
    ('GSE125449', S/'GSE125449/HCC_TAM_counts_QC.h5ad', 'Overlaps GSE151530; not an independent recurrence vote', 'UMI'),
    ('GSE140228_Droplet_HCC', S/'GSE140228/Droplet/macrophage_counts_QC.h5ad', 'HCC-only; CC excluded; shared donors with Smartseq2', 'UMI'),
    ('GSE140228_Smartseq2_HCC', S/'GSE140228/Smartseq2/macrophage_counts_QC.h5ad', 'HCC-only; read counts and partly shared donors', 'read'),
    ('GSE146115', S/'GSE146115/macrophage_counts_QC.h5ad', 'Low-support marker-inferred macrophages; C1 read counts', 'read')]


def main():
    common = pd.Index(pd.read_csv(R/'common_genes.csv').gene)
    summaries, provenance = [], []
    for cohort, path, role, count_type in SPECS:
        a = ad.read_h5ad(path)
        if 'Histology' in a.obs:
            a = a[a.obs.Histology.astype(str).eq('HCC')].copy()
        labels = symbols(a)
        measurable = common[common.isin(labels)]
        raw = remap(a.X, labels, measurable)
        ids = a.obs.cell_key.astype(str).to_numpy() if 'cell_key' in a.obs else a.obs_names.astype(str).to_numpy()
        donor_col = next((c for c in ['Donor', 'paper_sample', 'Sample'] if c in a.obs), None)
        donor = a.obs[donor_col].astype(str).to_numpy()
        obs = pd.DataFrame({'cell_id': ids, 'donor': donor, 'harmonized_group': 'Macrophage',
                            'PGAM5_counts': a.obs.PGAM5_counts.to_numpy(),
                            'total_counts': a.obs.total_counts.to_numpy()})
        assert np.array_equal(raw[:, measurable.get_loc('PGAM5')].toarray().ravel(), obs.PGAM5_counts)
        assert obs.cell_id.is_unique
        _, result = de(raw, obs.total_counts.to_numpy(float), obs, measurable, cohort)
        n1, n0 = int(obs.PGAM5_counts.gt(0).sum()), int(obs.PGAM5_counts.eq(0).sum())
        a1 = np.rint(result.positive_detection.to_numpy()*n1).astype(int)
        a0 = np.rint(result.undetected_detection.to_numpy()*n0).astype(int)
        cache = {}
        p = []
        for x, y in zip(a1, a0):
            key = (int(x), int(y))
            if key not in cache:
                cache[key] = fisher_exact([[x, n1-x], [y, n0-y]], alternative='two-sided').pvalue
            p.append(cache[key])
        result['detection_Fisher_p'] = p
        result['detection_Fisher_FDR'] = multipletests(p, method='fdr_bh')[1]
        result.to_csv(R/f'{cohort}_tie_corrected_DE.csv.gz', index=False)
        obs['cohort'] = cohort
        obs['count_type'] = count_type
        obs.to_csv(R/f'{cohort}_metadata.csv.gz', index=False)
        sparse.save_npz(R/f'{cohort}_raw_measured.npz', raw)
        pd.DataFrame({'gene': measurable}).to_csv(R/f'{cohort}_measured_genes.csv', index=False)
        summaries.append({'cohort': cohort, 'role': role, 'macrophages': len(obs), 'PGAM5_detected': n1,
            'donor_labels': len(set(donor)), 'count_type': count_type,
            'measured_primary_common_genes': len(measurable),
            'strict_up_excluding_PGAM5': int((result.strict_upregulated & result.gene.ne('PGAM5')).sum()),
            'strict_up_also_Fisher_FDR05': int((result.strict_upregulated & result.gene.ne('PGAM5') & result.detection_Fisher_FDR.lt(.05)).sum())})
        provenance.append({'cohort': cohort, 'source': str(path), 'bytes': path.stat().st_size, 'sha256': digest(path),
                           'role': role, 'count_type': count_type})
        print(summaries[-1], flush=True)
    pd.DataFrame(summaries).to_csv(R/'supplementary_dataset_summary.csv', index=False)
    (R/'supplementary_source_provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')


if __name__ == '__main__':
    main()
