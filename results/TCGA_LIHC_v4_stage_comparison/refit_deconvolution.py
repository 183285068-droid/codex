"""Reapply the frozen v4 reference to the saved, audited TCGA counts.

Raw counts here are only the 778 fitting genes. The CP10k denominator is
the full 23,368-gene assigned-count total from the prior RNA selection.
Never normalize these counts by the 778-gene subset sum.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.optimize import nnls, lsq_linear

R = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('--input-root', type=Path, default=R/'inputs')
args = parser.parse_args()
I = args.input_root
ref = pd.read_csv(I/'reference_850genes_17components.tsv', sep='\t', index_col=0)
sc = pd.read_csv(I/'reference_row_scales.csv', index_col='gene')
raw = pd.read_csv(I/'bulk_reference_gene_counts.csv.gz', index_col='patient')
previous = pd.read_csv(I/'previous_patient_scores.csv').set_index('patient')
previous_weights = pd.read_csv(I/'previous_component_estimates.csv', index_col='patient')
previous_cp10k = pd.read_csv(I/'previous_bulk_CP10k.csv.gz', index_col='patient')
assert ref.shape == (850,17) and ref.index.is_unique and ref.columns.is_unique
assert ref.index.equals(sc.index) and (sc.row_scale > 0).all()
assert raw.index.is_unique and raw.shape == (371,778)
assert raw.index.equals(previous.index) and previous.index.equals(previous_weights.index)
assert previous_weights.columns.equals(ref.columns)
genes = ref.index[ref.index.isin(raw.columns)]
assert len(genes) == 778 and len(ref.index.difference(genes)) == 72
assert set(raw.columns) == set(genes) and 'PGAM5' in genes
assert np.isfinite(raw.to_numpy()).all() and (raw.to_numpy() >= 0).all()
assert (previous.assigned_gene_counts > 0).all()
query = raw.loc[:,genes].div(previous.assigned_gene_counts, axis=0)*10000
cp10k_delta = float(np.max(np.abs(query.to_numpy()-previous_cp10k.loc[query.index,genes].to_numpy())))
assert cp10k_delta < 1e-9
scale = sc.loc[genes,'row_scale'].to_numpy()
A = ref.loc[genes].to_numpy()/scale[:,None]
target = ref.columns.get_loc('TAM_PGAM5_detected')
mac = [ref.columns.get_loc(g) for g in ['TAM_PGAM5_detected','TAM_PGAM5_undetected','Cycling_TAM_PGAM5_undetected']]
weights = []
diagnostics = []
independent = []
for j, (patient, row) in enumerate(query.iterrows()):
    b = row.to_numpy()/scale
    beta, residual = nnls(A, b, maxiter=5000)
    assert np.isfinite(beta).all() and beta.sum() > 0
    w = beta/beta.sum()
    weights.append(w)
    diagnostics.append({'patient':patient,'coefficient_sum_before_normalizing':float(beta.sum()),'weighted_residual':float(residual),'weighted_relative_residual':float(residual/np.linalg.norm(b)),'target_coefficient_before_normalizing':float(beta[target]),'measured_reference_genes':len(genes)})
    if j % 20 == 0:
        alternative = lsq_linear(A,b,bounds=(0,np.inf),method='bvls',tol=1e-10,max_iter=1000)
        assert alternative.success
        aw = alternative.x/alternative.x.sum()
        delta = float(np.max(np.abs(aw-w)))
        assert delta < 1e-5
        independent.append({'patient':patient,'max_normalized_coefficient_error':delta,'NNLS_target':float(w[target]),'independent_bounded_solver_target':float(aw[target])})
weights = pd.DataFrame(weights,index=query.index,columns=ref.columns)
assert np.allclose(weights.sum(axis=1),1,rtol=0,atol=1e-12)
previous_delta = float(np.max(np.abs(weights.to_numpy()-previous_weights.to_numpy())))
assert previous_delta < 1e-10
scores = previous.reset_index().copy()
scores['exploratory_target_fraction'] = weights.iloc[:,target].to_numpy()
scores['exploratory_total_macrophage_fraction'] = weights.iloc[:,mac].sum(axis=1).to_numpy()
scores['exploratory_target_within_macrophage_fraction'] = scores.exploratory_target_fraction/scores.exploratory_total_macrophage_fraction.replace(0,np.nan)
scores['group'] = np.where(scores.exploratory_target_fraction > 0,'High','Low')
assert scores.group.eq(previous.group.to_numpy()).all()
assert (scores.exploratory_target_fraction > 0).sum() == 108
assert scores.exploratory_target_fraction.eq(0).sum() == 263
scores.to_csv(R/'patient_scores.csv',index=False)
weights.to_csv(R/'all_reference_component_estimates.csv')
pd.DataFrame(diagnostics).to_csv(R/'sample_fit_diagnostics.csv',index=False)
pd.DataFrame(independent).to_csv(R/'independent_bounded_solver_checks.csv',index=False)
pd.DataFrame({'gene':ref.index,'bulk_available':ref.index.isin(genes)}).to_csv(R/'reference_bulk_gene_coverage.csv',index=False)
pd.DataFrame({'gene':ref.index.difference(genes)}).to_csv(R/'missing_reference_genes.csv',index=False)
hashes = {p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(I.iterdir()) if p.is_file()}
result = {'status':'PASS','reference_source_commit':'a7bf04f34197b3b888a59fc74bbbbfb187c99769','prior_TCGA_application_commit':'c02cfdbf761acd913223c71588a756d75451aafe','reference_training':'GSE151530 plus GSE149614; exact v4 reference retained','reference_genes_original':850,'reference_components':17,'bulk_matched_genes':778,'missing_reference_genes':72,'feature_coverage':778/850,'CP10k_denominator':'Frozen per-patient assigned counts across all23368measured genes, not fitting-gene subset sum','patients':371,'nonzero_target_estimates':108,'zero_target_estimates':263,'target':'TAM_PGAM5_detected','solver':'scipy NNLS with frozen training row_scale, all17components jointly fitted, beta normalized by its sum','max_saved_CP10k_difference':cp10k_delta,'max_saved_17component_estimate_difference':previous_delta,'independent_bounded_solver_checks':len(independent),'input_hashes':hashes,'limitations':'Frozen reference FAILED external biological validation and bulk gene coverage is below the original95%gate. Exploratory model mixture coefficients, not calibrated tissue cell fractions.'}
(R/'deconvolution_refit.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
