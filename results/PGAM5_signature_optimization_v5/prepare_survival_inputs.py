"""Freeze an outcome-blind grouping protocol; preserve survival source provenance."""
from pathlib import Path
import json
import shutil
import hashlib
import numpy as np
import pandas as pd
import pyreadr

R = Path(__file__).resolve().parent
OUT = R/'survival_exploratory'
SOURCE = Path('/workspace/scratch/TCGA_LIHC_candidate30_survival')


def main():
    OUT.mkdir(exist_ok=True)
    coefficients = pd.read_csv(R/'EXPLORATORY_TCGA_coefficients_NOT_CELL_FRACTIONS.csv')
    x = coefficients.exploratory_PGAM5_RNA_macrophage_coefficient
    assert len(coefficients) == 371 and coefficients.patient.is_unique and np.isfinite(x).all()
    protocol = {'population':'371 outcome-blind selected TCGA-LIHC primary-tumor patients',
        'predictor':'v5 robust_FCLS exploratory PGAM5 RNA-detected macrophage coefficient; reference failed validation, NOT cellular infiltration',
        'grouping':'High > global 371-patient coefficient median; Low <= median; same grouping for both endpoints, no optimal outcome cutpoint',
        'global_median':float(x.median()), 'global_IQR':float(x.quantile(.75)-x.quantile(.25)),
        'endpoints':['OS','PFI'], 'PFS_available':False,
        'endpoint_rule':'finite time >0 days and event exactly 0 or 1; no imputations; independent patient-level join',
        'primary':'Two-sided unweighted log-rank for OS and PFI; BH across these 2 tests. KM Greenwood log-log 95% CI, censor markers and risk table',
        'Cox':'High/low unadjusted descriptive HR; continuous coefficient per full-cohort IQR unadjusted; high/low and continuous models adjusted for age/10 years, male and pathological stage II, III-IV vs I. Efron ties; complete cases only',
        'secondary_multiplicity':'BH over 6 target-term Cox tests: continuous unadjusted, high/low adjusted, continuous adjusted, across both endpoints. Unadjusted high/low Cox accompanies primary log-rank effect estimation',
        'PH':'Rank-time proportional-hazards diagnostic for each covariate; report all p without selecting methods for significant outcomes',
        'clinical_source':'Same UCSCXenaShiny mirror used in previous tasks. Age at initial pathologic diagnosis, gender, exact pathological AJCC stage. III-IV pooled; missing/discrepancy stages excluded in adjusted models',
        'scope':'Exploratory coefficient-survival association only. Neither prognostic significance nor adjustment validates cell abundance, PGAM5 specificity or causality; no independent survival cohort',
        'posthoc_addendum':{'trigger':'Added after original model diagnostics showed PFI stage covariates violated PH; not a prespecified primary test',
            'analysis':'Continuous coefficient per original global IQR with age/sex adjustment and stage-stratified baseline hazards; both OS and PFI retained',
            'multiplicity':'BH over 2 posthoc target tests; original primary and secondary families unchanged'}}
    assert protocol['global_IQR'] > 0
    (OUT/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    survival_source = SOURCE/'TCGA_CDR_mirror_survival.csv.gz'
    clinical_source = SOURCE/'UCSCXenaShiny_source/data/tcga_clinical.rda'
    shutil.copy2(survival_source,OUT/survival_source.name)
    raw = pyreadr.read_r(clinical_source)['tcga_clinical']
    raw = raw[raw.type.eq('LIHC') & raw.patient.isin(coefficients.patient)].copy()
    fields = ['sample','patient','age_at_initial_pathologic_diagnosis','gender','ajcc_pathologic_tumor_stage']
    raw[fields].to_csv(OUT/'clinical_all_sample_aliases.csv',index=False)
    primary = raw[raw['sample'].str[13:15].eq('01')][fields].copy()
    assert primary.patient.is_unique
    primary.to_csv(OUT/'TCGA_LIHC_clinical_primary.csv',index=False)
    provenance = {'survival_repository':'https://github.com/openbiox/UCSCXenaShiny',
        'mirror_commit':'e012e4a3e77702863dc639bfc8ee23c8e48e3231',
        'survival_source':'data/tcga_surv.rda, Xena Toil TCGA-CDR collection',
        'clinical_source':'data/tcga_clinical.rda, Xena Toil clinical collection',
        'CDR_publication':'https://doi.org/10.1016/j.cell.2018.02.052',
        'input_hashes':{str(p):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
                        for p in [survival_source,clinical_source,R/'EXPLORATORY_TCGA_coefficients_NOT_CELL_FRACTIONS.csv']}}
    (OUT/'source_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(json.dumps({'median':protocol['global_median'],'IQR':protocol['global_IQR'],
                      'High':int(x.gt(x.median()).sum()),'Low':int(x.le(x.median()).sum())},indent=2))


if __name__ == '__main__':
    main()
