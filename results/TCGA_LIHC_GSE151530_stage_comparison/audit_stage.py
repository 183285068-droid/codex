"""Independent reconstruction of clinical joins, stage pooling and statistics."""
from pathlib import Path
import argparse
import gzip
import json
import numpy as np
import pandas as pd
import pyreadr
from scipy.stats import rankdata, kruskal, mannwhitneyu, spearmanr, chi2_contingency
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportion_confint

R = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('--estimate', type=Path, default=Path('/workspace/codex/results/GSE151530_only_TCGA_PGAM5_deconvolution/TCGA_LIHC_PGAM5_macrophage_estimates.csv'))
parser.add_argument('--clinical', type=Path, default=Path('/workspace/scratch/TCGA_LIHC_candidate30_survival/UCSCXenaShiny_source/data/tcga_clinical.rda'))
parser.add_argument('--geo', type=Path, default=Path('/workspace/scratch/TCGA_LIHC_candidate30_survival/GSE62944_clinical.txt.gz'))
args = parser.parse_args()
checks = []

def check(name, condition):
    assert condition, name
    checks.append(name)

def close(a, b):
    return np.allclose(a, b, atol=1e-12, rtol=1e-10, equal_nan=True)

original = pd.read_csv(args.estimate).set_index('patient').sort_index()
joined = pd.read_csv(R/'all371_estimates_and_stage_sources.csv').set_index('patient').sort_index()
check('Original 371 estimates and RNA aliquots unchanged', len(joined) == 371 and joined.index.is_unique and all(close(original[c],joined[c]) if pd.api.types.is_numeric_dtype(original[c]) else original[c].equals(joined[c]) for c in original.columns))
clinical = pyreadr.read_r(args.clinical)['tcga_clinical']
clinical = clinical[(clinical.type == 'LIHC') & (clinical['sample'].str.slice(13,15) == '01') & clinical.patient.isin(original.index)].set_index('patient').sort_index()
check('Exactly one primary LIHC clinical sample per patient', len(clinical) == 371 and clinical.index.is_unique)
check('Primary clinical sample matches selected RNA aliquot', clinical['sample'].equals(original.RNA_sample.str.slice(0,15).rename('sample')))
check('All primary raw stages faithfully copied', clinical.ajcc_pathologic_tumor_stage.fillna('MISSING').equals(joined.primary_stage_raw.fillna('MISSING').rename('ajcc_pathologic_tumor_stage')))

allowed = {f'Stage {major}{sub}':major for major in ['I','II','III','IV'] for sub in ['', 'A','B','C']}
for prefix in ['primary','GEO2015']:
    parsed = joined[prefix+'_stage_raw'].map(allowed)
    check(prefix+' original major stages independently reconstructed', parsed.fillna('MISSING').equals(joined[prefix+'_stage'].fillna('MISSING').rename(prefix+'_stage_raw')))
    merged = parsed.map({'I':'I','II':'II','III':'III-IV','IV':'III-IV'})
    check(prefix+' III and IV merged correctly', merged.fillna('MISSING').equals(joined[prefix+'_analysis_stage'].fillna('MISSING').rename(prefix+'_stage_raw')))

with gzip.open(args.geo, 'rt') as f:
    names = f.readline().strip('\r\n').split('\t')[3:]
    stage_values = None
    for line in f:
        fields = line.rstrip('\r\n').split('\t')
        if fields[0] == 'ajcc_pathologic_tumor_stage':
            stage_values = dict(zip(names, fields[3:]))
            break
geo_expected = original.RNA_sample.map(stage_values).replace('NA',np.nan)
check('2015 sensitivity raw stage matched exact selected RNA barcode', geo_expected.fillna('MISSING').equals(joined.GEO2015_stage_raw.fillna('MISSING').rename('RNA_sample')))
check('Only one disagreement where both major stages valid', int(joined.valid_both_but_different_major.sum()) == 1 and joined.index[joined.valid_both_but_different_major].tolist() == ['TCGA-FV-A2QQ'])

summary = pd.read_csv(R/'stage_descriptive_statistics.csv')
stats = pd.read_csv(R/'stage_overall_and_secondary_tests.csv')
pairs = pd.read_csv(R/'stage_pairwise_tests.csv')
for source, column, sizes in [('primary_Toil_clinical','primary_analysis_stage',[171,86,90]),('sensitivity_GEO2015','GEO2015_analysis_stage',[165,82,87])]:
    expected = joined[joined[column].notna()].copy()
    patients = pd.read_csv(R/(source+'_analysis_patients.csv')).set_index('patient').sort_index()
    check(source+' exact inclusion and all zero coefficients retained', patients.index.equals(expected.index) and close(patients.exploratory_target_fraction,expected.exploratory_target_fraction))
    check(source+' unique patients and correct stage grouping', patients.index.is_unique and patients.stage.equals(expected[column].rename('stage')))
    groups = [patients.loc[patients.stage == stage, 'exploratory_target_fraction'].to_numpy() for stage in ['I','II','III-IV']]
    check(source+' three group counts', list(map(len,groups)) == sizes)
    s = stats[stats.source == source].set_index('test')
    kw = kruskal(*groups)
    check(source+' Kruskal-Wallis H and asymptotic p', close([kw.statistic,kw.pvalue],s.loc['Kruskal_Wallis_all_stages',['statistic','asymptotic_p']].astype(float)))
    code = patients.stage.map({'I':1,'II':2,'III-IV':3}).to_numpy()
    y = patients.exploratory_target_fraction.to_numpy()
    sp = spearmanr(code,y)
    check(source+' ordered-stage rho and p', close([sp.statistic,sp.pvalue],s.loc['Spearman_ordered_stage',['statistic','asymptotic_p']].astype(float)))
    table = np.array([[(g>0).sum(),(g==0).sum()] for g in groups])
    chi = chi2_contingency(table, correction=False)
    check(source+' 3x2 nonzero chi-square statistic and p', close([chi.statistic,chi.pvalue],s.loc['Chi_square_nonzero_estimate_by_stage',['statistic','asymptotic_p']].astype(float)))
    for stage,g in zip(['I','II','III-IV'],groups):
        z = summary[(summary.source == source)&(summary.stage == stage)].iloc[0]
        lo,hi = proportion_confint(int((g>0).sum()),len(g),method='wilson')
        check(source+' '+stage+' descriptive counts, distribution and Wilson CI', close([len(g),(g>0).sum(),(g==0).sum(),g.mean()*100,np.median(g)*100,np.quantile(g,.25)*100,np.quantile(g,.75)*100,lo*100,hi*100],z[['patients','nonzero_estimates','zero_estimates','mean_relative_estimate_percent','median_relative_estimate_percent','Q1_relative_estimate_percent','Q3_relative_estimate_percent','nonzero_Wilson95CI_lower_percent','nonzero_Wilson95CI_upper_percent']].astype(float)))
    ps = pairs[pairs.source == source]
    check(source+' exactly three pairwise comparisons', len(ps) == 3)
    for row in ps.itertuples():
        ga = patients.loc[patients.stage == row.stage_A,'exploratory_target_fraction'].to_numpy()
        gb = patients.loc[patients.stage == row.stage_B,'exploratory_target_fraction'].to_numpy()
        u = rankdata(np.concatenate([ga,gb]))[:len(ga)].sum()-len(ga)*(len(ga)+1)/2
        m = mannwhitneyu(ga,gb,alternative='two-sided',method='asymptotic',use_continuity=False)
        check(source+' '+row.stage_A+'/'+row.stage_B+' independent ranks, U, effect size and p', close([u,m.pvalue,2*u/(len(ga)*len(gb))-1],[row.U,row.asymptotic_p,row.rank_biserial_A_vs_B]))
    check(source+' BH correction across three pairwise permutation p-values', close(multipletests(ps.permutation_p,method='fdr_bh')[1],ps.pairwise_BH_q_3))
    secondary = s.loc[['Spearman_ordered_stage','Chi_square_nonzero_estimate_by_stage']]
    check(source+' secondary family contains two tests with BH correction', close(multipletests(secondary.permutation_p,method='fdr_bh')[1],secondary.secondary_BH_q))
for name,table in [('overall',stats),('pairwise',pairs)]:
    check(name+' Monte Carlo p equals (extreme+1)/(20000+1)', table.permutations.eq(20000).all() and close((table.extreme_permutations+1)/(table.permutations+1),table.permutation_p))
check('Primary exclusions exactly 22 missing plus 2 stage discrepancies', joined.primary_stage_raw.isna().sum() == 22 and joined.primary_stage_raw.eq('[Discrepancy]').sum() == 2 and joined.primary_stage.isna().sum() == 24)
result = {'status':'PASS','checks_passed':len(checks),'checks':checks,'scope':'Clinical joins, unchanged NNLS coefficients, user-requested stage pooling, inclusion, descriptive statistics, rank statistics, Monte Carlo p formula and BH families. Coding audit does not validate biological reference.'}
(R/'audit.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
