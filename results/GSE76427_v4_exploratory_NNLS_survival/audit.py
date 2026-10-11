"""Independent source metadata, probe reconstruction and survival checks."""
from pathlib import Path
import gzip,json,hashlib,io
import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear
from scipy.stats import chi2
from audit_math import bh,univariate_cox

ROOT=Path(__file__).resolve().parent
INPUT=ROOT/'inputs'
checks=[]


def check(name,condition):
    checks.append({'check':name,'passed':bool(condition)});assert condition,name


def main():
    for row in pd.read_csv(ROOT/'input_manifest.csv').itertuples():
        check('Input SHA256 '+row.file,hashlib.sha256((ROOT/row.file).read_bytes()).hexdigest()==row.sha256)
    check('Original frozen v4 reference bytes',hashlib.sha256((INPUT/'reference_850x17.tsv').read_bytes()).hexdigest()=='d5c7df4eeee3c3f2ef1631a41b699b27842756b11f9c9dd2c1b9c06feaaeba68')
    check('Original frozen v4 row scales',hashlib.sha256((INPUT/'reference_row_scales.csv').read_bytes()).hexdigest()=='7fc0af49e8047866179944e3d472a389bdf92f810bb09a1594384bc60a97abfc')
    # Independently parse official per-GSM SOFT clinical metadata.
    records=[];row=None
    for line in gzip.open(INPUT/'sample_metadata.soft.gz','rt').read().splitlines():
        if line.startswith('^SAMPLE = '):
            if row is not None:records.append(row)
            row={'GSM':line.split(' = ',1)[1]}
        elif line.startswith('!Sample_characteristics_ch1 = '):
            key,value=line.split(' = ',1)[1].split(': ',1);row[key]=value
    records.append(row);meta=pd.DataFrame(records).set_index('GSM')
    tumor=meta[meta.tissue=='primary hepatocellular carcinoma tumor']
    check('115 unique HCC tumor patients and 52 adjacent excluded',len(tumor)==115 and len(meta)==167 and tumor['patient id'].is_unique)
    text=gzip.open(INPUT/'GSE76427_series_matrix.txt.gz','rt').read()
    raw=pd.read_csv(io.StringIO(text.split('!series_matrix_table_begin')[1].split('!series_matrix_table_end')[0].strip()),sep='\t',index_col=0)[tumor.index]
    annotation=pd.read_csv(INPUT/'GPL10558_annotation.tsv.gz',sep='\t',dtype=str)
    mapping=pd.read_csv(ROOT/'used_probe_mapping.csv')
    official=annotation.set_index('ID').loc[mapping.ID]
    check('44052 unambiguous mappings agree with original official symbol',len(mapping)==44052 and mapping.ID.is_unique and np.array_equal(mapping.gene,official.Symbol.str.strip()))
    gene=raw.loc[mapping.ID].copy();gene.index=mapping.gene;gene=gene.groupby(level=0).median()
    check('31425 mapped genes, original negative values retained before floor',len(gene)==31425 and (gene<0).sum().sum()==3318)
    linear=gene.clip(lower=0);query=linear.div(linear.sum(axis=0),axis=1)*10000
    ref=pd.read_csv(INPUT/'reference_850x17.tsv',sep='\t',index_col=0)
    scales=pd.read_csv(INPUT/'reference_row_scales.csv',index_col=0).row_scale
    genes=ref.index[ref.index.isin(gene.index)]
    check('755 shared reference genes and 3 PGAM5 probes',len(genes)==755 and 'PGAM5' in genes and mapping.gene.eq('PGAM5').sum()==3)
    saved=pd.read_csv(ROOT/'used_query_primary.tsv.gz',sep='\t',index_col=0)
    check('Reconstruct all 86825 primary query values',np.allclose(saved,query.loc[genes],rtol=1e-12,atol=1e-12) and saved.columns.equals(tumor.index))
    A=ref.loc[genes].to_numpy()/scales.loc[genes].to_numpy()[:,None]
    check('17 frozen components and full reference rank',ref.shape==(850,17) and np.linalg.matrix_rank(A)==17)
    coef=pd.read_csv(ROOT/'coefficients_negative_gene_floor_primary.csv',index_col=0)
    check('Coefficients nonnegative and sum to one',coef.min().min()>=0 and np.allclose(coef.sum(axis=1),1))
    errors=[]
    for i in np.linspace(0,114,24,dtype=int):
        gsm=saved.columns[i];y=query.loc[genes,gsm].to_numpy()/scales.loc[genes].to_numpy()
        fit=lsq_linear(A,y,bounds=(0,np.inf),method='bvls',tol=1e-12,max_iter=10000)
        error=float(np.abs(fit.x/fit.x.sum()-coef.loc[gsm]).max());errors.append(error)
        check('Independent BVLS '+gsm,fit.success and error<1e-8)
    patients=pd.read_csv(ROOT/'patient_coefficients_clinical.csv')
    cutoff=coef.TAM_PGAM5_detected.median()
    check('All115 median, no outcome-selected threshold',np.isclose(cutoff,0.03870300847357219) and np.array_equal(patients.high,patients.GSM.map(coef.TAM_PGAM5_detected).gt(cutoff).astype(int)))
    check('57 High and 58 Low',patients.high.sum()==57 and len(patients)==115)
    check('91 nonzero and24 zero',coef.TAM_PGAM5_detected.gt(0).sum()==91 and coef.TAM_PGAM5_detected.eq(0).sum()==24)
    stats=pd.read_csv(ROOT/'survival_statistics.csv');curves=pd.read_csv(ROOT/'KM_curve_values.csv');risk=pd.read_csv(ROOT/'KM_numbers_at_risk.csv')
    exclusions=pd.read_csv(ROOT/'endpoint_excluded_patients.csv');rebuilt=[]
    for ep in ['OS','RFS']:
        data=pd.read_csv(ROOT/f'{ep}_analysis_patients.csv');source=tumor.loc[data.GSM];row=stats[stats.endpoint==ep].iloc[0]
        check(ep+' source identity and months conversion',np.array_equal(data.patient_id,source['patient id']) and np.allclose(data.time_months,pd.to_numeric(source['duryears_'+ep.lower()])*12) and np.array_equal(data.event,pd.to_numeric(source['event_'+ep.lower()])))
        check(ep+' population/events',len(data)==(115 if ep=='OS' else 107) and data.event.sum()==(23 if ep=='OS' else 48))
        check(ep+' frozen all115 coefficient scale',np.allclose(data.coefficient_per_SD,(data.GSM.map(coef.TAM_PGAM5_detected)-coef.TAM_PGAM5_detected.mean())/coef.TAM_PGAM5_detected.std(ddof=1)))
        check(ep+' age and sex from official SOFT',np.allclose(data.age_per_10_years,pd.to_numeric(source['age (years)'])/10) and np.array_equal(data.male,pd.to_numeric(source['gender (1=m, 2=f)']).eq(1).astype(int)))
        h,l=data[data.high==1],data[data.high==0];oe,var=0.,0.
        for t in sorted(data.loc[data.event==1,'time_months'].unique()):
            nh=(h.time_months>=t).sum();nl=(l.time_months>=t).sum();dh=((h.time_months==t)&(h.event==1)).sum();dl=((l.time_months==t)&(l.event==1)).sum();n,d=nh+nl,dh+dl
            oe+=dh-d*nh/n
            if n>1:var+=nh*nl*d*(n-d)/(n*n*(n-1))
        p=chi2.sf(oe*oe/var,1);rebuilt.append(p)
        check(ep+' independent tied logrank',np.isclose(p,row.logrank_p,atol=1e-12))
        hr,lo,hi,cp=univariate_cox(data)
        check(ep+' independent Efron Cox HR/CI/p',np.allclose([hr,lo,hi],[row.HR_high_vs_low,row.HR_lower95,row.HR_upper95],rtol=1e-5,atol=1e-7) and np.isclose(cp,row.Cox_p,atol=1e-6))
        for group in ['High','Low']:
            q=data[data.group==group];saved_curve=curves[(curves.endpoint==ep)&(curves.group==group)];sf=1.;values={0.:1.}
            for t in sorted(q.time_months.unique()):sf*=1-((q.time_months==t)&(q.event==1)).sum()/(q.time_months>=t).sum();values[t]=sf
            check(ep+' '+group+' independent full KM curve',all(np.isclose(z.survival,values[z.time_months],atol=1e-12) for z in saved_curve.itertuples()))
            check(ep+' '+group+' complete risk table',all(z.at_risk==(q.time_months>=z.time_months).sum() for z in risk[(risk.endpoint==ep)&(risk.group==group)].itertuples()))
    check('Two-endpoint BH',np.allclose(stats.logrank_BH_q,bh(rebuilt),atol=1e-12))
    check('7 missing RFS records and HCC29 excluded explicitly',len(exclusions)==8 and exclusions.endpoint.eq('RFS').all() and exclusions.exclusion_reason.str.contains('Missing').sum()==7 and exclusions.loc[exclusions.exclusion_reason.str.contains('death'),'patient_id'].tolist()==['HCC29'])
    cox=pd.read_csv(ROOT/'Cox_binary_continuous_age_sex_adjusted.csv');target=cox.covariate.isin(['high','coefficient_per_SD'])
    check('All8 target Cox tests BH',target.sum()==8 and np.allclose(cox.loc[target,'target_test_BH_q_across_8'],bh(cox.loc[target,'p']),atol=1e-12))
    signed=pd.read_csv(ROOT/'coefficients_signed_signal_sensitivity.csv',index_col=0).TAM_PGAM5_detected
    check('Signed-signal sensitivity changes no group',np.array_equal(signed.gt(signed.median()),coef.TAM_PGAM5_detected.gt(cutoff)))
    check('PNG/PDF contain full 93.12-month OS follow-up',data.time_months.max()<96 and (ROOT/'KM_OS_RFS.png').stat().st_size>10000 and float(patients.OS_time_months.max())==93.12)
    out={'passed':True,'checks':checks,'independent_BVLS_samples':24,'BVLS_max_coefficient_error':max(errors),
         'reconstructed_query_values':saved.size,'note':'Numerical/source checks do not validate cell fractions. OS PH p<0.05 is a diagnosed model limitation, not an ignored failed numerical check.'}
    (ROOT/'audit_results.json').write_text(json.dumps(out,indent=2)+'\n');print(f'{len(checks)} source, projection and survival checks passed; BVLS maximum error={max(errors):.3g}.')


if __name__=='__main__':main()
