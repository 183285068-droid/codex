"""GSE76427 HCC tumor projection and source OS/RFS survival analysis."""
from pathlib import Path
import csv,gzip,io,json,importlib.metadata as md
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from statsmodels.stats.multitest import multipletests
import survival_tools as survival

ROOT=Path(__file__).resolve().parent
INPUT=ROOT/'inputs'
TARGET='TAM_PGAM5_detected'


def write(df,name):
    df.to_csv(ROOT/name,index=False,lineterminator='\n')


def load():
    text=gzip.open(INPUT/'GSE76427_series_matrix.txt.gz','rt').read()
    raw=pd.read_csv(io.StringIO(text.split('!series_matrix_table_begin')[1].split('!series_matrix_table_end')[0].strip()),sep='\t',index_col=0)
    rows=[next(csv.reader([l],delimiter='\t')) for l in text.splitlines() if l.startswith('!Sample')]
    fields={r[0]:r[1:] for r in rows if r[0]!='!Sample_characteristics_ch1'}
    meta=pd.DataFrame({'GSM':fields['!Sample_geo_accession'],'sample_title':fields['!Sample_title'],
                       'source_tissue':fields['!Sample_source_name_ch1'],'platform':fields['!Sample_platform_id']})
    for r in rows:
        if r[0]=='!Sample_characteristics_ch1':
            k=next(v for v in r[1:] if ':' in v).split(':',1)[0].strip()
            meta[k]=[v.split(':',1)[1].strip() if ':' in v else '' for v in r[1:]]
    assert raw.shape==(47322,167) and raw.columns.tolist()==meta.GSM.tolist()
    assert np.isfinite(raw.to_numpy()).all() and raw.index.is_unique and raw.columns.is_unique
    assert meta.platform.eq('GPL10558').all()
    tumor=meta[meta.tissue.eq('primary hepatocellular carcinoma tumor')].copy()
    assert len(tumor)==115 and tumor['patient id'].is_unique
    write(meta,'all167_sample_metadata.csv')
    write(meta[~meta.GSM.isin(tumor.GSM)],'excluded52_adjacent_tissues.csv')
    return raw[tumor.GSM],tumor


def project(raw):
    annotation=pd.read_csv(INPUT/'GPL10558_annotation.tsv.gz',sep='\t',dtype=str)
    symbol=annotation.Symbol.fillna('').str.strip()
    usable=annotation[annotation.ID.isin(raw.index)&symbol.ne('')&symbol.ne('---')&~symbol.str.contains(r'///|;|\|',regex=True)].copy()
    usable['gene']=symbol.loc[usable.index]
    assert usable.ID.is_unique
    gene=raw.loc[usable.ID].copy();gene.index=usable.gene.to_numpy()
    gene=gene.groupby(level=0,sort=True).median()
    ref=pd.read_csv(INPUT/'reference_850x17.tsv',sep='\t',index_col=0)
    scales=pd.read_csv(INPUT/'reference_row_scales.csv',index_col=0).row_scale
    genes=ref.index[ref.index.isin(gene.index)]
    assert ref.shape==(850,17) and len(genes)==755 and 'PGAM5' in genes
    A=ref.loc[genes].to_numpy()/scales.loc[genes].to_numpy()[:,None]
    assert np.linalg.matrix_rank(A)==17
    variants,diagnostics={},[]
    for variant in ['negative_gene_floor_primary','signed_signal_sensitivity']:
        linear=gene.clip(lower=0) if variant=='negative_gene_floor_primary' else gene.copy()
        denominator=linear.sum(axis=0)
        assert denominator.gt(0).all()
        query=linear.div(denominator,axis=1)*10000
        coefficients=[]
        for sample in query.columns:
            y=query.loc[genes,sample].to_numpy()/scales.loc[genes].to_numpy()
            beta,residual=nnls(A,y,maxiter=10000)
            assert beta.sum()>0
            gradient=A.T@(A@beta-y)
            norm=max(1.,np.linalg.norm(A)*np.linalg.norm(y))
            assert gradient.min()/norm>-1e-9
            assert np.max(np.abs(gradient*beta))/max(1.,np.linalg.norm(y)**2)<1e-9
            coefficients.append(beta/beta.sum())
            diagnostics.append({'GSM':sample,'variant':variant,'relative_fit_residual':residual/np.linalg.norm(y),
                                'raw_beta_sum':beta.sum(),'normalization_denominator':denominator[sample]})
        coef=pd.DataFrame(coefficients,index=query.columns,columns=ref.columns)
        variants[variant]=coef
        coef.rename_axis('GSM').to_csv(ROOT/f'coefficients_{variant}.csv')
        if variant=='negative_gene_floor_primary':query.loc[genes].to_csv(ROOT/'used_query_primary.tsv.gz',sep='\t')
    write(pd.DataFrame(diagnostics),'fit_diagnostics.csv')
    write(usable[['ID','gene','Entrez_Gene_ID','RefSeq_ID']],'used_probe_mapping.csv')
    write(pd.DataFrame({'gene':ref.index,'measured':ref.index.isin(genes)}),'reference_gene_coverage.csv')
    ref.loc[genes].to_csv(ROOT/'used_reference.tsv',sep='\t');scales.loc[genes].to_csv(ROOT/'used_reference_row_scales.csv')
    summary={'tumor_samples':115,'excluded_adjacent':52,'raw_probes':len(raw),'usable_probes':len(usable),
             'mapped_genes':len(gene),'covered_reference_genes':len(genes),'PGAM5_measured':True,
             'PGAM5_probes':usable.loc[usable.gene.eq('PGAM5'),'ID'].tolist(),
             'raw_signal_min':float(raw.min().min()),'raw_signal_max':float(raw.max().max()),
             'negative_probe_values':int(raw.lt(0).sum().sum()),'negative_gene_values':int(gene.lt(0).sum().sum()),
             'negative_reference_values':int(gene.loc[genes].lt(0).sum().sum()),'reference_query_values':int(gene.loc[genes].size),
             'weighted_reference_rank':17,'weighted_reference_condition':float(np.linalg.cond(A)),
             'median_relative_fit_residual':float(pd.DataFrame(diagnostics).query("variant=='negative_gene_floor_primary'").relative_fit_residual.median()),
             'cellular_abundance_validated':False}
    return variants,summary


def endpoint_frame(patients,ep):
    time=patients[f'{ep}_time_months'];event=patients[f'{ep}_event']
    reason=pd.Series('',index=patients.index)
    reason.loc[~np.isfinite(time)|time.le(0)|~event.isin([0,1])]='Missing or invalid time/event'
    if ep=='RFS':reason.loc[reason.eq('')&patients.OS_event.eq(1)&time.gt(patients.OS_time_months)]='RFS follow-up exceeds known death'
    data=patients.loc[reason.eq('')].copy();data['time_months']=time.loc[data.index];data['event']=event.loc[data.index].astype(int)
    ex=patients.loc[reason.ne(''),['GSM','patient_id','group',f'{ep}_time_months',f'{ep}_event']].copy()
    ex['endpoint'],ex['exclusion_reason']=ep,reason.loc[ex.index]
    return data,ex


def main():
    raw,tumor=load();variants,summary=project(raw)
    primary=variants['negative_gene_floor_primary'][TARGET]
    patients=tumor.rename(columns={'patient id':'patient_id'}).copy()
    patients[TARGET]=patients.GSM.map(primary)
    cutoff=float(primary.median());patients['group']=np.where(patients[TARGET].gt(cutoff),'High','Low')
    patients['high']=patients.group.eq('High').astype(int)
    patients['coefficient_per_SD']=(patients[TARGET]-primary.mean())/primary.std(ddof=1)
    patients['age_per_10_years']=pd.to_numeric(patients['age (years)'])/10
    patients['male']=pd.to_numeric(patients['gender (1=m, 2=f)']).eq(1).astype(int)
    for ep in ['OS','RFS']:
        suffix=ep.lower()
        patients[f'{ep}_event']=pd.to_numeric(patients[f'event_{suffix}'],errors='coerce')
        patients[f'{ep}_time_months']=pd.to_numeric(patients[f'duryears_{suffix}'],errors='coerce')*12
    write(patients,'patient_coefficients_clinical.csv')
    frames,stats,cox,exclusions,risks={},[],[],[],[]
    for ep in ['OS','RFS']:
        data,ex=endpoint_frame(patients,ep);frames[ep]=data
        assert len(data)==(115 if ep=='OS' else 107)
        row,risk=survival.compare(data,ep,'main_negative_floor',cutoff)
        stats.append(row);risks.extend(risk);cox.extend(survival.models(data,ep));exclusions.append(ex)
        write(data,f'{ep}_analysis_patients.csv')
    stats=pd.DataFrame(stats);stats['logrank_BH_q']=multipletests(stats.logrank_p,method='fdr_bh')[1]
    survival.CUTOFF=cutoff;survival.plot(frames,stats)
    write(stats,'survival_statistics.csv')
    cox=pd.DataFrame(cox);target=cox.covariate.isin(['high','coefficient_per_SD']);assert target.sum()==8
    cox.loc[target,'target_test_BH_q_across_8']=multipletests(cox.loc[target,'p'],method='fdr_bh')[1]
    write(cox,'Cox_binary_continuous_age_sex_adjusted.csv')
    write(pd.concat(exclusions,ignore_index=True),'endpoint_excluded_patients.csv')
    sensitivities=[]
    alternative=patients.copy();score=variants['signed_signal_sensitivity'][TARGET];alt_cutoff=float(score.median())
    alternative['high']=alternative.GSM.map(score).gt(alt_cutoff).astype(int);alternative['group']=np.where(alternative.high,'High','Low')
    for ep in ['OS','RFS']:
        data,_=endpoint_frame(alternative,ep);row,risk=survival.compare(data,ep,'signed_signal_sensitivity',alt_cutoff)
        sensitivities.append(row);risks.extend(risk)
    # Source-quality sensitivity: keep the discordant RFS censoring time exactly as published.
    data=patients.loc[patients.RFS_event.isin([0,1])&patients.RFS_time_months.gt(0)].copy()
    data['event']=data.RFS_event.astype(int);data['time_months']=data.RFS_time_months
    row,risk=survival.compare(data,'RFS','all108_source_RFS_including_beyond_death',cutoff)
    sensitivities.append(row);risks.extend(risk)
    write(pd.DataFrame(sensitivities),'survival_sensitivity.csv');write(pd.DataFrame(risks),'independent_logrank_risk_sets.csv')
    summary.update({'median_cutoff_all115':cutoff,'high_all115':int(patients.high.sum()),'low_all115':int((patients.high==0).sum()),
                    'nonzero':int((primary>0).sum()),'zero':int((primary==0).sum()),'OS_patients':115,'RFS_patients':107,
                    'RFS_missing':int(patients.RFS_event.isna().sum()),'RFS_beyond_death_excluded':['HCC29'],
                    'signed_signal_group_changes':int((patients.high!=alternative.high).sum()),
                    'source_OS_definition':'From resection to death from any cause, per original paper; GEO duryears_os is in years.',
                    'source_RFS_definition':'GEO event_rfs/duryears_rfs; original paper focused on OS, does not provide a detailed RFS composite-event definition.',
                    'packages':{n:md.version(n) for n in ['numpy','pandas','scipy','matplotlib','lifelines','statsmodels']}})
    (ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2));print(stats.to_string(index=False));print(pd.DataFrame(sensitivities).to_string(index=False))


if __name__=='__main__':main()
