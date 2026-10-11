"""Independent reconstruction of the input projection, clinical joins, and statistics."""
from pathlib import Path
import io,gzip,json,hashlib
import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear
from scipy.stats import chi2

R=Path(__file__).resolve().parent;I=R/'inputs'
checks=[]
def check(name,condition,detail=None):
    assert condition,name
    checks.append({'check':name,'passed':True,'detail':detail})

def main():
    reference=pd.read_csv(I/'reference_850x17.tsv',sep='\t',index_col=0)
    scales=pd.read_csv(I/'reference_row_scales.csv',index_col=0).row_scale
    used=pd.read_csv(R/'used_reference.tsv',sep='\t',index_col=0)
    genes=used.index.tolist()
    check('frozen_reference_projection',np.array_equal(used.to_numpy(),reference.loc[genes].to_numpy()))
    original=R.parent/'TCGA_LIHC_PGAM5_exploratory_deconvolution_survival/source_exploratory_reference.tsv'
    if original.exists():check('original_v4_reference_byte_identity',original.read_bytes()==(I/'reference_850x17.tsv').read_bytes())
    check('common_569_genes',len(genes)==569)
    check('all_17_components',reference.shape==(850,17) and used.shape==(569,17))
    check('PGAM5_omitted_not_imputed','PGAM5' not in genes and 'PGAM5' in reference.index)
    query=pd.read_csv(R/'primary_array_query_normalized_intensity.tsv.gz',sep='\t',index_col=0)
    scores=pd.read_csv(R/'patient_coefficients_clinical.csv')
    query=query[scores.Affy_GSM]
    mapped_sets=[]
    maximum_error=0.
    for gpl in ['GPL3921','GPL571']:
        text=gzip.open(I/f'matrix_{gpl}.txt.gz','rt').read()
        raw=pd.read_csv(io.StringIO(text.split('!series_matrix_table_begin')[1].split('!series_matrix_table_end')[0].strip()),sep='\t',index_col=0)
        annotation=pd.read_csv(I/f'{gpl}_annotation.tsv.gz',sep='\t',dtype=str)
        mapping={}
        for _,row in annotation.iterrows():
            symbol=row['Gene Symbol']
            if pd.notna(symbol):
                symbol=symbol.strip()
                if symbol and symbol!='---' and '///' not in symbol and row.ID in raw.index:
                    mapping.setdefault(symbol,[]).append(row.ID)
        mapped_sets.append(set(mapping))
        check(gpl+'_no_PGAM5_symbol','PGAM5' not in mapping)
        check(gpl+'_no_PGAM5_Entrez_or_transcript',not annotation.to_string(index=False).find('192111')>=0 and not annotation.to_string(index=False).find('NM_138575')>=0)
        patients=scores.loc[scores.platform.eq(gpl),'Affy_GSM'].tolist()
        denominator=np.zeros(len(patients))
        summaries={}
        # Direct numpy probe medians, independent of analysis pandas groupby.
        for gene,probes in mapping.items():
            linear=np.exp2(np.median(raw.loc[probes,patients].to_numpy(),axis=0))
            denominator+=linear
            if gene in genes:summaries[gene]=linear
        denominator_saved=scores.set_index('Affy_GSM').loc[patients,'linear_gene_denominator'].to_numpy()
        check(gpl+'_full_mapped_gene_denominator',np.allclose(denominator,denominator_saved,rtol=1e-12))
        rebuilt=np.array([summaries[g] for g in genes])/denominator[None,:]*1e4
        observed=query.loc[genes,patients].to_numpy()
        err=float(np.max(np.abs(observed-rebuilt)));maximum_error=max(maximum_error,err)
        check(gpl+'_all_gene_sample_values',np.allclose(observed,rebuilt,rtol=1e-12,atol=1e-12),{'values':int(rebuilt.size),'max_abs_error':err})
    check('exact_common_reference_gene_order',genes==[g for g in reference.index if all(g in s for s in mapped_sets)])
    clinical=pd.read_csv(I/'GSE14520_Extra_Supplement.txt.gz',sep='\t',na_values=['.',''])
    tumors=clinical[clinical['Tissue Type'].eq('Tumor')].set_index('Affy_GSM')
    check('clinical_tumor_patients_247_unique',len(scores)==247 and scores.Affy_GSM.is_unique and scores.ID.is_unique and set(scores.Affy_GSM)==set(tumors.index))
    for field in ['ID','Survival status','Survival months','Recurr status','Recurr months','Gender','Age','TNM staging']:
        observed=scores.set_index('Affy_GSM')[field].reindex(tumors.index)
        check('clinical_'+field,observed.fillna('MISSING').astype(str).equals(tumors[field].fillna('MISSING').astype(str)))
    cutoff=scores.target_coefficient.median()
    check('outcome_blind_global_median_zero',cutoff==0.)
    check('group_assignment_all_patients',np.array_equal(scores.group,np.where(scores.target_coefficient>cutoff,'High','Low')))
    check('retained_zero_coefficients_156',int(scores.target_coefficient.eq(0).sum())==156)
    components=pd.read_csv(R/'all_component_coefficients.csv').set_index('Affy_GSM')
    check('component_nonnegative_sum_one',(components.to_numpy()>=0).all() and np.allclose(components.sum(axis=1),1,atol=1e-12))
    check('target_component_correct',np.allclose(scores.target_coefficient,components.loc[scores.Affy_GSM,'TAM_PGAM5_detected'],atol=1e-12))
    A=used.to_numpy()/scales.loc[genes].to_numpy()[:,None]
    probes=[]
    # Evenly sample both platforms, plus extrema and zero/nonzero cases.
    for gpl in ['GPL3921','GPL571']:
        sub=scores[scores.platform.eq(gpl)].sort_values('target_coefficient')
        probes+=sub.iloc[np.unique(np.linspace(0,len(sub)-1,12).astype(int))].Affy_GSM.tolist()
    solver_rows=[]
    for sample in dict.fromkeys(probes):
        y=query.loc[genes,sample].to_numpy()/scales.loc[genes].to_numpy()
        sol=lsq_linear(A,y,bounds=(0,np.inf),method='bvls',tol=1e-12,max_iter=10000)
        estimate=sol.x/sol.x.sum()
        error=float(np.max(np.abs(estimate-components.loc[sample].to_numpy())))
        check('independent_BVLS_'+sample,sol.success and error<1e-8,{'max_component_abs_error':error})
        # Scalar query normalization must not change sum-normalized NNLS estimates.
        direct=lsq_linear(A,y*37,bounds=(0,np.inf),method='bvls',tol=1e-12,max_iter=10000)
        check('sample_scaling_invariance_'+sample,np.allclose(estimate,direct.x/direct.x.sum(),atol=1e-8))
        solver_rows.append({'Affy_GSM':sample,'max_component_error':error})
    pd.DataFrame(solver_rows).to_csv(R/'independent_bounded_solver_checks.csv',index=False)
    stats=pd.read_csv(R/'survival_statistics.csv')
    curve=pd.read_csv(R/'KM_curve_values.csv')
    risk=pd.read_csv(R/'KM_numbers_at_risk.csv')
    for endpoint,tfield,efield in [('OS','Survival months','Survival status'),('DFS','Recurr months','Recurr status')]:
        times=pd.to_numeric(scores[tfield],errors='coerce');events=pd.to_numeric(scores[efield],errors='coerce')
        valid=np.isfinite(times)&times.gt(0)&events.isin([0,1])
        data=scores.loc[valid].copy();data['t']=times[valid];data['e']=events[valid]
        saved=pd.read_csv(R/f'{endpoint}_analysis_patients.csv')
        check(endpoint+'_valid_patient_set',set(saved.Affy_GSM)==set(data.Affy_GSM) and len(saved)==242)
        check(endpoint+'_group_counts',int(saved.group.eq('High').sum())==89 and int(saved.group.eq('Low').sum())==153)
        high=data[data.group.eq('High')];low=data[data.group.eq('Low')]
        obs,ex,var=0.,0.,0.
        for time in np.unique(data.loc[data.e.eq(1),'t']):
            n1=int(high.t.ge(time).sum());n0=int(low.t.ge(time).sum());n=n1+n0
            d1=int((high.t.eq(time)&high.e.eq(1)).sum());d0=int((low.t.eq(time)&low.e.eq(1)).sum());d=d1+d0
            obs+=d1;ex+=d*n1/n
            if n>1:var+=d*(n-d)*n1*n0/(n*n*(n-1))
        p=chi2.sf((obs-ex)**2/var,1)
        row=stats[stats.endpoint.eq(endpoint)].iloc[0]
        check(endpoint+'_independent_logrank',np.isclose(p,row.logrank_p,atol=1e-12),{'p':float(p)})
        for group in ['High','Low']:
            sub=data[data.group.eq(group)]
            values=curve[(curve.endpoint==endpoint)&(curve.group==group)]
            survival=1.
            for _,c in values.iterrows():
                n=sub.t.ge(c.time_months).sum();d=(sub.t.eq(c.time_months)&sub.e.eq(1)).sum()
                if n:survival*=1-d/n
                assert np.isclose(survival,c.survival,atol=1e-12)
                assert c.CI95_lower<=c.survival+1e-12 and c.CI95_upper+1e-12>=c.survival
            check(endpoint+'_'+group+'_KM_steps_and_CI_bounds',True)
            check(endpoint+'_'+group+'_risk_counts',all(int(sub.t.ge(r.time_months).sum())==r.at_risk for _,r in risk[(risk.endpoint==endpoint)&(risk.group==group)].iterrows()))
    ps=stats.logrank_p.to_numpy();order=np.argsort(ps);q=np.empty(2)
    q[order]=np.minimum.accumulate((ps[order]*2/np.arange(1,3))[::-1])[::-1]
    check('two_endpoint_BH_adjustment',np.allclose(q,stats.logrank_BH_q,atol=1e-12))
    summary={'checks_passed':len(checks),'checks_failed':0,'independent_solver_patients':len(solver_rows),
             'all_projected_input_values_checked':int(query.size),'max_query_abs_error':maximum_error,
             'checks':checks,'scope':'Implementation audit only. Original reference failed biological/abundance validation; no PGAM5 probe is identified on supplied platform annotations.'}
    (R/'audit_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='checks'},indent=2))

if __name__=='__main__':main()
