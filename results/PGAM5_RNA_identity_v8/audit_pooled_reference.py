"""Independent checks of pooled selection, held mixtures, and fixed-weight QP."""
from pathlib import Path
import argparse
import json
import tempfile
import subprocess
import sys
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import linprog
from scipy.stats import pearsonr
from audit_results import manual_two_sided_mwu

R=Path(__file__).resolve().parent
P=json.loads((R/'protocol.json').read_text());S=Path(P['source_root'])
checks=[]


def check(name,value):
    assert bool(value),name
    checks.append(name)


def kkt(A,b,c,w,pg,target):
    gradient=A.T@(w*(A@c-b))/len(b)+.001*c
    budget=A[pg]*target;limit=b[pg]
    assert abs(c.sum()-1)<1e-7 and c.min()>-1e-8
    assert budget@c<=limit+1e-7
    tolerance=3e-5
    upper=[];rhs=[]
    for j in range(len(c)):
        # gradient+lambda+nu*budget>=0, equality for active coordinates.
        upper.append([-1.,-budget[j]]);rhs.append(gradient[j]+tolerance)
        if c[j]>1e-7:
            upper.append([1.,budget[j]]);rhs.append(-gradient[j]+tolerance)
    nubound=(0,None) if limit-budget@c<1e-7 else (0,0)
    result=linprog([0,0],A_ub=upper,b_ub=rhs,bounds=[(None,None),nubound],method='highs')
    return result.success


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--saved-only',action='store_true');args=parser.parse_args()
    versions=['continuous_rank','log2FC_1','log2FC_0p5'];units=['equalized_cell_fraction','library_RNA_contribution']
    split=pd.read_csv(R/'pooled_patient_sample_split.csv')
    check('each patient/sample label assigned once',not split.duplicated(['dataset','donor_or_sample']).any())
    check('GSE202642 independence remains unconfirmed',not split.loc[split.dataset.eq('GSE202642'),'independent_patient_confirmed'].any())
    members=pd.read_csv(R/'pooled_macrophage_DE_members.csv.gz')
    mapped=members.merge(split[['dataset','donor_or_sample','split']],left_on=['dataset','donor'],right_on=['dataset','donor_or_sample'],validate='many_to_one',suffixes=('','_expected'))
    check('pooled membership follows whole patient/sample split',mapped.split.eq(mapped.split_expected).all())
    check('all30533 macrophages retained',len(members)==30533)
    check('all1413 PGAM5 detections retained',members.PGAM5_counts.gt(0).sum()==1413)
    panel=set(pd.read_csv(R/'gene_roles.csv').gene)
    common=pd.Index(pd.read_csv(S/'common_genes.csv').gene) if not args.saved_only else None
    if not args.saved_only:
        genes=['PGAM5','DHFR','MKI67','NUSAP1','CD68','CSF1R','LPL']
        arrays=[];obsrows=[]
        for cohort in P['primary_cohorts']:
            obs=pd.read_csv(S/f'{cohort}_metadata.csv.gz');mac=obs.harmonized_group.eq('Macrophage')
            raw=sparse.load_npz(S/f'{cohort}_raw_common.npz')
            x=raw[mac.to_numpy()][:,common.get_indexer(genes)].toarray().astype(float)
            o=obs[mac].reset_index(drop=True)
            check(cohort+' direct source pooled membership',np.array_equal(o.cell_id,members.loc[members.dataset.eq(cohort),'cell_id']))
            # Match the declared multiplication by the full-library CP10k factor;
            # regrouping floating operations can spuriously break exact rank ties.
            arrays.append(x*(10000/o.total_counts.to_numpy())[:,None]);obsrows.append(o)
        values=np.concatenate(arrays);source=pd.concat(obsrows,ignore_index=True)
        check('direct raw source PG matches pooled membership',np.array_equal(source.PGAM5_counts,members.PGAM5_counts))
        np.savez_compressed(R/'actual_pooled_DE_audit_inputs.npz',values=values,PGAM5_counts=source.PGAM5_counts.to_numpy(),
            training=members.split.eq('training').to_numpy(),genes=np.array(genes))
    case=np.load(R/'actual_pooled_DE_audit_inputs.npz')
    for stage in ['POOLED_INTERNAL_TRAINING','POOLED_ALL_TRAINING']:
        folder=R/stage;de=pd.read_csv(folder/'pooled_DE.csv.gz')
        mask=case['training'] if stage=='POOLED_INTERNAL_TRAINING' else np.ones(len(case['PGAM5_counts']),dtype=bool)
        x=case['values'][mask];positive=case['PGAM5_counts'][mask]>0
        check(stage+' cell counts',de.RNA_detected_cells.eq(positive.sum()).all() and de.RNA_undetected_cells.eq((~positive).sum()).all())
        pv=de.Mann_Whitney_p.to_numpy();order=np.argsort(pv);q=np.empty(len(pv))
        q[order]=np.minimum(1,np.minimum.accumulate((pv[order]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1])
        check(stage+' independent full-universe BH',np.allclose(q,de.BH_FDR,atol=1e-12,rtol=1e-10))
        for j,gene in enumerate(case['genes']):
            row=de.set_index('gene').loc[gene]
            p=manual_two_sided_mwu(np.log1p(x[positive,j]),np.log1p(x[~positive,j]))
            fc=np.log2((x[positive,j].mean()+.001)/(x[~positive,j].mean()+.001))
            check(stage+' manual tie-corrected p '+gene,np.isclose(p,row.Mann_Whitney_p,atol=1e-10,rtol=1e-6))
            check(stage+' direct pooled log2FC '+gene,np.isclose(fc,row.log2FC,atol=1e-6))
        for level in versions:
            root=folder/level
            for unit in units:
                model=np.load(root/f'{unit}_model.npz')
                matrix=pd.read_csv(root/f'EXPLORATORY_reference_{unit}.tsv',sep='\t',index_col=0)
                expected=np.maximum(np.sqrt(np.mean(matrix.to_numpy()**2,axis=1)),.05)
                check(stage+level+unit+' reference scaling',np.allclose(expected,model['scale']) and np.allclose(matrix.to_numpy()/expected[:,None],model['A']))
                check(stage+level+unit+' PGAM5 retained','PGAM5' in model['genes'])
                check(stage+level+unit+' lineage/competition panel retained',(panel & set(de.gene))<=set(model['genes']))
                table=pd.read_csv(root/f'{unit}_gene_evidence.csv.gz')
                score=np.maximum(table.log2FC,0)*np.maximum(table.competitive_log2_ratio,0)*table.detected_group_detection
                check(stage+level+unit+' continuous rank independently computed',np.allclose(score,table.continuous_ranking_score))
                background=table.gene.eq('PGAM5')|table.gene.str.startswith(('MT-','RPL','RPS'))
                base=~background & ~table.cycle_list_member & ~table.stress_list_member & table.detected_group_detection.ge(.1)
                if level=='continuous_rank':eligible=base & score.gt(0)
                else:eligible=base & table.BH_FDR.lt(.05)&table.log2FC.ge(1 if level=='log2FC_1' else .5)&table.competitive_log2_ratio.ge(1)
                ix=np.flatnonzero(eligible);ix=ix[np.argsort(-score.to_numpy()[ix],kind='stable')[:50]]
                check(stage+level+unit+' training-only state feature ranking',np.array_equal(ix,np.flatnonzero(table.selected_candidate_state_feature)[np.argsort(-score.to_numpy()[table.selected_candidate_state_feature.to_numpy()],kind='stable')]))
                selected=pd.read_csv(root/f'{unit}_candidate_state_genes.csv')
                check(stage+level+unit+' selected gene list exact',np.array_equal(selected.gene,table.gene.iloc[ix]))
                definition=json.loads((root/f'{unit}_definition.json').read_text())
                check(stage+level+unit+' no protein/cross-cohort gate',definition['protein_validation_required'] is False and definition['cross_cohort_replication_required'] is False)
                if level=='continuous_rank':check(stage+unit+' no FDR/effect-size cutoff',definition['FDR_strict_max'] is None and definition['within_macrophage_log2FC_min'] is None)
    predictions=pd.read_csv(R/'pooled_internal_mixture_predictions.csv.gz')
    summary=pd.read_csv(R/'pooled_internal_validation_summary.csv')
    gates=pd.read_csv(R/'pooled_internal_validation_gates.csv')
    held=set(split.loc[split.split.eq('internal_validation'),'donor_or_sample'])
    check('all queries exclude training patient/sample labels',predictions.donor_or_sample.isin(held).all())
    check('prediction keys unique',not predictions.duplicated(['threshold_version','unit','mixture_id']).any())
    for (level,unit),d in predictions.groupby(['threshold_version','unit']):
        standard=d[d.scenario.eq('standard')];zero=d[d.truth.eq(0)]
        row=summary[summary.threshold_version.eq(level)&summary.unit.eq(unit)].iloc[0]
        check(level+unit+' independent MAE',np.isclose(np.abs(standard.truth-standard.predicted).mean()*100,row.standard_MAE_pp))
        r=np.corrcoef(standard.truth,standard.predicted)[0,1]
        check(level+unit+' independent correlation',np.isclose(r,row.standard_Pearson_r))
        check(level+unit+' independent zero-target percentile',np.isclose(np.percentile(zero.predicted,95)*100,row.all_zero_P95_percent))
        f=R/'POOLED_INTERNAL_TRAINING'/level
        model=np.load(f/f'{unit}_model.npz');queries=np.load(f/'actual_held_sample_mixtures.npz')[unit]
        meta=pd.read_csv(f/'held_sample_mixture_metadata.csv');meta=meta[meta.unit.eq(unit)].reset_index(drop=True)
        check(level+unit+' query truth tied to original saved mixture labels',np.array_equal(d.set_index('mixture_id').loc[meta.mixture_id,'truth'],meta.truth))
        if not args.saved_only:
            source=np.load(S/'all_mixture_expression.npz')['values']
            expected=source[meta.mixture_id.to_numpy(dtype=int)][:,common.get_indexer(model['genes'])]
            check(level+unit+' actual original query expression unchanged',np.array_equal(expected,queries))
        pg=int(np.flatnonzero(model['genes']=='PGAM5')[0])
        records=pd.read_csv(f/f'{unit}_fit_cases.csv')
        for record in records.itertuples():
            z=np.load(f/f'fit_case_{record.case}.npz')
            check(record.case+' real stored query',np.allclose(z['b'],queries[record.query_row]/model['scale']))
            actual=d[d.mixture_id.eq(record.mixture_id)].iloc[0]
            check(record.case+' saved coefficient matches actual outcome',np.isclose(z['c'][model['target_mask']].sum(),actual.predicted))
            check(record.case+' independent final convex-QP KKT',kkt(model['A'],z['b'],z['c'],z['w'],pg,model['target_mask']))
    with tempfile.TemporaryDirectory(dir=R) as temp:
        temp=Path(temp);f=R/'POOLED_INTERNAL_TRAINING/continuous_rank'
        for unit in units:
            model=np.load(f/f'{unit}_model.npz');query=np.load(f/'actual_held_sample_mixtures.npz')[unit][:2]
            pd.DataFrame(query.T,index=model['genes'],columns=['actual_mixture_0','actual_mixture_1']).to_csv(temp/f'{unit}.tsv',sep='\t')
            call=subprocess.run([sys.executable,str(R/'fit_reference.py'),'--reference-dir',str(f),'--unit',unit,
                '--expression-cp10k',str(temp/f'{unit}.tsv'),'--output-csv',str(temp/f'{unit}.csv')],capture_output=True,text=True)
            check(unit+' portable CLI executed',call.returncode==0)
            result=pd.read_csv(temp/f'{unit}.csv')
            original=predictions[predictions.threshold_version.eq('continuous_rank')&predictions.unit.eq(unit)].iloc[:2]
            check(unit+' actual portable predictions reproduce',np.allclose(result.RNA_detected_macrophage_coefficient,original.predicted,atol=1e-8))
            check(unit+' portable interpretation remains exploratory',not result.usable_as_validated_TCGA_cell_abundance.any())
    result={'status':'PASS','checks':len(checks),'check_names':checks,'scope':'saved actual inputs and independent arithmetic' if args.saved_only else 'original source mixtures, raw pooled contrasts and independent arithmetic',
        'KKT_scope':'Last fixed-weight convex subproblem only; not global nonconvex IRLS optimum or biological validation',
        'protein_and_cross_cohort_validation_required':False,'scientific_TCGA_cell_abundance_validation':'FAILED'}
    (R/('pooled_saved_audit.json' if args.saved_only else 'pooled_source_audit.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='check_names'},indent=2))


if __name__=='__main__':main()
