"""Pooled five-cohort RNA reference; no cross-cohort replication requirement."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import argparse
import json
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import mannwhitneyu,pearsonr
from statsmodels.stats.multitest import multipletests
from guard_solver import fit

R=Path(__file__).resolve().parent
P=json.loads((R/'protocol.json').read_text())
S=Path(P['source_root'])
C=P['primary_cohorts']
G=pd.Index(pd.read_csv(S/'common_genes.csv').gene)
CYCLE=set((S/'cell_cycle_genes.txt').read_text().split())
PANELS=pd.read_csv(R/'gene_roles.csv')
STRESS=set(PANELS[PANELS.role.eq('stress_diagnostic')].gene)


def de(norm,pg,folder):
    positive=pg>0; negative=~positive
    mean1=np.asarray(norm[positive].mean(0)).ravel(); mean0=np.asarray(norm[negative].mean(0)).ravel()
    freq1=np.asarray((norm[positive]>0).mean(0)).ravel()
    freq0=np.asarray((norm[negative]>0).mean(0)).ravel()
    p=np.ones(len(G))
    for start in range(0,len(G),128):
        end=min(len(G),start+128)
        x=np.log1p(norm[positive,start:end].toarray())
        y=np.log1p(norm[negative,start:end].toarray())
        p[start:end]=mannwhitneyu(x,y,axis=0,alternative='two-sided',method='asymptotic',use_continuity=True).pvalue
        if start%2048==0: print('Pooled DE',folder.name,start,'of',len(G),flush=True)
    p[(mean1==0)&(mean0==0)]=1
    assert np.isfinite(p).all()
    q=multipletests(p,method='fdr_bh')[1]
    fc=np.log2((mean1+.001)/(mean0+.001))
    result=pd.DataFrame({'gene':G,'RNA_detected_cells':int(positive.sum()),'RNA_undetected_cells':int(negative.sum()),
        'mean_CP10k_detected':mean1,'mean_CP10k_undetected':mean0,
        'detected_group_detection':freq1,'undetected_group_detection':freq0,
        'Mann_Whitney_p':p,'BH_FDR':q,'log2FC':fc,'strict_upregulated':(q<.05)&(fc>=1)&(freq1>=.1)})
    result.to_csv(folder/'pooled_DE.csv.gz',index=False)
    return result


def reference(index,stats,de_table,folder,unit,minimum_log2FC):
    support=index.groupby('group').agg(cells=('cells','sum'),donor_labels=('donor','nunique'))
    active=set(support[(support.cells>=50)&(support.donor_labels>=2)].index)
    assert 'TAM_PGAM5_detected' in active and 'TAM_PGAM5_undetected' in active
    def mapped(g):
        if g in active:return g
        if g.startswith('Cycling_') and g[8:] in active:return g[8:]
        return 'Unknown_other'
    ind=index.copy();ind['mapped_group']=ind.group.map(mapped)
    groups=sorted(ind.mapped_group.unique())
    refs=[];freqs=[];support_rows=[]
    for group in groups:
        cohort_profiles=[];cohort_freq=[];donors=0
        for cohort,part in ind[ind.mapped_group.eq(group)].groupby('dataset',sort=True):
            donorprofiles=[]
            for donor,d in part.groupby('donor',sort=True):
                v=stats[d.index].astype(float).sum(0)
                profile=v[0]/d.cells.sum() if unit=='equalized_cell_fraction' else v[1]/d.total_counts.sum()*10000
                donorprofiles.append(profile);donors+=1
            cohort_profiles.append(np.stack(donorprofiles).mean(0))
            cohort_freq.append(stats[part.index,2].astype(float).sum(0)/part.cells.sum())
        refs.append(np.stack(cohort_profiles).mean(0));freqs.append(np.stack(cohort_freq).mean(0))
        support_rows.append({'group':group,'cells':int(ind[ind.mapped_group.eq(group)].cells.sum()),
            'donor_or_sample_labels':donors,'cohorts':len(cohort_profiles)})
    full=np.column_stack(refs);frequency=np.column_stack(freqs)
    targetmask=np.array(['PGAM5_detected' in x for x in groups])
    macmask=np.array([x.startswith('TAM_') or x.startswith('Cycling_TAM_') for x in groups])
    # Aggregate target profiles within each donor/cohort including cycling target cells.
    targetind=ind[ind.group.str.contains('PGAM5_detected',regex=False)]
    cohorts=[]
    for _,part in targetind.groupby('dataset',sort=True):
        ds=[]
        for _,d in part.groupby('donor',sort=True):
            v=stats[d.index].astype(float).sum(0)
            ds.append(v[0]/d.cells.sum() if unit=='equalized_cell_fraction' else v[1]/d.total_counts.sum()*10000)
        cohorts.append(np.stack(ds).mean(0))
    targetmean=np.stack(cohorts).mean(0)
    maxother=full[:,~macmask].max(1)
    specificity=np.log2((targetmean+.05)/(maxother+.05))
    table=de_table.copy()
    table['cycle_list_member']=[x in CYCLE for x in G]
    table['stress_list_member']=[x in STRESS for x in G]
    table['balanced_target_mean_CP10k']=targetmean
    table['maximum_nonmacrophage_mean_CP10k']=maxother
    table['competitive_log2_ratio']=specificity
    table['strongest_nonmacrophage_competitor']=np.array(groups)[np.flatnonzero(~macmask)[full[:,~macmask].argmax(1)]]
    background=np.array([x=='PGAM5' or x.startswith(('MT-','RPL','RPS')) for x in G])
    score=np.maximum(table.log2FC.to_numpy(),0)*np.maximum(specificity,0)*table.detected_group_detection.to_numpy()
    if minimum_log2FC is None:
        table['within_macrophage_DE_eligible']=table.detected_group_detection>=.1
        eligible=table.within_macrophage_DE_eligible.to_numpy() & ~background & ~table.cycle_list_member.to_numpy() & ~table.stress_list_member.to_numpy() & (score>0)
    else:
        table['within_macrophage_DE_eligible']=(table.BH_FDR<.05)&(table.log2FC>=minimum_log2FC)&(table.detected_group_detection>=.1)
        eligible=table.within_macrophage_DE_eligible.to_numpy() & ~background & ~table.cycle_list_member.to_numpy() & ~table.stress_list_member.to_numpy() & (specificity>=1)
    table['continuous_ranking_score']=score
    states=np.flatnonzero(eligible); states=states[np.argsort(-score[states],kind='stable')[:50]]
    table['selected_candidate_state_feature']=False;table.loc[states,'selected_candidate_state_feature']=True
    table.to_csv(folder/f'{unit}_gene_evidence.csv.gz',index=False)
    table.loc[states].to_csv(folder/f'{unit}_candidate_state_genes.csv',index=False)
    panelgenes=set(PANELS.gene)|{'PGAM5'}
    selected=set(np.flatnonzero(G.isin(panelgenes)))|set(states)
    broadrows=[]
    for j,group in enumerate(groups):
        competitor=np.delete(full,j,axis=1).max(1)
        contrast=np.log2((full[:,j]+.05)/(competitor+.05))
        rank=contrast*np.minimum(full[:,j],1)
        choices=np.flatnonzero((frequency[:,j]>=.1)&(full[:,j]>=.05)&~np.array(G=='PGAM5'))
        choices=choices[np.argsort(-rank[choices],kind='stable')[:15]]
        selected.update(choices)
        broadrows.extend({'gene':G[i],'group':group,'reference_discrimination_log2_ratio':contrast[i]} for i in choices)
    ix=np.array(sorted(selected));scale=np.maximum(np.sqrt(np.mean(full[ix]**2,axis=1)),.05)
    A=full[ix]/scale[:,None]
    pd.DataFrame(full[ix],index=G[ix],columns=groups).to_csv(folder/f'EXPLORATORY_reference_{unit}.tsv',sep='\t',index_label='gene')
    pd.DataFrame({'gene':G[ix],'scale':scale}).to_csv(folder/f'{unit}_row_scales.csv',index=False)
    pd.DataFrame(support_rows).to_csv(folder/f'{unit}_group_support.csv',index=False)
    pd.DataFrame(broadrows).to_csv(folder/f'{unit}_broad_gene_roles.csv',index=False)
    np.savez_compressed(folder/f'{unit}_model.npz',A=A,scale=scale,genes=np.asarray(G[ix],dtype=str),groups=np.array(groups),
        target_mask=targetmask,macrophage_mask=macmask)
    definition={'unit':unit,'reference_genes':len(ix),'reference_components':len(groups),
        'within_macrophage_log2FC_min':minimum_log2FC,'FDR_strict_max':.05 if minimum_log2FC is not None else None,
        'candidate_state_genes':G[states].tolist(),'candidate_state_genes_count':len(states),
        'PGAM5_anchor_included':True,'all_measured_lineage_panel_genes_included':True,
        'RNA_target':'Existing macrophage and PGAM5 raw count>0; cycling target retained',
        'cross_cohort_replication_required':False,'protein_validation_required':False,
        'candidate_gene_list_is_validated_signature':False,'usable_as_validated_TCGA_cell_abundance':False}
    (folder/f'{unit}_definition.json').write_text(json.dumps(definition,indent=2)+'\n')
    return {'genes':len(ix),'groups':len(groups),'state_genes':len(states)}


def evaluate(task):
    level,unit=task
    folder=R/'POOLED_INTERNAL_TRAINING'/level
    model=np.load(folder/f'{unit}_model.npz');A=model['A'];scale=model['scale']
    query=np.load(folder/'actual_held_sample_mixtures.npz')[unit]
    meta=pd.read_csv(folder/'held_sample_mixture_metadata.csv')
    meta=meta[meta.unit.eq(unit)].reset_index(drop=True)
    assert len(meta)==len(query)
    pg=int(np.flatnonzero(model['genes']=='PGAM5')[0])
    target=model['target_mask'];mac=model['macrophage_mask']
    rows=[];case_rows=[]
    for j,row in enumerate(meta.itertuples()):
        b=query[j].astype(float)/scale
        c,w,res=fit(A,b,True,pg,target)
        rows.append({'threshold_version':level,'dataset':row.dataset,'donor_or_sample':row.donor,'unit':unit,'mixture_id':row.mixture_id,
            'scenario':row.scenario,'replicate':row.replicate,'nominal_target_cell_fraction':row.nominal_target_cell_fraction,
            'truth':row.truth,'predicted':float(c[target].sum()),'truth_macrophage':row.macrophage_cell_fraction if unit=='equalized_cell_fraction' else row.macrophage_RNA_contribution,
            'predicted_macrophage':float(c[mac].sum()),'relative_residual':res,
            'observed_PGAM5_scaled':b[pg],'target_PGAM5_scaled':float((A[pg]*target)@c),
            'usable_as_validated_cell_abundance':False})
        if (row.scenario=='standard' and row.replicate==0 and row.nominal_target_cell_fraction in [0,.05] and len(case_rows)<8) or j==len(meta)-1:
            key=f'{unit}_{j}'
            np.savez_compressed(folder/f'fit_case_{key}.npz',b=b,c=c,w=w,query_row=j,unit=unit,pg_index=pg)
            case_rows.append({'case':key,'mixture_id':row.mixture_id,'unit':unit,'query_row':j})
        if j%100==0: print('Internal holdout',unit,j,'of',len(meta),flush=True)
    pd.DataFrame(case_rows).to_csv(folder/f'{unit}_fit_cases.csv',index=False)
    return pd.DataFrame(rows)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--reuse-existing-DE',action='store_true',help='Reuse unchanged numerical DE after a serialization-only restart in the same fixed-input analysis')
    args=parser.parse_args()
    rng=np.random.default_rng(P['seed']);splitrows=[];data=[];obsall=[];statparts=[];indexparts=[]
    for cohort in C:
        obs=pd.read_csv(S/f'{cohort}_metadata.csv.gz');labels=sorted(obs.donor.unique())
        held=set(rng.choice(labels,max(1,round(.3*len(labels))),replace=False))
        for label in labels:
            splitrows.append({'dataset':cohort,'donor_or_sample':label,'split':'internal_validation' if label in held else 'training',
                'independent_patient_confirmed':cohort!='GSE202642','fresh_untouched_external_patient':False})
        mac=obs.harmonized_group.eq('Macrophage').to_numpy()
        raw=sparse.load_npz(S/f'{cohort}_raw_common.npz')[mac].astype(float)
        mo=obs[mac].reset_index(drop=True);mo['split']=np.where(mo.donor.isin(held),'internal_validation','training')
        data.append(raw.multiply(10000/mo.total_counts.to_numpy()[:,None]).tocsr());obsall.append(mo)
        ind=pd.read_csv(S/f'{cohort}_group_stats_index.csv');ind['split']=np.where(ind.donor.isin(held),'internal_validation','training')
        indexparts.append(ind);statparts.append(np.load(S/f'{cohort}_group_statistics.npz')['values'])
    pd.DataFrame(splitrows).to_csv(R/'pooled_patient_sample_split.csv',index=False)
    norm=sparse.vstack(data,format='csr');obs=pd.concat(obsall,ignore_index=True)
    obs[['dataset','cell_id','donor','PGAM5_counts','total_counts','split']].to_csv(R/'pooled_macrophage_DE_members.csv.gz',index=False)
    index=pd.concat(indexparts,ignore_index=True);stats=np.concatenate(statparts)
    units=['equalized_cell_fraction','library_RNA_contribution']
    dimensions=[]
    for name,mask in [('POOLED_INTERNAL_TRAINING',obs.split.eq('training').to_numpy()),('POOLED_ALL_TRAINING',np.ones(len(obs),dtype=bool))]:
        folder=R/name;folder.mkdir(exist_ok=True)
        if args.reuse_existing_DE and (folder/'pooled_DE.csv.gz').exists():
            table=pd.read_csv(folder/'pooled_DE.csv.gz')
            assert np.array_equal(table.gene,G)
            assert table.RNA_detected_cells.eq(int((obs.PGAM5_counts.to_numpy()[mask]>0).sum())).all()
            assert table.RNA_undetected_cells.eq(int((obs.PGAM5_counts.to_numpy()[mask]==0).sum())).all()
            print('Reusing unchanged pooled numerical DE after NPZ text-dtype correction:',name,flush=True)
        else:
            table=de(norm[mask],obs.PGAM5_counts.to_numpy()[mask],folder)
        ix=np.flatnonzero(index.split.eq('training').to_numpy()) if name=='POOLED_INTERNAL_TRAINING' else np.arange(len(index))
        ind=index.iloc[ix].reset_index(drop=True);v=stats[ix]
        for level,minfc in [('log2FC_1',1.),('log2FC_0p5',.5),('continuous_rank',None)]:
            output=folder/level;output.mkdir(exist_ok=True)
            for unit in units:
                info=reference(ind,v,table,output,unit,minfc)
                dimensions.append({'model':name,'threshold_version':level,'unit':unit,'macrophage_DE_cells':int(mask.sum()),
                    'PGAM5_RNA_detected_DE_cells':int(obs.PGAM5_counts.to_numpy()[mask].astype(bool).sum()),**info})
        if name=='POOLED_INTERNAL_TRAINING':
            held=pd.DataFrame(splitrows);held=set(held.loc[held.split.eq('internal_validation'),'donor_or_sample'])
            meta=pd.read_csv(S/'mixture_index.csv');meta=meta[meta.donor.isin(held)].copy()
            expression=np.load(S/'all_mixture_expression.npz')['values']
            for level in ['log2FC_1','log2FC_0p5','continuous_rank']:
                output=folder/level
                meta.to_csv(output/'held_sample_mixture_metadata.csv',index=False)
                queries={}
                for unit in units:
                    genes=np.load(output/f'{unit}_model.npz')['genes']
                    rows=meta.loc[meta.unit.eq(unit),'mixture_id'].to_numpy(dtype=int)
                    queries[unit]=expression[rows][:,G.get_indexer(genes)]
                np.savez_compressed(output/'actual_held_sample_mixtures.npz',**queries)
    pd.DataFrame(dimensions).to_csv(R/'pooled_model_dimensions.csv',index=False)
    tasks=[(level,unit) for level in ['log2FC_1','log2FC_0p5','continuous_rank'] for unit in units]
    with ProcessPoolExecutor(max_workers=2) as pool:predictions=pd.concat(list(pool.map(evaluate,tasks)),ignore_index=True)
    predictions.to_csv(R/'pooled_internal_mixture_predictions.csv.gz',index=False)
    rows=[];gates=[];doses=[]
    for (level,unit),d in predictions.groupby(['threshold_version','unit'],sort=True):
        standard=d[d.scenario.eq('standard')];zero=d[d.truth.eq(0)]
        mae=float(np.abs(standard.truth-standard.predicted).mean()*100)
        correlation=float(pearsonr(standard.truth,standard.predicted).statistic) if standard.predicted.nunique()>1 else np.nan
        p95=float(zero.predicted.quantile(.95)*100)
        patients=standard.loc[~standard.dataset.eq('GSE202642'),'donor_or_sample'].nunique()
        rows.append({'threshold_version':level,'unit':unit,'mixture_rows':len(d),'standard_rows':len(standard),
            'standard_patient_or_sample_labels':standard.donor_or_sample.nunique(),
            'standard_confirmed_patient_labels_excluding_GSE202642':patients,
            'standard_MAE_pp':mae,'standard_Pearson_r':correlation,'all_zero_P95_percent':p95,
            'all_zero_max_percent':zero.predicted.max()*100,'physical_mixtures':False,
            'fresh_external_validation':False,'cross_cohort_replication_required':False,
            'usable_as_validated_TCGA_cell_abundance':False})
        for name,result in [('MAE_le1pp',mae<=1),('r_ge0.7',correlation>=.7),('zero_P95_le0.5percent',p95<=.5),('internal_confirmed_patients_ge3',patients>=3)]:
            gates.append({'threshold_version':level,'unit':unit,'gate':name,'passed':bool(result)})
        for dose,part in standard.groupby('nominal_target_cell_fraction'):
            doses.append({'threshold_version':level,'unit':unit,'nominal_cell_fraction':dose,'mixtures':len(part),
                'patient_or_sample_labels':part.donor_or_sample.nunique(),'median_truth':part.truth.median(),
                'median_prediction':part.predicted.median(),'prediction_q25':part.predicted.quantile(.25),
                'prediction_q75':part.predicted.quantile(.75),'fraction_estimate_above0.5percent':part.predicted.gt(.005).mean()})
    pd.DataFrame(rows).to_csv(R/'pooled_internal_validation_summary.csv',index=False)
    pd.DataFrame(gates).to_csv(R/'pooled_internal_validation_gates.csv',index=False)
    pd.DataFrame(doses).to_csv(R/'pooled_internal_dose_recovery.csv',index=False)
    state=json.loads((R/'validation_result.json').read_text())
    state.update(cross_cohort_replication_required=False,independent_RNA_identity_required_for_operational_label=False,
        protein_validation_required=False,pooled_internal_mixture_gates_passed=all(x['passed'] for x in gates if x['threshold_version']=='continuous_rank'),
        interpretation='Pooled RNA-detected macrophage reference; no cross-cohort requirement; computational internal validation and RNA-vs-cell/bulk-platform limits separately reported',
        pooled_state_gene_counts={x['threshold_version']+'_'+x['unit']:x['state_genes'] for x in dimensions if x['model']=='POOLED_ALL_TRAINING'})
    (R/'validation_result.json').write_text(json.dumps(state,indent=2)+'\n')
    print(pd.DataFrame(dimensions).to_string(index=False),flush=True)
    print(pd.DataFrame(rows).to_string(index=False),flush=True)


if __name__=='__main__':main()
