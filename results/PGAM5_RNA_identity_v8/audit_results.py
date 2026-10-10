"""Independent arithmetic, source identity, import, and RNA-validator checks."""
from pathlib import Path
import argparse
import json
import hashlib
import tempfile
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import norm,rankdata
import anndata as ad
from attach_annotations import attach,FIELDS
from validate_RNA_inputs import validate,SAMPLE_FIELDS,IDENTITY_FIELDS,MIXTURE_FIELDS

R=Path(__file__).resolve().parent
P=json.loads((R/'protocol.json').read_text())
S=Path(P['source_root'])
checks=[]


def check(name,value):
    assert bool(value),name
    checks.append(name)


def manual_two_sided_mwu(x,y):
    n1,n0=len(x),len(y)
    together=np.concatenate([x,y]); N=n1+n0
    ranks=rankdata(together,method='average')
    U=ranks[:n1].sum()-n1*(n1+1)/2
    tie=np.unique(together,return_counts=True)[1].astype(float)
    sd=np.sqrt(n1*n0/12*((N+1)-np.sum(tie**3-tie)/(N*(N-1))))
    return min(1.,float(2*norm.sf((abs(U-n1*n0/2)-.5)/sd))) if sd>0 else 1.


def saved_audit():
    cover=pd.read_csv(R/'annotation_coverage.csv').set_index('dataset')
    for c in P['primary_cohorts']+P['supplementary']:
        a=pd.read_csv(R/f'{c}_RNA_annotation.csv.gz')
        mac=a.harmonized_group.eq('Macrophage'); det=mac & a.PGAM5_counts.gt(0)
        check(c+' IDs unique',a.cell_id.is_unique)
        check(c+' zero remains undetected',a.loc[a.PGAM5_counts.eq(0),'PGAM5_RNA_status'].eq('RNA_undetected').all())
        check(c+' positive remains detected',a.loc[a.PGAM5_counts.gt(0),'PGAM5_RNA_status'].eq('RNA_detected').all())
        check(c+' target follows lineage and raw RNA',np.array_equal(a.PGAM5_macrophage_annotation.eq('Macrophage_PGAM5_RNA_detected'),det))
        check(c+' counts reproduced',len(a)==cover.loc[c,'cells'] and mac.sum()==cover.loc[c,'macrophages'] and det.sum()==cover.loc[c,'PGAM5_RNA_detected_macrophages'])
        check(c+' corroboration summary',int((det&a.macrophage_RNA_marker_corroboration).sum())==cover.loc[c,'marker_corroborated_PGAM5_RNA_detected_macrophages'])
        check(c+' state unresolved',a.loc[mac,'PGAM5_related_state'].eq('Unresolved_no_validated_RNA_state').all())
        check(c+' validation never promoted',not a.usable_for_validated_TCGA_PGAM5_cell_abundance.any() and not a.RNA_state_validation_passed.any())
    primary=cover.loc[P['primary_cohorts']]
    check('primary200068 cells',primary.cells.sum()==200068)
    check('primary30533 macrophages',primary.macrophages.sum()==30533)
    check('primary1413 detected',primary.PGAM5_RNA_detected_macrophages.sum()==1413)
    screen=pd.read_csv(R/'RNA_state_screen_by_fold.csv')
    for _,row in screen.iterrows():
        t=pd.read_csv(R/row.held_cohort/'training_gene_evidence.csv.gz')
        train=row.training_cohorts.split(',')
        check(row.held_cohort+' training boundary',row.held_cohort not in train or row.held_cohort=='ALL_TRAINING')
        votes=t[[x+'_strict_up' for x in train]].to_numpy().sum(1)
        check(row.held_cohort+' votes independently counted',np.array_equal(votes,t.strict_cohort_votes))
        specificity=np.log2((t.training_target_mean_CP10k+.05)/(t.training_max_competitor_mean_CP10k+.05))
        check(row.held_cohort+' ratios independently calculated',np.allclose(specificity,t.training_competitive_log2_ratio))
        base=~t.gene.str.startswith(('MT-','RPL','RPS')) & t.gene.ne('PGAM5')
        recurrent=base & (votes>=row.required_cohort_votes) & t.donor_positive_direction_fraction.ge(.7)
        eligible=recurrent & ~t.cycle_list_member & ~t.stress_list_member & specificity.ge(1)
        check(row.held_cohort+' eligibility independently recomputed',np.array_equal(eligible,t.RNA_state_feature_eligible))
        check(row.held_cohort+' insufficient panel gate',eligible.sum()==row.competitive_specific_candidates and not row.feature_count_gate_passed)
    x=pd.read_csv(R/'actual_cell_marker_counts.csv').to_numpy()
    genes=pd.Index(pd.read_csv(R/'actual_cell_measured_genes.csv').gene)
    a=pd.read_csv(R/'actual_cell_annotations.csv')
    check('actual60 PGAM5 counts',np.array_equal(x[:,genes.get_loc('PGAM5')],a.PGAM5_counts))
    core=['CD68','CD163','CSF1R','C1QA','C1QB','C1QC','MERTK','MSR1']
    broad=['LST1','TYROBP','FCER1G','AIF1','SPI1','CTSS']
    check('actual60 core detection',(x[:,genes.get_indexer(core)]>0).sum(1).tolist()==a.macrophage_core_detected_genes.tolist())
    check('actual60 broad detection',(x[:,genes.get_indexer(broad)]>0).sum(1).tolist()==a.broad_myeloid_detected_genes.tolist())
    with tempfile.TemporaryDirectory(dir=R) as tmp:
        tmp=Path(tmp); source=R/'import_example/actual_60_HCC_cells.h5ad'
        original=ad.read_h5ad(source)
        attached=attach(source,R/'GSE149614_RNA_annotation.csv.gz',tmp/'attached.h5ad','cell_id')
        check('import original X preserved',(original.X!=attached.X).nnz==0)
        check('import object row order preserved',np.array_equal(original.obs_names,attached.obs_names))
        truth=a.set_index('cell_id').loc[original.obs.cell_id]
        for field in FIELDS:
            left=attached.obs['pgam5_v8_'+field].to_numpy();right=truth[field].to_numpy()
            equal=np.allclose(left,right,equal_nan=True) if pd.api.types.is_numeric_dtype(truth[field].dtype) else np.array_equal(left,right)
            check('import field '+field,equal)
        bad=original.copy(); bad.obs.loc[bad.obs_names[0],'cell_id']='unmatched_cell'
        bad.write_h5ad(tmp/'unmatched.h5ad')
        failed=False
        try: attach(tmp/'unmatched.h5ad',R/'GSE149614_RNA_annotation.csv.gz',tmp/'should_not_exist.h5ad','cell_id')
        except AssertionError: failed=True
        check('unmatched stops without output',failed and not (tmp/'should_not_exist.h5ad').exists())
        replay=R/'computational_replay'
        samples=pd.read_csv(replay/'samples.tsv',sep='\t')
        identity=pd.read_csv(replay/'RNA_identity.tsv',sep='\t')
        mix=pd.read_csv(replay/'mixtures.tsv',sep='\t')
        result=validate(samples,identity,mix,tmp/'replay',computational_replay=True)
        check('actual replay never becomes real validation',not result['mixture_gates_passed'] and result['status']=='COMPUTATIONAL_REPLAY_NOT_INDEPENDENT')
        check('no invented independent RNA identity',result['eligible_independent_RNA_identity_patients']==0)
        check('historical samples not called new validation',result['eligible_validation_patients']==0)
        zeros=mix.copy(); zeros['estimated_target_fraction']=0
        result=validate(samples,identity,zeros,tmp/'zeros',computational_replay=True)
        check('all-zero estimator fails recovery',not result['numerical_mixture_gates_passed'])
        corrupt=mix.copy(); corrupt.loc[0,'total_cells']+=1
        failed=False
        try: validate(samples,identity,corrupt,tmp/'badcounts',computational_replay=True)
        except AssertionError: failed=True
        check('inconsistent mixture cell total rejected',failed)
        leaked=samples.copy(); extra=leaked.iloc[[0]].copy(); extra['sample_id']='extra_sample_same_patient'; extra['split']='discovery'
        failed=False
        try: validate(pd.concat([leaked,extra]),identity,mix,tmp/'leak',computational_replay=True)
        except AssertionError: failed=True
        check('patient split overlap rejected',failed)
        missing=validate(pd.DataFrame(columns=SAMPLE_FIELDS),pd.DataFrame(columns=IDENTITY_FIELDS),pd.DataFrame(columns=MIXTURE_FIELDS),tmp/'missing')
        check('empty RNA templates mean NOT_RUN',missing['status']=='NOT_RUN_REQUIRED_RNA_INPUTS_MISSING' and not missing['mixture_gates_passed'])
    status=json.loads((R/'validation_result.json').read_text())
    check('protein not an acceptance criterion',status['protein_validation_required'] is False)
    check('unsupported TCGA claim suppressed',status['usable_for_validated_TCGA_PGAM5_cell_abundance'] is False)


def source_audit():
    genes=pd.Index(pd.read_csv(S/'common_genes.csv').gene)
    selected=['PGAM5','DHFR','MKI67','NUSAP1','TOP2A','CENPK','TYMS']
    profiles={}
    for c in P['primary_cohorts']+P['supplementary']:
        obs=pd.read_csv(S/f'{c}_metadata.csv.gz')
        raw=sparse.load_npz(S/f'{c}_raw_common.npz' if c in P['primary_cohorts'] else S/f'{c}_raw_measured.npz')
        measured=genes if c in P['primary_cohorts'] else pd.Index(pd.read_csv(S/f'{c}_measured_genes.csv').gene)
        a=pd.read_csv(R/f'{c}_RNA_annotation.csv.gz')
        check(c+' source IDs/order preserved',np.array_equal(obs.cell_id,a.cell_id))
        check(c+' source lineage unchanged',np.array_equal(obs.harmonized_group,a.harmonized_group))
        check(c+' source zero/cycling retained',np.array_equal(obs.PGAM5_counts,a.PGAM5_counts) and len(obs)==len(a))
        count=raw[:,measured.get_loc('PGAM5')].toarray().ravel()
        check(c+' actual full-source PGAM5',np.array_equal(count,a.PGAM5_counts))
        core=[g for g in ['CD68','CD163','CSF1R','C1QA','C1QB','C1QC','MERTK','MSR1'] if g in measured]
        check(c+' actual full-source core markers',np.array_equal(np.asarray((raw[:,measured.get_indexer(core)]>0).sum(1)).ravel(),a.macrophage_core_detected_genes))
        if c not in P['primary_cohorts']: continue
        d=pd.read_csv(S/f'{c}_tie_corrected_DE.csv.gz')
        # BH computed independently using sorted p and suffix minima.
        p=d.Mann_Whitney_two_sided_p.to_numpy(); order=np.argsort(p); q=np.empty(len(p))
        q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
        check(c+' full-universe independent BH',np.allclose(q,d.BH_FDR,atol=1e-12,rtol=1e-10))
        mac=obs.harmonized_group.eq('Macrophage').to_numpy(); target=mac & (count>0); other=mac & (count==0)
        small=raw[:,genes.get_indexer(selected)].toarray().astype(float)*10000/obs.total_counts.to_numpy()[:,None]
        for j,gene in enumerate(selected):
            expected=d.set_index('gene').loc[gene]
            pv=manual_two_sided_mwu(np.log1p(small[target,j]),np.log1p(small[other,j]))
            check(c+' manual tie-corrected p '+gene,np.isclose(pv,expected.Mann_Whitney_two_sided_p,atol=1e-10,rtol=1e-6))
            fc=np.log2((small[target,j].mean()+.001)/(small[other,j].mean()+.001))
            check(c+' source fold-change '+gene,np.isclose(fc,expected.log2FC_linear_CP10k,atol=1e-6))
        profiles[c]={}
        broad=obs.fine_group.str.replace(r'^Cycling_','',regex=True)
        broad[broad.str.startswith('TAM_')]='Macrophage'
        for group in sorted(broad.unique()):
            if group=='Macrophage': continue
            part=obs[broad.eq(group)]
            means=[small[d.index].mean(0) for _,d in part.groupby('donor')]
            profiles[c][group]={'values':np.stack(means).mean(0),'cells':len(part),'donors':len(means)}
        part=obs[target]
        means=[small[d.index].mean(0) for _,d in part.groupby('donor')]
        profiles[c]['Target_RNA_detected_macrophage']={'values':np.stack(means).mean(0),'cells':len(part),'donors':len(means)}
        print('Source verified',c,flush=True)
    for held in P['primary_cohorts']+['ALL_TRAINING']:
        training=[c for c in P['primary_cohorts'] if c!=held]
        t=pd.read_csv(R/held/'training_gene_evidence.csv.gz').set_index('gene').loc[selected]
        target=np.stack([profiles[c]['Target_RNA_detected_macrophage']['values'] for c in training]).mean(0)
        groups=sorted(set().union(*[set(profiles[c]) for c in training])-{'Target_RNA_detected_macrophage'})
        competitors=[]
        for group in groups:
            available=[profiles[c][group] for c in training if group in profiles[c]]
            if sum(x['cells'] for x in available)>=50 and sum(x['donors'] for x in available)>=2:
                competitors.append(np.stack([x['values'] for x in available]).mean(0))
        maximum=np.stack(competitors).max(0)
        check(held+' raw-source balanced target profiles',np.allclose(target,t.training_target_mean_CP10k,rtol=1e-5,atol=1e-6))
        check(held+' raw-source competitor maxima',np.allclose(maximum,t.training_max_competitor_mean_CP10k,rtol=1e-5,atol=1e-6))
    provenance=json.loads((R/'source_provenance.json').read_text())
    check('protocol hash',hashlib.sha256((R/'protocol.json').read_bytes()).hexdigest()==provenance['protocol_sha256'])


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--saved-only',action='store_true'); args=parser.parse_args()
    saved_audit()
    if not args.saved_only: source_audit()
    output={'status':'PASS','checks':len(checks),'check_names':checks,'scope':'saved results, actual import and validator' if args.saved_only else 'source/independent arithmetic, actual import and validator',
        'scientific_RNA_state_validation':'FAILED','TCGA_cellular_abundance_validation':'FAILED'}
    (R/('saved_audit.json' if args.saved_only else 'source_audit.json')).write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k!='check_names'},indent=2))


if __name__=='__main__': main()
