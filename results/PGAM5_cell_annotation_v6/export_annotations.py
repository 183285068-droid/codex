"""Separate observed RNA positivity from rejected inferred programs in every export."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import sparse
from discover_programs import project,state_masks,association

R = Path(__file__).resolve().parent
P = json.loads((R/'protocol.json').read_text())
S = Path(P['source_root'])


def annotate(obs,program):
    result = obs.copy()
    macrophage = result.harmonized_group.eq('Macrophage')
    result['cell_type'] = result.harmonized_group
    result['PGAM5_RNA_status'] = np.where(result.PGAM5_counts.gt(0),'RNA_detected','RNA_undetected')
    result['PGAM5_macrophage_annotation'] = np.where(macrophage,
        np.where(result.PGAM5_counts.gt(0),'Macrophage_PGAM5_RNA_detected','Macrophage_PGAM5_RNA_undetected'),'Not_macrophage')
    fields = ['dominant_program','mapping_margin','candidate_program_score','candidate_state_member']
    if 'bootstrap_membership_agreement' in program: fields += ['bootstrap_membership_agreement']
    source = program.set_index('cell_id')
    assert source.index.is_unique
    for field in fields:
        result[field] = result.cell_id.map(source[field])
    result['PGAM5_related_state_annotation'] = 'Not_applicable'
    result.loc[macrophage,'PGAM5_related_state_annotation'] = 'Other_or_no_candidate_program'
    result.loc[macrophage & result.candidate_state_member.eq(True),'PGAM5_related_state_annotation'] = 'Exploratory_REJECTED_PGAM5_program_candidate'
    result['annotation_level'] = np.where(macrophage,'Observed_RNA_status_plus_exploratory_program','Existing_lineage_identity')
    result['usable_for_validated_TCGA_PGAM5_cell_abundance'] = False
    assert result.cell_id.is_unique
    assert result.loc[macrophage,'candidate_state_member'].notna().all()
    return result


def main():
    definition = json.loads((R/'ALL_TRAINING/definition.json').read_text())
    assert not definition['training_selection_eligible'], 'Review wording: a supported training state would still need external validation'
    summaries,macrophages = [],[]
    for cohort in P['primary_cohorts']:
        obs = pd.read_csv(S/f'{cohort}_metadata.csv.gz')
        program = pd.read_csv(R/f'{cohort}_macrophage_annotations.csv.gz')
        result = annotate(obs,program)
        result.to_csv(R/f'{cohort}_all_cell_annotation.csv.gz',index=False)
        mac = result[result.cell_type.eq('Macrophage')].copy()
        macrophages.append(mac)
        summaries.append({'dataset':cohort,'scope':'primary_HCC','QC_cells':len(result),'macrophages':len(mac),
            'PGAM5_RNA_detected_macrophages':int(mac.PGAM5_counts.gt(0).sum()),
            'candidate_program_cells':int(mac.candidate_state_member.sum()),
            'PGAM5_RNA_detected_in_candidate':int((mac.candidate_state_member&mac.PGAM5_counts.gt(0)).sum()),
            'annotation_usable_for_validated_TCGA_abundance':False})
    h = pd.read_csv(R/'ALL_TRAINING/program_loadings.tsv',sep='\t',index_col=0).to_numpy().T
    features = pd.read_csv(R/'ALL_TRAINING/features.csv')
    cohort_rows,donor_rows = [],[]
    for cohort in ['GSE125449','GSE140228_Droplet_HCC','GSE140228_Smartseq2_HCC','GSE146115']:
        obs = pd.read_csv(S/f'{cohort}_metadata.csv.gz')
        measured = pd.Index(pd.read_csv(S/f'{cohort}_measured_genes.csv').gene)
        shared = features.gene.isin(measured).to_numpy()
        coverage = shared.mean()
        assert coverage>=.95, 'Insufficient program features for supplementary projection'
        raw = sparse.load_npz(S/f'{cohort}_raw_measured.npz')
        cycle_markers = [g for g in ['MKI67','TOP2A','UBE2C','CENPF','CDK1','BIRC5'] if g in measured]
        obs['cycling_flag'] = np.asarray((raw[:,measured.get_indexer(cycle_markers)]>0).sum(1)).ravel()>=2
        ix = measured.get_indexer(features.loc[shared,'gene'])
        x = raw[:,ix].toarray().astype(float)/obs.total_counts.to_numpy()[:,None]*10000
        x = np.minimum(np.log1p(x)/features.loc[shared,'SD'].to_numpy(),10)
        w = project(x,h[:,shared])
        masks,winner,margin = state_masks(w,np.asarray(definition['activation_thresholds']))
        obs['dataset'] = cohort
        a,d = association(obs,masks,'ALL_TRAINING','supplementary_projection')
        cohort_rows.append(a)
        donor_rows.append(d)
        obs['dominant_program'],obs['mapping_margin'] = winner,margin
        obs['candidate_program_score'] = w[:,definition['target_program']]
        obs['candidate_state_member'] = masks[:,definition['target_program']]
        obs['program_feature_coverage'] = coverage
        result = annotate(obs,obs)
        result.to_csv(R/f'{cohort}_macrophage_annotation.csv.gz',index=False)
        summaries.append({'dataset':cohort,'scope':'supplementary_not_independent_primary_vote','QC_cells':len(result),
            'macrophages':len(result),'PGAM5_RNA_detected_macrophages':int(result.PGAM5_counts.gt(0).sum()),
            'candidate_program_cells':int(result.candidate_state_member.sum()),
            'PGAM5_RNA_detected_in_candidate':int((result.candidate_state_member&result.PGAM5_counts.gt(0)).sum()),
            'program_feature_coverage':coverage,'annotation_usable_for_validated_TCGA_abundance':False})
    pd.concat(cohort_rows,ignore_index=True).to_csv(R/'supplementary_program_association_by_cohort.csv',index=False)
    pd.concat(donor_rows,ignore_index=True).to_csv(R/'supplementary_program_association_by_donor.csv',index=False)
    pd.concat(macrophages,ignore_index=True).to_csv(R/'ALL_PRIMARY_macrophage_annotation.csv.gz',index=False)
    pd.DataFrame(summaries).to_csv(R/'annotation_coverage.csv',index=False)
    dictionary = {'cell_id':'Exact existing cell identifier; join by ID, not row order',
        'cell_type':'Existing conservative lineage identity',
        'PGAM5_counts':'Symbol-summed observed raw RNA counts',
        'PGAM5_RNA_status':'RNA_detected if raw>0; zero does not mean protein negative',
        'PGAM5_macrophage_annotation':'Literal macrophage RNA status, separate from program',
        'dominant_program':'Algorithmic NMF index, no biological subtype name',
        'candidate_state_member':'Fixed dominant/activation/margin gate of rejected exploratory candidate',
        'PGAM5_related_state_annotation':'Rejected exploratory program label; cannot be called confirmed PGAM5-specific state',
        'bootstrap_membership_agreement':'Agreement in 2 conditional within-donor bootstrap fits; not a posterior probability',
        'usable_for_validated_TCGA_PGAM5_cell_abundance':'False: RNA state/program and bulk reference not validated for cellular abundance'}
    (R/'annotation_data_dictionary.json').write_text(json.dumps(dictionary,indent=2)+'\n')
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__=='__main__':
    main()
