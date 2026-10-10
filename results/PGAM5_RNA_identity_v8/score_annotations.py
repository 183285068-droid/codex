"""Add descriptive RNA module scores without changing raw-detection labels."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import sparse

R=Path(__file__).resolve().parent
P=json.loads((R/'protocol.json').read_text());S=Path(P['source_root'])


def main():
    versions=['continuous_rank','log2FC_1','log2FC_0p5']
    definitions={v:json.loads((R/'POOLED_ALL_TRAINING'/v/'equalized_cell_fraction_definition.json').read_text()) for v in versions}
    for c in P['primary_cohorts']+P['supplementary']:
        primary=c in P['primary_cohorts']
        genes=pd.Index(pd.read_csv(S/'common_genes.csv' if primary else S/f'{c}_measured_genes.csv').gene)
        raw=sparse.load_npz(S/f'{c}_raw_common.npz' if primary else S/f'{c}_raw_measured.npz')
        a=pd.read_csv(R/f'{c}_RNA_annotation.csv.gz')
        before=a.PGAM5_macrophage_annotation.copy()
        for v in versions:
            selected=definitions[v]['candidate_state_genes']
            measured=[g for g in selected if g in genes]
            a['candidate_module_'+v+'_feature_coverage']=len(measured)/len(selected) if selected else np.nan
            if measured and len(measured)==len(selected):
                x=raw[:,genes.get_indexer(measured)].toarray().astype(float)*10000/a.total_counts.to_numpy()[:,None]
                a['candidate_module_'+v+'_score']=np.log1p(x).mean(1)
            else:
                a['candidate_module_'+v+'_score']=np.nan
        assert np.array_equal(before,a.PGAM5_macrophage_annotation)
        a.to_csv(R/f'{c}_RNA_annotation.csv.gz',index=False)
        print('Descriptive candidate module scores:',c,flush=True)
    frozen={'RNA_target_definition':'Existing macrophage identity AND symbol-summed PGAM5 raw count>0',
        'RNA_undetected_definition':'Existing macrophage AND raw count=0; never relabel from a module score',
        'macrophage_identity_source':'Existing independent-of-PGAM5 lineage annotations, supplemented by descriptive multigene RNA corroboration',
        'candidate_module':'Mean log1p(CP10k) of selected whole-pool exploratory genes; full feature coverage required. No binary threshold is defined and scores do not establish a stable subtype',
        'primary_reference':'POOLED_ALL_TRAINING/continuous_rank',
        'primary_features_selected_without_FDR_or_log2FC_hard_cutoff':True,
        'protein_validation_required':False,'cross_cohort_replication_required':False,
        'candidate_scores_change_observed_RNA_label':False,'validated_TCGA_cell_abundance':False}
    (R/'frozen_RNA_annotation_definition.json').write_text(json.dumps(frozen,indent=2)+'\n')
    example=pd.read_csv(R/'actual_cell_annotations.csv')
    full=pd.read_csv(R/'GSE149614_RNA_annotation.csv.gz').set_index('cell_id')
    full.loc[example.cell_id].reset_index().to_csv(R/'actual_cell_annotations.csv',index=False)
    actual=pd.read_csv(S/'GSE149614_metadata.csv.gz')
    ids=pd.Index(actual.cell_id).get_indexer(example.cell_id)
    genes=pd.Index(pd.read_csv(S/'common_genes.csv').gene)
    raw=sparse.load_npz(S/'GSE149614_raw_common.npz')
    symbols=set(pd.read_csv(R/'gene_roles.csv').gene)
    for definition in definitions.values():symbols.update(definition['candidate_state_genes'])
    measured=genes[genes.isin(symbols)]
    pd.DataFrame(raw[ids][:,genes.get_indexer(measured)].toarray(),columns=measured).to_csv(R/'actual_cell_marker_counts.csv',index=False)
    pd.DataFrame({'gene':measured}).to_csv(R/'actual_cell_measured_genes.csv',index=False)


if __name__=='__main__':main()
