"""Header-only RNA templates and replay of real saved computational inputs."""
from pathlib import Path
import csv
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from validate_RNA_inputs import SAMPLE_FIELDS,IDENTITY_FIELDS,MIXTURE_FIELDS,validate

R=Path(__file__).resolve().parent
P=json.loads((R/'protocol.json').read_text())
S=Path(P['source_root'])


def main():
    for file,fields in [('RNA_sample_template.tsv',SAMPLE_FIELDS),('RNA_identity_template.tsv',IDENTITY_FIELDS),
        ('RNA_physical_mixture_template.tsv',MIXTURE_FIELDS)]:
        pd.DataFrame(columns=fields).to_csv(R/file,sep='\t',index=False)
    validate(pd.DataFrame(columns=SAMPLE_FIELDS),pd.DataFrame(columns=IDENTITY_FIELDS),
        pd.DataFrame(columns=MIXTURE_FIELDS),R/'pending_experimental_inputs')
    saved=pd.read_csv(Path(P['v7_root'])/'all_predictions.csv.gz')
    old=saved[saved.method.eq('prototypes_target_only_budget') & saved.unit.eq('equalized_cell_fraction')].copy()
    meta=pd.read_csv(S/'mixture_index.csv')
    meta=meta[meta.unit.eq('equalized_cell_fraction')].copy()
    joined=old.merge(meta[['mixture_id','macrophage_cell_fraction']],on='mixture_id',validate='one_to_one')
    samples=[]
    for donor in sorted(joined.donor.unique()):
        samples.append(dict(sample_id=donor,patient_id=donor,split='validation',
            independent_patient_confirmed=not donor.startswith('GSE202642:'),previously_used_for_model=True,
            RNA_platform='historical_scRNA_computational',RNA_expression_unit='idealized_equal_cell_CP10k',
            QC_status='PASS',frozen_RNA_definition_id='v6_RNA_detected_macrophage',independent_RNA_evidence_id=''))
    samples=pd.DataFrame(samples,columns=SAMPLE_FIELDS)
    mixtures=[]
    for _,row in joined.iterrows():
        target=int(round(row.truth*2000)); macro=int(round(row.macrophage_cell_fraction*2000))
        # Remaining classes are not measured by this replay; retain as other_cells.
        mixtures.append(dict(mixture_id=int(row.mixture_id),sample_id=row.donor,
            preparation_replicate=row.replicate,mixture_kind=('standard' if row.scenario=='standard' else
                'dose_sensitivity' if row.scenario.startswith('macrophage_fraction_') else 'zero_target_competitor'),
            target_cells=target,other_macrophage_cells=macro-target,monocyte_cells=0,DC_cells=0,tumor_cells=0,
            other_cells=2000-macro,total_cells=2000,actual_target_fraction=row.truth,
            estimated_target_fraction=row.predicted,physical_mixture=False,
            bulk_platform='in_silico_equalized_scRNA',bulk_expression_unit='idealized_equal_cell_CP10k',
            reference_expression_unit='mean_cell_CP10k',unit_calibration_evidence_id='',QC_status='PASS'))
    mix=pd.DataFrame(mixtures,columns=MIXTURE_FIELDS)
    example=pd.read_csv(R/'actual_cell_annotations.csv')
    # Original PGAM5 observations only; no invented independent RNA measurements.
    identity=pd.DataFrame({'cell_id':example.cell_id,'sample_id':example.donor,
        'assigned_lineage':example.harmonized_group,'PGAM5_raw_counts':example.PGAM5_counts,
        'independent_RNA_truth_available':False,'independent_RNA_truth_macrophage':False,
        'independent_RNA_truth_PGAM5_detected':False,'RNA_identity_method':'','QC_status':'PASS'})
    # Some example cells' donors have no standard mixtures; include their real metadata.
    for donor in sorted(set(identity.sample_id)-set(samples.sample_id)):
        samples.loc[len(samples)]=[donor,donor,'validation',True,True,'historical_scRNA_computational',
            'idealized_equal_cell_CP10k','PASS','v6_RNA_detected_macrophage','']
    replay=R/'computational_replay'; replay.mkdir(exist_ok=True)
    samples.to_csv(replay/'samples.tsv',sep='\t',index=False,quoting=csv.QUOTE_ALL)
    identity.to_csv(replay/'RNA_identity.tsv',sep='\t',index=False)
    mix.to_csv(replay/'mixtures.tsv',sep='\t',index=False)
    validate(samples,identity,mix,replay,computational_replay=True)
    (replay/'README.md').write_text('These are historical computational mixtures, not physical experiments or new patients. All samples are marked previously used. Independent RNA truth is unavailable. Nonmacrophage subtype counts are not supplied by the replay and are grouped into other_cells; zero entries are bookkeeping placeholders, not measurements of subtype absence. Cell counts sum to2000 and preserve original target/total macrophage counts. No protein measurements are involved.\n')
    actual=pd.read_csv(S/'GSE149614_metadata.csv.gz')
    ids=pd.Index(actual.cell_id).get_indexer(example.cell_id)
    assert (ids>=0).all()
    raw=sparse.load_npz(S/'GSE149614_raw_common.npz')[ids]
    genes=pd.read_csv(S/'common_genes.csv').gene.to_numpy()
    # Reverse order and arbitrary obs names exercise exact cell_id mapping.
    order=np.arange(len(example))[::-1]
    a=ad.AnnData(raw[order],obs=pd.DataFrame({'cell_id':example.cell_id.to_numpy()[order]},
        index=[f'original_object_row_{i}' for i in range(len(order))]),var=pd.DataFrame(index=genes))
    path=R/'import_example'; path.mkdir(exist_ok=True)
    a.write_h5ad(path/'actual_60_HCC_cells.h5ad',compression='gzip')
    print('RNA templates empty; saved computational mixtures replayed; actual60-cell import example prepared.')


if __name__=='__main__': main()
