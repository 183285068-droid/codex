"""Attach exact v8 RNA annotations, preserving counts/order and existing fields."""
from pathlib import Path
import argparse
import anndata as ad
import pandas as pd

FIELDS=['macrophage_identity','PGAM5_counts','PGAM5_RNA_status','PGAM5_macrophage_annotation',
    'macrophage_core_detected_genes','broad_myeloid_detected_genes','macrophage_RNA_marker_corroboration',
    'other_lineage_RNA_review_flag','RNA_lineage_evidence_level','PGAM5_related_state',
    'RNA_state_validation_passed','usable_for_validated_TCGA_PGAM5_cell_abundance']
FIELDS += ['candidate_module_'+v+'_'+field for v in ['continuous_rank','log2FC_1','log2FC_0p5'] for field in ['score','feature_coverage']]


def attach(input_path,annotation_path,output_path,cell_id_column=None):
    assert Path(input_path).resolve()!=Path(output_path).resolve(), 'Preserve original object'
    assert not Path(output_path).exists(), 'Output already exists; choose a new output path'
    a=ad.read_h5ad(input_path)
    ids=a.obs[cell_id_column].astype(str).to_numpy() if cell_id_column else a.obs_names.astype(str).to_numpy()
    assert pd.Index(ids).is_unique, 'Duplicate object cell IDs'
    table=pd.read_csv(annotation_path).set_index('cell_id')
    assert table.index.is_unique, 'Duplicate annotation cell IDs'
    assert pd.Index(ids).isin(table.index).all(), 'Unmatched exact cell IDs; use the matching cohort and HCC scope'
    assert not any('pgam5_v8_'+f in a.obs for f in FIELDS), 'v8 fields already exist'
    aligned=table.loc[ids]
    for field in FIELDS:
        values=aligned[field].to_numpy()
        a.obs['pgam5_v8_'+field]=pd.Categorical(values) if aligned[field].dtype=='object' else values
    a.write_h5ad(output_path,compression='gzip')
    return a


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--input-h5ad',type=Path,required=True)
    p.add_argument('--annotation-csv',type=Path,required=True)
    p.add_argument('--output-h5ad',type=Path,required=True)
    p.add_argument('--cell-id-column')
    a=p.parse_args()
    result=attach(a.input_h5ad,a.annotation_csv,a.output_h5ad,a.cell_id_column)
    print(f'Attached {len(result)} exact RNA annotations; stable state and TCGA cellular abundance remain unvalidated.')


if __name__=='__main__': main()
