"""Safely attach exact cell-ID annotations to AnnData; never relabel unmatched cells."""
from pathlib import Path
import argparse
import anndata as ad
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-h5ad',type=Path,required=True)
    parser.add_argument('--annotation-csv',type=Path,required=True)
    parser.add_argument('--output-h5ad',type=Path,required=True)
    parser.add_argument('--cell-id-column',help='Existing obs column containing exact cell IDs; default obs_names')
    args = parser.parse_args()
    assert args.input_h5ad.resolve()!=args.output_h5ad.resolve(), 'Preserve the original object'
    a = ad.read_h5ad(args.input_h5ad)
    ids = a.obs[args.cell_id_column].astype(str) if args.cell_id_column else pd.Series(a.obs_names.astype(str),index=a.obs_names)
    assert ids.is_unique, 'Duplicate object cell IDs'
    table = pd.read_csv(args.annotation_csv).set_index('cell_id')
    assert table.index.is_unique, 'Duplicate annotation cell IDs'
    missing = ids[~ids.isin(table.index)]
    assert not len(missing), f'{len(missing)} object cells lack exact annotations; use the matching HCC scope and IDs'
    table = table.loc[ids.to_numpy()]
    fields = ['cell_type','PGAM5_RNA_status','PGAM5_macrophage_annotation','PGAM5_related_state_annotation',
              'dominant_program','candidate_program_score','candidate_state_member','usable_for_validated_TCGA_PGAM5_cell_abundance']
    for field in fields:
        values = table[field].to_numpy()
        if field in ['cell_type','PGAM5_RNA_status','PGAM5_macrophage_annotation','PGAM5_related_state_annotation']:
            a.obs['pgam5_v6_'+field] = pd.Categorical(values)
        elif field=='candidate_state_member':
            a.obs['pgam5_v6_'+field] = pd.Categorical(pd.Series(values).map({True:'True',False:'False'}).fillna('Not_applicable'))
        else:
            a.obs['pgam5_v6_'+field] = values
    a.write_h5ad(args.output_h5ad,compression='gzip')
    print(f'Attached {len(a)} exact cell annotations; PGAM5 cellular abundance validation remains false.')


if __name__=='__main__':
    main()
