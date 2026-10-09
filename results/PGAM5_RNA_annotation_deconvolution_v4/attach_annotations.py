"""Attach exported RNA-status annotation by exact cell_id; never silently reorder mismatched IDs."""
import argparse
import anndata as ad
import pandas as pd
p=argparse.ArgumentParser();p.add_argument('--input-h5ad',required=True);p.add_argument('--annotation-csv',required=True);p.add_argument('--output-h5ad',required=True);args=p.parse_args()
if args.input_h5ad==args.output_h5ad:raise ValueError('Output must differ from input to preserve original data')
a=ad.read_h5ad(args.input_h5ad);o=pd.read_csv(args.annotation_csv,dtype={'cell_id':str});assert o.cell_id.is_unique and a.obs_names.is_unique
missing=a.obs_names.difference(o.cell_id)
if len(missing):raise ValueError(f'{len(missing)} input cells have no exact annotation match; check sample-prefixed cell IDs before attaching')
o=o.set_index('cell_id').reindex(a.obs_names)
for col in ['cell_type','PGAM5_RNA_status','PGAM5_macrophage_annotation','cycling_flag']:a.obs[col]=o[col].to_numpy()
ad.settings.allow_write_nullable_strings=True;a.write_h5ad(args.output_h5ad,compression='gzip');print(f'Attached checked annotation to {len(a)} cells')
