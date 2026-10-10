"""Fit frozen exploratory RNA reference to measured, compatible CP10k inputs."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from guard_solver import fit


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--reference-dir',type=Path,required=True)
    p.add_argument('--unit',choices=['equalized_cell_fraction','library_RNA_contribution'],required=True)
    p.add_argument('--expression-cp10k',type=Path,required=True,help='TSV: gene rows, sample columns; CP10k using the full measured library denominator')
    p.add_argument('--output-csv',type=Path,required=True)
    args=p.parse_args()
    assert not args.output_csv.exists(), 'Preserve existing outputs'
    model=np.load(args.reference_dir/f'{args.unit}_model.npz')
    x=pd.read_csv(args.expression_cp10k,sep='\t',index_col=0)
    assert x.index.is_unique and x.columns.is_unique,'Duplicate genes/samples'
    assert pd.Index(model['genes']).isin(x.index).all(),'Reference genes missing; do not fill unmeasured genes with zeros'
    values=x.loc[model['genes']].to_numpy(float)
    assert np.isfinite(values).all() and (values>=0).all(),'Invalid expression'
    pg=int(np.flatnonzero(model['genes']=='PGAM5')[0]);rows=[]
    for j,sample in enumerate(x.columns):
        c,w,res=fit(model['A'],values[:,j]/model['scale'],True,pg,model['target_mask'])
        row={'sample':sample,'RNA_detected_macrophage_coefficient':float(c[model['target_mask']].sum()),
            'total_macrophage_coefficient':float(c[model['macrophage_mask']].sum()),'relative_residual':res,
            'unit':args.unit,'usable_as_validated_TCGA_cell_abundance':False}
        row.update({'component_'+g:float(v) for g,v in zip(model['groups'],c)})
        rows.append(row)
    pd.DataFrame(rows).to_csv(args.output_csv,index=False)
    print(f'Fitted {len(rows)} samples. Outputs are exploratory coefficients; platform and cell-RNA calibration remain unvalidated.')


if __name__=='__main__': main()
