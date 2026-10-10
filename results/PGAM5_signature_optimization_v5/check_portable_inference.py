"""Check that the exported model reproduces actual TCGA fits without large sources."""
from pathlib import Path
import json
import subprocess
import sys
import numpy as np
import pandas as pd

R = Path(__file__).resolve().parent


def main():
    bulk = pd.read_csv(R/'TCGA_used_reference_gene_CP10k.csv.gz', index_col=0).iloc[:, :3]
    bulk.to_csv(R/'portable_example_CP10k.tsv', sep='\t')
    subprocess.run([sys.executable, str(R/'fit_saved_reference.py'), '--bulk', str(R/'portable_example_CP10k.tsv'),
                    '--unit', 'CP10k', '--output', str(R/'portable_example_coefficients.csv')], check=True)
    actual = pd.read_csv(R/'EXPLORATORY_TCGA_all_component_coefficients.csv', index_col=0).loc[bulk.columns]
    saved = pd.read_csv(R/'portable_example_coefficients.csv').set_index('sample').loc[bulk.columns]
    difference = float(abs(actual.to_numpy()-saved[actual.columns].to_numpy()).max())
    assert difference < 1e-7, difference
    assert not saved.usable_as_cellular_infiltration.any()
    assert not saved.reference_validation_passed.any()
    result = {'status':'PASS', 'patients':3, 'max_coefficient_difference':difference,
              'input_unit':'CP10k, original full-library denominator preserved before feature subset',
              'biological_validation_passed':False}
    (R/'portable_inference_check.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
