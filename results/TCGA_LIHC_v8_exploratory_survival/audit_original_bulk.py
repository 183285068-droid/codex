"""Optional original-matrix audit; large GEO source is not included in the ZIP."""
from pathlib import Path
import argparse,gzip,json
import numpy as np
import pandas as pd

R = Path(__file__).resolve().parent

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source',type=Path,default=Path('/workspace/scratch/TCGA_LIHC_candidate30_survival/TCGA24_tumor_featurecounts.txt.gz'))
    args = p.parse_args()
    stored = pd.read_csv(R/'TCGA_reference_union_raw_counts.csv.gz',index_col=0)
    meta = pd.read_csv(R/'TCGA_primary_aliquot_selection.csv').set_index('patient').loc[stored.columns]
    totals = np.zeros(len(meta),dtype=np.int64)
    found,allgenes,rows = set(),set(),0
    with gzip.open(args.source,'rt') as f:
        header = f.readline().rstrip('\r\n').split('\t')[1:]
        indices = [header.index(s) for s in meta['sample']]
        for line in f:
            fields = line.rstrip().split('\t'); gene = fields[0]
            assert len(fields) == len(header)+1
            assert gene not in allgenes
            allgenes.add(gene);rows += 1
            values = np.array([int(fields[j+1]) for j in indices],dtype=np.int64)
            assert (values >= 0).all()
            totals += values
            if gene in stored.index:
                assert np.array_equal(values,stored.loc[gene].to_numpy()),gene
                found.add(gene)
    assert found == set(stored.index)
    assert np.array_equal(totals,meta.full_library_total_counts.to_numpy())
    coverage = {}
    for unit in ['library_RNA_contribution','equalized_cell_fraction']:
        model = np.load(R/'reference'/unit/f'{unit}_model.npz')
        missing = pd.read_csv(R/unit/'unmeasured_reference_genes.csv').gene
        assert set(model['genes'])-allgenes == set(missing)
        coverage[unit] = {'original_genes':len(model['genes']),'exact_matched_genes':len(model['genes'])-len(missing)}
    result = {'status':'PASS','original_matrix_gene_rows':rows,'original_matrix_sample_columns':len(header),
              'patients':len(meta),'checked_count_rows':len(found),'all_371_full_library_denominators_checked':True,
              'missing_genes_checked_against_full_original_matrix':True,'reference_coverage':coverage}
    (R/'original_bulk_source_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__ == '__main__': main()
