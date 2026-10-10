"""Reuse outcome-blind 371 primary-tumor aliquots and full-library raw denominators."""
from pathlib import Path
import gzip
import json
import hashlib
import numpy as np
import pandas as pd

R = Path(__file__).resolve().parent
SOURCE = Path('/workspace/scratch/TCGA_LIHC_candidate30_survival')


def main():
    selection = pd.read_csv(SOURCE/'RNA_aliquot_selection.csv')
    keep = selection[selection.selected].sort_values('patient').copy()
    assert len(keep) == 371 and keep.patient.is_unique
    assert keep['sample'].str.split('-').str[3].str.startswith('01').all()
    genes = set(pd.read_csv(R/'common_genes.csv').gene)
    values, names = [], []
    library = np.zeros(len(keep))
    measured_genes = 0
    source = SOURCE/'TCGA24_tumor_featurecounts.txt.gz'
    with gzip.open(source, 'rt') as f:
        header = f.readline().rstrip().split('\t')[1:]
        index = pd.Index(header).get_indexer(keep['sample'])
        assert (index >= 0).all()
        for line in f:
            gene, separator, counts = line.rstrip().partition('\t')
            assert separator
            row = np.fromstring(counts, sep='\t')
            assert len(row) == len(header)
            row = row[index]
            assert np.isfinite(row).all() and (row >= 0).all() and (row == np.floor(row)).all()
            library += row
            measured_genes += 1
            if gene in genes:
                values.append(row.astype(np.int64))
                names.append(gene)
    assert len(set(names)) == len(names)
    assert np.allclose(library, keep.assigned_gene_counts, atol=1e-6, rtol=0)
    np.savez_compressed(R/'TCGA_primary_common_gene_counts.npz',
                        counts=np.stack(values), genes=np.array(names), patients=keep.patient.to_numpy(dtype=str))
    keep['full_library_total_counts'] = library
    keep.to_csv(R/'TCGA_primary_aliquot_selection.csv', index=False)
    h = hashlib.sha256()
    with source.open('rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''):
            h.update(chunk)
    provenance = {'patients': len(keep), 'primary_tumor_sample_type': '01',
        'source': str(source), 'source_SHA256': h.hexdigest(),
        'source_bytes': source.stat().st_size, 'all_measured_bulk_genes': measured_genes,
        'measured_common_reference_universe': len(names),
        'selection': 'Reuse previous outcome-blind aliquot selection, one selected primary tumor per patient',
        'bulk_unit': 'GSE62944 FeatureCounts read counts; full-library CPM/CP10k. Not TPM.',
        'new_TPM_access': 'Public GDC and Xena endpoint checks blocked by egress proxy HTTP403; no TPM fabricated or inferred using gene-span length',
        'limitation': 'Bulk read counts and scRNA UMI have uncalibrated transcript-length/platform effects; coefficients cannot be validated cellular abundance'}
    (R/'TCGA_input_provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    print(json.dumps(provenance, indent=2))


if __name__ == '__main__':
    main()
