"""Exercise ID-based AnnData attachment on actual source cells, including failure."""
from pathlib import Path
import json
import argparse
import subprocess
import sys
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

R = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--saved-example',action='store_true',help='Use the packaged 60 actual source cells without large caches')
    args = parser.parse_args()
    cohort = 'GSE149614'
    out = R/'import_example'
    out.mkdir(exist_ok=True)
    before,after = out/'source60cells.h5ad',out/'annotated60cells.h5ad'
    if args.saved_example:
        a = ad.read_h5ad(before)
        small = a.obs[['cell_id']].copy()
    else:
        source = Path(json.loads((R/'protocol.json').read_text())['source_root'])
        obs = pd.read_csv(source/f'{cohort}_metadata.csv.gz')
        raw = sparse.load_npz(source/f'{cohort}_raw_common.npz')
        genes = pd.Index(pd.read_csv(source/'common_genes.csv').gene)
        mac = obs.harmonized_group.eq('Macrophage').to_numpy()
        positive = obs.PGAM5_counts.gt(0).to_numpy()
        rows = np.r_[np.flatnonzero(~mac)[:20],np.flatnonzero(mac&positive)[:20],np.flatnonzero(mac&~positive)[:20]]
        small = obs.iloc[rows][['cell_id']].copy()
        small.index = ['arbitrary_object_row_'+str(j) for j in range(len(small))]
        selected_genes = ['PGAM5','CD68','C1QA']
        a = ad.AnnData(X=raw[rows][:,genes.get_indexer(selected_genes)].copy(),obs=small,
                       var=pd.DataFrame(index=selected_genes))
        a.write_h5ad(before)
    command = [sys.executable,str(R/'attach_annotations.py'),'--input-h5ad',str(before),
        '--annotation-csv',str(R/f'{cohort}_all_cell_annotation.csv.gz'),'--output-h5ad',str(after),
        '--cell-id-column','cell_id']
    subprocess.run(command,check=True)
    result = ad.read_h5ad(after)
    expected = pd.read_csv(R/f'{cohort}_all_cell_annotation.csv.gz').set_index('cell_id').loc[small.cell_id]
    assert np.array_equal(a.X.toarray(),result.X.toarray())
    assert result.obs_names.equals(a.obs_names)
    assert np.array_equal(result.obs.pgam5_v6_PGAM5_RNA_status.astype(str),expected.PGAM5_RNA_status)
    assert np.array_equal(result.obs.pgam5_v6_cell_type.astype(str),expected.cell_type)
    expected_membership = expected.candidate_state_member.map({True:'True',False:'False'}).fillna('Not_applicable')
    assert np.array_equal(result.obs.pgam5_v6_candidate_state_member.astype(str),expected_membership)
    assert not result.obs.pgam5_v6_usable_for_validated_TCGA_PGAM5_cell_abundance.any()
    a.obs.loc[a.obs.index[0],'cell_id'] = 'INTENTIONALLY_UNMATCHED'
    unmatched = out/'unmatched_ID_test.h5ad'
    a.write_h5ad(unmatched)
    bad_command = command.copy()
    bad_command[bad_command.index('--input-h5ad')+1] = str(unmatched)
    bad_output = out/'should_not_exist.h5ad'
    bad_command[bad_command.index('--output-h5ad')+1] = str(bad_output)
    bad = subprocess.run(bad_command,capture_output=True,text=True)
    assert bad.returncode!=0 and not bad_output.exists()
    checks = {'status':'PASS','actual_source_cells':len(a),
        'RNA_status_and_lineage_exact':True,'candidate_membership_exact':True,
        'original_counts_and_order_unchanged':True,
        'unmatched_ID_rejected':True,'Seurat_R_import_executed':False,
        'interpretation':'ID-based attachment verified; biological PGAM5 subtype and abundance remain unvalidated'}
    (R/'import_check.json').write_text(json.dumps(checks,indent=2)+'\n')
    print(json.dumps(checks,indent=2))


if __name__=='__main__':
    main()
