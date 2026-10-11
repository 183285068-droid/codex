"""Independent checks against downloaded sparse Matrix Market counts."""
from pathlib import Path
import gzip
import hashlib
import json
import os
import anndata as ad
import numpy as np
import pandas as pd
import scipy.io

ROOT=Path(__file__).resolve().parent
WORK=Path(os.environ.get('GSE283406_WORK',ROOT/'work'))
PROGRAMS=json.loads((ROOT/'programs.json').read_text())
checks=[]
def check(name,valid,detail=''):
    checks.append(dict(check=name,passed=bool(valid),detail=detail))
    if not valid:raise AssertionError(name)

manifest=json.loads((ROOT/'source_manifest.json').read_text())
for item in manifest:
    p=WORK/'inputs'/item['name']
    check('bytes_'+p.name,p.stat().st_size==item['bytes'])
    check('sha256_'+p.name,hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256'])
check('missing_coordinates_recorded',not any('position' in x['name'] or 'scalefactor' in x['name'] for x in manifest))
summary=pd.read_csv(ROOT/'sample_QC_summary.csv').set_index('gsm')
raw_barcodes={}
coverage=pd.read_csv(ROOT/'program_gene_coverage.csv')
for gsm,section in [('GSM8661789','J259B'),('GSM8661790','J259X'),('GSM8661791','J259B'),('GSM8661792','J259X')]:
    prefix=WORK/'inputs'/f'{gsm}_{section}'
    f=pd.read_csv(str(prefix)+'.features.tsv.gz',header=None,sep='\t')
    b=pd.read_csv(str(prefix)+'.barcodes.tsv.gz',header=None)[0].astype(str).to_numpy()
    raw_barcodes[gsm]=set(b)
    with gzip.open(str(prefix)+'.matrix.mtx.gz','rb') as handle:x=scipy.io.mmread(handle,spmatrix=True).T.tocsr()
    x.eliminate_zeros()
    total=np.asarray(x.sum(1)).ravel();n=x.getnnz(1)
    qc=(total>=500)&(n>=200)
    pg=x[:,np.flatnonzero(f[1].eq('PGAM5'))[0]].toarray().ravel()
    allobs=pd.read_csv(ROOT/f'{gsm}_all_spot_QC.csv').set_index('barcode').loc[b]
    check(gsm+'_raw_total_UMI',np.array_equal(total,allobs.total_counts))
    check(gsm+'_raw_n_genes',np.array_equal(n,allobs.n_genes))
    check(gsm+'_raw_PGAM5',np.array_equal(pg,allobs.PGAM5_raw))
    check(gsm+'_QC_gate',np.array_equal(qc,allobs.pass_qc))
    check(gsm+'_unique_ensembl_ids',f[0].is_unique)
    check(gsm+'_duplicate_symbols_count',int(f[1].duplicated().sum())==3)
    obs=pd.read_csv(ROOT/f'{gsm}_QC_expression_annotations.csv').set_index('barcode')
    check(gsm+'_QC_barcode_alignment',list(obs.index)==list(b[qc]))
    expected_pg=np.log1p(1e4*pg[qc]/total[qc])
    check(gsm+'_full_library_PGAM5_normalization',np.allclose(expected_pg,obs.PGAM5_log1p_CP10k,atol=2e-6))
    check(gsm+'_zero_retention',np.array_equal(pg[qc]==0,obs.PGAM5_log1p_CP10k.eq(0)))
    check(gsm+'_summary_count',len(obs)==int(summary.loc[gsm,'qc_spots']))
    for program,genes in PROGRAMS.items():
        ids=[np.flatnonzero(f[1].eq(g))[0] for g in genes if f[1].eq(g).any()]
        counts=x[qc][:,ids].toarray()
        expected=np.log1p(1e4*counts/total[qc,None]).mean(1)
        check(gsm+'_'+program+'_score',np.allclose(expected,obs[program],atol=2e-6))
        check(gsm+'_'+program+'_detection',np.array_equal((counts>0).sum(1),obs[program+'_n_detected']))
        c=coverage[coverage.gsm.eq(gsm)&coverage.program.eq(program)].iloc[0]
        check(gsm+'_'+program+'_coverage',len(ids)==c.present and len(genes)==c.requested)
    mac=(obs.Macrophage_core.ge(obs.Macrophage_core.quantile(.75))&obs.Macrophage_core_n_detected.ge(3)&obs.Pan_myeloid.ge(obs.Pan_myeloid.median()))
    t=obs.T_lineage.ge(obs.T_lineage.quantile(.75))&obs.T_lineage_n_detected.ge(2)
    check(gsm+'_Mac_RNA_gate',np.array_equal(mac,obs.Mac_RNA_enriched))
    check(gsm+'_T_RNA_gate',np.array_equal(t,obs.T_RNA_enriched))
    for k in ['Activation_early','Exhaustion_core','Effector_cytotoxic']:
        expected=t&obs[k].ge(obs.loc[t,k].quantile(.75))&obs[k].gt(0)&obs[k+'_n_detected'].ge(2)
        check(gsm+'_'+k+'_state',np.array_equal(expected,obs[k+'_supported']))
    if gsm in ['GSM8661791','GSM8661792']:
        a=ad.read_h5ad(ROOT/f'{section}_AVITI_expression_marker_subset.h5ad')
        check(gsm+'_no_fake_spatial_coordinates','spatial' not in a.obsm and not bool(a.uns['coordinates_available']))
        check(gsm+'_UMAP_finite',np.isfinite(a.obsm['X_umap']).all())
        check(gsm+'_portable_barcode_alignment',list(a.obs_names)==list(obs.index))
        for gene in a.var_names:
            symbol=a.var.loc[gene,'original_gene_symbol'];geneid=a.var.loc[gene,'gene_id']
            i=np.flatnonzero(f[0].eq(geneid))[0]
            counts=x[qc,i].toarray().ravel()
            check(gsm+'_'+gene+'_portable_counts',np.array_equal(counts,a[:,[gene]].layers['counts'].toarray().ravel()))
            expected=np.log1p(1e4*counts/total[qc])
            check(gsm+'_'+gene+'_portable_normalization',np.allclose(expected,a[:,[gene]].X.toarray().ravel(),atol=2e-6))

for section,ng,ag in [('J259B','GSM8661789','GSM8661791'),('J259X','GSM8661790','GSM8661792')]:
    check(section+'_same_raw_barcode_set',raw_barcodes[ng]==raw_barcodes[ag])
    for res in ['hires','lowres']:
        n=gzip.decompress((WORK/'inputs'/f'{ng}_{section}.tissue_{res}_image.png.gz').read_bytes())
        a=gzip.decompress((WORK/'inputs'/f'{ag}_{section}.tissue_{res}_image.png.gz').read_bytes())
        check(section+'_'+res+'_same_original_image',n==a)
    pairs=pd.read_csv(ROOT/f'{section}_paired_platform_expression.csv')
    check(section+'_unique_paired_barcodes',pairs.barcode.is_unique)
    n=pd.read_csv(ROOT/f'{ng}_QC_expression_annotations.csv').set_index('barcode')
    a=pd.read_csv(ROOT/f'{ag}_QC_expression_annotations.csv').set_index('barcode')
    check(section+'_shared_QC_set',set(pairs.barcode)==set(n.index)&set(a.index))
    for platform,data in [('NextSeq',n),('AVITI',a)]:
        check(section+'_'+platform+'_paired_expression',np.allclose(pairs[platform+'_PGAM5_log1p_CP10k'],data.loc[pairs.barcode,'PGAM5_log1p_CP10k']))
evidence=pd.read_csv(ROOT/'expression_cluster_naming_evidence.csv')
check('all_naming_markers_enriched',evidence.cluster_mean_log1p_CP10k.gt(evidence.other_mean_log1p_CP10k).all())
check('all_naming_programs_above_mean',evidence.program_mean.gt(evidence.program_all_spot_mean).all())
pd.DataFrame(checks).to_csv(ROOT/'audit_checks.csv',index=False)
(ROOT/'validation.json').write_text(json.dumps({'passed':len(checks),'failed':0,'physical_spatial_analysis':'blocked_missing_coordinates'},indent=2)+'\n')
print(f'Passed {len(checks)} independent raw-count, annotation and technical-pair checks.')
