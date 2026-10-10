"""Independent audits from portable raw counts and saved cell-level outputs."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import adjusted_rand_score
from sklearn.decomposition import PCA

R=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--full-source',type=Path);args=ap.parse_args()
    rows=[]
    def check(label,value):
        assert value,label;rows.append({'check':label,'passed':bool(value)})
    obs=pd.read_csv(R/'cell_annotations.csv.gz')
    genes=pd.read_csv(R/'portable_feature_genes.csv').gene.to_numpy()
    raw=sparse.load_npz(R/'portable_raw_feature_counts.npz')
    check('12,135 unique cells, seven user-confirmed patients',len(obs)==12135 and obs.cell_id.is_unique and obs.patient_id.nunique()==7)
    check('patient unit uses current user definition',np.array_equal(obs.patient_id,obs.sample_name) and obs.patient_identity_source.eq('User_confirmed_seven_independent_patients').all())
    check('portable raw matrix identifiers and integer counts',raw.shape==(len(obs),len(genes)) and len(set(genes))==len(genes) and (raw.data>=0).all() and np.array_equal(raw.data,np.round(raw.data)))
    pg=int(np.flatnonzero(genes=='PGAM5')[0]);pgraw=raw[:,pg].toarray().ravel()
    pos=pgraw>0
    check('PGAM5 RNA labels reconstructed from raw counts',np.array_equal(pos,obs.PGAM5_detected) and np.array_equal(pgraw,obs.PGAM5_counts) and int(pos.sum())==868)
    if (R/'RNA_subgroup_cell_annotation.csv.gz').exists():
        annotation=pd.read_csv(R/'RNA_subgroup_cell_annotation.csv.gz')
        check('exported annotation preserves all cell and patient IDs',np.array_equal(annotation.cell_id,obs.cell_id) and np.array_equal(annotation.patient_id,obs.patient_id))
        candidate=obs.primary_cluster.astype(str).eq('7').to_numpy()
        check('candidate subgroup uses PGAM5-free cluster C7',np.array_equal(annotation.RNA_subgroup_annotation.eq('Cycling_macrophage_PGAM5_enriched_candidate'),candidate) and int(candidate.sum())==403)
        check('RNA detected subset of candidate has79 cells',np.array_equal(annotation.PGAM5_detected_within_candidate,candidate&pos) and int((candidate&pos).sum())==79)
    norm=raw.multiply((1e4/obs.total_counts.to_numpy())[:,None]).tocsr()
    check('PGAM5 CP10k full-library normalization',np.allclose(norm[:,pg].toarray().ravel(),obs.PGAM5_CP10k,rtol=1e-6,atol=1e-7))
    cycle_genes=['MKI67','TOP2A','UBE2C','CENPF','CDK1','BIRC5']
    ix=[int(np.flatnonzero(genes==g)[0]) for g in cycle_genes]
    cy=np.asarray((raw[:,ix]>0).sum(axis=1)).ravel()>=2
    check('PGAM5-independent cycling definition',np.array_equal(cy,obs.cycling_flag) and int(cy.sum())==649 and int((cy&pos).sum())==132)
    pca=np.load(R/'primary_PCA_model.npz')
    check('PGAM5 excluded from all 2,000 fitted features','PGAM5' not in pca['genes'] and len(pca['genes'])==2000)
    check('cycle genes retained in fitted features',any(g in pca['genes'] for g in cycle_genes))
    feature_index=[int(np.flatnonzero(genes==g)[0]) for g in pca['genes']]
    x=norm[:,feature_index].toarray().astype(np.float32);x=np.log1p(x)
    check('PCA feature means and SD',np.allclose(x.mean(0),pca['mean'],atol=2e-6) and np.allclose(np.maximum(x.std(0,ddof=1),.1),pca['SD'],atol=2e-6))
    computed=(np.clip((x-pca['mean'])/pca['SD'],-10,10)-pca['pca_mean'])@pca['components'].T
    # sklearn randomized PCA fit_transform uses truncated U*S whereas transform
    # computes X*V; stochastic SVD approximation may leave small discrepancies.
    difference=np.linalg.norm(computed-pca['PC'])/max(np.linalg.norm(pca['PC']),1e-12)
    check('PCA component orthogonality',np.allclose(pca['components']@pca['components'].T,np.eye(30),atol=3e-6))
    # U*S from a randomized range approximation is not algebraically identical to
    # X*V. Reexecute the documented fit_transform instead of imposing a false
    # exact-projection identity. Retain the measured projection discrepancy.
    scaled=np.clip((x-pca['mean'])/pca['SD'],-10,10)
    refitted=PCA(n_components=30,svd_solver='randomized',random_state=0).fit_transform(scaled)
    check('actual randomized-PCA fit_transform reproduced from raw counts',np.allclose(refitted,pca['PC'],atol=1e-4,rtol=1e-5))
    a=pd.read_csv(R/'cluster_PGAM5_association.csv',dtype={'cluster':str})
    by_patient=pd.read_csv(R/'cluster_PGAM5_by_sample.csv',dtype={'cluster':str})
    for geometry,column in [('primary','primary_cluster'),('primary_resolution_0.4','primary_resolution_0.4'),('primary_resolution_1.2','primary_resolution_1.2')]:
        labels=obs[column].astype(str).to_numpy();table=a[a.geometry.eq(geometry)]
        check(geometry+' clusters partition all cells',int(table.cells.sum())==len(obs) and int(table.PGAM5_detected.sum())==int(pos.sum()))
        pvalues=[]
        for _,r in table.iterrows():
            mask=labels==r.cluster;n=int(mask.sum());p=int((mask&pos).sum());n0=len(obs)-n;p0=int((~mask&pos).sum())
            check(geometry+' C'+r.cluster+' detection and capture',n==r.cells and p==r.PGAM5_detected and np.isclose(p/n,r.PGAM5_detection_fraction) and np.isclose(p/pos.sum(),r.capture_of_all_PGAM5_detected))
            represented=obs.loc[mask,'patient_id'].nunique()
            check(geometry+' C'+r.cluster+' nonempty patient representation',represented==r.represented_samples)
            pvalue=fisher_exact([[p,n-p],[p0,n0-p0]],alternative='two-sided').pvalue;pvalues.append(pvalue)
            check(geometry+' C'+r.cluster+' Fisher p',np.isclose(pvalue,r.Fisher_cell_level_p,rtol=1e-10,atol=0))
            s=by_patient[by_patient.geometry.eq(geometry)&by_patient.cluster.eq(r.cluster)]
            for _,v in s.iterrows():
                inside=mask&obs.patient_id.eq(v.sample_name).to_numpy();outside=~mask&obs.patient_id.eq(v.sample_name).to_numpy()
                check(geometry+' C'+r.cluster+' '+v.sample_name+' raw count cells',int(inside.sum())==v.cluster_cells and int((inside&pos).sum())==v.cluster_PGAM5_detected)
                if inside.sum() and outside.sum():
                    m1=float(norm[inside,pg].mean());m0=float(norm[outside,pg].mean())
                    check(geometry+' C'+r.cluster+' '+v.sample_name+' zero-inclusive means',np.allclose([m1,m0],[v.PGAM5_CP10k_cluster_mean,v.PGAM5_CP10k_other_mean],rtol=1e-6,atol=1e-7))
        check(geometry+' BH family',np.allclose(multipletests(pvalues,method='fdr_bh')[1],table.Fisher_cell_level_BH_q,rtol=1e-10,atol=0))
    stable=pd.read_csv(R/'cluster_stability.csv',dtype={'cluster':str,'matched_cluster':str})
    for _,r in stable.iterrows():
        if r.method=='graph_seed':
            ref=obs.primary_cluster.astype(str).to_numpy();alt=obs['primary_seed'+str(int(r.replicate))+'_cluster'].astype(str).to_numpy()
        else:
            subset=pd.read_csv(R/f'primary_subsample_{int(r.replicate)}_labels.csv.gz',dtype={'cluster':str})
            ref=obs.set_index('cell_id').loc[subset.cell_id,'primary_cluster'].astype(str).to_numpy();alt=subset.cluster.to_numpy()
            counts=subset.merge(obs[['cell_id','patient_id']],on='cell_id').groupby('patient_id').size()
            expected=obs.groupby('patient_id').size().mul(.8).apply(np.floor).astype(int)
            check('subsampling retains80pct per patient seed'+str(int(r.replicate)),counts.equals(expected))
        inter=int(((ref==r.cluster)&(alt==r.matched_cluster)).sum());union=int(((ref==r.cluster)|(alt==r.matched_cluster)).sum())
        check('Jaccard '+r.method+' '+str(r.replicate)+' C'+r.cluster,np.isclose(inter/union,r.Jaccard,atol=1e-12))
    W=sparse.load_npz(R/'primary_neighbor_graph.npz').tocsr();d=np.asarray(W.sum(axis=1)).ravel()
    observed=float(np.mean((W@pos.astype(float)/d)[pos]))
    graphstats=pd.read_csv(R/'PGAM5_graph_concentration.csv').iloc[0]
    check('graph concentration uses actual PGAM5 labels',np.isclose(observed,graphstats.observed_positive_neighbor_fraction,atol=1e-12))
    null=pd.read_csv(R/'PGAM5_graph_permutation_null.csv').positive_weighted_neighbor_fraction.to_numpy()
    check('500-permutation empirical p with plus-one correction',len(null)==500 and np.isclose((1+int((null>=observed).sum()))/501,graphstats.one_sided_permutation_p))
    marker=pd.read_csv(R/'cluster_marker_directions_by_patient.csv',dtype={'cluster':str})
    for _,r in marker.iterrows():
        j=int(np.flatnonzero(genes==r.gene)[0]);inside=obs.patient_id.eq(r.patient_id).to_numpy()&obs.primary_cluster.astype(str).eq(r.cluster).to_numpy();outside=obs.patient_id.eq(r.patient_id).to_numpy()&~obs.primary_cluster.astype(str).eq(r.cluster).to_numpy()
        if inside.sum() and outside.sum():
            ratio=np.log2((float(norm[inside,j].mean())+.05)/(float(norm[outside,j].mean())+.05))
            check('marker direction '+r.cluster+' '+r.gene+' '+r.patient_id,np.isclose(ratio,r.descriptive_log2_ratio_pseudocount0p05,atol=1e-10))
    dots=pd.read_csv(R/'marker_dotplot_values.csv',dtype={'cluster':str});log=norm.copy();log.data=np.log1p(log.data)
    for _,r in dots.iterrows():
        mask=obs.primary_cluster.astype(str).eq(r.cluster).to_numpy();j=int(np.flatnonzero(genes==r.gene)[0])
        check('dotplot '+r.cluster+' '+r.gene,np.allclose([float((raw[mask,j]>0).mean()),float(log[mask,j].mean())],[r.detection_fraction,r.mean_log1p_CP10k],atol=1e-6))
    full_source_checked=False
    if args.full_source:
        import anndata as ad
        original=ad.read_h5ad(args.full_source);original=original[obs.cell_id.to_numpy()]
        check('original full raw matrix matches portable counts',(original[:,genes].X.tocsr()!=raw).nnz==0)
        check('all full-cell normalization denominators',np.array_equal(np.asarray(original.X.sum(axis=1)).ravel(),obs.total_counts))
        full_source_checked=True
    pd.DataFrame(rows).to_csv(R/'independent_audit_checks.csv',index=False)
    result={'status':'PASS','checks':len(rows),'all_passed':True,'original_full_count_matrix_checked':full_source_checked,
            'PCA_loading_relative_reconstruction_error':float(difference),
            'audit_scope':'Raw labels, full-library normalization, actual cluster/patient counts, Fisher/BH, resampling Jaccard, graph concentration/permutation p, marker directions and plotted values. Does not prove biological identity or calibrated abundance.'}
    (R/('original_source_audit.json' if full_source_checked else 'portable_audit.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
