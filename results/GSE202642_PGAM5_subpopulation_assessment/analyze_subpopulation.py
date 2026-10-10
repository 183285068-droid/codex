"""Retrospective within-GSE202642 RNA subgroup assessment; no PGAM5 in geometry."""
from pathlib import Path
import argparse,hashlib,json,random,shutil
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import igraph as ig
from scipy import sparse
from scipy.stats import fisher_exact,mannwhitneyu,spearmanr
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score,silhouette_samples
from statsmodels.stats.multitest import multipletests

R = Path(__file__).resolve().parent
DEFAULT_SOURCE = Path('/workspace/scratch/PGAM5_reproducible_macrophage_states/GSE202642_marker_consistent_macrophages.h5ad')
CYCLE_ANCHORS = ['MKI67','TOP2A','UBE2C','CENPF','CDK1','BIRC5']
PANELS = {
    'Macrophage':['C1QA','C1QB','C1QC','CD68','CSF1R','MERTK','MSR1','CD163'],
    'Cycling':CYCLE_ANCHORS+['TK1','RRM2','PCLAF','LMNB1'],
    'Inflammatory':['IL1B','CXCL8','TNF','NFKBIA','FCN1','S100A8','S100A9'],
    'TAM_lipid':['TREM2','APOE','LPL','FABP5','LGALS3','SPP1'],
    'Antigen_presentation':['HLA-DRA','HLA-DRB1','HLA-DPA1','CD74','CTSS'],
    'DC_competition':['CD1C','FCER1A','CLEC10A','CLEC9A','LILRA4'],
    'Stress':['FOS','JUN','JUNB','DUSP1','ATF3','HSPA1A','HSPA1B','HSP90AA1']}
BACKGROUND = set('ALB APOA1 APOA2 APOC3 TTR RBP4 APOH AHSG FGB FGA FGG ORM1 ORM2 PRAP1 TF GSTA1 CES1 CD3D CD3E TRAC CD79A MS4A1 COL1A1 COL1A2 DCN PECAM1 VWF'.split())

def hash_file(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for c in iter(lambda:f.read(1024**2),b''):h.update(c)
    return {'path':str(p),'bytes':p.stat().st_size,'sha256':h.hexdigest()}

def graph_from_pc(pc,obs):
    b=ad.AnnData(np.zeros((len(pc),1),dtype=np.float32),obs=obs.copy())
    b.obsm['X_pca']=pc
    sc.pp.neighbors(b,n_neighbors=20,n_pcs=30,use_rep='X_pca',transformer='sklearn',random_state=0)
    return b

def cluster(b,resolution,seed,key):
    sc.tl.leiden(b,resolution=resolution,random_state=seed,key_added=key,
                 flavor='igraph',directed=False,n_iterations=4)
    return b.obs[key].astype(str).to_numpy()

def match_clusters(ref,alt):
    rows=[]
    for group in sorted(set(ref),key=int):
        mask=ref==group
        best=None
        for other in sorted(set(alt),key=int):
            hit=alt==other; inter=int((mask&hit).sum());union=int((mask|hit).sum())
            row={'cluster':group,'matched_cluster':other,'Jaccard':inter/union,
                 'original_cells':int(mask.sum()),'matched_cells':int(hit.sum()),'intersection_cells':inter}
            if best is None or row['Jaccard']>best['Jaccard']:best=row
        rows.append(best)
    return rows

def association(obs,label,geometry):
    rows,sample_rows=[],[]
    pos=obs.PGAM5_detected.to_numpy(bool)
    for group in sorted(set(label),key=int):
        inside=label==group;n=int(inside.sum());p=int((inside&pos).sum())
        n0=len(obs)-n;p0=int((~inside&pos).sum())
        odds,pval=fisher_exact([[p,n-p],[p0,n0-p0]],alternative='two-sided')
        directions=[];expression_ratios=[]
        for sample in sorted(obs.sample_name.unique()):
            s=obs.sample_name.to_numpy()==sample;ii=inside&s;oo=~inside&s
            ns,no=int(ii.sum()),int(oo.sum());ps,po=int((ii&pos).sum()),int((oo&pos).sum())
            eligible=ns>=20 and no>=20 and int((s&pos).sum())>=3
            rate1,rate0=ps/ns if ns else np.nan,po/no if no else np.nan
            mean1=float(obs.loc[ii,'PGAM5_CP10k'].mean());mean0=float(obs.loc[oo,'PGAM5_CP10k'].mean())
            ratio=(mean1+.05)/(mean0+.05) if ns and no else np.nan
            if eligible:directions.append(rate1>rate0);expression_ratios.append(ratio)
            sample_rows.append({'geometry':geometry,'cluster':group,'sample_name':sample,
                'cluster_cells':ns,'other_cells':no,'cluster_PGAM5_detected':ps,'other_PGAM5_detected':po,
                'cluster_detection_rate':rate1,'other_detection_rate':rate0,
                'PGAM5_CP10k_cluster_mean':mean1,'PGAM5_CP10k_other_mean':mean0,
                'CP10k_mean_ratio_pseudocount_0p05':ratio,'evaluable_sample':eligible})
        samples=obs.loc[inside,'sample_name'].value_counts()
        samples=samples[samples.gt(0)]
        rows.append({'geometry':geometry,'cluster':group,'cells':n,'PGAM5_detected':p,
            'PGAM5_detection_fraction':p/n,'capture_of_all_PGAM5_detected':p/int(pos.sum()),
            'other_detection_fraction':p0/n0,'pooled_detection_rate_ratio':(p/n)/(p0/n0) if p0 else np.nan,
            'Fisher_cell_level_odds_ratio':odds,'Fisher_cell_level_p':pval,
            'cluster_PGAM5_mean_CP10k':float(obs.loc[inside,'PGAM5_CP10k'].mean()),
            'other_PGAM5_mean_CP10k':float(obs.loc[~inside,'PGAM5_CP10k'].mean()),
            'cycling_fraction':float(obs.loc[inside,'cycling_flag'].mean()),
            'represented_samples':len(samples),'samples_with_ge20_cells':int(samples.ge(20).sum()),
            'largest_sample_fraction':float(samples.max()/n),'evaluable_samples':len(directions),
            'fraction_evaluable_samples_detection_enriched':float(np.mean(directions)) if directions else np.nan,
            'median_sample_CP10k_ratio':float(np.median(expression_ratios)) if expression_ratios else np.nan})
    tab=pd.DataFrame(rows);tab['Fisher_cell_level_BH_q']=multipletests(tab.Fisher_cell_level_p,method='fdr_bh')[1]
    return tab,pd.DataFrame(sample_rows)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,default=DEFAULT_SOURCE);args=ap.parse_args()
    assert not (R/'result_summary.json').exists(),'Preserve completed analysis; use a fresh directory to rerun.'
    cycle_source=Path('/workspace/scratch/PGAM5_signature_optimization_v5/cell_cycle_provenance.json')
    cycle_data=json.loads(cycle_source.read_text());cycles=set(cycle_data['S']+cycle_data['G2M']+cycle_data['extension'])
    shutil.copy2(cycle_source,R/'cell_cycle_provenance.json')
    protocol={
        'question':'Can all PGAM5 RNA-detected macrophages in GSE202642 be treated as an independently supported macrophage subtype?',
        'scope':'GSE202642 HCC tumor libraries suffix5-11 only. Existing v8 12,135 marker-consistent macrophages; seven sample labels treated as seven independent patients according to current user confirmation. No other cohort, protein, or TCGA analysis.',
        'history':'Retrospective assessment. Earlier analyses already suggested cycling association; no claim of prospective blinding.',
        'RNA_label':'PGAM5 raw count>0 means detected, raw count0 means undetected. Zero counts and cycling cells retained.',
        'primary_geometry':'Full-cell CP10k then log1p; remove PGAM5, MT/RPL/RPS/Ig/HBA/HBB and explicit nonmacrophage-background genes from features only. Seurat HVG selection batch_key=sample_name,2000genes; per-gene z-scale SDfloor0.1,clip[-10,10],PCA30 randomized seed0;20NN exact sklearn,Scanpy UMAP fuzzy graph;igraph Leiden resolution0.8,4iterations,seed0.',
        'cycling_scope':'Cycle genes and cycling cells retained in the primary clustering. No requirement that a subgroup remain separate after cycle-program removal. Cycling markers characterize biology and do not disqualify a subgroup.',
        'technical_robustness':'Same graph seeds1/2, resolutions0.4/1.2, plus80%-of-each-sample cell subsampling seeds11/12/13 with new scaling/PCA/graph. Fixed features for subsample checks. Report best-match Jaccard and ARI; no result-driven choice of primary resolution.',
        'sample_repetition':'Report each cluster across seven user-confirmed independent patient labels. Association evaluable if cluster>=20 and other>=20cells and totalPGAM5detected>=3.',
        'statistics':'Two-sided Fisher per cluster with BH within geometry, descriptive cell-level inference; within-patient-label permutation of graph PGAM5-positive-neighbor concentration,500permutations; pooled arithmetic-mean CP10k expression contrasts and detection differences are descriptive marker profiles, post-clustering and not independent subtype validation. No Welch or paired gene tests performed.',
        'not_performed':'No depth matching/regression, no patient-paired DE, no cells removed by PGAM5 count threshold, no ambient-RNA/doublet correction, no external cohort/protein tests.',
        'decision':'Assess independent cluster structure, concentration/capture of PGAM5 detections, non-PGAM5 marker programs, sample representation and robustness jointly. Detection enrichment alone does not establish a subtype; no universal raw-count/FDR/log2FC threshold is claimed.',
        'patient_identity_source':'Current user explicitly confirmed that the seven labels represent seven independent patients; this is the governing definition for the analysis.',
        'abundance_scope':'Even a reproducible within-dataset RNA state does not by itself validate TCGA cell-fraction deconvolution.'}
    (R/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    a=ad.read_h5ad(args.source);raw=a.X.tocsr().astype(np.int32)
    obs=a.obs.copy();obs.index.name='cell_id_index';obs['cell_id']=a.obs_names.astype(str)
    obs=obs.rename(columns={'donor_id':'upstream_donor_id','donor_kind':'upstream_donor_kind'})
    obs['patient_id']=obs.sample_name.astype(str)
    obs['patient_identity_source']='User_confirmed_seven_independent_patients'
    genes=np.array(a.var_names,dtype=str)
    assert len(obs)==12135 and obs.index.is_unique and a.var_names.is_unique
    assert (raw.data>=0).all() and np.array_equal(raw.sum(axis=1).A.ravel(),obs.total_counts.to_numpy())
    pg=int(np.flatnonzero(genes=='PGAM5')[0]);obs['PGAM5_counts']=raw[:,pg].toarray().ravel()
    obs['PGAM5_detected']=obs.PGAM5_counts.gt(0)
    norm=raw.multiply((1e4/obs.total_counts.to_numpy())[:,None]).tocsr().astype(np.float32)
    log=norm.copy();log.data=np.log1p(log.data)
    obs['PGAM5_CP10k']=norm[:,pg].toarray().ravel()
    cy_ix=[int(np.flatnonzero(genes==g)[0]) for g in CYCLE_ANCHORS]
    obs['cycling_flag']=np.asarray((raw[:,cy_ix]>0).sum(axis=1)).ravel()>=2
    for panel,members in PANELS.items():
        ids=np.flatnonzero(np.isin(genes,members));obs[panel+'_module']=np.asarray(log[:,ids].mean(axis=1)).ravel()
    obs.to_csv(R/'input_cell_metadata.csv.gz')
    (R/'source_provenance.json').write_text(json.dumps({'input':hash_file(args.source),'cell_cycle_source':hash_file(cycle_source),
        'GEO':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202642',
        'RNA_identity_source':'https://github.com/183285068-droid/codex/blob/main/results/PGAM5_RNA_identity_v8/GSE202642_RNA_annotation.csv.gz',
        'annotation_limit':'Existing marker-inferred C1Q-rich macrophages; not author labels or exhaustive all-macrophage coverage.'},indent=2)+'\n')
    print('INPUT',len(obs),'macrophages',int(obs.PGAM5_detected.sum()),'PGAM5detected',flush=True)
    bases,stability,global_stability,associations,sample_associations,all_hvg=[],[],[],[],[],[]
    graphs={};primary_labels={};markers=[];module_profiles=[];model_parameters={};features=set()
    for geometry in ['primary']:
        eligible=np.array([g!='PGAM5' and g not in BACKGROUND and not g.startswith(('MT-','RPL','RPS','IGH','IGK','IGL','HBA','HBB')) for g in genes])
        ids=np.flatnonzero(eligible)
        c=ad.AnnData(log[:,ids],obs=obs[['sample_name']].copy(),var=pd.DataFrame(index=genes[ids]))
        sc.pp.highly_variable_genes(c,n_top_genes=2000,flavor='seurat',batch_key='sample_name')
        hv=c.var.copy();hv['gene']=hv.index;hv['geometry']=geometry;all_hvg.append(hv)
        selected=ids[c.var.highly_variable.to_numpy()]
        assert len(selected)==2000 and 'PGAM5' not in genes[selected]
        x=log[:,selected].toarray();mu=x.mean(axis=0);sd=np.maximum(x.std(axis=0,ddof=1),.1)
        scaled=np.clip((x-mu)/sd,-10,10)
        pca=PCA(n_components=30,svd_solver='randomized',random_state=0);pc=pca.fit_transform(scaled)
        b=graph_from_pc(pc,obs);lab=cluster(b,.8,0,'state')
        primary_labels[geometry]=lab;graphs[geometry]=b
        obs[geometry+'_cluster']=lab
        for j in range(30):obs[geometry+'_PC'+str(j+1)]=pc[:,j]
        pd.DataFrame(pca.components_.T,index=genes[selected],columns=['PC'+str(j+1) for j in range(30)]).to_csv(R/(geometry+'_PCA_loadings.csv.gz'))
        np.savez_compressed(R/(geometry+'_PCA_model.npz'),genes=genes[selected],mean=mu,SD=sd,
                            pca_mean=pca.mean_,components=pca.components_,cell_ids=obs.cell_id.to_numpy(dtype=str),PC=pc)
        sparse.save_npz(R/(geometry+'_neighbor_graph.npz'),b.obsp['connectivities'])
        features.update(genes[selected])
        at,st=association(obs,lab,geometry);associations.append(at);sample_associations.append(st)
        sil=silhouette_samples(pc,obs.PGAM5_detected.to_numpy(int))
        silhouette=pd.DataFrame({'cell_id':obs.cell_id,'geometry':geometry,'PGAM5_detected':obs.PGAM5_detected,'silhouette_PGAM5_RNA_label':sil})
        silhouette.to_csv(R/(geometry+'_PGAM5_label_silhouette.csv.gz'),index=False)
        print('GEOMETRY',geometry,'clusters',len(set(lab)),'maxPGdetection',at.PGAM5_detection_fraction.max(),'maxPGcapture',at.capture_of_all_PGAM5_detected.max(),flush=True)
        for seed in [1,2]:
            alt=cluster(b,.8,seed,'seed'+str(seed))
            obs[geometry+'_seed'+str(seed)+'_cluster']=alt
            for row in match_clusters(lab,alt):stability.append({'geometry':geometry,'method':'graph_seed','replicate':seed,**row})
            global_stability.append({'geometry':geometry,'method':'graph_seed','replicate':seed,'ARI':adjusted_rand_score(lab,alt),'clusters':len(set(alt))})
        for res in [.4,1.2]:
            alt=cluster(b,res,0,'resolution'+str(res));obs[geometry+'_resolution_'+str(res)]=alt
            alt_at,alt_st=association(obs,alt,geometry+'_resolution_'+str(res))
            associations.append(alt_at);sample_associations.append(alt_st)
        for seed in [11,12,13]:
            rng=np.random.default_rng(seed);chosen=[]
            for sample in sorted(obs.sample_name.unique()):
                index=np.flatnonzero(obs.sample_name.to_numpy()==sample)
                chosen.extend(rng.choice(index,int(np.floor(.8*len(index))),replace=False))
            chosen=np.array(sorted(chosen));xs=x[chosen]
            mus=xs.mean(0);sds=np.maximum(xs.std(axis=0,ddof=1),.1)
            qs=PCA(n_components=30,svd_solver='randomized',random_state=seed).fit_transform(np.clip((xs-mus)/sds,-10,10))
            sub=graph_from_pc(qs,obs.iloc[chosen]);alt=cluster(sub,.8,seed,'substate')
            pd.DataFrame({'cell_id':obs.iloc[chosen].cell_id,'cluster':alt}).to_csv(R/f'{geometry}_subsample_{seed}_labels.csv.gz',index=False)
            for row in match_clusters(lab[chosen],alt):stability.append({'geometry':geometry,'method':'80pct_sample_subsampling','replicate':seed,**row})
            global_stability.append({'geometry':geometry,'method':'80pct_sample_subsampling','replicate':seed,'ARI':adjusted_rand_score(lab[chosen],alt),'clusters':len(set(alt))})
            print('ROBUSTNESS',geometry,seed,'ARI',global_stability[-1]['ARI'],flush=True)
        # All-gene descriptive markers use normalized arithmetic means and detection,
        # and are not independent confirmatory tests after feature-based clustering.
        for group in sorted(set(lab),key=int):
            mask=lab==group
            mean=np.asarray(norm[mask].mean(axis=0)).ravel();other=np.asarray(norm[~mask].mean(axis=0)).ravel()
            det=np.asarray((raw[mask]>0).mean(axis=0)).ravel();od=np.asarray((raw[~mask]>0).mean(axis=0)).ravel()
            contrast=np.log2((mean+.05)/(other+.05));score=contrast*np.maximum(det-od,0)
            allowed=np.array([g!='PGAM5' and g not in BACKGROUND and not g.startswith(('MT-','RPL','RPS','IGH','IGK','IGL','HBA','HBB')) for g in genes])&(det>=.1)&(contrast>0)
            order=np.flatnonzero(allowed);order=order[np.argsort(-score[order],kind='stable')[:30]]
            for rank,j in enumerate(order,1):
                markers.append({'geometry':geometry,'cluster':group,'rank':rank,'gene':genes[j],
                    'cluster_mean_CP10k':float(mean[j]),'other_mean_CP10k':float(other[j]),
                    'cluster_detection_fraction':float(det[j]),'other_detection_fraction':float(od[j]),
                    'descriptive_log2FC_CP10k_pseudocount0p05':float(contrast[j]),
                    'cycling_list_gene':genes[j] in cycles,'stress_list_gene':genes[j] in PANELS['Stress']})
                features.add(genes[j])
            for panel in PANELS:
                module_profiles.append({'geometry':geometry,'cluster':group,'panel':panel,
                                        'mean_log1p_CP10k':float(obs.loc[mask,panel+'_module'].mean())})
        model_parameters[geometry]={'HVGs':2000,'PCs':30,'clusters':len(set(lab)),
            'mean_PGAM5_label_silhouette':float(sil.mean()),
            'mean_silhouette_PGAM5_detected':float(sil[obs.PGAM5_detected].mean()),
            'mean_silhouette_PGAM5_undetected':float(sil[~obs.PGAM5_detected].mean()),
            'positive_label_vs_cluster_ARI':float(adjusted_rand_score(lab,obs.PGAM5_detected))}
    obs.to_csv(R/'cell_annotations.csv.gz')
    pd.concat(all_hvg).to_csv(R/'HVG_selection_statistics.csv.gz',index=False)
    association_table=pd.concat(associations,ignore_index=True);association_table.to_csv(R/'cluster_PGAM5_association.csv',index=False)
    pd.concat(sample_associations).to_csv(R/'cluster_PGAM5_by_sample.csv',index=False)
    pd.DataFrame(stability).to_csv(R/'cluster_stability.csv',index=False)
    pd.DataFrame(global_stability).to_csv(R/'global_clustering_stability.csv',index=False)
    pd.DataFrame(markers).to_csv(R/'cluster_descriptive_markers.csv',index=False)
    pd.DataFrame(module_profiles).to_csv(R/'cluster_module_profiles.csv',index=False)
    # Fixed PGAM5-free cycling definition checked separately from unsupervised clusters.
    cy=obs.cycling_flag.to_numpy(bool);pos=obs.PGAM5_detected.to_numpy(bool)
    cy_n=int(cy.sum());cy_p=int((cy&pos).sum());rest_n=len(obs)-cy_n;rest_p=int((~cy&pos).sum())
    odds,cy_pval=fisher_exact([[cy_p,cy_n-cy_p],[rest_p,rest_n-rest_p]])
    cycling_samples=[]
    for sample in sorted(obs.sample_name.unique()):
        s=obs.sample_name.to_numpy()==sample;inside=s&cy;outside=s&~cy
        mean1=float(obs.loc[inside,'PGAM5_CP10k'].mean());mean0=float(obs.loc[outside,'PGAM5_CP10k'].mean())
        cycling_samples.append({'sample_name':sample,'cycling_cells':int(inside.sum()),'other_cells':int(outside.sum()),
            'cycling_PGAM5_detected':int((inside&pos).sum()),'other_PGAM5_detected':int((outside&pos).sum()),
            'cycling_mean_PGAM5_CP10k':mean1,'other_mean_PGAM5_CP10k':mean0,
            'CP10k_ratio_pseudocount0p05':(mean1+.05)/(mean0+.05)})
    pd.DataFrame(cycling_samples).to_csv(R/'fixed_cycling_state_by_sample.csv',index=False)
    cycling_result={'genes':CYCLE_ANCHORS,'rule':'raw detection of >=2of6 genes; PGAM5 not used',
        'cycling_cells':cy_n,'cycling_PGAM5_detected':cy_p,'cycling_PGAM5_detection_fraction':cy_p/cy_n,
        'other_cells':rest_n,'other_PGAM5_detected':rest_p,'other_detection_fraction':rest_p/rest_n,
        'detection_rate_ratio':(cy_p/cy_n)/(rest_p/rest_n),'capture_of_all_PGAM5_detected':cy_p/int(pos.sum()),
        'Fisher_cell_level_p':float(cy_pval),'Fisher_odds_ratio':float(odds),
        'samples_normalized_expression_enriched':sum(row['CP10k_ratio_pseudocount0p05']>1 for row in cycling_samples)}
    # Label-permutation conditions on library membership, not depth; it is descriptive.
    neighbor_stats=[];nulls=[]
    rng=np.random.default_rng(202642)
    for geometry,b in graphs.items():
        W=b.obsp['connectivities'].tocsr();degree=np.asarray(W.sum(axis=1)).ravel()
        def stat(v):return float(np.mean((W@v/np.maximum(degree,1e-12))[v.astype(bool)]))
        observed=stat(pos.astype(float));values=[]
        for replicate in range(500):
            shuffled=pos.copy()
            for sample in sorted(obs.sample_name.unique()):
                index=np.flatnonzero(obs.sample_name.to_numpy()==sample);shuffled[index]=rng.permutation(shuffled[index])
            value=stat(shuffled.astype(float));values.append(value);nulls.append({'geometry':geometry,'replicate':replicate,'positive_weighted_neighbor_fraction':value})
        neighbor_stats.append({'geometry':geometry,'observed_positive_neighbor_fraction':observed,
            'sample_conditioned_null_mean':float(np.mean(values)),'null95_lower':float(np.quantile(values,.025)),
            'null95_upper':float(np.quantile(values,.975)),
            'one_sided_permutation_p':(1+int(np.sum(np.array(values)>=observed)))/501,'permutations':500})
    pd.DataFrame(neighbor_stats).to_csv(R/'PGAM5_graph_concentration.csv',index=False)
    pd.DataFrame(nulls).to_csv(R/'PGAM5_graph_permutation_null.csv',index=False)
    # Portable inputs retain all fitted feature sets and every exported marker/panel gene.
    features.update(['PGAM5']+[g for panel in PANELS.values() for g in panel if g in genes])
    ix=np.flatnonzero(np.isin(genes,sorted(features)))
    sparse.save_npz(R/'portable_raw_feature_counts.npz',raw[:,ix])
    pd.DataFrame({'gene':genes[ix]}).to_csv(R/'portable_feature_genes.csv',index=False)
    pd.DataFrame({'gene':genes,'whole_cohort_detection_fraction':np.asarray((raw>0).mean(0)).ravel(),
                  'whole_cohort_mean_CP10k':np.asarray(norm.mean(0)).ravel()}).to_csv(R/'full_gene_source_summary.csv.gz',index=False)
    dist=obs.PGAM5_counts.value_counts().sort_index();dist.rename_axis('raw_PGAM5_count').rename('cells').to_csv(R/'PGAM5_raw_count_distribution.csv')
    pd.DataFrame([{'sample_name':s,'macrophages':len(t),'PGAM5_detected':int(t.PGAM5_detected.sum()),
       'detection_fraction':float(t.PGAM5_detected.mean()),'cycling_cells':int(t.cycling_flag.sum())}
       for s,t in obs.groupby('sample_name',observed=True)]).to_csv(R/'sample_summary.csv',index=False)
    depth={'PGAM5_detected_median_total_counts':float(obs.loc[pos,'total_counts'].median()),
           'PGAM5_undetected_median_total_counts':float(obs.loc[~pos,'total_counts'].median()),
           'depth_adjusted':False,'statement':'Descriptive only; user previously requested no depth adjustment.'}
    result={'dataset':'GSE202642','HCC_independent_patients':obs.patient_id.nunique(),
        'patient_identity_source':'User confirmed seven labels correspond to seven independent patients',
        'macrophages':len(obs),'PGAM5_RNA_detected':int(pos.sum()),'PGAM5_RNA_undetected':int((~pos).sum()),
        'geometries':model_parameters,'fixed_cycling_state':cycling_result,'depth_diagnostic':depth,
        'portable_input_genes':len(ix),'biological_assessment':'Pending joint interpretation of cluster concentration, multigene markers, patient representation and reproducibility; cycling is not an exclusion criterion.',
        'validated_TCGA_cell_abundance_reference_obtained':False,
        'decision_note':'Interpret all exported cluster concentration, marker, sample, and stability evidence in README; negative binary-label structure is not proof that no PGAM5-dependent biology exists.'}
    (R/'result_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
