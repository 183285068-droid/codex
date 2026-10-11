"""Apply reviewed RNA labels and add fixed no-CTLA4/TIGIT exhaustion sensitivity panel."""
import json
import numpy as np,pandas as pd,anndata as ad
from spatial_core import ROOT,WORK,PROGRAMS
rules=json.loads((ROOT/'annotation_rules.json').read_text());records=[]
for assay,path in [('singlecell',ROOT/'singlecell_marker_subset.h5ad'),('spatial',ROOT/'NYU_LIHC1_Vis_spatial_marker_subset.h5ad')]:
    a=ad.read_h5ad(path);mapping=rules[assay]
    assert set(a.obs.expression_cluster)==set(mapping)
    labels={c:c+': '+v[0]+' | '+v[1] for c,v in mapping.items()};a.obs['RNA_annotation']=a.obs.expression_cluster.map(labels).astype(str)
    for p,genes in PROGRAMS.items():
        if p not in a.obs:
            present=[g for g in genes if g in a.var_names];a.obs[p]=a[:,present].X.toarray().mean(1);a.obs[p+'_n_detected']=(a[:,present].layers['counts'].toarray()>0).sum(1)
    if assay=='singlecell':
        a.obs['Mac_RNA_cluster_confirmed']=a.obs.Mac_RNA_supported&a.obs.expression_cluster.isin(['C3','C6'])
    a.write_h5ad(path,compression='gzip')
    a.obs.to_csv(ROOT/('singlecell_annotations.csv' if assay=='singlecell' else 'NYU_LIHC1_Vis_spot_annotations.csv'),float_format='%.17g')
    profile=ROOT/('singlecell_cluster_profiles.csv' if assay=='singlecell' else 'NYU_LIHC1_Vis_cluster_profiles.csv');d=pd.read_csv(profile)
    d['RNA_annotation']=d.cluster.map(labels)
    for p in PROGRAMS:d[p]=d.cluster.map(a.obs.groupby('expression_cluster',observed=True)[p].mean())
    d.to_csv(profile,index=False)
    for c,v in mapping.items():records.append(dict(assay=assay,cluster=c,marker_pair=v[0],RNA_annotation=v[1],n_observations=int(a.obs.expression_cluster.eq(c).sum()),provisional=True))
    full=WORK/('singlecell_full_processed.h5ad' if assay=='singlecell' else 'NYU_LIHC1_Vis_full_processed.h5ad')
    if full.exists():
        b=ad.read_h5ad(full);b.obs=a.obs.copy();b.write_h5ad(full,compression='gzip')
    coverage=[]
    for p,genes in PROGRAMS.items():
        present=[g for g in genes if g in a.var.original_gene_symbol.to_numpy()];coverage.append(dict(program=p,requested=len(genes),present=len(present),present_genes=';'.join(present),missing_genes=';'.join(g for g in genes if g not in present)))
    pd.DataFrame(coverage).to_csv(ROOT/(assay+'_program_gene_coverage_final.csv'),index=False)
pd.DataFrame(records).to_csv(ROOT/'RNA_group_annotations.csv',index=False)
