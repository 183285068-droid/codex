"""Spatial RNA pattern statistics and descriptive macrophage/T candidate edges."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import false_discovery_control
from spatial_core import ROOT, SAMPLES, PROGRAMS, SEED

PAIRS=pd.read_csv(ROOT/'selected_LR_pairs.csv')

def moran(a,item):
    z=a.obs.PGAM5_log1p_CP10k.to_numpy(dtype=float)
    z=z-z.mean()
    graph=a.obsp['spatial_connectivities'].astype(float)
    coo=graph.tocoo();i,j=coo.row,coo.col
    scale=len(a)/(graph.sum()*np.dot(z,z))
    observed=float(scale*np.dot(z,graph@z))
    bins=pd.qcut(a.obs.total_counts,5,labels=False,duplicates='drop').astype(int).to_numpy()*3+pd.qcut(a.obs.pct_mt,3,labels=False,duplicates='drop').astype(int).to_numpy()
    row=[]
    rng=np.random.default_rng(SEED+int(item['gsm'][3:]))
    for method in ['unrestricted','UMI_MT_stratified']:
        groups=np.zeros(len(a),dtype=int) if method=='unrestricted' else bins
        unique=np.unique(groups);ids=[np.flatnonzero(groups==g) for g in unique]
        means=np.zeros(len(a));correction=np.zeros(len(i))
        for index in ids:
            means[index]=z[index].mean()
            if len(index)>1:
                same=(groups[i]==groups[index[0]])&(groups[j]==groups[index[0]])
                correction[same]=-np.var(z[index])/(len(index)-1)
        expected=float(scale*(np.dot(means,graph@means)+correction.sum()))
        null=[];examples=[]
        for start in range(0,4999,128):
            b=min(128,4999-start)
            permutations=np.tile(np.arange(len(a)),(b,1))
            for index in ids:permutations[:,index]=rng.permuted(permutations[:,index],axis=1)
            x=z[permutations]
            values=scale*np.sum(x*(graph@x.T).T,axis=1)
            null.extend(values.tolist())
            if start==0:examples=permutations[:3].copy()
        null=np.asarray(null)
        p=(1+np.sum(np.abs(null-expected)>=np.abs(observed-expected)-1e-12))/5000
        row.append(dict(**item,method=method,n_spots=len(a),directed_edges=graph.nnz,isolated_spots=int(a.obs.spatial_degree.eq(0).sum()),
            Moran_I=observed,conditional_expected_I=expected,null_mean_I=null.mean(),permutation_p=p,permutations=4999,
            strata=len(ids),permutable_spots=sum(len(index) for index in ids if len(index)>1)))
        np.savez_compressed(ROOT/f'{item["sample"]}_{method}_Moran_audit.npz',z=z,groups=groups,null_I=null,
                            example_permutations=examples,observed=observed,expected=expected,scale=scale)
    return row

def communication(a,item):
    graph=a.obsp['spatial_connectivities'].tocoo()
    mac=a.obs.Mac_RNA_enriched.to_numpy();t=a.obs.T_RNA_enriched.to_numpy()
    choose=mac[graph.row]&t[graph.col]
    sources,targets=graph.row[choose],graph.col[choose]
    edges=pd.DataFrame(dict(source_barcode=a.obs_names[sources],target_barcode=a.obs_names[targets],
        source_PGAM5_group=a.obs.PGAM5_mac_group.iloc[sources].to_numpy(),
        source_PGAM5_raw=a.obs.PGAM5_raw.iloc[sources].to_numpy(),
        source_PGAM5_log1p_CP10k=a.obs.PGAM5_log1p_CP10k.iloc[sources].to_numpy(),
        target_activation_supported=a.obs.Activation_early_supported.iloc[targets].to_numpy(),
        target_exhaustion_supported=a.obs.Exhaustion_core_supported.iloc[targets].to_numpy()))
    expression={g:np.expm1(a[:,[g]].X.toarray().ravel()).astype(float) for g in a.var_names}
    columns=[]
    for pair in PAIRS.itertuples():
        receptor=pair.receptor_complex.split('_')
        name=pair.ligand+'--'+pair.receptor_complex
        assert pair.ligand in expression and all(g in expression for g in receptor), name
        lig=expression[pair.ligand][sources]
        rec=np.min(np.vstack([expression[g][targets] for g in receptor]),axis=0)
        edges[name]=np.sqrt(lig*rec)
        columns.append(name)
    edges.to_csv(ROOT/f'{item["sample"]}_mac_T_candidate_edges.csv.gz',index=False)
    summaries=[]
    for group in ['All_mac_sources','PGAM5_high','Other_detected','Undetected']:
        subset=edges if group=='All_mac_sources' else edges[edges.source_PGAM5_group.eq(group)]
        for state in ['All_T_targets','Activation_supported','Exhaustion_supported']:
            data=subset if state=='All_T_targets' else subset[subset['target_'+('activation_supported' if state=='Activation_supported' else 'exhaustion_supported')]]
            for pair in columns:
                supported=data[pair].gt(0)
                summaries.append(dict(**item,source_group=group,target_state=state,pair=pair,
                    eligible_edges=len(data),unique_sources=data.source_barcode.nunique(),unique_targets=data.target_barcode.nunique(),
                    supported_edges=int(supported.sum()),supported_unique_sources=data.loc[supported,'source_barcode'].nunique(),
                    supported_unique_targets=data.loc[supported,'target_barcode'].nunique(),
                    edge_codetection_pct=100*supported.mean() if len(data) else np.nan))
    source_summary=[]
    for name,mask in [('All_mac_RNA',mac),('Mac_with_T_neighbor',mac&a.obs.T_neighbor_count.gt(0)),
                      ('Mac_with_Tex_neighbor',mac&a.obs.T_exhaustion_neighbor_count.gt(0)),
                      ('Mac_with_activation_neighbor',mac&a.obs.T_activation_neighbor_count.gt(0))]:
        d=a.obs.loc[mask]
        source_summary.append(dict(**item,region=name,n_source_spots=len(d),PGAM5_detected=int(d.PGAM5_raw.gt(0).sum()),
            PGAM5_detection_pct=100*d.PGAM5_raw.gt(0).mean(),PGAM5_mean_log1p_CP10k=d.PGAM5_log1p_CP10k.mean()))
    return summaries,source_summary

if __name__=='__main__':
    rows=[];lr=[];regions=[]
    for item in SAMPLES:
        a=ad.read_h5ad(ROOT/f'{item["sample"]}_spatial_marker_subset.h5ad')
        rows.extend(moran(a,item))
        l,r=communication(a,item);lr.extend(l);regions.extend(r)
        print('Spatial statistics:',item['sample'],flush=True)
    table=pd.DataFrame(rows);table['BH_q_2_tests']=false_discovery_control(table.permutation_p.to_numpy())
    table.to_csv(ROOT/'PGAM5_spatial_autocorrelation.csv',index=False)
    pd.DataFrame(lr).to_csv(ROOT/'candidate_LR_codetection_summary.csv',index=False)
    pd.DataFrame(regions).to_csv(ROOT/'PGAM5_T_neighbor_source_summary.csv',index=False)
