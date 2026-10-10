"""Extract only PGAM5>0 AND MYO19>0 macrophages from frozen prior annotations."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
import anndata as ad

R=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--cohorts',type=Path,default=R/'cohorts.json')
args=parser.parse_args()
specs=json.loads(args.cohorts.read_text())
summaries=[];cells=[];provenance=[]

def first_column(obs,names,default):
    for name in names:
        if name in obs:return obs[name].astype(str).to_numpy()
    return np.asarray(default,dtype=str)

for spec in specs:
    path=Path(spec['path']);a=ad.read_h5ad(path,backed='r')
    mask=a.obs.harmonized_group.astype(str).eq('Macrophage').to_numpy() if 'harmonized_group' in a.obs else np.ones(a.n_obs,bool)
    assert mask.sum()==spec['expected_macrophages']
    labels=a.var.gene.astype(str).to_numpy() if 'gene' in a.var else a.var_names.astype(str).to_numpy()
    locations={g:np.flatnonzero(labels==g) for g in ['PGAM5','MYO19']}
    assert all(len(v)>0 for v in locations.values())
    indices=np.unique(np.concatenate(list(locations.values())))
    b=a[mask,indices].to_memory()
    x=b.X.toarray() if hasattr(b.X,'toarray') else np.asarray(b.X)
    counts=np.column_stack([x[:,labels[indices]==g].sum(axis=1) for g in ['PGAM5','MYO19']]).astype(float)
    assert np.isfinite(counts).all() and (counts>=0).all() and (counts==np.floor(counts)).all()
    assert np.array_equal(counts[:,0],b.obs.PGAM5_counts.to_numpy())
    total=b.obs.total_counts.to_numpy(dtype=float)
    assert (total>0).all()
    ids=first_column(b.obs,['cell_id','cell_key'],b.obs_names.astype(str))
    assert len(set(ids))==len(ids)
    donor=first_column(b.obs,['donor_id','Donor','patient','patient_proxy','paper_sample','Sample'],b.obs_names.astype(str))
    sample=first_column(b.obs,['Sample','sample','sample_name','paper_sample'],donor)
    histology=first_column(b.obs,['Histology','diagnosis','tissue'],['HCC_scope_from_prior_analysis']*len(ids))
    if 'annotation_path' in spec:
        anno=pd.read_csv(spec['annotation_path'])
        anno=anno[anno.harmonized_group.eq('Macrophage')].set_index('cell_id')
        assert anno.index.is_unique and set(anno.index)==set(ids)
        assert np.array_equal(anno.loc[ids,'PGAM5_counts'].to_numpy(),counts[:,0])
    co=(counts[:,0]>0)&(counts[:,1]>0)
    assert int(co.sum())==spec['expected_codetected']
    d=pd.DataFrame({'cohort':spec['cohort'],'dataset':spec['dataset'],'platform':spec['platform'],'count_type':spec['count_type'],'cell_id':ids[co],'source_obs_row':np.flatnonzero(mask)[co],'source_obs_name':b.obs_names.astype(str).to_numpy()[co],'donor_label':donor[co],'sample_label':sample[co],'histology':histology[co],'total_counts':total[co],'PGAM5_raw_count':counts[co,0].astype(int),'MYO19_raw_count':counts[co,1].astype(int)})
    for g in ['PGAM5','MYO19']:
        d[g+'_CP10k']=d[g+'_raw_count']/d.total_counts*10000
        d[g+'_log1p_CP10k']=np.log1p(d[g+'_CP10k'])
    cells.append(d)
    summaries.append({'cohort':spec['cohort'],'dataset':spec['dataset'],'platform':spec['platform'],'count_type':spec['count_type'],'scope':spec['scope'],'QC_macrophages':len(b),'PGAM5_detected':int((counts[:,0]>0).sum()),'MYO19_detected':int((counts[:,1]>0).sum()),'both_detected':int(co.sum()),'both_detected_percent_of_macrophages':float(co.mean()*100),'co_detected_donor_labels':len(set(donor[co])),'co_detected_sample_libraries':len(set(sample[co])),'co_detected_histology':json.dumps(d.histology.value_counts().to_dict()),'status':'eligible_n>=3' if co.sum()>=3 else 'not_estimable_n<3'})
    inputs=[path]+([Path(spec['annotation_path'])] if 'annotation_path' in spec else [])
    provenance.append({'cohort':spec['cohort'],'source_matrix_shape':list(a.shape),'macrophage_obs_rows':int(mask.sum()),'target_gene_column_indices':{g:locations[g].tolist() for g in locations},'target_gene_column_names':{g:a.var_names[locations[g]].astype(str).tolist() for g in locations},'input_hashes':{str(p):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in inputs}})
    a.file.close()
    print(spec['cohort'],'macrophages',len(b),'both detected',int(co.sum()),flush=True)
allcells=pd.concat(cells,ignore_index=True)
assert allcells[['PGAM5_raw_count','MYO19_raw_count']].gt(0).all().all()
assert not allcells[['cohort','cell_id']].duplicated().any()
allcells.to_csv(R/'codetected_macrophage_expression.csv.gz',index=False)
pd.DataFrame(summaries).to_csv(R/'cohort_detection_summary.csv',index=False)
(R/'source_provenance.json').write_text(json.dumps({'cohorts':provenance,'selected_cells':len(allcells),'selection':'Strict PGAM5 raw symbol-summed count>0 AND MYO19 raw symbol-summed count>0, among prior QC macrophages','normalization':'Preserved full-library total_counts denominator; CP10k and natural log1p. Read counts remain read counts, not UMI','latest_annotations':'Five harmonized/marker-consistent v4 cohorts; original prior author/marker labels for125449,140228,146115','excluded_dataset':'GSE154906 was cancelled and is not restarted','no_pooled_independent_test':'GSE125449 overlaps GSE151530; GSE140228 platforms share donors; cohorts/platforms kept separate','source_expression_sha256':hashlib.sha256((R/'codetected_macrophage_expression.csv.gz').read_bytes()).hexdigest()},indent=2))
