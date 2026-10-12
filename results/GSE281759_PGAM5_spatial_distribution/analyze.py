"""Random-label spatial distribution of fixed PGAM5-high macrophage-RNA spots."""
from pathlib import Path
from collections import deque
import json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import connected_components
from scipy.stats import false_discovery_control
R=Path(__file__).resolve().parent;BASE=R.parent/'GSE281759_spatial_HCC_PGAM5';SAMPLES=json.loads((BASE/'sample_scope.json').read_text());B=4999;SEED=281759602
METRICS=['high_high_undirected_edges','median_nearest_high_spacing','largest_component_fraction','regional_heterogeneity','boundary_proxy_median_difference']
DIRECTIONS=['greater','less','greater','greater','two_sided_equal_tail']
def boundary_proxy(q):
 q=q[q.in_tissue.eq(1)];coords=q[['array_row','array_col']].to_numpy(int);lookup={tuple(v):i for i,v in enumerate(coords)};neighbors=[]
 for row,col in coords:neighbors.append([lookup[(row+dr,col+dc)] for dr,dc in [(0,2),(0,-2),(1,1),(1,-1),(-1,1),(-1,-1)] if (row+dr,col+dc) in lookup])
 degree=np.array([len(v) for v in neighbors]);distance=np.full(len(q),-1,int);queue=deque(np.flatnonzero(degree<6));distance[degree<6]=0
 while queue:
  i=queue.popleft()
  for j in neighbors[i]:
   if distance[j]<0:distance[j]=distance[i]+1;queue.append(j)
 assert (distance>=0).all()
 return pd.Series(distance,index=q.barcode),degree
def values(h,W,D,tiles,boundary):
 ix=np.flatnonzero(h);sub=W[ix][:,ix];join=float(sub.nnz/2);n=len(ix)
 nn=float(np.median(D[np.ix_(ix,ix)].min(1))) if n>1 else np.nan
 num,labels=connected_components(sub,directed=False);largest=float(np.bincount(labels).max()/n) if n else np.nan
 size=np.bincount(tiles,minlength=9);count=np.bincount(tiles[ix],minlength=9);p=n/len(h);valid=size>0;hetero=float(np.sum((count[valid]-size[valid]*p)**2/(size[valid]*p*(1-p))))
 edge_delta=float(np.median(boundary[h])-np.median(boundary[~h]));return np.array([join,nn,largest,hetero,edge_delta])
def run():
 config=dict(seed=SEED,permutations=B,metrics=dict(zip(METRICS,DIRECTIONS)),primary_null='label shuffle within global source UMI half x macrophage-core half; detected Mac-RNA candidates only',sensitivity_null='same quality strata x full-section 3x3 spatial bins',primary_tests=80,joint_tests=160,coordinate_units='nearest grid spacing; array x=col/2, y=sqrt(3)*row/2',boundary='shortest lattice hops to degree<6 boundary of pre-QC matrix-available in_tissue coordinates; not histological tumor margin')
 (R/'analysis_config.json').write_text(json.dumps(config,indent=2));rows=[];summ=[];refs=[];region_rows=[];checks=[]
 for item in SAMPLES:
  name=item['sample'];path=BASE/f'{name}_spatial_marker_subset.h5ad';qpath=BASE/f'{name}_all_spot_QC.csv'
  for p in [path,qpath]:refs.append(dict(path=str(p.relative_to(R.parent)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
  a=ad.read_h5ad(path);o=a.obs.copy();xy=o[['array_col','array_row']].to_numpy(float);xy[:,1]*=np.sqrt(3);xy/=2
  parent=o.Mac_RNA_enriched.to_numpy(bool)&o.PGAM5_raw.to_numpy().astype(bool);index=np.flatnonzero(parent);high=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy();h=high[index];n1=int(h.sum());n0=len(h)-n1
  W=a.obsp['spatial_connectivities'][index][:,index].tocsr();checkedges=W-W.T;assert checkedges.nnz==0 and W.diagonal().sum()==0;checks.append(name+' symmetric no-self graph')
  D=cdist(xy[index],xy[index]);np.fill_diagonal(D,np.inf)
  q=pd.read_csv(qpath);proxy,degree=boundary_proxy(q);bd=proxy.loc[o.index[index]].to_numpy();o['matrix_mask_boundary_hops']=proxy.loc[o.index].to_numpy()
  arr=o[['array_col','array_row']].to_numpy();bins=np.minimum(((arr-arr.min(0))/(np.ptp(arr,axis=0)+1)*3).astype(int),2);tile=bins[:,0]+3*bins[:,1];tiles=tile[index];depth=o.total_counts.to_numpy()[index];core=o.Macrophage_core.to_numpy()[index];quality=(depth>np.median(depth)).astype(int)+2*(core>np.median(core)).astype(int)
  observed=values(h,W,D,tiles,bd);full_sub=a.obsp['spatial_connectivities'][high][:,high];_,labels=connected_components(full_sub,directed=False);sizes=np.bincount(labels);nn=D[np.ix_(np.flatnonzero(h),np.flatnonzero(h))].min(1);nearhigh=np.asarray(a.obsp['spatial_connectivities']@high.astype(int)).ravel()>0
  summ.append(dict(**item,qc_spots=len(o),matrix_preQC_in_tissue=int(q.in_tissue.eq(1).sum()),matrix_boundary_spots=int((degree<6).sum()),mac_spots=int(o.Mac_RNA_enriched.sum()),n_detected=len(index),n_high=n1,n_other_detected=n0,high_pct_qc=100*n1/len(o),high_pct_mac=100*n1/o.Mac_RNA_enriched.sum(),high_with_high_neighbor=int(nearhigh[high].sum()),high_high_edges=observed[0],high_components=len(sizes),singleton_high_spots=int((sizes==1).sum()),largest_high_component=int(sizes.max()),largest_component_fraction=observed[2],median_nearest_high_spacing=observed[1],median_high_boundary_hops=float(np.median(bd[h])),median_other_boundary_hops=float(np.median(bd[~h])),boundary_delta=observed[4]))
  for t in range(9):region_rows.append(dict(**item,tile=t,n_qc=int((tile==t).sum()),n_mac=int((o.Mac_RNA_enriched.to_numpy(bool)&(tile==t)).sum()),n_detected=int((tiles==t).sum()),n_high=int((h&(tiles==t)).sum())))
  o['high_component_id']=-1;o.loc[high,'high_component_id']=labels;o['nearest_high_spacing']=np.nan;o.loc[high,'nearest_high_spacing']=nn;o['source_candidate_detected']=parent;o['spatial_bin']=tile;o.to_csv(R/f'{name}_annotations.csv.gz',float_format='%.17g')
  rng=np.random.default_rng(SEED+int(item['gsm'][3:]))
  for null_kind,codes in [('quality_stratified',quality),('spatial_quality_stratified',quality+4*tiles)]:
   groups=[np.flatnonzero(codes==v) for v in np.unique(codes)];null=np.empty((B,len(METRICS)));examples=[]
   for b in range(B):
    hp=h.copy()
    for g in groups:hp[g]=rng.permutation(hp[g])
    null[b]=values(hp,W,D,tiles,bd)
    if b<3:examples.append(hp.copy())
   np.savez_compressed(R/f'{name}_{null_kind}_audit.npz',barcodes=o.index[index].to_numpy(str),high=h,groups=codes,observed=observed,null=null,example_labels=np.asarray(examples),tiles=tiles,boundary=bd)
   for j,k in enumerate(METRICS):
    usable=n1>=5 and n0>=10 and np.ptp(null[:,j])>1e-10
    if usable:
     upper=(1+(null[:,j]>=observed[j]-1e-10).sum())/(B+1);lower=(1+(null[:,j]<=observed[j]+1e-10).sum())/(B+1);p=upper if DIRECTIONS[j]=='greater' else lower if DIRECTIONS[j]=='less' else min(1,2*min(upper,lower))
    else:p=np.nan
    rows.append(dict(**item,null_kind=null_kind,metric=k,direction=DIRECTIONS[j],n_high=n1,n_other_detected=n0,observed=observed[j],null_mean=float(null[:,j].mean()),null_q025=float(np.quantile(null[:,j],.025)),null_q975=float(np.quantile(null[:,j],.975)),permutation_p=p,status='tested' if usable else 'not_testable',reason='' if usable else 'small groups or constant conditional null'))
  print(name,'high',n1,'joins',observed[0],'NN',round(observed[1],2),'boundary delta',observed[4],flush=True)
 st=pd.DataFrame(rows);st['joint_q_160']=false_discovery_control(st.permutation_p.fillna(1).to_numpy());st['primary_q_80']=np.nan;sel=st.null_kind.eq('quality_stratified');st.loc[sel,'primary_q_80']=false_discovery_control(st.loc[sel,'permutation_p'].fillna(1).to_numpy());st['primary_only_q_75']=np.nan;sel2=sel&st.tumor_status.eq('primary');st.loc[sel2,'primary_only_q_75']=false_discovery_control(st.loc[sel2,'permutation_p'].fillna(1).to_numpy());st.to_csv(R/'spatial_distribution_statistics.csv',index=False);sd=pd.DataFrame(summ);sd.to_csv(R/'section_summary.csv',index=False);pd.DataFrame(region_rows).to_csv(R/'spatial_bin_counts.csv',index=False)
 p=sd.groupby('patient',as_index=False).agg(sections=('sample','size'),qc_spots=('qc_spots','sum'),mac_spots=('mac_spots','sum'),n_high=('n_high','sum'),n_other_detected=('n_other_detected','sum'),high_with_high_neighbor=('high_with_high_neighbor','sum'),high_high_edges=('high_high_edges','sum'),singleton_high_spots=('singleton_high_spots','sum'),median_section_NN=('median_nearest_high_spacing','median'),median_section_boundary_delta=('boundary_delta','median'));p['high_pct_mac']=100*p.n_high/p.mac_spots;p['high_neighbor_pct']=100*p.high_with_high_neighbor/p.n_high;p.to_csv(R/'patient_descriptive.csv',index=False)
 (R/'reference_manifest.json').write_text(json.dumps(refs,indent=2));(R/'audit_results.json').write_text(json.dumps(dict(passed=len(checks),checks=checks),indent=2));print('SIGNIFICANT',st[st.primary_q_80<.05].to_string(index=False));print(p.to_string(index=False))
if __name__=='__main__':run()
