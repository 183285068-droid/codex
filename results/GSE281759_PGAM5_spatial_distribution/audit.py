"""Independent KD-tree lattice, Dijkstra boundary, DFS components and null p checks."""
import json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.stats import false_discovery_control
from analyze import R,BASE,SAMPLES,METRICS,DIRECTIONS
checks=[]
def check(name,v):
 assert v,name
 checks.append(name)
def coords(o):
 xy=o[['array_col','array_row']].to_numpy(float);xy[:,0]/=2;xy[:,1]*=np.sqrt(3)/2;return xy
def independent(h,xy,pairs,tiles,bd):
 ix=np.flatnonzero(h);selected=set(ix);chosen=[(int(i),int(j)) for i,j in pairs if i in selected and j in selected];adj={i:[] for i in ix}
 for i,j in chosen:adj[i].append(j);adj[j].append(i)
 remaining=set(ix);sizes=[]
 while remaining:
  start=remaining.pop();stack=[start];n=1
  while stack:
   node=stack.pop()
   for nxt in adj[node]:
    if nxt in remaining:remaining.remove(nxt);stack.append(nxt);n+=1
  sizes.append(n)
 nn=cKDTree(xy[ix]).query(xy[ix],k=2)[0][:,1];n1=len(ix);n0=len(h)-n1;chi=0
 for tile in np.unique(tiles):
  count=int((tiles==tile).sum());high=int((h&(tiles==tile)).sum());e=count*n1/len(h);chi+=(high-e)**2/e+((count-high)-count*n0/len(h))**2/(count*n0/len(h))
 return np.array([len(chosen),np.median(nn),max(sizes)/n1,chi,np.median(bd[h])-np.median(bd[~h])]),sizes
def run():
 stats=pd.read_csv(R/'spatial_distribution_statistics.csv');summary=pd.read_csv(R/'section_summary.csv')
 for v in json.loads((R/'reference_manifest.json').read_text()):check(v['path']+' SHA',hashlib.sha256((R.parent/v['path']).read_bytes()).hexdigest()==v['sha256'])
 for item in SAMPLES:
  name=item['sample'];a=ad.read_h5ad(BASE/f'{name}_spatial_marker_subset.h5ad');o=a.obs;xy=coords(o);q=pd.read_csv(BASE/f'{name}_all_spot_QC.csv').query('in_tissue==1');pre=coords(q);pair=cKDTree(pre).query_pairs(1.000001,output_type='ndarray');rows=np.r_[pair[:,0],pair[:,1]];cols=np.r_[pair[:,1],pair[:,0]];G=coo_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(q),len(q))).tocsr();boundary=np.asarray(G.sum(1)).ravel()<6;dist=dijkstra(G,directed=False,indices=np.flatnonzero(boundary),min_only=True);proxy=pd.Series(dist,index=q.barcode);annotation=pd.read_csv(R/f'{name}_annotations.csv.gz',index_col=0).loc[o.index];check(name+' Dijkstra boundary proxy',np.array_equal(proxy.loc[o.index].to_numpy(),annotation.matrix_mask_boundary_hops))
  whole_pairs=cKDTree(xy).query_pairs(1.000001,output_type='ndarray');e=a.obsp['spatial_connectivities'].tocoo();expected=set(map(tuple,whole_pairs))|set(map(tuple,whole_pairs[:,::-1]));check(name+' independent QC lattice',expected==set(zip(e.row,e.col)))
  mac=(o.Macrophage_core>=o.Macrophage_core.quantile(.75))&(o.Macrophage_core_n_detected>=3)&(o.Pan_myeloid>=o.Pan_myeloid.median());positive=o.PGAM5_raw.gt(0);high=(mac&positive&o.PGAM5_log1p_CP10k.ge(o.loc[mac&positive,'PGAM5_log1p_CP10k'].quantile(.75))).to_numpy();check(name+' source high gate',np.array_equal(high,o.PGAM5_mac_group.eq('PGAM5_high')))
  for kind in ['quality_stratified','spatial_quality_stratified']:
   z=np.load(R/f'{name}_{kind}_audit.npz');ix=o.index.get_indexer(z['barcodes']);check(name+'/'+kind+' parent all detected Mac-RNA',set(z['barcodes'])==set(o.index[mac&positive]));check(name+'/'+kind+' labels',np.array_equal(z['high'],high[ix]));local=cKDTree(xy[ix]).query_pairs(1.000001,output_type='ndarray');obs,sizes=independent(z['high'],xy[ix],local,z['tiles'],z['boundary']);check(name+'/'+kind+' independent five metrics',np.allclose(obs,z['observed'],atol=1e-10));check(name+'/'+kind+' boundary vector',np.array_equal(proxy.loc[o.index[ix]],z['boundary']))
   for j,h in enumerate(z['example_labels']):
    for group in np.unique(z['groups']):check(name+'/'+kind+f'/example{j}/stratum{group}',h[z['groups']==group].sum()==z['high'][z['groups']==group].sum())
    metric,_=independent(h,xy[ix],local,z['tiles'],z['boundary']);check(name+'/'+kind+f'/example{j} five null metrics',np.allclose(metric,z['null'][j],atol=1e-10))
   records=stats[(stats['sample']==name)&(stats.null_kind==kind)].set_index('metric')
   for j,k in enumerate(METRICS):
    row=records.loc[k];null=z['null'][:,j];check(name+'/'+kind+'/'+k+' null quantiles',np.allclose([row.null_q025,row.null_q975],np.quantile(null,[.025,.975])))
    if row.status=='tested':
     u=(1+(null>=obs[j]-1e-10).sum())/5000;l=(1+(null<=obs[j]+1e-10).sum())/5000;p=u if DIRECTIONS[j]=='greater' else l if DIRECTIONS[j]=='less' else min(1,2*min(u,l));check(name+'/'+kind+'/'+k+' permutation p',abs(row.permutation_p-p)<1e-12)
   sr=summary[summary['sample']==name].iloc[0];check(name+' summary singletons/components',sr.singleton_high_spots==sum(v==1 for v in sizes) and sr.largest_high_component==max(sizes))
 for col,mask in [('joint_q_160',np.ones(len(stats),bool)),('primary_q_80',stats.null_kind.eq('quality_stratified')),('primary_only_q_75',stats.null_kind.eq('quality_stratified')&stats.tumor_status.eq('primary'))]:check(col+' BH',np.allclose(false_discovery_control(stats.loc[mask,'permutation_p'].fillna(1).to_numpy()),stats.loc[mask,col]))
 (R/'independent_audit_results.json').write_text(json.dumps(dict(passed=len(checks),checks=checks),indent=2));print('PASS',len(checks))
if __name__=='__main__':run()
