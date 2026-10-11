"""Independent raw-input reconstruction and mathematical verification."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE245908_mpl')
import h5py,numpy as np,pandas as pd,scipy.sparse as sp,scanpy as sc,json,hashlib,gzip
from scipy.spatial import cKDTree
from scipy.stats import rankdata
ROOT=Path(__file__).resolve().parent
CHECKS=[]
def check(name,ok,detail=''):
 CHECKS.append({'check':name,'pass':bool(ok),'detail':str(detail)});
 if not ok:raise AssertionError((name,detail))
def close(name,x,y,atol=1e-6):check(name,np.allclose(x,y,atol=atol,rtol=1e-6,equal_nan=True),np.nanmax(np.abs(np.asarray(x)-np.asarray(y))))
def BH(p):
 p=np.asarray(p);order=np.argsort(p);adj=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];out=np.zeros(len(p));out[order]=np.minimum(adj,1);return out

def residual(y,z):
 q,_=np.linalg.qr(np.column_stack([np.ones(len(y)),z]));v=rankdata(np.round(np.asarray(y,dtype=float),12));v-=q@(q.T@v);return v-v.mean()
def cor(x,y):return float(x@y/np.linalg.norm(x)/np.linalg.norm(y))
if __name__=='__main__':
 manifest=json.loads((ROOT/'source_manifest.json').read_text());programs=json.loads((ROOT/'programs.json').read_text());association=pd.read_csv(ROOT/'PGAM5_spatial_associations.csv');moran_table=pd.read_csv(ROOT/'spatial_autocorrelation.csv');dot=pd.read_csv(ROOT/'region_marker_dotplot_values.csv')
 for rec in manifest:
  p=ROOT/rec['path'];check('source_sha_'+p.name,hashlib.sha256(p.read_bytes()).hexdigest()==rec['sha256']);check('source_size_'+p.name,p.stat().st_size==rec['bytes'])
 for s,gsm in [('CHC20','GSM7850822'),('CHC23','GSM7850823')]:
  prefix=ROOT/'inputs'/s/f'{gsm}_{s}_'
  with h5py.File(str(prefix)+'filtered_feature_bc_matrix.h5','r') as f:
   m=f['matrix'];shape=m['shape'][:];X=sp.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=shape).T.tocsr();names=m['features']['name'].asstr()[:];ids=m['features']['id'].asstr()[:];bars=m['barcodes'].asstr()[:]
  totals=np.asarray(X.sum(1)).ravel();ng=np.asarray((X>0).sum(1)).ravel();mt=np.asarray(X[:,np.array([v.startswith('MT-') for v in names])].sum(1)).ravel()/totals*100
  pos=pd.read_csv(str(prefix)+'tissue_positions_list.csv.gz',header=None,index_col=0,names=['barcode','in_tissue','array_row','array_col','pixel_row','pixel_col']).loc[bars]
  qc=pd.read_csv(ROOT/f'{s}_all_spot_QC.csv',index_col=0);check(s+'_unique_barcodes',len(set(bars))==len(bars));check(s+'_all_tissue',pos.in_tissue.eq(1).all());check(s+'_same_QC_order',list(qc.index)==list(bars));close(s+'_UMI',qc.total_counts,totals);close(s+'_genes',qc.n_genes,ng);close(s+'_MT',qc.pct_mt,mt,atol=1e-5)
  pgidx=np.flatnonzero(names=='PGAM5');check(s+'_one_PGAM5_gene',len(pgidx)==1);raw=X[:,pgidx].toarray().ravel();close(s+'_raw_PGAM5_all',qc.PGAM5_raw,raw)
  keep=(totals>=500)&(ng>=200)&(mt<=25);check(s+'_QC_mask',np.array_equal(keep,qc.pass_qc));goodbars=bars[keep];obs=pd.read_csv(ROOT/f'{s}_spot_annotations.csv',index_col=0);check(s+'_retained_IDs',list(obs.index)==list(goodbars));check(s+'_all_named_regions',obs.region_name.notna().all())
  a=sc.read_h5ad(ROOT/f'{s}_spatial_marker_subset.h5ad');check(s+'_portable_order',list(a.obs_names)==list(goodbars));geneidx={gene:i for i,gene in enumerate(ids)};idx=[geneidx[gene] for gene in a.var.gene_ids]
  xr=X[keep][:,idx];close(s+'_portable_raw_counts',a.layers['counts'].toarray(),xr.toarray(),atol=0);logx=np.log1p(xr.toarray()/totals[keep,None]*10000);close(s+'_portable_normalized',a.X.toarray(),logx)
  expectedpg=np.log1p(raw[keep]/totals[keep]*10000);close(s+'_normalized_PGAM5',obs.PGAM5_log1p_CP10k,expectedpg)
  close(s+'_original_pixel_coords',a.obsm['spatial'],pos.loc[goodbars,['pixel_col','pixel_row']].to_numpy(),atol=0)
  low=a.uns['spatial'][s]['images']['lowres'];sf=a.uns['spatial'][s]['scalefactors']['tissue_lowres_scalef'];xy=a.obsm['spatial']*sf;check(s+'_scaled_points_within_image',((xy[:,0]>=0)&(xy[:,0]<low.shape[1])&(xy[:,1]>=0)&(xy[:,1]<low.shape[0])).all())
  for name,genes in programs.items():
   gi=[np.flatnonzero(names==v)[0] for v in genes];score=np.log1p(X[keep][:,gi].toarray()/totals[keep,None]*10000).mean(1);close(s+'_program_'+name,obs[name],score)
  # Hex first-ring geometry independently via Euclidean KD tree, not index offsets.
  coords=np.column_stack([obs.array_col,np.sqrt(3)*obs.array_row]);pairs=cKDTree(coords).query_pairs(2.001);ed=pd.read_csv(ROOT/f'{s}_spatial_edges.csv.gz');ix={v:i for i,v in enumerate(obs.index)};saved={(min(ix[r.barcode_from],ix[r.barcode_to]),max(ix[r.barcode_from],ix[r.barcode_to])) for r in ed.itertuples()};check(s+'_geometric_graph',saved==pairs);check(s+'_directed_symmetric_edges',len(ed)==2*len(pairs))
  row=[ix[v] for v in ed.barcode_from];col=[ix[v] for v in ed.barcode_to];W=sp.csr_matrix((np.ones(len(ed)),(row,col)),shape=(len(obs),len(obs)));degree=np.asarray(W.sum(1)).ravel();close(s+'_degrees',obs.spatial_degree,degree)
  connected=degree>0;b=obs[connected];W=W[connected][:,connected];deg=np.asarray(W.sum(1)).ravel();cov=np.column_stack([rankdata(np.log1p(b.total_counts)),rankdata(b.n_genes),rankdata(b.pct_mt)]);px=residual(b.PGAM5_log1p_CP10k.to_numpy(),cov)
  for r in association[association['sample']==s].itertuples():
   y=b[r.program].to_numpy()
   if r.relationship=='adjacent_spot':y=(W@y)/deg
   rx=rankdata(np.round(b.PGAM5_log1p_CP10k,12));ry=rankdata(np.round(y,12));close(s+'_spearman_'+r.program+'_'+r.relationship,cor(rx-rx.mean(),ry-ry.mean()),r.raw_Spearman_rho,atol=1e-11)
   close(s+'_partial_'+r.program+'_'+r.relationship,cor(px,residual(y,cov)),r.QC_adjusted_partial_rank_rho,atol=1e-10)
  for k in programs:
   yr=(W@b[k].to_numpy())/deg;close(s+'_neighbor_map_'+k,obs.loc[connected,'neighbor_'+k],yr,atol=1e-10)
  for r in moran_table[moran_table['sample']==s].itertuples():
   x=b[r.feature].to_numpy();x=x-x.mean();I=len(x)/W.sum()*(x@(W@x))/(x@x);close(s+'_Moran_'+r.feature,I,r.Moran_I,atol=1e-11)
  coeff=np.load(ROOT/f'{s}_spectral_null_coefficients.npz');rng=np.random.default_rng(int(coeff['seed']));signs=rng.choice([-1.,1.],size=(int(coeff['n_null']),int(coeff['n_spots']))).astype(np.float32);products=coeff['products'];den=coeff['denominators'];null=signs@products.T/den[None,:]
  for i,label in enumerate(coeff['labels']):
   k,rel=label.split('|');r=association[(association['sample']==s)&(association.program==k)&(association.relationship==rel)].iloc[0];observed=products[i].sum()/den[i];close(s+'_spectral_observed_'+label,observed,r.QC_adjusted_partial_rank_rho,atol=1e-11);p=(1+(np.abs(null[:,i])>=abs(observed)).sum())/(1+len(null));close(s+'_spectral_p_'+label,p,r.spatial_spectral_two_sided_p,atol=1e-11)
  profiles=pd.read_csv(ROOT/f'{s}_cluster_profiles.csv')
  for r in profiles.itertuples():
   dd=obs[obs.cluster==r.cluster];check(s+'_cluster_n_'+r.cluster,len(dd)==r.n_spots);close(s+'_cluster_pgmean_'+r.cluster,dd.PGAM5_log1p_CP10k.mean(),r.PGAM5_mean_log1p_CP10k);check(s+'_cluster_detected_'+r.cluster,(dd.PGAM5_raw>0).sum()==r.PGAM5_detected_spots)
  # Every plotted dot statistic independently calculated from the source matrix.
  for r in dot[dot['sample']==s].itertuples():
   gi=np.flatnonzero(names==r.gene)[0];mask=obs.cluster.eq(r.cluster).to_numpy();counts=X[keep][:,gi].toarray().ravel();z=np.log1p(counts/totals[keep]*10000);close(s+'_dot_mean_'+r.cluster+'_'+r.gene,z[mask].mean(),r.mean_log1p_CP10k);close(s+'_dot_detection_'+r.cluster+'_'+r.gene,(counts[mask]>0).mean()*100,r.detection_pct)
 close('association_BH_family48',BH(association.spatial_spectral_two_sided_p),association.spatial_BH_q,atol=1e-12);check('48_association_tests',len(association)==48);close('Moran_BH_family26',BH(moran_table.two_sided_label_permutation_p),moran_table.Moran_BH_q,atol=1e-12)
 preservation=pd.read_csv(ROOT/'spectral_null_numerical_checks.csv');check('null_preserves_variance',preservation.norm_relative_error.max()<1e-11);check('null_preserves_MoranI',preservation.moran_absolute_error.max()<1e-11)
 check('all_region_markers_enriched',(pd.read_csv(ROOT/'region_naming_marker_evidence.csv').mean_log1p_CP10k_in_region>pd.read_csv(ROOT/'region_naming_marker_evidence.csv').mean_log1p_CP10k_outside).all())
 pd.DataFrame(CHECKS).to_csv(ROOT/'independent_audit_checks.csv',index=False);print(f'Passed {len(CHECKS)} independent checks')
