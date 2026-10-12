"""TLS-compatible RNA patch exploration. This does NOT identify pathologically confirmed TLS."""
from pathlib import Path
import json,hashlib,shutil
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.stats import false_discovery_control
ROOT=Path(__file__).resolve().parent;REF=ROOT.parent/'GSE281759_spatial_HCC_PGAM5'
SAMPLES=json.loads((ROOT/'sample_scope.json').read_text());SEED=281759104;B=4999
TIERS=['TB_aggregate','TB_CXCL13_MHCII','Strict_separate_TB_CXCL13']
COMPS=['high_vs_other_detected','high_vs_all_other_mac','nonTB_high_vs_other_detected']
ENDPOINTS=['within_patch','adjacent_outside_patch']
MARKERS={'TLS_chemokines':['CCL19','CCL21','CXCL13','LTB','LTA'],
 'B':['MS4A1','CD79A','CD79B','CD19','CD22'],'T':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2'],
 'Tfh':['CXCR5','BCL6','ICOS','PDCD1','CD40LG'],'FDC':['CR2','FCER2','FDCSP','CXCL13'],
 'HEV':['CHST4','ACKR1','SELE','CCL21'],'MHCII':['HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','CD74']}
def dense(x):return np.asarray(x.toarray() if sparse.issparse(x) else x)
def annotate(a):
 o=a.obs;G=a.obsp['spatial_connectivities'].astype(bool).astype(np.int8);R=(G+G@G+sparse.eye(len(a),format='csr')).astype(bool).astype(np.int8)
 def raw(g):return dense(a[:,[g]].layers['counts']).ravel()
 t=o.T_RNA_enriched.to_numpy(bool);bs=dense(a[:,MARKERS['B']].X).mean(1);bd=(dense(a[:,MARKERS['B']].layers['counts'])>0).sum(1)
 b=(bs>=np.quantile(bs,.75))&(bs>0)&(bd>=2)&((raw('MS4A1')>0)|(raw('CD19')>0));bo=b&~t;to=t&~b;cx=raw('CXCL13')>0
 mh=(dense(a[:,MARKERS['MHCII']].layers['counts'])>0).sum(1)>=3
 nb=np.asarray(R@b.astype(int)).ravel();nt=np.asarray(R@t.astype(int)).ravel();nu=np.asarray(R@(b|t).astype(int)).ravel();nbo=np.asarray(R@bo.astype(int)).ravel();nto=np.asarray(R@to.astype(int)).ravel();ncx=np.asarray(R@cx.astype(int)).ravel();nmh=np.asarray(R@mh.astype(int)).ravel()
 basic=b&(nb>=2)&(nt>=3)&(nu>=5)&(nbo>=1)&(nto>=1)
 # Require a direct contact between distinct B-only/T-only RNA spots INSIDE the local patch.
 contact=np.zeros(len(a),bool)
 for i in np.flatnonzero(basic):
  nodes=R.indices[R.indptr[i]:R.indptr[i+1]];sub=G[nodes][:,nodes].tocoo();contact[i]=bool(np.any(bo[nodes[sub.row]]&to[nodes[sub.col]]))
 basic&=contact
 centers={TIERS[0]:basic,TIERS[1]:basic&(ncx>=1)&(nmh>=1),TIERS[2]:basic&(nbo>=2)&(nto>=3)&(ncx>=2)&(nmh>=1)}
 o['gate_T']=t;o['gate_B']=b;o['B_only']=bo;o['T_only']=to;o['CXCL13_detected']=cx;o['MHCII_3_detected']=mh
 for name,v in [('B_count',nb),('T_count',nt),('union_count',nu),('B_only_count',nbo),('T_only_count',nto),('CXCL13_count',ncx),('MHCII_count',nmh)]:o['patch_'+name]=v
 for k,c in centers.items():
  region=np.asarray(R@c.astype(int)).ravel()>0;adj=(np.asarray(G@region.astype(int)).ravel()>0)&~region
  o['center_'+k]=c;o['region_'+k]=region;o['adjacent_'+k]=adj
 return G,R,centers
def groups(o,ix):
 xy=o[['array_col','array_row']].to_numpy();tile=np.minimum(((xy[ix]-xy.min(0))/(np.ptp(xy,axis=0)+1)*3).astype(int),2);code=tile[:,0]+3*tile[:,1];out=[]
 for c in np.unique(code):
  z=np.flatnonzero(code==c);d=o.total_counts.to_numpy()[ix[z]]>np.median(o.total_counts.to_numpy()[ix[z]])
  for v in [False,True]:
   zz=z[d==v]
   if not len(zz):continue
   m=o.Macrophage_core.to_numpy()[ix[zz]]>np.median(o.Macrophage_core.to_numpy()[ix[zz]])
   for u in [False,True]:
    p=zz[m==u]
    if len(p):out.append(p)
 return out
def evaluate(a,item,tier,endpoint,comp,rng):
 o=a.obs;mac=o.Mac_RNA_enriched.to_numpy(bool)&o.spatial_degree.gt(0).to_numpy();high=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy();eligible=mac.copy();region=o['region_'+tier].to_numpy(bool)
 if comp!='high_vs_all_other_mac':eligible&=o.PGAM5_raw.gt(0).to_numpy()
 if comp.startswith('nonTB'):eligible&=~o.gate_T.to_numpy(bool)&~o.gate_B.to_numpy(bool)
 if endpoint=='adjacent_outside_patch':eligible&=~region
 ix=np.flatnonzero(eligible);h=high[ix].astype(float);y=(region if endpoint=='within_patch' else o['adjacent_'+tier].to_numpy(bool))[ix].astype(float);n1=int(h.sum());n0=len(h)-n1
 base=dict(**item,tier=tier,endpoint=endpoint,comparison=comp,n_high=n1,n_comparator=n0,n_centers=int(o['center_'+tier].sum()),n_patch_spots=int(region.sum()),high_supported=int(y[h>0].sum()),comparator_supported=int(y[h==0].sum()),high_pct=100*y[h>0].mean() if n1 else np.nan,comparator_pct=100*y[h==0].mean() if n0 else np.nan,status='not_testable',reason='',p=np.nan,delta_pp=np.nan)
 if n1<5 or n0<10:base['reason']='high<5 or comparator<10';return base
 strata=groups(o,ix);H=np.tile(h,(B,1));code=np.zeros(len(ix),int)
 for c,g in enumerate(strata):H[:,g]=rng.permuted(H[:,g],axis=1);code[g]=c
 obs=100*(h@y/n1-(1-h)@y/n0);null=100*(H@y/n1-(1-H)@y/n0);base['delta_pp']=obs;base['null_mean_pp']=float(null.mean())
 if np.ptp(null)<1e-10:base['reason']='no candidate / constant conditional outcome'
 else:base.update(status='tested',p=min(1.,2*min((1+(null>=obs-1e-10).sum())/(B+1),(1+(null<=obs+1e-10).sum())/(B+1))))
 np.savez_compressed(ROOT/f'{item["sample"]}_{tier}_{endpoint}_{comp}_audit.npz',barcodes=a.obs_names[ix].to_numpy(str),high=h,outcome=y,groups=code,assignments_first10=H[:10].astype(np.uint8),observed=obs,null=null)
 return base
def run():
 (ROOT/'inputs').mkdir(exist_ok=True)
 config=dict(seed=SEED,permutations=B,tiers=TIERS,comparisons=COMPS,endpoints=ENDPOINTS,planned_tests=288,primary_tests=32,primary='CXCL13/MHCII-supported patches, high vs other detected, both endpoints',candidate='B-gated center; radius-2 hex neighborhood incl center; >=2 B spots, >=3 T spots, >=5 unique B|T spots, >=1 B-only and T-only; direct B-only/T-only contact inside patch',CXCL13_support='>=1 CXCL13-detected spot and >=1 MHCII spot with >=3 of 5 markers detected',strict='>=2 B-only, >=3 T-only, >=2 CXCL13-detected spots; MHCII support',region='union of radius-2 patches around qualifying centers; not measured TLS boundaries',adjacent='first-ring adjacency to region, only outside-region sources; self-edges absent',thresholds_frozen_before_results=True,pathological_TLS_confirmed=False,minimum_groups=[5,10],null='high labels within 3x3 spatial tile x UMI halves x macrophage-core halves; +1 two-sided equal-tail p',limits=['missing canonical TLS/FDC/HEV genes','no TLS histology/IHC annotations','mixed spots cannot attribute PGAM5 to macrophages','two-ring cutoff is exploratory, not validated TLS classifier'])
 (ROOT/'analysis_config.json').write_text(json.dumps(config,indent=2)+'\n');rows=[];desc=[];summary=[];coverage=[];components=[];manifest=[]
 for item in SAMPLES:
  sample=item['sample'];name=sample+'_spatial_marker_subset.h5ad';dst=ROOT/'inputs'/name
  if not dst.exists():shutil.copyfile(REF/name,dst)
  manifest.append(dict(path='inputs/'+name,bytes=dst.stat().st_size,sha256=hashlib.sha256(dst.read_bytes()).hexdigest()));a=ad.read_h5ad(dst);G,R,centers=annotate(a);o=a.obs;high=o.PGAM5_mac_group.eq('PGAM5_high').to_numpy();mac=o.Mac_RNA_enriched.to_numpy(bool)&o.spatial_degree.gt(0).to_numpy()
  for panel,genes in MARKERS.items():
   for g in genes:coverage.append(dict(**item,panel=panel,gene=g,present=g in a.var_names,detected_spots=int((dense(a[:,[g]].layers['counts'])>0).sum()) if g in a.var_names else np.nan))
  rng=np.random.default_rng(SEED+int(item['gsm'][3:]))
  for tier in TIERS:
   region=o['region_'+tier].to_numpy(bool);adj=o['adjacent_'+tier].to_numpy(bool);ix=np.flatnonzero(region);ncomp,labels=connected_components(G[ix][:,ix],directed=False) if len(ix) else (0,np.array([],int))
   summary.append(dict(**item,tier=tier,qc_spots=len(a),B_spots=int(o.gate_B.sum()),T_spots=int(o.gate_T.sum()),centers=int(centers[tier].sum()),patch_spots=int(region.sum()),patch_components=ncomp,high_all=int(high.sum()),high_within=int((high&region).sum()),high_adjacent_outside=int((high&adj).sum())))
   for c in range(ncomp):
    nodes=ix[labels==c];components.append(dict(**item,tier=tier,component_id=c,n_spots=len(nodes),n_B=int(o.gate_B.to_numpy()[nodes].sum()),n_T=int(o.gate_T.to_numpy()[nodes].sum()),n_CXCL13=int(o.CXCL13_detected.to_numpy()[nodes].sum()),n_high=int(high[nodes].sum()),center_x=float(o.pixel_col.to_numpy()[nodes].mean()),center_y=float(o.pixel_row.to_numpy()[nodes].mean())))
   for group in ['PGAM5_high','Other_detected','Undetected']:
    mask=mac&o.PGAM5_mac_group.eq(group).to_numpy();desc.append(dict(**item,tier=tier,source_group=group,n_sources=int(mask.sum()),within_patch=int((mask&region).sum()),adjacent_outside_patch=int((mask&adj).sum())))
   for endpoint in ENDPOINTS:
    for comp in COMPS:rows.append(evaluate(a,item,tier,endpoint,comp,rng))
  o.to_csv(ROOT/f'{sample}_TLS_proxy_annotations.csv.gz',float_format='%.17g');print(sample,'centers',[int(centers[k].sum()) for k in TIERS],flush=True)
 stat=pd.DataFrame(rows);stat['joint_q_288']=false_discovery_control(stat.p.fillna(1).to_numpy());primary=stat.tier.eq(TIERS[1])&stat.comparison.eq(COMPS[0]);stat['primary_q_32']=np.nan;stat.loc[primary,'primary_q_32']=false_discovery_control(stat.loc[primary,'p'].fillna(1).to_numpy());stat.to_csv(ROOT/'spatial_statistics.csv',index=False)
 pd.DataFrame(summary).to_csv(ROOT/'section_candidate_summary.csv',index=False);pd.DataFrame(desc).to_csv(ROOT/'source_candidate_counts.csv',index=False);pd.DataFrame(coverage).to_csv(ROOT/'TLS_marker_coverage.csv',index=False);pd.DataFrame(components).to_csv(ROOT/'candidate_patch_components.csv',index=False);(ROOT/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
if __name__=='__main__':run()
