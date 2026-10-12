"""Independent target reconstruction, lattice and permutation/FDR verification."""
import json
import numpy as np,pandas as pd,anndata as ad
from scipy.spatial import cKDTree
from scipy.stats import false_discovery_control
from analyze import R,BASE,K,PANELS,SAMPLES,dense
checks=[]
def check(name,v):
 assert v,name
 checks.append(name)
def run():
 d=pd.read_csv(R/'neighborhood_counts.csv');st=pd.read_csv(R/'spatial_statistics.csv')
 for item in SAMPLES:
  sample=item['sample'];a=ad.read_h5ad(R/f'{sample}_supplemental_T_markers.h5ad');o=a.obs;raw=dense(a.layers['counts']);x=dense(a.X)
  check(sample+' full-library supplemental normalization',np.allclose(x,np.log1p(raw/o.total_counts.to_numpy()[:,None]*10000),rtol=1e-6,atol=1e-6))
  xy=o[['array_row','array_col']].to_numpy(float);xy[:,0]*=np.sqrt(3);edges=cKDTree(xy).query_pairs(2.000001,output_type='ndarray');edge_set=set(map(tuple,edges))|set(map(tuple,edges[:,::-1]));e=a.obsp['spatial_connectivities'].tocoo();check(sample+' independent lattice reconstruction',edge_set==set(zip(e.row,e.col)))
  t=o.T_RNA_enriched.to_numpy(bool)
  def detected(g):return raw[:,a.var_names.get_loc(g)]>0
  for k,gs in PANELS.items():
   ix=a.var_names.get_indexer(gs);score=x[:,ix].mean(1);n=(raw[:,ix]>0).sum(1);minimum=3 if k in ['Naive_like','Central_memory_like'] else 2;gate=t&(score>=np.quantile(score[t],.75))&(score>0)&(n>=minimum)
   if k=='CD4_related':gate&=detected('CD4')
   if k=='CD8_supported':gate&=detected('CD8A')&detected('CD8B')
   if k=='Treg_supported':gate&=detected('FOXP3')
   if k=='Naive_like':gate&=(detected('CCR7')|detected('SELL'))&(detected('TCF7')|detected('LEF1'))
   if k=='Central_memory_like':gate&=detected('CCR7')&detected('IL7R')
   if k=='Effector_memory_like':gate&=detected('GZMK')&(detected('CCL5')|detected('IL7R'))
   check(sample+'/'+k+' independently reconstructed gate',np.array_equal(gate,o['gate_'+k]))
  for kind in st.comparison.unique():
   p=R/f'{sample}_{kind}_audit.npz'
   if not p.exists():continue
   z=np.load(p);h=z['high'];y=z['outcomes'];idx=o.index.get_indexer(z['barcodes']);check(sample+'/'+kind+' high labels',np.array_equal(h,o.PGAM5_mac_group.eq('PGAM5_high').to_numpy()[idx]))
   outcomes=np.column_stack([np.asarray(a.obsp['spatial_connectivities']@o['gate_'+k].to_numpy(int)).ravel()>0 for k in K]);check(sample+'/'+kind+' saved outcomes from geometry',np.array_equal(y,outcomes[idx]))
   n1=h.sum();n0=len(h)-n1;obs=100*(h@y/n1-(1-h)@y/n0);check(sample+'/'+kind+' observed statistic',np.allclose(obs,z['observed']))
   rows=st[(st['sample']==sample)&(st.comparison==kind)].set_index('target')
   for j,k in enumerate(K):
    row=rows.loc[k];null=z['null'][:,j]
    if row.status=='tested':
     pvalue=min(1,2*min((1+(null>=obs[j]-1e-10).sum())/5000,(1+(null<=obs[j]+1e-10).sum())/5000));check(sample+'/'+kind+'/'+k+' equal-tail permutation p',abs(pvalue-row.p)<1e-12)
 for field,mask in [('joint_q_288',np.ones(len(st),bool)),('primary_q_96',st.comparison.eq('high_vs_other_detected')),('primary_only_q_90',st.comparison.eq('high_vs_other_detected')&st.tumor_status.eq('primary'))]:
  check(field+' BH reconstruction',np.allclose(false_discovery_control(st.loc[mask,'p'].fillna(1).to_numpy()),st.loc[mask,field]))
 (R/'independent_audit_results.json').write_text(json.dumps(dict(passed=len(checks),checks=checks),indent=2));print('PASS',len(checks))
if __name__=='__main__':run()
