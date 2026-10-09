from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.optimize import nnls
R=Path(__file__).resolve().parent;S=Path('/workspace/scratch');ref=pd.read_csv(R/'EXPLORATORY_REFERENCE_NOT_VALIDATED.tsv',sep='\t',index_col=0);scales=pd.read_csv(R/'reference_row_scales.csv').row_scale.to_numpy();rows=[];coverage=[];target=ref.columns.get_loc('TAM_PGAM5_detected')
for n in ['GSE189903','GSE242889','GSE202642']:
 path=S/'PGAM5_myeloid_reference_v2'/(n+'_harmonized_counts_QC.h5ad') if n=='GSE189903' else S/n/('all_tumor_counts_QC.h5ad' if n=='GSE242889' else 'all_library_counts_QC.h5ad');a=ad.read_h5ad(path);o=pd.read_csv(R/(n+'_all_cell_annotation.csv.gz'));idx=a.obs_names.get_indexer(o.cell_id);assert (idx>=0).all();sym=a.var.gene.astype(str) if 'gene' in a.var else pd.Series(a.var_names);fidx=ref.index.get_indexer(sym);present=fidx>=0;mer=sparse.csr_matrix((np.ones(present.sum()),(np.flatnonzero(present),fidx[present])),shape=(len(sym),len(ref)));x=(a.X[idx]@mer).tocsr();norm=x.multiply(1e4/o.total_counts.to_numpy()[:,None]).tocsr();measurable=ref.index.isin(sym);A=ref.to_numpy()[measurable]/scales[measurable,None]
 for d,t in o.groupby('donor'):
  pos=t.PGAM5_RNA_status.eq('detected').sum();neg=t.PGAM5_RNA_status.eq('undetected').sum();tum=t.cell_type.eq('Tumor_hepatocyte').sum();tnk=t.cell_type.eq('T_NK').sum();coverage.append({'dataset':n,'donor':d,'target_cells':int(pos),'undetected_macrophages':int(neg),'tumor_epithelial_proxy_cells':int(tum),'T_NK_cells':int(tnk),'complete_mixture_eligible':bool(pos>=5 and neg>=20 and tum>=20 and tnk>=20),'note':'No inferred epithelial group does not establish biological absence; conservative annotation limits coverage'})
  pools={g:t.index[t.cell_type.eq(g)].to_numpy() for g in ['Monocyte_enriched','DC_enriched','pDC','Mixed_APC','T_NK','Tumor_hepatocyte','Unknown_other']};pools['PGAM5_undetected_macrophages']=t.index[t.PGAM5_RNA_status.eq('undetected')].to_numpy();pools['PGAM5_detected_nonmacrophages']=t.index[t.cell_type.ne('Macrophage')&t.PGAM5_counts.gt(0)].to_numpy()
  for group,ix in pools.items():
   if len(ix)<20:continue
   for unit in ['equal_RNA_cell_fraction','pooled_count_RNA_fraction']:
    y=np.asarray(norm[ix].mean(0)).ravel() if unit.startswith('equal') else np.asarray(x[ix].sum(0)).ravel()/o.total_counts.iloc[ix].sum()*1e4;b=y[measurable]/scales[measurable]
    for solver in ['NNLS','nonnegative_ridge_alpha0.01']:
     AA=A if solver=='NNLS' else np.vstack([A,np.sqrt(.01)*np.eye(A.shape[1])]);bb=b if solver=='NNLS' else np.r_[b,np.zeros(A.shape[1])];coef,_=nnls(AA,bb,maxiter=2000);w=coef/max(coef.sum(),1e-12);rows.append({'dataset':n,'donor':d,'pure_pool':group,'cells':len(ix),'solver':solver,'unit':unit,'true_target_fraction':0.,'predicted_target_fraction':float(w[target]),'interpretation':'Supplemental pressure test, not a substitute for complete mixed-tumor gates; unknown pools have uncertain true lineage'})
 del a,x,norm
pd.DataFrame(rows).to_csv(R/'all_donor_pure_competitor_tests.csv',index=False);pd.DataFrame(coverage).to_csv(R/'all_external_donor_test_coverage.csv',index=False);print('Supplemental pure competitor tests',len(rows),flush=True)
