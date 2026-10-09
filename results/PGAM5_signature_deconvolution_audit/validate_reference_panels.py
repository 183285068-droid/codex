from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy.optimize import nnls
from scipy.stats import pearsonr
R=Path(__file__).resolve().parent
cons=pd.read_csv(R/'primary_consensus_at_least3.csv').gene.tolist();anchors=['C1QA','C1QB','C1QC','CSF1R','CD68','TYROBP','FCER1G','LST1','AIF1','CD163','MSR1','MRC1','SPP1']
maptypes={'T cells':'T','T-cell':'T','B cells':'B','B-cell':'B','Malignant cells':'Malignant','Malignant cell':'Malignant','CAFs':'CAF','CAF':'CAF','TECs':'Endothelial','TEC':'Endothelial','unclassified':'Unclassified','TAM_PGAM5_detected':'TAM_detected','TAM_PGAM5_undetected':'TAM_undetected'}
data={}
for name in ['GSE151530','GSE189903']:
 z=np.load(R/(name+'_selected_gene_cells.npz'));o=pd.read_csv(R/(name+'_QC_cell_metadata.csv.gz'));o['group']=o.reference_group.map(maptypes);assert o.group.notna().all();data[name]=(z['expression_cp10k'],z['genes'].tolist(),o)
groups=['TAM_detected','TAM_undetected','Malignant','T','B','CAF','Endothelial','Unclassified']
panels={'state21_only':cons,'lineage13_plus_state21':anchors+cons,'lineage13_plus_state21_plus_PGAM5':anchors+cons+['PGAM5']}
x,genes,o=data['GSE151530'];fullref=np.column_stack([x[o.group.eq(g)].mean(0) for g in groups]);pd.DataFrame(fullref,index=genes,columns=groups).to_csv(R/'EXPLORATORY_GSE151530_reference_CP10K_NOT_VALIDATED.tsv',sep='\t')
rng=np.random.default_rng(42);results=[];conditions=[]
# Reference generation leaves a donor out; candidate discovery was pooled, so this remains INTERNAL diagnostic validation.
# Pseudo-mixtures average CP10K-normalized cells, assuming equal RNA output per cell. They are not actual TCGA bulk.
for source,(xx,gg,oo) in data.items():
 assert gg==genes
 for donor in sorted(oo.donor.unique()):
  test=oo.donor.eq(donor).to_numpy()
  if (test&oo.group.eq('TAM_detected').to_numpy()).sum()<5:continue
  available={g:np.flatnonzero(test&oo.group.eq(g).to_numpy()) for g in groups}
  if any(len(v)==0 for v in available.values()):continue
  train=~o.donor.eq(donor).to_numpy() if source=='GSE151530' else np.ones(len(o),dtype=bool)
  ref=np.column_stack([x[train&o.group.eq(g).to_numpy()].mean(0) for g in groups]);assert np.isfinite(ref).all()
  for panel,chosen in panels.items():
   ix=[genes.index(k) for k in chosen];A=ref[ix];scale=np.sqrt(np.mean(A*A,axis=1));A=A/np.maximum(scale[:,None],1e-8)
   conditions.append({'test_dataset':source,'donor':donor,'panel':panel,'matrix_condition_number':float(np.linalg.cond(A)),'TAM_column_cosine':float(np.dot(A[:,0],A[:,1])/(np.linalg.norm(A[:,0])*np.linalg.norm(A[:,1])))})
   rng=np.random.default_rng(42+sum((source+str(donor)).encode())) # identical mixtures for each panel
   for f in [0,.005,.01,.02,.05]:
    weights=np.array([f,.25-f,.50,.15,.03,.025,.025,.02]);assert np.isclose(weights.sum(),1)
    for rep in range(10):
     # 2,000-cell equal-RNA pseudo-bulk; 0.5% is 10 positive cells.
     counts=np.rint(weights*2000).astype(int);y=np.zeros(len(genes))
     for g,n in zip(groups,counts):
      if n:y+=xx[rng.choice(available[g],size=n,replace=True)].sum(0)/2000
     beta,res=nnls(A,y[ix]/np.maximum(scale,1e-8));pred=beta[0]/max(beta.sum(),1e-12)
     results.append({'test_dataset':source,'donor':donor,'panel':panel,'true_equal_RNA_cell_fraction':f,'estimated_fraction':pred,'replicate':rep,'nnls_residual':res})
r=pd.DataFrame(results);assert len(r)>0;r.to_csv(R/'internal_equal_RNA_pseudobulk_predictions.csv',index=False)
pd.DataFrame(conditions).to_csv(R/'reference_identifiability_diagnostics.csv',index=False)
summary=[]
for (source,panel),t in r.groupby(['test_dataset','panel']):
 zero=t[t.true_equal_RNA_cell_fraction.eq(0)];nonzero=t[t.true_equal_RNA_cell_fraction.gt(0)];err=t.estimated_fraction-t.true_equal_RNA_cell_fraction
 row={'test_dataset':source,'panel':panel,'heldout_donors':t.donor.nunique(),'mixtures':len(t),'MAE_percentage_points':float(abs(err).mean()*100),'zero_true_median_estimated_percent':float(zero.estimated_fraction.median()*100),'zero_true_p95_estimated_percent':float(zero.estimated_fraction.quantile(.95)*100),'pearson_r':float(pearsonr(t.true_equal_RNA_cell_fraction,t.estimated_fraction).statistic),'status':'exploratory internal diagnostics, not TCGA validation'}
 summary.append(row)
s=pd.DataFrame(summary);s.to_csv(R/'internal_pseudobulk_validation_summary.csv',index=False);print(s.round(3).to_string(index=False));print('donors',r.groupby('test_dataset').donor.unique().to_dict())
