from pathlib import Path
import gzip,json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy.io import mmread
r=Path('/workspace/scratch/GSE283406_restored');reports=[]; markers=['PGAM5','CD3D','CD3E','CD3G','TRAC','CD8A','CD8B','CD4','FOXP3','IKZF2','IL2RA','CTLA4','TIGIT','PDCD1','LAG3','HAVCR2','TOX','CXCL13','CD69','TNFRSF9','ICOS','CCR7','SELL','TCF7','LEF1','IL7R','GZMK','ITGAE','ZNF683','PRF1','GZMB','CCL5']
for gsm in ['GSM8661790','GSM8661791','GSM8661792']:
 paths={k:next(Path('/workspace/attachments').glob('*/'+gsm+'_*.'+v+'.gz')) for k,v in [('matrix','matrix.mtx'),('features','features.tsv'),('barcodes','barcodes.tsv')]}
 f=pd.read_csv(paths['features'],sep='\t',header=None);b=pd.read_csv(paths['barcodes'],sep='\t',header=None)[0];x=mmread(paths['matrix']).tocsr();assert x.shape==(len(f),len(b));assert b.is_unique;assert (x.data>=0).all() and np.all(x.data==np.floor(x.data))
 d={'sample':gsm,'genes':len(f),'barcodes':len(b),'nonzero':int(x.nnz),'missing_markers':[g for g in markers if g not in set(f[1])],'files':{k:{'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for k,p in paths.items()}}
 if gsm=='GSM8661792':
  a=ad.AnnData(X=x.T.tocsr(),obs=pd.DataFrame(index=b.tolist()),var=pd.DataFrame({'gene_id':f[0].tolist()},index=f[1].tolist()));a.var_names_make_unique();a.write_h5ad(r/'J259X_AVITI_full_panel_counts.h5ad',compression='gzip')
  old=ad.read_h5ad('/workspace/codex/results/GSE283406_PGAM5_Tcell_expression_candidates/inputs/J259X_AVITI_expression_marker_subset.h5ad');sub=a[old.obs_names,old.var_names];d['previous_subset_counts_match']=bool((sub.X!=old.layers['counts']).nnz==0);assert d['previous_subset_counts_match']
  qc=(np.asarray(a.X.sum(1)).ravel()>=500)&(a.X.getnnz(axis=1)>=200)
  d['QC_barcodes']=int(qc.sum());d['marker_detection_QC']={g:int((a.X[qc,a.var_names.get_loc(g)]>0).sum()) for g in markers if g in a.var_names}
 reports.append(d);print(json.dumps(d,ensure_ascii=False))
(r/'upload_validation.json').write_text(json.dumps(reports,indent=2,ensure_ascii=False))
