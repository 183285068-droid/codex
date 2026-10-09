from pathlib import Path
import gzip,json,hashlib,sys
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
from scipy.io import mmread
from scipy.sparse import csr_matrix
from threadpoolctl import threadpool_limits
R=Path(__file__).resolve().parent
summaries={}
platforms=sys.argv[1:] or ['Droplet','Smartseq2']
assert set(platforms).issubset({'Droplet','Smartseq2'})
for tech in platforms:
 print('Preparing',tech,flush=True)
 pref=R/tech;pref.mkdir(exist_ok=True)
 if tech=='Droplet':
  o=pd.read_csv(R/'GSE140228_UMI_counts_Droplet_cellinfo.tsv.gz',sep='\t')
  g=pd.read_csv(R/'GSE140228_UMI_counts_Droplet_genes.tsv.gz',sep='\t')
  b=pd.read_csv(R/'GSE140228_UMI_counts_Droplet_barcodes.tsv.gz',header=None)[0]
  assert o.Barcode.equals(b)
 else:
  o=pd.read_csv(R/'GSE140228_cell_info_Smartseq2.tsv.gz',sep='\t')
  g=pd.read_csv(R/'GSE140228_gene_info_Smartseq2.tsv.gz',sep='\t')
 assert o.Barcode.is_unique and g.ENSEMBL.is_unique
 o.index=o.Barcode.astype(str)
 tumor=o.Tissue.eq('Tumor');macro=o.celltype_sub.str.startswith('Mφ-')
 selected=o[tumor & macro].copy()
 o[tumor].groupby(['Sample','Donor','Histology','Tissue_sub','celltype_sub'],observed=True).size().rename('cells').reset_index().to_csv(pref/'tumor_sample_celltypes.csv',index=False)
 if tech=='Droplet':
  with threadpool_limits(limits=4):
   with gzip.open(R/'GSE140228_UMI_counts_Droplet.mtx.gz','rb') as f:x=mmread(f).tocsr()
  assert x.shape==(len(g),len(o)),(x.shape,len(g),len(o))
  x=x[:,np.flatnonzero((tumor&macro).to_numpy())].T.tocsr().astype(np.int32)
 else:
  with gzip.open(R/'GSE140228_read_counts_Smartseq2.csv.gz','rt') as f:header=f.readline().rstrip().split(',')
  assert set(header[1:])==set(o.Barcode) and len(header)-1==len(o)
  c=pd.read_csv(R/'GSE140228_read_counts_Smartseq2.csv.gz',usecols=['gene']+selected.Barcode.tolist(),index_col='gene')
  assert len(c)==len(g)
  actual=c.index.to_numpy();plain=g.SYMBOL.to_numpy();qualified=(g.SYMBOL+'_'+g.ENSEMBL).to_numpy()
  assert np.all((actual==plain)|(actual==qualified)), 'Gene row identity mismatch'
  g['count_row_label']=actual
  x=csr_matrix(c[selected.Barcode].to_numpy(dtype=np.int32).T)
  del c
 assert (x.data>0).all() and x.shape==(len(selected),len(g))
 g.index=g.ENSEMBL;g['gene']=g.SYMBOL
 a=ad.AnnData(x,obs=selected,var=g)
 pg=np.flatnonzero(g.gene.eq('PGAM5'));assert len(pg)==1
 a.obs['total_counts']=np.asarray(a.X.sum(1)).ravel();a.obs['n_genes']=np.diff(a.X.indptr)
 a.obs['pct_mt']=np.asarray(a.X[:,g.gene.str.startswith('MT-').to_numpy()].sum(1)).ravel()/a.obs.total_counts.to_numpy()*100
 a.obs['PGAM5_counts']=a.X[:,pg[0]].toarray().ravel();a.obs['PGAM5_detected']=a.obs.PGAM5_counts.gt(0)
 a.obs.to_csv(pref/'macrophage_metadata_before_QC.csv.gz')
 t=a[a.obs.n_genes.ge(500).to_numpy()&a.obs.pct_mt.lt(20).to_numpy()].copy()
 t.obs.to_csv(pref/'macrophage_analysis_metadata.csv.gz');t.write_h5ad(pref/'macrophage_counts_QC.h5ad',compression='gzip')
 t.obs.groupby(['Sample','Donor','Histology'],observed=True).agg(macrophages=('Barcode','size'),PGAM5_positive=('PGAM5_detected','sum'),rate=('PGAM5_detected','mean')).to_csv(pref/'PGAM5_by_sample.csv')
 pos=t.obs.PGAM5_detected.to_numpy()
 assert min(pos.sum(),(~pos).sum())>=2,(tech,pos.sum(),(~pos).sum())
 print(tech,'tumor',tumor.sum(),'raw macrophages',len(a),'QC',len(t),'PGAM5+',pos.sum(),'zero',(~pos).sum(),flush=True)
 fp=np.asarray((t.X[pos]>0).mean(0)).ravel();fn=np.asarray((t.X[~pos]>0).mean(0)).ravel()
 t.obs['PGAM5_group']=pd.Categorical(np.where(pos,'detected','undetected'))
 sc.pp.normalize_total(t,target_sum=1e4);mp=np.asarray(t.X[pos].mean(0)).ravel();mn=np.asarray(t.X[~pos].mean(0)).ravel();sc.pp.log1p(t)
 sc.tl.rank_genes_groups(t,groupby='PGAM5_group',groups=['detected'],reference='undetected',method='wilcoxon',tie_correct=True,corr_method='benjamini-hochberg',use_raw=False)
 r=sc.get.rank_genes_groups_df(t,group='detected').set_index('names').rename(columns={'scores':'Wilcoxon_score','logfoldchanges':'log2FC','pvals':'p_value','pvals_adj':'FDR'})
 r['gene']=t.var.loc[r.index,'gene'].values
 aux=pd.DataFrame({'fraction_PGAM5_positive':fp,'fraction_PGAM5_undetected':fn,'mean_10k_detected':mp,'mean_10k_undetected':mn},index=t.var_names);r=r.join(aux);r['defining_gene']=r.gene.eq('PGAM5');r.to_csv(pref/'full_DE.csv')
 mask=r.FDR.lt(.05)&r.log2FC.abs().ge(1)&r[['fraction_PGAM5_positive','fraction_PGAM5_undetected']].max(axis=1).ge(.1)
 passed=r[mask].sort_values(['FDR','log2FC'],ascending=[True,False]);passed.to_csv(pref/'DE_including_PGAM5.csv')
 candidates=passed[~passed.defining_gene];up=candidates[candidates.log2FC.ge(1)];down=candidates[candidates.log2FC.le(-1)].sort_values(['FDR','log2FC']);up.to_csv(pref/'upregulated.csv');down.to_csv(pref/'downregulated.csv');up.head(20).to_csv(pref/'top20_upregulated.csv')
 qa=[]
 for gene in ['C1QA','C1QB','C1QC','CD68','CSF1R','LST1','TYROBP','SPP1','MMP9','ALB','APOA1','TTR','EPCAM','KRT19','CD3D','TRAC','MS4A1','FCN1','S100A8']:
  idx=np.flatnonzero(g.gene.eq(gene))
  if not len(idx):continue
  v=t.X[:,idx].toarray().sum(axis=1)
  for label,sel in [('PGAM5_detected',pos),('PGAM5_undetected',~pos)]:qa.append({'gene':gene,'group':label,'cells':int(sel.sum()),'fraction':float((v[sel]>0).mean()),'mean_log1p_10k':float(v[sel].mean())})
 pd.DataFrame(qa).to_csv(pref/'marker_QA.csv',index=False)
 s={'platform':tech,'tumor_samples':int(o.loc[tumor,'Sample'].nunique()),'tumor_donors':int(o.loc[tumor,'Donor'].nunique()),'tumor_cells':int(tumor.sum()),'tumor_cell_histology':o[tumor].Histology.value_counts().to_dict(),'tumor_sample_histology':o[tumor][['Sample','Histology']].drop_duplicates().Histology.value_counts().to_dict(),'raw_macrophages':len(a),'QC_macrophages':len(t),'PGAM5_positive':int(pos.sum()),'PGAM5_undetected':int((~pos).sum()),'QC_macrophage_histology':t.obs.Histology.value_counts().to_dict(),'PGAM5_positive_histology':t.obs[pos].Histology.value_counts().to_dict(),'samples_with_QC_macrophages':int(t.obs.Sample.nunique()),'samples_with_positive_macrophages':int(t.obs.loc[pos,'Sample'].nunique()),'tested_genes':len(r),'upregulated_excluding_PGAM5':len(up),'downregulated_excluding_PGAM5':len(down),'top20_up':up.gene.head(20).tolist(),'top10_down':down.gene.head(10).tolist(),'criteria':'FDR <0.05, |Scanpy approximate log2FC| >=1, >=10% detection in either group; PGAM5 excluded from candidate lists','method':'Pooled-cell Wilcoxon, tie correction, BH FDR across all genes; no patient pairing, depth matching, sample or depth covariates','QC':'n_genes>=500, pct_mt<20; normalization to 10000 then log1p','annotation':'Original celltype_sub labels starting with Mφ-, all six macrophage subtypes retained','count_type':'UMI' if tech=='Droplet' else 'read counts'}
 (pref/'analysis_summary.json').write_text(json.dumps(s,indent=2,ensure_ascii=False));summaries[tech]=s;print(json.dumps(s,indent=2,ensure_ascii=False),flush=True)
 del x,a,t
summaries={tech:json.loads((R/tech/'analysis_summary.json').read_text()) for tech in ['Droplet','Smartseq2']}
(R/'analysis_summary.json').write_text(json.dumps(summaries,indent=2,ensure_ascii=False))
