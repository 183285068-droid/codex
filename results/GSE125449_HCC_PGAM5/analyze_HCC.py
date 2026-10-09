from pathlib import Path
import gzip,json,hashlib
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
from scipy.io import mmread
from threadpoolctl import threadpool_limits
R=Path('/workspace/scratch/GSE125449');records=[];rec=None
with gzip.open(R/'GSE125449_family.soft.gz','rt') as f:
 for line in f:
  line=line.strip()
  if line.startswith('^SAMPLE = '):
   if rec:records.append(rec)
   rec={'GSM':line.split(' = ')[1]}
  elif rec is not None and line.startswith('!Sample_title = '):rec['Sample']=line.split(' = ')[1]
  elif rec is not None and line.startswith('!Sample_characteristics_ch1 = cancer type: '):rec['diagnosis']=line.split('cancer type: ')[1]
if rec:records.append(rec)
geo=pd.DataFrame(records);assert geo.diagnosis.eq('Hepatocellular carcinoma').sum()==9
sets=[];gene_tables=[];all_obs=[]
for k in [1,2]:
 o=pd.read_csv(R/f'GSE125449_Set{k}_samples.txt.gz',sep='\t').rename(columns={'Cell Barcode':'barcode'})
 b=pd.read_csv(R/f'GSE125449_Set{k}_barcodes.tsv.gz',header=None)[0];g=pd.read_csv(R/f'GSE125449_Set{k}_genes.tsv.gz',sep='\t',header=None,names=['ensembl','gene'])
 assert o.barcode.equals(b) and g.ensembl.is_unique
 o=o.merge(geo,on='Sample',how='left',validate='many_to_one',sort=False);assert o.barcode.equals(b) and o.diagnosis.notna().all()
 o['set']=str(k);o['paper_sample']=np.where(o.diagnosis.eq('Hepatocellular carcinoma'),'H','C')+o.Sample.str.extract(r'LCP(\d+)$')[0]
 o['cell_key']=o.paper_sample+'_'+o.barcode.str.replace(r'-\d+$','',regex=True)
 o.index='set'+str(k)+'_'+o.Sample+'_'+o.barcode
 with threadpool_limits(limits=4):
  with gzip.open(R/f'GSE125449_Set{k}_matrix.mtx.gz','rb') as f:x=mmread(f).tocsr().T.tocsr().astype(np.int32)
 assert x.shape==(len(o),len(g)) and (x.data>0).all()
 g.index=g.ensembl;gene_tables.append(g)
 a=ad.AnnData(x,obs=o,var=g);sets.append(a);all_obs.append(o)
 print('Set',k,x.shape,flush=True)
all_o=pd.concat(all_obs);all_o.to_csv(R/'all_cell_metadata.csv.gz')
all_o[['Sample','paper_sample','set','GSM','diagnosis']].drop_duplicates().to_csv(R/'sample_diagnosis_mapping.csv',index=False)
a=ad.concat(sets,join='outer',merge='first');genes=pd.concat(gene_tables).drop_duplicates('ensembl').set_index('ensembl',drop=False);a.var=genes.reindex(a.var_names)
assert a.var.gene.notna().all() and a.obs_names.is_unique
hcc=a.obs.diagnosis.eq('Hepatocellular carcinoma').to_numpy();hc=a[hcc].copy();hc.obs.groupby(['Sample','paper_sample','Type'],observed=True).size().rename('cells').reset_index().to_csv(R/'HCC_celltypes_by_sample.csv',index=False)
pg=np.flatnonzero(hc.var.gene.eq('PGAM5'));assert len(pg)==1
hc.obs['total_counts']=np.asarray(hc.X.sum(1)).ravel();hc.obs['n_genes']=np.diff(hc.X.indptr)
mt=hc.var.gene.str.startswith('MT-').to_numpy();hc.obs['pct_mt']=np.asarray(hc.X[:,mt].sum(1)).ravel()/hc.obs.total_counts.to_numpy()*100
hc.obs['PGAM5_counts']=hc.X[:,pg[0]].toarray().ravel();hc.obs['PGAM5_detected']=hc.obs.PGAM5_counts.gt(0)
raw=hc[hc.obs.Type.eq('TAM').to_numpy()].copy();raw.obs.to_csv(R/'HCC_TAM_metadata_before_QC.csv.gz')
t=raw[raw.obs.n_genes.ge(500).to_numpy()&raw.obs.pct_mt.lt(20).to_numpy()].copy();t.write_h5ad(R/'HCC_TAM_counts_QC.h5ad',compression='gzip')
t.obs.to_csv(R/'HCC_TAM_analysis_metadata.csv.gz');t.obs.groupby(['paper_sample','Sample'],observed=True).agg(TAMs=('barcode','size'),PGAM5_positive=('PGAM5_detected','sum'),rate=('PGAM5_detected','mean')).to_csv(R/'HCC_TAMs_by_sample.csv')
pos=t.obs.PGAM5_counts.to_numpy()>0
print('HCC',len(hc),'raw TAM',len(raw),'QC TAM',len(t),'positive',pos.sum(),flush=True)
assert pos.any() and (~pos).any()
fp=np.asarray((t.X[pos]>0).mean(0)).ravel();fn=np.asarray((t.X[~pos]>0).mean(0)).ravel()
t.obs['PGAM5_group']=pd.Categorical(np.where(pos,'detected','undetected'));sc.pp.normalize_total(t,target_sum=1e4)
mp=np.asarray(t.X[pos].mean(0)).ravel();mn=np.asarray(t.X[~pos].mean(0)).ravel();sc.pp.log1p(t)
sc.tl.rank_genes_groups(t,groupby='PGAM5_group',groups=['detected'],reference='undetected',method='wilcoxon',tie_correct=True,corr_method='benjamini-hochberg',use_raw=False)
r=sc.get.rank_genes_groups_df(t,group='detected').set_index('names').rename(columns={'scores':'Wilcoxon_score','logfoldchanges':'log2FC','pvals':'p_value','pvals_adj':'FDR'});r['gene']=t.var.loc[r.index,'gene'].values
aux=pd.DataFrame({'fraction_PGAM5_positive':fp,'fraction_PGAM5_undetected':fn,'mean_10k_detected':mp,'mean_10k_undetected':mn},index=t.var_names);r=r.join(aux);r['defining_gene']=r.gene.eq('PGAM5');r.to_csv(R/'GSE125449_HCC_PGAM5_full_DE.csv')
mask=r.FDR.lt(.05)&r.log2FC.abs().ge(1)&r[['fraction_PGAM5_positive','fraction_PGAM5_undetected']].max(axis=1).ge(.1)
passed=r[mask].sort_values(['FDR','log2FC'],ascending=[True,False]);passed.to_csv(R/'GSE125449_HCC_PGAM5_DE_including_PGAM5.csv')
candidates=passed[~passed.defining_gene];up=candidates[candidates.log2FC.ge(1)];down=candidates[candidates.log2FC.le(-1)].sort_values(['FDR','log2FC']);up.to_csv(R/'GSE125449_HCC_PGAM5_upregulated.csv');down.to_csv(R/'GSE125449_HCC_PGAM5_downregulated.csv');up.head(20).to_csv(R/'GSE125449_HCC_PGAM5_top20_upregulated.csv')
# Match sample identifiers and barcode cores to assess sample/cell reuse, not independent replication.
old=pd.read_csv('/workspace/scratch/GSE151530/all_cell_metadata.csv.gz');old['cell_key']=old.Sample.str.replace(r'[a-z]+$','',regex=True)+'_'+old.Cell.str.replace(r'-\d+$','',regex=True)
old_keys=set(old.cell_key);overlap=t.obs.copy();overlap['cell_key_found_in_GSE151530']=overlap.cell_key.isin(old_keys)
overlap[['Sample','paper_sample','barcode','cell_key','PGAM5_counts','cell_key_found_in_GSE151530']].to_csv(R/'GSE151530_sample_barcode_overlap_check.csv',index=False)
summary={'HCC_samples':int(hc.obs.Sample.nunique()),'HCC_sample_ids':sorted(hc.obs.paper_sample.unique().tolist()),'HCC_cells':len(hc),'raw_HCC_TAMs':len(raw),'QC_HCC_TAMs':len(t),'QC_samples_with_TAMs':int(t.obs.Sample.nunique()),'PGAM5_positive_cells':int(pos.sum()),'PGAM5_undetected_cells':int((~pos).sum()),'samples_with_PGAM5_positive_TAMs':int(t.obs.loc[pos,'Sample'].nunique()),'tested_genes':len(r),'upregulated_excluding_PGAM5':len(up),'downregulated_excluding_PGAM5':len(down),'top20_up':up.gene.head(20).tolist(),'top10_down':down.gene.head(10).tolist(),'cells_sample_barcode_matching_GSE151530':int(overlap.cell_key_found_in_GSE151530.sum()),'criteria':'FDR <0.05, |Scanpy approximate log2FC| >=1, >=10% detection in either group; PGAM5 excluded from candidate lists','method':'Pooled-cell Wilcoxon, tie correction, BH FDR; no patient pairing, depth matching or sample/depth covariates','QC':'n_genes>=500, pct_mt<20; 10k total-count normalization then log1p','cell_annotation':'Original TAM labels'}
(R/'analysis_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
# Check macrophage identity and hepatic/nonmyeloid signals without removing cells post hoc.
qa=[]
for gene in ['C1QA','C1QB','C1QC','CD68','CSF1R','LST1','TYROBP','SPP1','MMP9','ALB','APOA1','TTR','EPCAM','KRT19','CD3D','TRAC','MS4A1']:
 idx=np.flatnonzero(t.var.gene.eq(gene))
 if not len(idx):continue
 v=t.X[:,idx].toarray().sum(axis=1)
 for label,sel in [('PGAM5_detected',pos),('PGAM5_undetected',~pos)]:qa.append({'gene':gene,'group':label,'cells':int(sel.sum()),'fraction':float((v[sel]>0).mean()),'mean_log1p_10k':float(v[sel].mean())})
pd.DataFrame(qa).to_csv(R/'TAM_marker_QA.csv',index=False)
