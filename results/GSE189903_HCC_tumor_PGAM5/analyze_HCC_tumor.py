from pathlib import Path
import gzip,json,sys
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
from scipy.io import mmread
from threadpoolctl import threadpool_limits
R=Path(__file__).resolve().parent
records=[];rec=None
with gzip.open(R/'GSE189903_family.soft.gz','rt') as f:
 for line in f:
  line=line.strip()
  if line.startswith('^SAMPLE = '):
   if rec:records.append(rec)
   rec={'GSM':line.split(' = ')[1]}
  elif rec is not None and line.startswith('!Sample_title = '):
   rec['GEO_title']=line.split(' = ')[1];rec['S_ID']=line.split()[-1]
  elif rec is not None and line.startswith('!Sample_characteristics_ch1 = cancer type: '):
   rec['diagnosis']=line.split('cancer type: ')[1]
  elif rec is not None and line.startswith('!Sample_characteristics_ch1 = tissue type: '):
   rec['tissue']=line.split('tissue type: ')[1]
if rec:records.append(rec)
geo=pd.DataFrame(records);assert len(geo)==34 and geo.S_ID.is_unique
o=pd.read_csv(R/'GSE189903_Info.txt.gz',sep='\t')
b=pd.read_csv(R/'GSE189903_barcodes.tsv.gz',header=None)[0]
g=pd.read_csv(R/'GSE189903_genes.tsv.gz',sep='\t',header=None,names=['ensembl','gene'])
assert o.Cell.equals(b) and o.Cell.is_unique and g.ensembl.is_unique and g.gene.notna().all()
o=o.merge(geo,on='S_ID',how='left',validate='many_to_one',sort=False);assert o.Cell.equals(b) and o.diagnosis.notna().all()
o['patient']=o.Sample.str.extract(r'^(\d+[HC])')[0];assert o.patient.notna().all()
o.index=o.Cell
selected=o.diagnosis.eq('Hepatocellular carcinoma')&o.tissue.isin(['Tumor core','Tumor border'])
assert o.loc[selected,'Sample'].nunique()==16 and o.loc[selected,'patient'].nunique()==4
mapping=o[['S_ID','Sample','patient','GSM','diagnosis','tissue']].drop_duplicates()
mapping['included']=mapping.S_ID.isin(o.loc[selected,'S_ID'])
mapping.to_csv(R/'sample_diagnosis_tissue_mapping.csv',index=False)
o.to_csv(R/'all_cell_metadata.csv.gz')
o[selected].groupby(['S_ID','Sample','patient','tissue','Type'],observed=True).size().rename('cells').reset_index().to_csv(R/'HCC_tumor_celltypes_by_sample.csv',index=False)
tam=selected&o.Type.eq('TAM')
print('Source cells',len(o),'selected HCC tumor cells',int(selected.sum()),'raw TAM',int(tam.sum()),flush=True)
with threadpool_limits(limits=4):
 with gzip.open(R/'GSE189903_matrix.mtx.gz','rb') as f:x=mmread(f).tocsr()
assert x.shape==(len(g),len(o)),(x.shape,len(g),len(o))
x=x[:,np.flatnonzero(tam.to_numpy())].T.tocsr().astype(np.int32)
assert (x.data>0).all()
g.index=g.ensembl;a=ad.AnnData(x,obs=o[tam].copy(),var=g)
pg=np.flatnonzero(g.gene.eq('PGAM5'));assert len(pg)==1
a.obs['total_counts']=np.asarray(a.X.sum(1)).ravel();a.obs['n_genes']=np.diff(a.X.indptr)
a.obs['pct_mt']=np.asarray(a.X[:,g.gene.str.startswith('MT-').to_numpy()].sum(1)).ravel()/a.obs.total_counts.to_numpy()*100
a.obs['PGAM5_counts']=a.X[:,pg[0]].toarray().ravel();a.obs['PGAM5_detected']=a.obs.PGAM5_counts.gt(0)
a.obs.to_csv(R/'HCC_tumor_TAM_metadata_before_QC.csv.gz')
t=a[a.obs.n_genes.ge(500).to_numpy()&a.obs.pct_mt.lt(20).to_numpy()].copy()
t.obs.to_csv(R/'HCC_tumor_TAM_analysis_metadata.csv.gz');t.write_h5ad(R/'HCC_tumor_TAM_counts_QC.h5ad',compression='gzip')
t.obs.groupby(['S_ID','Sample','patient','tissue'],observed=True).agg(TAMs=('Cell','size'),PGAM5_positive=('PGAM5_detected','sum'),rate=('PGAM5_detected','mean')).to_csv(R/'PGAM5_by_sample.csv')
t.obs.groupby(['patient','tissue'],observed=True).agg(TAMs=('Cell','size'),PGAM5_positive=('PGAM5_detected','sum'),rate=('PGAM5_detected','mean')).to_csv(R/'PGAM5_by_patient_tissue.csv')
t.obs.groupby(['tissue'],observed=True).agg(TAMs=('Cell','size'),PGAM5_positive=('PGAM5_detected','sum'),rate=('PGAM5_detected','mean')).to_csv(R/'PGAM5_by_tissue.csv')
pos=t.obs.PGAM5_detected.to_numpy();assert min(pos.sum(),(~pos).sum())>=2
print('QC TAM',len(t),'PGAM5 detected',int(pos.sum()),'undetected',int((~pos).sum()),flush=True)
fp=np.asarray((t.X[pos]>0).mean(0)).ravel();fn=np.asarray((t.X[~pos]>0).mean(0)).ravel()
t.obs['PGAM5_group']=pd.Categorical(np.where(pos,'detected','undetected'))
sc.pp.normalize_total(t,target_sum=1e4);mp=np.asarray(t.X[pos].mean(0)).ravel();mn=np.asarray(t.X[~pos].mean(0)).ravel();sc.pp.log1p(t)
sc.tl.rank_genes_groups(t,groupby='PGAM5_group',groups=['detected'],reference='undetected',method='wilcoxon',tie_correct=True,corr_method='benjamini-hochberg',use_raw=False)
r=sc.get.rank_genes_groups_df(t,group='detected').set_index('names').rename(columns={'scores':'Wilcoxon_score','logfoldchanges':'log2FC','pvals':'p_value','pvals_adj':'FDR'})
r['gene']=t.var.loc[r.index,'gene'].values
aux=pd.DataFrame({'fraction_PGAM5_positive':fp,'fraction_PGAM5_undetected':fn,'mean_10k_detected':mp,'mean_10k_undetected':mn},index=t.var_names);r=r.join(aux);r['defining_gene']=r.gene.eq('PGAM5')
r.to_csv(R/'GSE189903_HCC_tumor_PGAM5_full_DE.csv')
mask=r.FDR.lt(.05)&r.log2FC.abs().ge(1)&r[['fraction_PGAM5_positive','fraction_PGAM5_undetected']].max(axis=1).ge(.1)
passed=r[mask].sort_values(['FDR','log2FC'],ascending=[True,False]);passed.to_csv(R/'GSE189903_HCC_tumor_PGAM5_DE_including_PGAM5.csv')
candidates=passed[~passed.defining_gene];up=candidates[candidates.log2FC.ge(1)];down=candidates[candidates.log2FC.le(-1)].sort_values(['FDR','log2FC'])
up.to_csv(R/'GSE189903_HCC_tumor_PGAM5_upregulated.csv');down.to_csv(R/'GSE189903_HCC_tumor_PGAM5_downregulated.csv');up.head(20).to_csv(R/'GSE189903_HCC_tumor_PGAM5_top20_upregulated.csv')
qa=[]
for gene in ['C1QA','C1QB','C1QC','CD68','CSF1R','LST1','TYROBP','SPP1','MMP9','ALB','APOA1','TTR','EPCAM','KRT19','CD3D','TRAC','MS4A1','FCN1','S100A8']:
 idx=np.flatnonzero(g.gene.eq(gene))
 if not len(idx):continue
 v=t.X[:,idx].toarray().sum(axis=1)
 for label,sel in [('PGAM5_detected',pos),('PGAM5_undetected',~pos)]:qa.append({'gene':gene,'group':label,'cells':int(sel.sum()),'fraction':float((v[sel]>0).mean()),'mean_log1p_10k':float(v[sel].mean())})
pd.DataFrame(qa).to_csv(R/'TAM_marker_QA.csv',index=False)
s={'HCC_tumor_samples':int(o.loc[selected,'Sample'].nunique()),'HCC_tumor_patients':int(o.loc[selected,'patient'].nunique()),'HCC_tumor_cells':int(selected.sum()),'HCC_tumor_sample_ids':o.loc[selected,'Sample'].drop_duplicates().tolist(),'HCC_selected_sample_tissues':mapping[mapping.included].tissue.value_counts().to_dict(),'raw_TAMs':len(a),'QC_TAMs':len(t),'PGAM5_positive_cells':int(pos.sum()),'PGAM5_undetected_cells':int((~pos).sum()),'QC_samples_with_TAMs':int(t.obs.Sample.nunique()),'samples_with_positive_TAMs':int(t.obs.loc[pos,'Sample'].nunique()),'patients_with_positive_TAMs':int(t.obs.loc[pos,'patient'].nunique()),'tested_genes':len(r),'upregulated_excluding_PGAM5':len(up),'downregulated_excluding_PGAM5':len(down),'top20_up':up.gene.head(20).tolist(),'top10_down':down.gene.head(10).tolist(),'scope':'Official Hepatocellular carcinoma diagnosis, Tumor core or Tumor border only; adjacent non-tumor and ICC excluded per user clarification','criteria':'FDR<0.05, |Scanpy approximate log2FC|>=1, >=10% detection in either group; PGAM5 excluded from candidates','method':'Pooled-cell Wilcoxon with tie correction, BH FDR across all genes; no patient pairing, depth matching or sample/depth covariates','QC':'n_genes>=500,pct_mt<20; normalize_total 10000 then log1p','cell_annotation':'Original author TAM labels retained'}
(R/'analysis_summary.json').write_text(json.dumps(s,indent=2));print(json.dumps(s,indent=2),flush=True)
