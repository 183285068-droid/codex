"""Within-slice Visium QC and expression-based clustering; spots are not cells."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE245908_mpl')
os.environ.setdefault('NUMBA_CACHE_DIR','/tmp/GSE245908_numba')
import numpy as np,pandas as pd,scanpy as sc,gzip,json,scipy.sparse as sp
from PIL import Image
from io import BytesIO
ROOT=Path(__file__).resolve().parent
WORK=Path(os.environ.get('GSE245908_WORK','/workspace/scratch/GSE245908'));WORK.mkdir(parents=True,exist_ok=True)
SAMPLES={'CHC20':'GSM7850822','CHC23':'GSM7850823'}
PROGRAMS={
 'Macrophage_core':['C1QA','C1QB','C1QC','CD68','CSF1R','MERTK','MSR1'],
 'Pan_myeloid':['LST1','TYROBP','FCER1G','AIF1','LYZ','CTSS'],
 'T_cell':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2'],
 'Cytotoxic_T_NK':['NKG7','GNLY','PRF1','GZMB','GZMA','CTSW'],
 'DC':['CD1C','FCER1A','CLEC10A','CLEC9A','XCR1','CST3'],
 'Monocyte_inflammatory':['S100A8','S100A9','FCN1','VCAN','IL1B'],
 'Fibroblast_matrix':['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA'],
 'Endothelial':['PECAM1','VWF','CDH5','KDR','EMCN'],
 'Hepatocyte':['ALB','APOA1','APOA2','TTR','ASGR1','CYP2E1'],
 'Epithelial':['EPCAM','KRT8','KRT18','KRT19'],
 'Cycling':['MKI67','TOP2A','UBE2C','CENPF','BIRC5','PCNA'],
 'THBS1_TREM1_myeloid':['THBS1','TREM1','CD163','TGFB1','IL13RA1'],
}
EXTRA=['PGAM5','MYO19','SPP1','MMP12','LPL','STAB1','MARCO','CD5L','SPIC','DNASE1L3','HLA-DRA','HLA-DRB1','CXCL8','CCL2','CCL8','CD74','CD14','CD79A','MS4A1','MZB1','JCHAIN','IGKC','IGLC2','FOS','JUN','HSPA1A','HSPA1B','MT1G','MT2A','CA9','VEGFA','GLUL','AFP','GPC3','SLC2A1','S100A10','FABP4','RGS5','ACTA2','TAGLN']
EXTRA+=['CYP2A6','ALDOB','RBP4','FGA','COMP','MALAT1','HMGB1','AGR2','IGFBP7','BHMT','GSTA1','CYP2E1','IGFBP2','MTRNR2L8','SAA1','CRP','ORM1']
def array(x):return np.asarray(x.toarray() if sp.issparse(x) else x)
def read(s,g):
 p=ROOT/'inputs'/s;prefix=f'{g}_{s}_'
 a=sc.read_10x_h5(p/(prefix+'filtered_feature_bc_matrix.h5'));a.var_names_make_unique()
 assert a.obs_names.is_unique and a.var['gene_ids'].is_unique
 c=pd.read_csv(p/(prefix+'tissue_positions_list.csv.gz'),header=None,names=['barcode','in_tissue','array_row','array_col','pixel_row','pixel_col']).set_index('barcode')
 assert c.index.is_unique and set(a.obs_names)<=set(c.index)
 a.obs=c.loc[a.obs_names].copy();assert a.obs.in_tissue.eq(1).all()
 a.obsm['spatial']=a.obs[['pixel_col','pixel_row']].to_numpy()
 sf=json.loads(gzip.decompress((p/(prefix+'scalefactors_json.json.gz')).read_bytes()))
 im=np.asarray(Image.open(BytesIO(gzip.decompress((p/(prefix+'tissue_hires_image.png.gz')).read_bytes()))))
 a.uns['spatial']={s:{'images':{'hires':im},'scalefactors':sf}}
 # Sum over all gene-expression features; preserve zeros and use original counts.
 a.obs['sample']=s;a.obs['total_counts']=np.asarray(a.X.sum(1)).ravel();a.obs['n_genes']=a.X.getnnz(axis=1)
 a.obs['pct_mt']=np.asarray(a.X[:,a.var_names.str.startswith('MT-')].sum(1)).ravel()/a.obs.total_counts*100
 a.obs['PGAM5_raw']=array(a.X[:,a.var_names=='PGAM5']).ravel().astype(int)
 a.obs['pass_qc']=(a.obs.total_counts>=500)&(a.obs.n_genes>=200)&(a.obs.pct_mt<=25)
 a.obs['pass_qc_relaxed_mt30']=(a.obs.total_counts>=500)&(a.obs.n_genes>=200)&(a.obs.pct_mt<=30)
 a.obs.to_csv(ROOT/f'{s}_all_spot_QC.csv',index_label='barcode')
 return a
if __name__=='__main__':
 summaries=[];coverage=[]
 (ROOT/'programs.json').write_text(json.dumps(PROGRAMS,indent=2))
 for s,g in SAMPLES.items():
  original=read(s,g);a=original[original.obs.pass_qc].copy();a.layers['counts']=a.X.copy()
  sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a);a.raw=a
  a.obs['PGAM5_log1p_CP10k']=array(a[:,['PGAM5']].X).ravel()
  for k,genes in PROGRAMS.items():
   found=[v for v in genes if v in a.var_names];coverage.append({'sample':s,'program':k,'n_requested':len(genes),'n_present':len(found),'genes_present':';'.join(found)})
   a.obs[k]=array(a[:,found].X).mean(1)
  h=a.copy();sc.pp.filter_genes(h,min_cells=10)
  # Exclude mitochondrial/ribosomal genes from HVG selection, preserve full expression for scoring.
  eligible=~(h.var_names.str.startswith('MT-')|h.var_names.str.match(r'^RP[SL]\d'))
  h=h[:,eligible].copy();sc.pp.highly_variable_genes(h,n_top_genes=2000,flavor='seurat');h=h[:,h.var.highly_variable].copy()
  h.var.to_csv(ROOT/f'{s}_HVGs.csv',index_label='gene')
  sc.pp.scale(h,max_value=10);sc.tl.pca(h,n_comps=40,random_state=245908);sc.pp.neighbors(h,n_neighbors=15,n_pcs=30,random_state=245908)
  sc.tl.leiden(h,resolution=.6,key_added='cluster',random_state=245908,flavor='igraph',n_iterations=2,directed=False)
  sc.tl.umap(h,random_state=245908)
  a.obs['cluster']=pd.Categorical(['S'+str(x) for x in h.obs.cluster]);a.obsm['X_umap']=h.obsm['X_umap'];a.obsm['X_pca']=h.obsm['X_pca']
  a.obsp['connectivities']=h.obsp['connectivities'];a.obsp['distances']=h.obsp['distances'];a.uns['neighbors']=h.uns['neighbors'];a.uns['pca']=h.uns['pca']
  sc.tl.rank_genes_groups(a,groupby='cluster',method='wilcoxon',tie_correct=True,use_raw=False,n_genes=100,pts=True)
  markers=sc.get.rank_genes_groups_df(a,group=None);markers.to_csv(ROOT/f'{s}_cluster_top100_markers.csv',index=False)
  rows=[]
  for c in a.obs.cluster.cat.categories:
   d=a[a.obs.cluster==c];r={'sample':s,'cluster':c,'n_spots':len(d),'PGAM5_detected_spots':int((d.obs.PGAM5_raw>0).sum()),'PGAM5_detection_pct':float((d.obs.PGAM5_raw>0).mean()*100),'PGAM5_mean_log1p_CP10k':float(d.obs.PGAM5_log1p_CP10k.mean()),'PGAM5_mean_raw':float(d.obs.PGAM5_raw.mean())}
   r.update({k:float(d.obs[k].mean()) for k in PROGRAMS});rows.append(r)
  pd.DataFrame(rows).to_csv(ROOT/f'{s}_cluster_profiles.csv',index=False)
  a.write_h5ad(WORK/f'{s}_full_processed.h5ad',compression='gzip')
  panel=list(dict.fromkeys(EXTRA+[x for v in PROGRAMS.values() for x in v]));panel=[x for x in panel if x in a.var_names]
  portable=a[:,panel].copy();portable.raw=None
  portable.uns['spatial'][s]['images']={'lowres':np.asarray(Image.open(BytesIO(gzip.decompress((ROOT/'inputs'/s/f'{g}_{s}_tissue_lowres_image.png.gz').read_bytes()))))}
  portable.uns.pop('rank_genes_groups',None)
  portable.write_h5ad(ROOT/f'{s}_spatial_marker_subset.h5ad',compression='gzip')
  a.obs.to_csv(ROOT/f'{s}_spot_annotations.csv',index_label='barcode',float_format='%.17g')
  r={'sample':s,'gsm':g,'input_spots':len(original),'qc_spots':len(a),'removed_spots':len(original)-len(a),'relaxed_qc_spots':int(original.obs.pass_qc_relaxed_mt30.sum()),'n_genes_matrix':a.n_vars,'n_clusters':a.obs.cluster.nunique(),'PGAM5_detected_spots':int((a.obs.PGAM5_raw>0).sum()),'PGAM5_detection_pct':float((a.obs.PGAM5_raw>0).mean()*100),'median_UMI':float(a.obs.total_counts.median()),'median_genes':float(a.obs.n_genes.median()),'median_pct_mt':float(a.obs.pct_mt.median())};summaries.append(r);print(r,flush=True)
 pd.DataFrame(summaries).to_csv(ROOT/'sample_QC_summary.csv',index=False);pd.DataFrame(coverage).to_csv(ROOT/'program_gene_coverage.csv',index=False)
