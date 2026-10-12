from pathlib import Path
import json,hashlib,itertools,os
import numpy as np,pandas as pd,h5py
from scipy import sparse
from scipy.io import mmread
from scipy.stats import false_discovery_control
R=Path(__file__).resolve().parent;INPUT=Path(os.environ.get('GSE238264_INPUT','/workspace/scratch/GSE238264/inputs'))
SAMPLES=['HCC1R','HCC2R','HCC3R','HCC4R','HCC5NR','HCC6NR','HCC7NR']
PANELS={'Macrophage':['C1QA','C1QB','C1QC','CD68','CSF1R','MERTK','MSR1'],'Myeloid':['LST1','TYROBP','FCER1G','AIF1','LYZ','CTSS'],'T_lineage':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2'],'Fibroblast':['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA'],'Exhaustion':['PDCD1','LAG3','TIGIT','CTLA4','TOX','CXCL13'],'Activation':['CD69','IL2RA','TNFRSF9','ICOS','CD40LG','IL2'],'Effector':['IFNG','PRF1','GZMB','GZMA','CCL5','TBX21'],'Naive_memory':['CCR7','SELL','TCF7','LEF1','IL7R'],'Treg':['FOXP3','IKZF2','IL2RA','CTLA4']}
PRIMARY=['PGAM5_detected_mac_fraction','high_neighbor_mac_fraction','high_neighbor_T_fraction','high_neighbor_fibroblast_score']
SECONDARY=['high_spot_density','high_neighbor_Exhaustion_score','high_neighbor_Activation_score','high_neighbor_Effector_score','high_neighbor_Naive_memory_score','high_neighbor_Treg_score']
def run():
 rows=[];coverage=[];manifest=[];checks=[]
 for sample in SAMPLES:
  root=INPUT/sample
  with h5py.File(root/'filtered_feature_bc_matrix.h5') as f:
   m=f['matrix'];x=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=m['shape'][:]);names=np.array([g.decode() for g in m['features']['name'][:]]);bc=[b.decode() for b in m['barcodes'][:]]
  assert np.all(x.data>=0) and np.all(x.data==np.floor(x.data));mtx=mmread(root/'filtered_feature_bc_matrix/matrix.mtx.gz').tocsc();assert (x!=mtx).nnz==0;fs=pd.read_csv(root/'filtered_feature_bc_matrix/features.tsv.gz',sep='\t',header=None);bs=pd.read_csv(root/'filtered_feature_bc_matrix/barcodes.tsv.gz',header=None)[0].tolist();assert fs[1].tolist()==names.tolist() and bs==bc
  p=pd.read_csv(root/'spatial/tissue_positions_list.csv',header=None,names=['barcode','in_tissue','array_row','array_col','pixel_row','pixel_col']).set_index('barcode');assert p.index.is_unique;pos=p.loc[bc].copy();umi=x.sum(0).A1;ng=(x>0).sum(0).A1;mt=x[np.char.startswith(names,'MT-')].sum(0).A1/umi*100;qc=(umi>=500)&(ng>=200)&(mt<=25)&pos.in_tissue.eq(1).to_numpy();o=pos.iloc[np.flatnonzero(qc)].copy();o['total_counts']=umi[qc];o['n_genes']=ng[qc];o['pct_mt']=mt[qc];n=len(o);xc=x[:,qc]
  genes=list(dict.fromkeys(['PGAM5','CD8A','CD8B','CD4']+[g for gs in PANELS.values() for g in gs]));z=np.column_stack([xc[names==g].sum(0).A1 for g in genes]);y=np.log1p(z/o.total_counts.to_numpy()[:,None]*1e4);o['PGAM5_raw']=z[:,0];o['PGAM5_log1p_CP10k']=y[:,0]
  for k,gs in PANELS.items():
   assert all(np.any(names==g) for g in gs),(sample,k);ix=[genes.index(g) for g in gs];o['score_'+k]=y[:,ix].mean(1);o['detected_'+k]=(z[:,ix]>0).sum(1);coverage.append(dict(sample=sample,panel=k,genes=';'.join(gs),n_present=len(gs)))
  lookup={tuple(v):i for i,v in enumerate(o[['array_row','array_col']].to_numpy())};edges=[]
  for (rr,cc),i in lookup.items():
   for dr,dc in [(0,2),(0,-2),(1,1),(1,-1),(-1,1),(-1,-1)]:
    j=lookup.get((rr+dr,cc+dc))
    if j is not None:edges.append((i,j))
  ij=np.array(edges);G=sparse.csr_matrix((np.ones(len(ij)),(ij[:,0],ij[:,1])),shape=(n,n));degree=G.sum(1).A1;mac=(o.score_Macrophage>=o.score_Macrophage.quantile(.75))&(o.detected_Macrophage>=3)&(o.score_Myeloid>=o.score_Myeloid.median())&(degree>0);t=(o.score_T_lineage>=o.score_T_lineage.quantile(.75))&(o.detected_T_lineage>=2)&(degree>0);det=mac&(o.PGAM5_raw>0);cut=o.loc[det,'PGAM5_log1p_CP10k'].quantile(.75);high=det&(o.PGAM5_log1p_CP10k>=cut);group=np.where(mac,np.where(high,'PGAM5_high',np.where(det,'Other_detected','Undetected')),'Outside_mac_gate');o['source_group']=group;o['Mac_RNA']=mac;o['T_RNA']=t;o['spatial_degree']=degree;nt=G@t.to_numpy(float)
  row=dict(sample=sample,response='non-responder' if sample.endswith('NR') else 'responder',QC_spots=n,input_spots=len(bc),mac_spots=int(mac.sum()),PGAM5_detected_mac_spots=int(det.sum()),PGAM5_high_spots=int(high.sum()),other_detected_spots=int((det&~high).sum()),T_spots=int(t.sum()),high_with_T_neighbors=int((high&(nt>0)).sum()),high_cutoff_logCP10k=float(cut),PGAM5_detected_mac_fraction=float(det.sum()/mac.sum()),high_spot_density=float(high.sum()/n))
  for name,vals in [('mac_fraction',np.divide(G@mac.to_numpy(float),degree,out=np.zeros(n),where=degree>0)),('T_fraction',np.divide(nt,degree,out=np.zeros(n),where=degree>0)),('fibroblast_score',np.divide(G@o.score_Fibroblast.to_numpy(),degree,out=np.zeros(n),where=degree>0))]:
   o['neighbor_'+name]=vals;row['high_neighbor_'+name]=float(vals[high].mean()) if high.sum()>=3 else np.nan;row['other_neighbor_'+name]=float(vals[det&~high].mean()) if (det&~high).sum()>=3 else np.nan
  for k in ['Exhaustion','Activation','Effector','Naive_memory','Treg']:
   vals=np.divide(G@(o['score_'+k].to_numpy()*t.to_numpy()),nt,out=np.full(n,np.nan),where=nt>0);o['neighbor_T_'+k]=vals;valid=high&(nt>0);row['high_neighbor_'+k+'_score']=float(vals[valid].mean()) if valid.sum()>=3 else np.nan
   mask=t&(o['score_'+k]>=o.loc[t,'score_'+k].quantile(.75))&(o['score_'+k]>0)&(o['detected_'+k]>=2)
   if k=='Treg':mask&=z[:,genes.index('FOXP3')]>0
   o['T_'+k+'_supported']=mask
  o.to_csv(R/f'{sample}_spot_annotations.csv.gz',index_label='barcode',float_format='%.17g');pd.DataFrame({'source_barcode':o.index[ij[:,0]],'target_barcode':o.index[ij[:,1]]}).to_csv(R/f'{sample}_spatial_edges.csv.gz',index=False)
  for file in sorted(root.rglob('*')):
   if file.is_file():manifest.append(dict(sample=sample,path=str(file.relative_to(INPUT)),bytes=file.stat().st_size,sha256=hashlib.sha256(file.read_bytes()).hexdigest()))
  checks.append(dict(sample=sample,matrix_formats_equal=True,features_barcodes_equal=True,integer_counts=True,unique_spot_coordinates=len(lookup)==n));rows.append(row);print(sample,row['QC_spots'],'mac',row['mac_spots'],'High',row['PGAM5_high_spots'],flush=True)
 d=pd.DataFrame(rows);d.to_csv(R/'patient_spatial_metrics.csv',index=False);tests=[];audit=[]
 for family,metrics in [('primary',PRIMARY),('secondary',SECONDARY)]:
  for metric in metrics:
   dd=d[d[metric].notna()].copy();v=dd[metric].to_numpy();resp=dd.response.eq('responder').to_numpy();a=int(resp.sum());b=int((~resp).sum());obs=v[resp].mean()-v[~resp].mean() if a and b else np.nan;r=dict(family=family,metric=metric,n_responders=a,n_nonresponders=b,mean_responder=v[resp].mean() if a else np.nan,mean_nonresponder=v[~resp].mean() if b else np.nan,difference=obs,p_exact=np.nan,status='not_testable: require all 7 patients')
   if len(dd)==7:
    null=[]
    for selected in itertools.combinations(range(7),4):
     h=np.zeros(7,bool);h[list(selected)]=True;null.append(v[h].mean()-v[~h].mean())
    r.update(p_exact=float((np.abs(null)>=abs(obs)-1e-12).mean()),status='tested');audit.extend([dict(metric=metric,assignment=i,delta=vv) for i,vv in enumerate(null)])
   tests.append(r)
 table=pd.DataFrame(tests);table['BH_q']=np.nan
 for family in ['primary','secondary']:
  use=table.family.eq(family);q=false_discovery_control(table.loc[use,'p_exact'].fillna(1).to_numpy());table.loc[use,'BH_q']=np.where(table.loc[use,'p_exact'].notna(),q,np.nan)
 table.to_csv(R/'response_association_tests.csv',index=False);pd.DataFrame(audit).to_csv(R/'exact_patient_permutation_audit.csv',index=False);pd.DataFrame(coverage).to_csv(R/'marker_coverage.csv',index=False);pd.DataFrame(checks).to_csv(R/'input_checks.csv',index=False);(R/'input_manifest.json').write_text(json.dumps(manifest,indent=2));(R/'analysis_config.json').write_text(json.dumps(dict(primary_metrics=PRIMARY,secondary_metrics=SECONDARY,panels=PANELS,QC='UMI>=500, detected_features>=200, MT<=25%, in_tissue=1',mac_gate='core Q75, >=3 markers, myeloid >=median, degree>0',high='Q75 of PGAM5 detected mac spots within each patient',minimum_high_sources_for_metrics=3,minimum_high_with_T_for_programs=3,patients=7,responders=4,nonresponders=3,patient_test='all 35 assignments, two-sided absolute mean difference',FDR='4 primary fixed tests; 6 secondary separately, missing p treated as 1 only for correction',timepoint='post-neoadjuvant resection',RNA_proxy=True),indent=2));print(table.to_string(index=False),flush=True)
if __name__=='__main__':run()
