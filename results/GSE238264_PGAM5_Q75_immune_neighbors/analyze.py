from pathlib import Path
import os,json,hashlib
import numpy as np,pandas as pd,h5py
from scipy import sparse
from scipy.stats import false_discovery_control
import frozen_spatial_methods as f
R=Path(__file__).resolve().parent;BASE=R.parent/'GSE238264_PGAM5_spatial_response';RAW=Path(os.environ.get('GSE238264_INPUT','/workspace/scratch/GSE238264/inputs'));f.ROOT=R
SAMPLES=['HCC1R','HCC2R','HCC3R','HCC4R','HCC5NR','HCC6NR','HCC7NR']
PANELS={'T':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2'],'B':['MS4A1','CD79A','CD79B','CD19','CD22'],'Plasma':['MZB1','JCHAIN','SDC1','XBP1','DERL3'],'NK':['NKG7','GNLY','KLRD1','KLRF1','NCR1'],'cDC1':['CLEC9A','XCR1','BATF3','IRF8','CADM1'],'cDC2':['CD1C','FCER1A','CLEC10A','CD1E'],'pDC':['CLEC4C','IL3RA','TCF4','GZMB'],'Monocyte':['FCN1','VCAN','S100A8','S100A9','LYZ'],'Neutrophil':['FCGR3B','CSF3R','CXCR2','S100A8','S100A9']}
ANCHORS={'B':['MS4A1','CD79A'],'Plasma':['MZB1','JCHAIN'],'NK':['GNLY','KLRD1','NCR1'],'cDC1':['CLEC9A','XCR1'],'cDC2':['CD1C','FCER1A'],'pDC':['CLEC4C','IL3RA'],'Monocyte':['FCN1','VCAN'],'Neutrophil':['FCGR3B','CSF3R','CXCR2']}
PAIRS=[('SPP1','CD44'),('LGALS9','HAVCR2'),('PVR','TIGIT'),('NECTIN2','TIGIT'),('CD274','PDCD1'),('CXCL9','CXCR3'),('CXCL10','CXCR3'),('CCL5','CCR5'),('CXCL12','CXCR4'),('CXCL13','CXCR5'),('TNFSF13B','TNFRSF13C'),('TNFSF13','TNFRSF17'),('TGFB1','TGFBR1_TGFBR2'),('IL10','IL10RA_IL10RB')]
GROUPS=['PGAM5_high','Other_detected','Undetected']
def run():
 stats=[];summ=[];lrs=[];cover=[];support=[];manifest=[];resource=pd.read_csv(R.parent/'GSE245908_PGAM5_Tcell_spatial_interactions/source_inputs/omni_resource.csv');selected=[]
 for l,r in PAIRS:
  hit=resource[(resource.source_genesymbol==l)&(resource.target_genesymbol==r)];assert len(hit);selected.append(dict(ligand=l,receptor=r,resources=';'.join(sorted(hit.resource.unique())),resource_snapshot='LIANA 1.10.0 omni_resource.csv'))
 pd.DataFrame(selected).to_csv(R/'selected_LR_pairs.csv',index=False)
 for si,sample in enumerate(SAMPLES):
  p=BASE/f'{sample}_spot_annotations.csv.gz';o=pd.read_csv(p,index_col=0);edges=pd.read_csv(BASE/f'{sample}_spatial_edges.csv.gz');ix={b:i for i,b in enumerate(o.index)};G=sparse.csr_matrix((np.ones(len(edges)),([ix[b] for b in edges.source_barcode],[ix[b] for b in edges.target_barcode])),shape=(len(o),len(o)));deg=G.sum(1).A1;n=len(o)
  hp=RAW/sample/'filtered_feature_bc_matrix.h5'
  with h5py.File(hp) as h:
   m=h['matrix'];raw=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=m['shape'][:]);names=np.array([g.decode() for g in m['features']['name'][:]]);bc=[b.decode() for b in m['barcodes'][:]]
  bi={b:i for i,b in enumerate(bc)};raw=raw[:,[bi[b] for b in o.index]];assert np.allclose(raw.sum(0).A1,o.total_counts)
  genes=list(dict.fromkeys([g for gs in PANELS.values() for g in gs]+[l for l,r in PAIRS]+[g for l,r in PAIRS for g in r.split('_')]));count=np.column_stack([raw[names==g].sum(0).A1 for g in genes]);cp=count/o.total_counts.to_numpy()[:,None]*1e4;x=np.log1p(cp)
  masks={};cutoffs={}
  for k,gs in PANELS.items():
   present=[g for g in gs if g in names];cover.append(dict(sample=sample,panel=k,present=';'.join(present),missing=';'.join(g for g in gs if g not in names)));indices=[genes.index(g) for g in gs];score=x[:,indices].mean(1);nd=(count[:,indices]>0).sum(1);cutoff=float(np.quantile(score,.75));mask=(score>=cutoff)&(score>0)&(nd>=2)&(deg>0)
   if k=='T':mask=o.T_RNA.to_numpy(bool)
   else:
    anchors=(count[:,[genes.index(g) for g in ANCHORS[k]]]>0)
    mask&=anchors.all(1) if k=='Plasma' else anchors.any(1)
   masks[k]=mask;cutoffs[k]=cutoff;o['score_'+k]=score;o['gate_'+k]=mask
  mac=o.Mac_RNA.to_numpy(bool);pos=mac&(o.PGAM5_raw.to_numpy()>0);high=o.source_group.eq('PGAM5_high').to_numpy();assert np.array_equal(high,pos&(o.PGAM5_log1p_CP10k>=o.loc[pos,'PGAM5_log1p_CP10k'].quantile(.75)).to_numpy());Y=np.column_stack([np.divide(G@mask.astype(float),deg,out=np.zeros(n),where=deg>0) for mask in masks.values()])
  rng=np.random.default_rng(23826475+si);ss,info=f.permutation_results(o,pos,high,Y,list(PANELS),'immune_coverage',rng,sample)
  if info['movable_sources']==0:
   for row in ss:row['spatial_conditional_p']=np.nan;row['status']='not_testable: no exchangeable source labels'
  stats+=ss;support.append(info)
  for j,(k,mask) in enumerate(masks.items()):
   nt=G@mask.astype(float);o['neighbor_fraction_'+k]=Y[:,j]
   for group in GROUPS:
    use=o.source_group.eq(group).to_numpy();summ.append(dict(sample=sample,response='non-responder' if sample.endswith('NR') else 'responder',subtype=k,n_target_spots=int(mask.sum()),source_group=group,n_sources=int(use.sum()),n_sources_with_target=int((use&(nt>0)).sum()),n_edges=int(nt[use].sum()),mean_neighbor_fraction=float(Y[use,j].mean())))
  ri,ci=G.nonzero();eligible=mac[ri]&np.column_stack(list(masks.values())).any(1)[ci];ri=ri[eligible];ci=ci[eligible];e=pd.DataFrame({'source_barcode':o.index[ri],'target_barcode':o.index[ci],'source_group':o.source_group.iloc[ri].to_numpy()})
  for k,mask in masks.items():e['target_'+k]=mask[ci]
  for l,r in PAIRS:
   rec=cp[:,[genes.index(g) for g in r.split('_')]].min(1);score=np.sqrt(cp[ri,genes.index(l)]*rec[ci]);e[l+'--'+r]=score
   for k,mask in masks.items():
    nt=G@mask.astype(float);selection=mask[ci];agg=np.bincount(ri[selection],weights=score[selection],minlength=n);source_mean=np.divide(agg,nt,out=np.full(n,np.nan),where=nt>0)
    for group in GROUPS:
     keep=e.source_group.eq(group).to_numpy()&selection;v=score[keep];use=o.source_group.eq(group).to_numpy()&(nt>0);lrs.append(dict(sample=sample,subtype=k,source_group=group,ligand=l,receptor=r,n_edges=len(v),n_sources=len(np.unique(ri[keep])),n_targets=len(np.unique(ci[keep])),codetected_edges=int((v>0).sum()),codetection_pct=100*(v>0).mean() if len(v) else np.nan,mean_source_RNA_compatibility=float(source_mean[use].mean()) if use.any() else np.nan))
  o.to_csv(R/f'{sample}_immune_spots.csv.gz',index_label='barcode',float_format='%.17g');e.to_csv(R/f'{sample}_directed_LR_edges.csv.gz',index=False)
  for path in [p,BASE/f'{sample}_spatial_edges.csv.gz',hp]:manifest.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),bytes=path.stat().st_size))
  print(sample,'High',high.sum(),'targets', {k:int(v.sum()) for k,v in masks.items()},flush=True)
 t=pd.DataFrame(stats);t['BH_q_63_planned']=false_discovery_control(t.spatial_conditional_p.fillna(1).to_numpy());t.loc[t.spatial_conditional_p.isna(),'BH_q_63_planned']=np.nan;assert len(t)==63;t.to_csv(R/'spatial_conditional_tests.csv',index=False);pd.DataFrame(summ).to_csv(R/'source_immune_neighborhoods.csv',index=False);pd.DataFrame(lrs).to_csv(R/'LR_by_immune_subtype.csv',index=False);pd.DataFrame(support).to_csv(R/'permutation_support.csv',index=False);pd.DataFrame(cover).to_csv(R/'marker_coverage.csv',index=False);(R/'input_manifest.json').write_text(json.dumps(manifest,indent=2));(R/'analysis_config.json').write_text(json.dumps(dict(PGAM5_quantile=.75,source_gates='original frozen baseline',panels=PANELS,anchors=ANCHORS,Plasma_requires_both_MZB1_JCHAIN=True,immune_labels_nonexclusive=True,minimum_high_other_for_test=3,permutations=4999,seed=23826475,FDR_family=63,LR='14 resource-supported candidates; descriptive only',space='hexagonal first ring; 3x3 tiles and per-tile depth median stratification',clinical_response='labels retained for context, no new response-group tests'),indent=2));print('Computed 63 fixed tests; significant q<0.05:',int(t.BH_q_63_planned.lt(.05).sum()))
if __name__=='__main__':run()
