from pathlib import Path
import os,json,itertools,hashlib,zipfile
import numpy as np,pandas as pd
from scipy import sparse
from scipy.stats import false_discovery_control
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
R=Path(__file__).resolve().parent;BASE=R.parent/'GSE238264_PGAM5_spatial_response';RAW=Path(os.environ.get('GSE238264_INPUT','/workspace/scratch/GSE238264/inputs'));QS=[.25,.5,.6,.75];METRICS=['mac_fraction','T_fraction','fibroblast_score'];rows=[];checks=[];manifest=[]
for sample in ['HCC1R','HCC2R','HCC3R','HCC4R','HCC5NR','HCC6NR','HCC7NR']:
 p=BASE/f'{sample}_spot_annotations.csv.gz';o=pd.read_csv(p,index_col=0);edge=pd.read_csv(BASE/f'{sample}_spatial_edges.csv.gz');ix={b:i for i,b in enumerate(o.index)};G=sparse.csr_matrix((np.ones(len(edge)),([ix[b] for b in edge.source_barcode],[ix[b] for b in edge.target_barcode])),shape=(len(o),len(o)));deg=G.sum(1).A1;mac=o.Mac_RNA.to_numpy();positive=mac&(o.PGAM5_raw.to_numpy()>0);fig,axs=plt.subplots(1,4,figsize=(17,5));root=RAW/sample/'spatial';im=np.array(Image.open(root/'tissue_hires_image.png').convert('RGB'));scale=json.loads((root/'scalefactors_json.json').read_text());xy=o[['pixel_col','pixel_row']].to_numpy()*scale['tissue_hires_scalef'];assert (xy>=0).all() and (xy[:,0]<im.shape[1]).all() and (xy[:,1]<im.shape[0]).all();group='non-responder' if sample.endswith('NR') else 'responder'
 for ax,q in zip(axs,QS):
  cutoff=o.loc[positive,'PGAM5_log1p_CP10k'].quantile(q);high=positive&(o.PGAM5_log1p_CP10k.to_numpy()>=cutoff);row=dict(sample=sample,response=group,quantile=q,threshold=cutoff,n_detected=int(positive.sum()),n_high=int(high.sum()),low_support=bool(high.sum()<3))
  for metric in METRICS:
   vals=o['neighbor_'+metric].to_numpy();row[metric]=float(vals[high].mean()) if high.any() else np.nan
  check=bool(np.array_equal(high,o.source_group.eq('PGAM5_high').to_numpy())) if q==.75 else bool(high.sum()>=o.source_group.eq('PGAM5_high').sum());checks.append(dict(sample=sample,quantile=q,passed=check));assert check
  rows.append(row);o['high_Q'+str(int(q*100))]=high;ax.imshow(im,origin='upper');ax.scatter(xy[high,0],xy[high,1],s=22,facecolors='none',edgecolors='#e65e19',lw=1);ax.set_axis_off();ax.set_title(f'Q{int(q*100)} | High={high.sum()} / detected={positive.sum()}');
 fig.suptitle(sample+' | '+group+' | PGAM5 threshold sensitivity over H&E');fig.tight_layout();fig.savefig(R/'figures'/f'{sample}_HE_thresholds.png',dpi=180,bbox_inches='tight');fig.savefig(R/'figures'/f'{sample}_HE_thresholds.pdf',bbox_inches='tight');plt.close(fig);o[['high_Q25','high_Q50','high_Q60','high_Q75']].to_csv(R/f'{sample}_threshold_labels.csv.gz',index_label='barcode');manifest.append(dict(path=str(p.relative_to(R.parent)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
d=pd.DataFrame(rows);d.to_csv(R/'slice_metrics_by_threshold.csv',index=False);tests=[];nullrows=[]
for q in QS:
 dd=d[d['quantile']==q]
 for metric in METRICS:
  valid=dd[dd[metric].notna()];v=valid[metric].to_numpy();h=valid.response.eq('responder').to_numpy();delta=v[h].mean()-v[~h].mean();row=dict(quantile=q,metric=metric,n_responders=int(h.sum()),n_nonresponders=int((~h).sum()),mean_responder=v[h].mean(),mean_nonresponder=v[~h].mean(),difference=delta,p_exact=np.nan)
  if len(valid)==7:
   null=[]
   for i,c in enumerate(itertools.combinations(range(7),4)):
    a=np.zeros(7,bool);a[list(c)]=True;z=v[a].mean()-v[~a].mean();null.append(z);nullrows.append(dict(quantile=q,metric=metric,assignment=i,delta=z))
   row['p_exact']=float((np.abs(null)>=abs(delta)-1e-12).mean());assert len(null)==35
  tests.append(row)
t=pd.DataFrame(tests);t['BH_q_12_tests']=false_discovery_control(t.p_exact.fillna(1).to_numpy());t.to_csv(R/'slice_response_tests.csv',index=False);pd.DataFrame(nullrows).to_csv(R/'exact_permutation_audit.csv',index=False);pd.DataFrame(checks).to_csv(R/'threshold_checks.csv',index=False);(R/'input_manifest.json').write_text(json.dumps(manifest,indent=2))
fig,axs=plt.subplots(4,3,figsize=(14,14))
for i,q in enumerate(QS):
 for j,metric in enumerate(METRICS):
  ax=axs[i,j];dd=d[d['quantile']==q];factor=100 if metric.endswith('fraction') else 1
  for k,g in enumerate(['responder','non-responder']):
   ss=dd[dd.response==g];xx=k+np.linspace(-.15,.15,len(ss));ax.scatter(xx,ss[metric]*factor,c='#2678a8' if k==0 else '#b95e27')
   for x,row in zip(xx,ss.itertuples()):ax.annotate(f'{row.sample} (n={row.n_high})',(x,getattr(row,metric)*factor),xytext=(3,14 if k==1 and row.sample=='HCC7NR' and getattr(row,metric)<.01 else 3),textcoords='offset points',fontsize=7)
  rr=t[(t['quantile']==q)&(t.metric==metric)].iloc[0];ax.set_title(f'Q{int(q*100)} | {metric}\np={rr.p_exact:.3f}; joint q={rr.BH_q_12_tests:.3f}',fontsize=9);ax.set_xticks([0,1],['Responder (4)','Non-responder (3)']);ax.set_xlim(-.4,1.7);ax.set_ylim(bottom=0);ax.set_ylabel('%' if factor==100 else 'Mean log1p(CP10k)')
fig.suptitle('Exploratory threshold sensitivity | each dot is one slice/patient\nSource counts shown; low-count slice estimates remain unstable');fig.tight_layout();fig.savefig(R/'figures/threshold_response_comparisons.png',dpi=180,bbox_inches='tight');fig.savefig(R/'figures/threshold_response_comparisons.pdf',bbox_inches='tight');plt.close(fig)
print(d[['sample','quantile','n_detected','n_high']].to_string(index=False));print(t.to_string(index=False))
