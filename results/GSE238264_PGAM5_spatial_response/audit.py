import json,hashlib,itertools
import numpy as np,pandas as pd,h5py
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import false_discovery_control
from analyze import R,INPUT,SAMPLES,PANELS,PRIMARY,SECONDARY
checks=[]
def check(name,ok):checks.append(dict(check=name,passed=bool(ok)));assert ok,name
metrics=pd.read_csv(R/'patient_spatial_metrics.csv');tests=pd.read_csv(R/'response_association_tests.csv');manifest=json.loads((R/'input_manifest.json').read_text())
for rec in manifest:check('SHA '+rec['path'],hashlib.sha256((INPUT/rec['path']).read_bytes()).hexdigest()==rec['sha256'])
for sample in SAMPLES:
 o=pd.read_csv(R/f'{sample}_spot_annotations.csv.gz',index_col=0);root=INPUT/sample
 with h5py.File(root/'filtered_feature_bc_matrix.h5') as f:
  m=f['matrix'];X=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=m['shape'][:]);names=np.array([g.decode() for g in m['features']['name'][:]]);bc=[b.decode() for b in m['barcodes'][:]]
 pos=pd.read_csv(root/'spatial/tissue_positions_list.csv',header=None,index_col=0).loc[bc];total=X.sum(0).A1;ng=(X>0).sum(0).A1;mt=X[np.char.startswith(names,'MT-')].sum(0).A1/total*100;qc=(total>=500)&(ng>=200)&(mt<=25)&pos[1].eq(1).to_numpy();check(sample+' QC',list(o.index)==list(np.array(bc)[qc]));check(sample+' coordinates',np.array_equal(pos.loc[o.index,[4,5]].to_numpy(),o[['pixel_row','pixel_col']].to_numpy()));X=X[:,qc];check(sample+' full-library UMI',np.allclose(X.sum(0).A1,o.total_counts));coords=np.column_stack([o.array_col*np.sqrt(3)/2,o.array_row*1.5]);pairs=cKDTree(coords).query_pairs(np.sqrt(3)+1e-8);expected={(i,j) for i,j in pairs}|{(j,i) for i,j in pairs};lookup={b:i for i,b in enumerate(o.index)};edge=pd.read_csv(R/f'{sample}_spatial_edges.csv.gz');observed={(lookup[a],lookup[b]) for a,b in zip(edge.source_barcode,edge.target_barcode)};check(sample+' geometric graph',observed==expected);ij=np.array(list(expected));G=sparse.csr_matrix((np.ones(len(ij)),(ij[:,0],ij[:,1])),shape=(len(o),len(o)));degree=G.sum(1).A1
 for k,gs in PANELS.items():
  z=np.column_stack([X[names==g].sum(0).A1 for g in gs]);v=np.log1p(z/o.total_counts.to_numpy()[:,None]*1e4);check(sample+k+' normalized program',np.allclose(v.mean(1),o['score_'+k]));check(sample+k+' detection',np.array_equal((z>0).sum(1),o['detected_'+k]))
 mac=(o.score_Macrophage>=o.score_Macrophage.quantile(.75))&(o.detected_Macrophage>=3)&(o.score_Myeloid>=o.score_Myeloid.median())&(degree>0);t=(o.score_T_lineage>=o.score_T_lineage.quantile(.75))&(o.detected_T_lineage>=2)&(degree>0);raw=X[names=='PGAM5'].sum(0).A1;check(sample+' raw PGAM5',np.array_equal(raw,o.PGAM5_raw));det=mac&(raw>0);cut=o.loc[det,'PGAM5_log1p_CP10k'].quantile(.75);high=det&(o.PGAM5_log1p_CP10k>=cut);check(sample+' source gates',np.array_equal(mac,o.Mac_RNA) and np.array_equal(t,o.T_RNA) and np.array_equal(high,o.source_group.eq('PGAM5_high')))
 row=metrics[metrics['sample']==sample].iloc[0];check(sample+' clinical response',row.response==('non-responder' if sample.endswith('NR') else 'responder'));check(sample+' patient counts',row.QC_spots==len(o) and row.mac_spots==mac.sum() and row.PGAM5_high_spots==high.sum());check(sample+' PGAM5 detection metric',np.isclose(row.PGAM5_detected_mac_fraction,det.sum()/mac.sum()))
 for name,vals in [('mac_fraction',np.divide(G@mac.to_numpy(float),degree,out=np.zeros(len(o)),where=degree>0)),('T_fraction',np.divide(G@t.to_numpy(float),degree,out=np.zeros(len(o)),where=degree>0)),('fibroblast_score',np.divide(G@o.score_Fibroblast.to_numpy(),degree,out=np.zeros(len(o)),where=degree>0))]:
  check(sample+name+' per source',np.allclose(vals,o['neighbor_'+name]));v=vals[high].mean() if high.sum()>=3 else np.nan;check(sample+name+' patient aggregation',np.isnan(row['high_neighbor_'+name]) if np.isnan(v) else np.isclose(v,row['high_neighbor_'+name]))
for row in tests.itertuples():
 dd=metrics[metrics[row.metric].notna()];v=dd[row.metric].to_numpy();h=dd.response.eq('responder').to_numpy()
 if row.status=='tested':
  null=[]
  for selected in itertools.combinations(range(7),4):
   mask=np.zeros(7,bool);mask[list(selected)]=True;null.append(v[mask].mean()-v[~mask].mean())
  delta=v[h].mean()-v[~h].mean();p=(np.abs(null)>=abs(delta)-1e-12).mean();check(row.metric+' exact 35-label test',len(null)==35 and np.isclose(p,row.p_exact) and np.isclose(delta,row.difference))
 else:check(row.metric+' missing p',np.isnan(row.p_exact))
for family in ['primary','secondary']:
 dd=tests[tests.family==family];q=false_discovery_control(dd.p_exact.fillna(1).to_numpy());ok=dd.p_exact.notna().to_numpy();check(family+' BH',np.allclose(q[ok],dd.BH_q.to_numpy()[ok]))
check('4 responder and 3 non-responder patients',len(metrics)==7 and metrics.response.eq('responder').sum()==4 and metrics.response.eq('non-responder').sum()==3)
pd.DataFrame(checks).to_csv(R/'audit_checks.csv',index=False);print('PASS',len(checks),'checks')
