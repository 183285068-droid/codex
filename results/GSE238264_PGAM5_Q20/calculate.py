from pathlib import Path
import os,json,itertools,hashlib,zipfile
os.environ.setdefault('MPLCONFIGDIR','/tmp/gse245908_mpl')
import numpy as np,pandas as pd
from scipy import sparse
from scipy.stats import false_discovery_control
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;B=R.parent/'GSE238264_PGAM5_spatial_response';S=R.parent/'GSE238264_PGAM5_threshold_sensitivity';samples=['HCC1R','HCC2R','HCC3R','HCC4R','HCC5NR','HCC6NR','HCC7NR'];rows=[];checks=[];manifest=[];metrics=['mac_fraction','T_fraction','fibroblast_score']
for sample in samples:
 p=B/f'{sample}_spot_annotations.csv.gz';o=pd.read_csv(p,index_col=0);ep=B/f'{sample}_spatial_edges.csv.gz';edge=pd.read_csv(ep);ix={b:i for i,b in enumerate(o.index)};G=sparse.csr_matrix((np.ones(len(edge)),([ix[b] for b in edge.source_barcode],[ix[b] for b in edge.target_barcode])),shape=(len(o),len(o)));deg=G.sum(1).A1;positive=o.Mac_RNA.to_numpy()&(o.PGAM5_raw.to_numpy()>0);cut=o.loc[positive,'PGAM5_log1p_CP10k'].quantile(.2);high=positive&(o.PGAM5_log1p_CP10k.to_numpy()>=cut);assert high.any()
 values={'T_fraction':np.divide(G@o.T_RNA.to_numpy(float),deg,out=np.zeros(len(o)),where=deg>0),'mac_fraction':np.divide(G@o.Mac_RNA.to_numpy(float),deg,out=np.zeros(len(o)),where=deg>0),'fibroblast_score':np.divide(G@o.score_Fibroblast.to_numpy(),deg,out=np.zeros(len(o)),where=deg>0)}
 row=dict(sample=sample,response='non-responder' if sample.endswith('NR') else 'responder',quantile=.2,threshold=float(cut),n_detected=int(positive.sum()),n_high=int(high.sum()),low_support=bool(high.sum()<3))
 for metric,v in values.items():
  assert np.allclose(v,o['neighbor_'+metric]);row[metric]=float(v[high].mean());checks.append(dict(sample=sample,check='independent spatial aggregation '+metric,passed=True))
 q25=positive&(o.PGAM5_log1p_CP10k.to_numpy()>=o.loc[positive,'PGAM5_log1p_CP10k'].quantile(.25));assert np.all(high[q25]);checks.append(dict(sample=sample,check='Q25 sources nested in Q20',passed=True));o[['PGAM5_raw','PGAM5_log1p_CP10k']].assign(Q20_source=high).to_csv(R/f'{sample}_Q20_labels.csv.gz',index_label='barcode');rows.append(row)
 for pp in [p,ep]:manifest.append(dict(path=str(pp.relative_to(R.parent)),sha256=hashlib.sha256(pp.read_bytes()).hexdigest()))
d=pd.DataFrame(rows);d.to_csv(R/'Q20_slice_metrics.csv',index=False);tests=[];audit=[];h=d.response.eq('responder').to_numpy();assert h.sum()==4 and (~h).sum()==3
for metric in metrics:
 v=d[metric].to_numpy();delta=float(v[h].mean()-v[~h].mean());null=[]
 for i,c in enumerate(itertools.combinations(range(7),4)):
  mask=np.zeros(7,bool);mask[list(c)]=True;z=float(v[mask].mean()-v[~mask].mean());null.append(z);audit.append(dict(metric=metric,assignment=i,delta=z))
 count=int((np.abs(null)>=abs(delta)-1e-12).sum());p=count/35;assert len(null)==35;tests.append(dict(quantile=.2,metric=metric,n_responders=4,n_nonresponders=3,mean_responder=v[h].mean(),mean_nonresponder=v[~h].mean(),difference=delta,p_exact=p,extreme_assignments=count))
t=pd.DataFrame(tests);old=pd.read_csv(S/'slice_response_tests.csv');family=pd.concat([t,old[['quantile','metric','n_responders','n_nonresponders','mean_responder','mean_nonresponder','difference','p_exact']]],ignore_index=True);assert len(family)==15;family['BH_q_15_tests']=false_discovery_control(family.p_exact.to_numpy());t['BH_q_15_tests']=family.BH_q_15_tests.iloc[:3].to_numpy();t.to_csv(R/'Q20_response_tests.csv',index=False);family.to_csv(R/'all_threshold_response_tests.csv',index=False);pd.DataFrame(audit).to_csv(R/'exact_permutation_audit.csv',index=False);pd.DataFrame(checks).to_csv(R/'calculation_checks.csv',index=False);(R/'input_manifest.json').write_text(json.dumps(manifest,indent=2))
# Independent p recalculation using excluded nonresponder triples.
for row in t.itertuples():
 v=d[row.metric].to_numpy();observed=v[h].mean()-v[~h].mean();counts=0
 for nr in itertools.combinations(range(7),3):
  select=np.ones(7,bool);select[list(nr)]=False;counts+=abs(v[select].mean()-v[~select].mean())>=abs(observed)-1e-12
 assert np.isclose(counts/35,row.p_exact)
fig,ax=plt.subplots(figsize=(8,5))
for i,g in enumerate(['responder','non-responder']):
 dd=d[d.response==g];xx=i+np.linspace(-.14,.14,len(dd));ax.scatter(xx,dd.T_fraction*100,color='#2678a8' if i==0 else '#b95e27',s=55)
 for x,row in zip(xx,dd.itertuples()):ax.annotate(f'{row.sample} (n={row.n_high})',(x,row.T_fraction*100),xytext=(4,14 if row.sample=='HCC7NR' else 4),textcoords='offset points',fontsize=9)
 ax.hlines(dd.T_fraction.mean()*100,i-.2,i+.2,color='#2678a8' if i==0 else '#b95e27',lw=2)
r=t[t.metric=='T_fraction'].iloc[0];ax.set_xlim(-.4,1.6);ax.set_ylim(0,max(d.T_fraction)*100+9);ax.set_xticks([0,1],['Responder (4 slices)','Non-responder (3 slices)']);ax.set_ylabel('Mean T-RNA neighbor coverage (%)');ax.set_title(f'PGAM5 Q20 | one dot per slice/patient\nExact p={r.p_exact:.5f}; 15-comparison BH q={r.BH_q_15_tests:.5f}');fig.tight_layout();fig.savefig(R/'figures/Q20_T_RNA_coverage.png',dpi=180);fig.savefig(R/'figures/Q20_T_RNA_coverage.pdf');plt.close(fig)
(R/'README.md').write_text(f'''# GSE238264：Q20切片邻近T-RNA覆盖

按用户指定Q20：每张切片PGAM5检出巨噬RNA源内的log1p(CP10k)第20百分位及以上，通常纳入检出源的约80%，零计数不纳入。不是全部巨噬点的第20百分位，也不是跨样本统一绝对表达阈值。

沿用原QC、巨噬及T门控和真实六角邻接；每源先计算T-RNA邻居数/全部保留邻居数，再切片内平均，最后组内切片等权平均。7张切片对应7位患者，4应答/3未应答；至少1个源点即可汇总，n<3仍标记低支持。

Q20应答组均值{100*r.mean_responder:.4f}%，未应答组{100*r.mean_nonresponder:.4f}%，差值{100*r.difference:.4f}个百分点。全部35种4/3分配的双侧绝对均值差检验，极端配置{int(r.extreme_assignments)}种，p={r.p_exact:.8f}。这是精确值，没有加Monte Carlo修正。追加Q20三指标与此前四阈值三指标共15项联合BH后，T覆盖q={r.BH_q_15_tests:.8f}；先前12项q不直接移用到Q20。

这些是用户追加阈值的探索性比较，不作选取最小p的确认性推断。差异是否显著以0.05阈值判断；治疗后混合空间RNA与病理应答的关联不能等同细胞机制、T细胞比例或治疗前预测。HCC2R原始PGAM5检出源只有3个，门槛放宽不消除稀疏性。

逐切片数值见Q20_slice_metrics.csv，参数见Q20_response_tests.csv，全部阈值校正见all_threshold_response_tests.csv。已从原始邻接边独立重算覆盖并验证缓存一致；所有Q25源纳入Q20；反向枚举三例未应答标签独立复算p。PNG/PDF见figures。源输入哈希已记录。

在完整仓库中使用原分析requirements.txt环境运行python calculate.py即可复现，不需要重新下载原始计数。代码依赖相邻原分析和阈值敏感性目录；结果ZIP不重复原始文件。
''');(R/'.gitignore').write_text('__pycache__/\n');z=R/'GSE238264_PGAM5_Q20_all_results.zip';files=[p for p in sorted(R.rglob('*')) if p.is_file() and p!=z and '__pycache__' not in p.parts]
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as f:
 for p in files:f.write(p,p.relative_to(R))
with zipfile.ZipFile(z) as f:assert f.testzip() is None
print(d.to_string(index=False));print(t.to_string(index=False));print('All checks passed; ZIP bytes',z.stat().st_size)
