from pathlib import Path
import json,hashlib,platform,importlib.metadata as md
import numpy as np,pandas as pd,anndata as ad
from scipy.optimize import lsq_linear,nnls
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;S=Path('/workspace/scratch');rng=np.random.default_rng(20261009);cells=[];counts=[];sources=[]
for n in ['GSE151530','GSE149614','GSE189903','GSE242889','GSE202642']:
 o=pd.read_csv(R/(n+'_all_cell_annotation.csv.gz'));path=S/'PGAM5_myeloid_reference_v2'/(n+'_harmonized_counts_QC.h5ad') if n in ['GSE151530','GSE149614','GSE189903'] else S/n/('all_tumor_counts_QC.h5ad' if n=='GSE242889' else 'all_library_counts_QC.h5ad');sources.append(path);a=ad.read_h5ad(path,backed='r');pgidx=np.flatnonzero((a.var.gene.astype(str).to_numpy() if 'gene' in a.var else a.var_names.to_numpy())=='PGAM5');assert len(pgidx)>0
 assert o.cell_id.is_unique and o.cell_id.isin(a.obs_names).all();idx=np.r_[rng.choice(np.flatnonzero(o.cell_type.eq('Macrophage')),10,replace=False),rng.choice(np.flatnonzero(~o.cell_type.eq('Macrophage')),10,replace=False)]
 for j in idx:
  z=o.iloc[j];row=a.obs_names.get_loc(z.cell_id);value=float(a.X[row,:][:,pgidx].sum());assert value==z.PGAM5_counts
  expected=('Macrophage_PGAM5_detected' if value>0 else 'Macrophage_PGAM5_undetected') if z.cell_type=='Macrophage' else z.cell_type;assert z.PGAM5_macrophage_annotation==expected;cells.append({'dataset':n,'cell_id':z.cell_id,'source_PGAM5_count':value,'stored_PGAM5_count':z.PGAM5_counts,'annotation':expected,'matched':True})
 mac=o.cell_type.eq('Macrophage');det=mac&o.PGAM5_counts.gt(0);assert o.loc[det,'PGAM5_RNA_status'].eq('detected').all() and o.loc[mac&~det,'PGAM5_RNA_status'].eq('undetected').all();counts.append({'dataset':n,'QC_HCC_cells':len(o),'macrophages':int(mac.sum()),'PGAM5_detected_macrophages':int(det.sum()),'PGAM5_undetected_macrophages':int((mac&~det).sum()),'donor_units':o.donor.nunique(),'donor_unit_caveat':'sample proxies' if n=='GSE202642' else 'previous author/inferred patient labels'});a.file.close()
pd.DataFrame(cells).to_csv(R/'independent_annotation_raw_count_checks.csv',index=False);ct=pd.DataFrame(counts);ct.to_csv(R/'annotation_dataset_summary.csv',index=False)
solverchecks=[]
for f in sorted(R.glob('solver_audit_case_*.npz')):
 z=np.load(f);A=z['A'];b=z['b'];solver=str(z['solver']);ti=int(z['target_index']);w=z['weights']
 if solver!='NNLS':b=np.r_[b,np.zeros(A.shape[1])];A=np.vstack([A,np.sqrt(.01)*np.eye(A.shape[1])])
 alt=lsq_linear(A,b,bounds=(0,np.inf),method='bvls',tol=1e-10,max_iter=1000);assert alt.success;weights=alt.x/max(alt.x.sum(),1e-12);err=float(abs(weights-w).max());assert err<1e-5;solverchecks.append({'file':f.name,'solver':solver,'max_normalized_coefficient_difference':err,'target_prediction':float(w[ti]),'independent_target_prediction':float(weights[ti])})
pd.DataFrame(solverchecks).to_csv(R/'independent_solver_checks.csv',index=False)
ref=pd.read_csv(R/'EXPLORATORY_REFERENCE_NOT_VALIDATED.tsv',sep='\t',index_col=0);scale=pd.read_csv(R/'reference_row_scales.csv');assert ref.index.tolist()==scale.gene.tolist();A=ref.to_numpy()/scale.row_scale.to_numpy()[:,None];exact=[]
for i in range(20):
 w=rng.dirichlet(np.ones(A.shape[1]));b=A@w;coef,res=nnls(A,b);err=float(abs(coef-w).max());assert err<1e-6;exact.append({'case':i,'max_error':err,'residual':res})
pd.DataFrame(exact).to_csv(R/'exact_reference_recovery_checks.csv',index=False)
pred=pd.read_csv(R/'validation_predictions.csv');su=pd.read_csv(R/'validation_summary.csv');summarychecks=[]
for _,z in su.iterrows():
 t=pred[(pred.dataset==z.dataset)&(pred.phase==z.phase)&(pred.solver==z.solver)&(pred.unit==z.unit)&(pred.scenario==z.scenario)];mae=float(abs(t.predicted-t.truth).mean()*100);assert np.isclose(mae,z.MAE_pp,atol=1e-12);assert t.donor.nunique()==z.donors;summarychecks.append({'dataset':z.dataset,'solver':z.solver,'unit':z.unit,'scenario':z.scenario,'recomputed_MAE_pp':mae,'matched':True})
pd.DataFrame(summarychecks).to_csv(R/'independent_summary_checks.csv',index=False);audit={'raw_count_cell_annotations_checked':len(cells),'independent_bounded_solver_cases':len(solverchecks),'exact_reference_mixtures_checked':len(exact),'validation_summary_rows_checked':len(summarychecks),'status':'PASSED within checked implementation scope; biological validity separate'};(R/'implementation_audit_summary.json').write_text(json.dumps(audit,indent=2))
fig,axs=plt.subplots(1,2,figsize=(11,4.8),constrained_layout=True)
for ax,unit in zip(axs,['equal_RNA_cell_fraction','pooled_count_RNA_fraction']):
 for n,c in zip(['GSE189903','GSE242889','GSE202642'],['#b35c33','#247ba0','#558b2f']):
  t=pred[(pred.dataset==n)&(pred.solver=='NNLS')&(pred.unit==unit)&(pred.scenario=='standard')];ax.scatter(t.truth*100,t.predicted*100,s=20,alpha=.65,label=n)
 ax.plot([0,12],[0,12],'k--',lw=1);ax.set_xlabel('Known target proportion (%)');ax.set_ylabel('Estimated target proportion (%)');ax.set_title('Equalized RNA per cell' if unit.startswith('equal') else 'Observed library RNA contribution');ax.legend(fontsize=8);ax.grid(alpha=.2)
fig.savefig(R/'external_target_recovery.png',dpi=180);fig.savefig(R/'external_target_recovery.pdf');plt.close(fig)
versions={'python':platform.python_version(),**{k:md.version(k) for k in ['numpy','pandas','anndata','scanpy','scipy','scikit-learn','statsmodels','matplotlib']}};(R/'environment_versions.json').write_text(json.dumps(versions,indent=2));(R/'requirements.txt').write_text('\n'.join(k+'=='+v for k,v in versions.items() if k!='python')+'\n')
sha={}
for p in sources:
 h=hashlib.sha256()
 with p.open('rb') as stream:
  for b in iter(lambda:stream.read(1024*1024),b''):h.update(b)
 sha[str(p)]={'bytes':p.stat().st_size,'sha256':h.hexdigest()}
(R/'raw_source_sha256.json').write_text(json.dumps(sha,indent=2))
res=json.loads((R/'result.json').read_text());lines=['# PGAM5 RNA检出巨噬细胞注释与反卷积验证','', '**已完成可复现的RNA检出状态注释；当前参考未通过外部验证，不能作为TCGA-LIHC的PGAM5⁺巨噬细胞丰度工具。**','', '## 注释定义与覆盖','', '目标标签`Macrophage_PGAM5_detected`=已有保守谱系标记确认的巨噬细胞且符号合并后的PGAM5原始RNA计数>0；对应未检出标签不是蛋白阴性。5个*_all_cell_annotation.csv.gz可按cell_id回填Seurat/AnnData。保留cell_type、PGAM5_RNA_status、PGAM5_macrophage_annotation和cycling_flag四个独立字段。', '', '| 数据集 | HCC QC细胞 | 巨噬细胞 | PGAM5检出巨噬细胞 | 未检出巨噬细胞 |','|---|---:|---:|---:|---:|']
for _,z in ct.iterrows():lines.append(f'| {z.dataset} | {z.QC_HCC_cells} | {z.macrophages} | {z.PGAM5_detected_macrophages} | {z.PGAM5_undetected_macrophages} |')
lines+=['','GSE202642限定末尾库编号5–11的7个HCC肿瘤样本，编号映射沿用已核实的上一轮来源记录；7个样本不等于已确认7名独立患者。GSE149614只纳入原发T样本，GSE189903仅HCC core/border，GSE242889含所有5例HCC且不按MVI筛选。取消的GSE154906未重启。','', '## 参考与锁定验证','', '沿用GSE151530+GSE149614联合训练的donor_equal_50参考：850基因、17组分，包含目标/未检出巨噬细胞、增殖组分、单核、DC/pDC、Mixed_APC及肿瘤/其他谱系与Unknown_other。低支持组会回退到Unknown_other；部分组跨患者支持有限。特征和配置来自之前发现集内验证，此次没有使用新增外部性能选基因或算法。固定NNLS与非负岭回归alpha=.01，后者不是CIBERSORTx/MuSiC/BayesPrism。所有队列此前被探索过，不是前瞻性盲法验证。','', '新增加242889/202642完整HCC细胞的外部伪bulk。其非目标群以已有大簇和PGAM5以外标记保守注释；未知细胞全部保留而未丢弃。未知比例高、恶性上皮没有CNV确认，竞争标签是推断参考而非金标准。目标巨噬细胞成员沿用上一轮PGAM5-free髓系身份筛选，不能用其他PGAM5表达细胞替代目标。每位合格患者/样本必须有>=5目标、>=20未检出巨噬细胞、>=20肿瘤/上皮代理及>=20 T/NK；不合格者仍保留注释并记录覆盖，但不能贡献完整剂量验证。','', 'GSE242889全部850特征可测；GSE202642缺少7个旧符号，对查询拟合只使用可测的843个，未把缺失特征当表达0。特征可测率要求>=95%；不根据表达/性能选择保留行。缺失名见*_missing_reference_genes.csv，未来TCGA也必须处理基因版本匹配。','', '固定2000细胞伪bulk，总巨噬细胞25%，目标0/0.5/1/2/5%，其余50%肿瘤/上皮代理、25%其他谱系；每情形重复2次。额外测试无目标时单核/DC/pDC/Mixed_APC、增殖及PGAM5检出非巨噬细胞富集。抽样重复不是独立患者。内部固定配置LODO是补充诊断，不替代之前嵌套NNLS验证。','', '预设外部门槛沿用之前参考验证：MAE<=1百分点、相关>=.7、零目标预测95分位<=.5%；每个主要外部队列>=3独立患者，单核/DC等关键挑战>=2可测试患者。缺失的挑战不算通过。','', '## 外部标准混合物结果','', '| 队列 | 算法 | 单位 | 可测试患者/样本 | MAE百分点 | Pearson r | 零目标预测95分位% |','|---|---|---|---:|---:|---:|---:|']
for _,z in su[(su.scenario=='standard')&su.dataset.ne('joint_discovery')].iterrows():lines.append(f'| {z.dataset} | {z.solver} | {z.unit} | {z.donors} | {z.MAE_pp:.3f} | {z.Pearson_r:.3f} | {z.zero_target_p95_percent:.3f} |')
lines+=['','![External target recovery](external_target_recovery.png)','', '两种固定求解器的主要外部综合门槛：'+json.dumps(res['solver_external_gates'])+'。完整失败项见predeclared_gate_checks.csv，不能只选择表现最好的队列/单位/算法。','', '## 丰度单位和TCGA用途','', 'equal_RNA_cell_fraction是每细胞CP10k等权平均的理想化细胞混合比例；pooled_count_RNA_fraction真值是观测文库计数贡献。后者不是已校准的真实细胞RNA含量，两者都不能直接当作临床组织细胞比例。没有测序深度匹配/协变量，没有患者配对DE，没有环境RNA/双细胞校正。','', '本轮没有生成TCGA丰度，也没有用既有30/31基因生存评分替代丰度。EXPLORATORY_REFERENCE_NOT_VALIDATED.tsv仅供研究复核；拟合时必须应用reference_row_scales.csv并使用同时拟合组分，不能取一个目标列单独打分。即使算法拟合门槛通过，仍需bulk平台和RNA/细胞丰度校准，以及独立PGAM5蛋白联合巨噬细胞标记的生物学确认。当前结果支持RNA检出注释，尚不支持PGAM5特异细胞亚群或TCGA细胞丰度。','', '## 复现与实现核验','', f'独立核验{len(cells)}个细胞的原始PGAM5计数/注释，{len(solverchecks)}个保存混合输入的有界求解器复算，20个精确参考恢复，以及{len(summarychecks)}行汇总复算均通过。这证明抽查范围内实现一致，不证明生物学注释准确。','', '安装requirements.txt（Python版本见environment_versions.json），调整源码默认/workspace/scratch路径。先运行analyze.py再audit_report.py。legacy_reference_core.py来自前一轮联合参考工作流，默认依赖/workspace/scratch/PGAM5_myeloid_reference_v2的validation_protocol.json和三个harmonized_counts_QC.h5ad；可通过PGAM5_REFERENCE_SOURCE_ROOT改路径。新增输入依赖242/202全细胞QC计数和簇标签，以及前一轮状态注释metadata。包内包含结果、代码和审计用小型混合输入；不包含大型h5ad原始矩阵。此前DE、生存结果均未修改。']
(R/'README.md').write_text('\n'.join(lines)+'\n');print(json.dumps(audit,indent=2),flush=True)
