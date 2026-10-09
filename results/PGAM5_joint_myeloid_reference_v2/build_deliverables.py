from pathlib import Path
import json,hashlib,platform
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent
names=['GSE151530','GSE149614','GSE189903'];counts=[];coverage=[];qa=[]
for n in names:
 o=pd.read_csv(R/(n+'_harmonized_cell_metadata.csv.gz'))
 for (donor,g),t in o.groupby(['donor_id','harmonized_group']):coverage.append({'dataset':n,'donor':donor,'group':g,'cells':len(t),'PGAM5_detected_macrophages':int(t.fine_group.eq('TAM_PGAM5_detected').sum())})
 s=json.loads((R/(n+'_annotation_summary.json')).read_text());counts.append({'dataset':n,'QC_cells':len(o),'target_detected':s['target_positive_cells'],'target_undetected':s['target_undetected_cells'],**s['groups']})
 q=pd.read_csv(R/(n+'_cluster_annotation_QA.csv'));q['dataset']=n;qa.append(q)
pd.DataFrame(counts).fillna(0).to_csv(R/'harmonized_cohort_counts.csv',index=False);pd.DataFrame(coverage).to_csv(R/'donor_lineage_coverage.csv',index=False)
markers=['C1QA','C1QC','CSF1R','TYROBP','FCN1','S100A9','CD1C','FCER1A','TCF4','IL3RA','GZMB']
heat=[];labels=[]
for q in qa:
 for g,t in q.groupby('assigned_group'):
  if g in ['Macrophage','Monocyte_enriched','DC_enriched','pDC','Mixed_APC']:
   heat.append([np.average(t['fraction_'+m],weights=t.cells) for m in markers]);labels.append(q.dataset.iloc[0]+': '+g)
fig,ax=plt.subplots(figsize=(10,7));im=ax.imshow(heat,vmin=0,vmax=1,cmap='viridis',aspect='auto');ax.set_xticks(range(len(markers)),markers,rotation=45,ha='right');ax.set_yticks(range(len(labels)),labels);fig.colorbar(im,ax=ax,label='RNA detection fraction');fig.tight_layout();fig.savefig(R/'lineage_marker_audit.png',dpi=180);plt.close(fig)
v=pd.read_csv(R/'validation_predictions.csv');su=pd.read_csv(R/'validation_summary.csv');result=json.loads((R/'validation_result.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(10,4))
for ax,phase in zip(axes,['outer_nested','external_locked']):
 t=v[v.phase.eq(phase)&v.unit.eq('equal_RNA_cell_fraction')&v.scenario.eq('standard')];ax.scatter(100*t.truth,100*t.predicted,s=12,alpha=.5);ax.plot([0,5],[0,5],color='black',ls='--');ax.set(xlabel='True target fraction (%)',ylabel='Estimated fraction (%)',title=phase)
fig.tight_layout();fig.savefig(R/'heldout_prediction.png',dpi=180);plt.close(fig)
lines=['# PGAM5 RNA检出巨噬细胞：补齐髓系竞争参考后的联合训练与验证','', '已完成三个数据集的统一注释、GSE151530与GSE149614联合训练、嵌套患者留出验证，以及锁定参考后GSE189903外部性能验证。原有差异表达结果未改动。','', '## 数据与注释','', '| 数据集 | QC细胞 | PGAM5检出巨噬细胞 | PGAM5未检出巨噬细胞 | 单核富集 | DC富集 | pDC | 混合APC |','|---|---:|---:|---:|---:|---:|---:|---:|']
for c in counts:lines.append('| '+' | '.join(str(c.get(k,0)) for k in ['dataset','QC_cells','target_detected','target_undetected','Monocyte_enriched','DC_enriched','pDC','Mixed_APC'])+' |')
lines+=['','采用独立数据集内聚类、PGAM5以外的标记基因进行髓系谱系判定。单核使用FCN1/S100A8/S100A9/TYROBP，DC富集使用CD1C/FCER1A，巨噬细胞使用C1QA/B/C、CSF1R、CD68及TYROBP。共表达DC和巨噬细胞标记的群保留为Mixed_APC，未强制归入一种谱系；未支持的亚型不推断为cDC1/cDC2。','', '注释改变了巨噬细胞成员范围，故细胞数与原先作者标签/旧推断分析不同。PGAM5原始count>0的定义不变；count=0仅为RNA未检出。TAM与不明细胞候选群重新判定，原有明确非候选细胞标签仍沿用。初次标记审计发现双谱系群，需要在任何联合模型训练前增加Mixed_APC规则；初版和修正版规则均保留。GSE189903的标记汇总参与注释合理性检查，且此前已被探索过，因此并非前瞻性盲法外部队列。外部表达未用于模型特征或配置选择。','', '## 参考与验证','',f"锁定配置：{result['selected_config']}；{result['selected_genes']}个基因、{result['reference_components']}个参考组分；内部{result['internal_donors']}位合格患者，外部{result['external_donors']}位。训练保留低支持群的Unknown_other回退，不用外部性能调整参考。",'', '| 验证 | 单位 | MAE（百分点） | Pearson r | 真值为0时预测95分位（%） |','|---|---|---:|---:|---:|']
for _,t in su[su.phase.isin(['outer_nested','external_locked'])&su.scenario.eq('standard')].iterrows():lines.append(f'| {t.phase} | {t.unit} | {t.MAE_percentage_points:.3f} | {t.Pearson_r:.3f} | {t.zero_target_p95_percent:.3f} |')
lines+=['','预设标准：MAE≤1个百分点、相关系数≥0.7、零目标预测95分位≤0.5%，同时要求外部单核/DC挑战各至少2位可测试患者。模拟每次抽样2000细胞，目标0–5%，总巨噬细胞25%；重复抽样不是独立生物学患者。CP10k细胞均值与原始计数池的RNA权重分别报告；后者使用观测文库大小，不等于校准的真实细胞RNA含量。','', '最终验证结论：'+('通过预设数值标准。' if result['validation_passed'] else '**未通过预设标准，不能据此报告TCGA-LIHC的PGAM5⁺巨噬细胞丰度。**'),'', 'GSE189903未出现符合规则的纯DC富集群，不能独立验证这类混淆风险；缺失不视为通过。单核、Mixed_APC及循环相关的零目标挑战详见validation_summary.csv与predeclared_gate_checks.csv。此次执行的是NNLS基线，不能将其结果推广为CIBERSORTx/MuSiC/BayesPrism必然失败。没有运行TCGA反卷积，也没有生成患者丰度数值。','', '## 文件与复现','', '- LOCKED_REFERENCE_CP10K_VALIDATION_REQUIRED.tsv：研究验证参考矩阵，行=基因、列=同时拟合组分；需按训练规则缩放，不能直接解释为临床验证签名。','- donor_lineage_coverage.csv、各数据集annotation_transitions与逐细胞压缩metadata：患者覆盖、新旧注释对应。','- validation_predictions/summary、predeclared_gate_checks：全部性能、敏感性分析和未通过项。','- exact_reference_recovery_checks、independent_bounded_solver_checks、implementation_audit_summary：实现复核。','- reannotate_myeloid.py、train_joint_reference.py、audit_implementation.py：核心脚本。','', '共同基因空间为15,923个符号，外部数据仅用于确认基因是否可测，特征排序不使用外部表达。非成对的患者留出是预测验证，未执行按患者配对的差异表达，也未增加测序深度协变量。原始矩阵及大型h5ad不随压缩包发布。151530/149614输入由此前参考构建流程生成；189903使用官方matrix/genes/metadata，限定HCC core+border。脚本采用本工作区原始输入路径，复现时需调整路径。Python环境版本见environment_versions.json。未做未经验证的环境RNA纠正或双细胞剔除。后续需独立含DC的HCC队列检验、跨平台与丰度单位校准；当前结果不足以提供可直接用于TCGA丰度的最终签名。']
audit=json.loads((R/'implementation_audit_summary.json').read_text()); lines += ['', f"实现复核：{audit['exact_reference_mixtures_checked']}组完全匹配参考的混合恢复通过，最大系数误差{audit['max_exact_coefficient_error']:.3g}；{audit['independent_bounded_solver_cases']}组独立有界求解器结果一致；{audit['external_predictions_exactly_reproduced']}条外部预测逐条精确复现。上述检查未发现数值求解或结果写入错误，但不能证明注释/参考具有生物学准确性。", '', '外部零目标竞争挑战（等RNA单位）的预测95分位：单核富集6.98%，Mixed_APC富集4.56%，循环富集2.86%，均超过0.5%标准。全部4位外部患者的补充纯单核均值测试中，有一位被误估为6.63%的目标细胞；即使补齐竞争组分，仍存在跨谱系混淆。', '', 'Mixed_APC在训练集中364细胞，其中361来自GSE149614:HCC06；pDC的86细胞中57来自GSE151530:H73。总细胞数并不代表充分的跨患者支持，折内稀少群会回退到Unknown_other。外部仅3H/4H达到目标阳性>=5的要求；1H/2H各4个阳性仍可用于补充纯竞争细胞测试，但不能替代预设的完整目标剂量验证。', '', 'pure_competitor_stress_tests.csv额外评估纯单核/DC/Mixed_APC等细胞均值；外部可用所有4位患者的纯竞争池，训练队列使用患者留出参考。这是补充压力测试，不取代既定混合肿瘤门槛。']
(R/'README.md').write_text('\n'.join(lines)+'\n')
import importlib.metadata as md
(R/'environment_versions.json').write_text(json.dumps({'python':platform.python_version(),**{p:md.version(p) for p in ['numpy','pandas','scipy','anndata','scanpy','igraph','matplotlib']}},indent=2))
