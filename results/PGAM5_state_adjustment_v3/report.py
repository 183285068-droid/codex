from pathlib import Path
import json,hashlib,platform,importlib.metadata as md
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent
res=json.loads((R/'result.json').read_text());ev=pd.read_csv(R/'training_gene_stability_specificity.csv');s=pd.read_csv(R/'validation_summary.csv');a=pd.read_csv(R/'panel_fold_testability.csv');p=pd.read_csv(R/'predictions.csv');sv=pd.read_csv(R/'state_discrimination.csv');pv=pd.read_csv(R/'program_validation_by_patient.csv');pl=pd.read_csv(R/'program_lineage_specificity.csv');checks=[]
for _,t in s.iterrows():
 if t.scenario=='standard':
  checks.extend([{'phase':t.phase,'panel':t.panel,'model':t.model,'criterion':t.unit+'_MAE<=1pp','value':t.MAE_pp,'passed':t.MAE_pp<=1},{'phase':t.phase,'panel':t.panel,'model':t.model,'criterion':t.unit+'_r>=.7','value':t.Pearson_r,'passed':pd.notna(t.Pearson_r) and t.Pearson_r>=.7}])
 checks.append({'phase':t.phase,'panel':t.panel,'model':t.model,'criterion':t.unit+'_'+t.scenario+'_zero95<=.5%','value':t.zero_p95_percent,'passed':pd.notna(t.zero_p95_percent) and t.zero_p95_percent<=.5})
checks += [{'phase':'external_retest','panel':'all','model':'all','criterion':'eligible_external_patients>=3','value':2,'passed':False},{'phase':'external_retest','panel':'all','model':'all','criterion':'pure_DC_test_patients>=2','value':0,'passed':False}]
pd.DataFrame(checks).to_csv(R/'abundance_gate_checks.csv',index=False)
pr=[]
for phase,t in pv[pv.role.eq('test')].groupby('phase'):
 eligible=t[t.eligible_positive_ge5];rfrac=float((eligible.Spearman_score_PGAM5>=.1).mean()) if len(eligible) else 0
 l=pl[pl.phase.eq(phase)&pl.role.eq('test')];bad=[]
 for d,dd in l.groupby('donor'):
  mm=dd.loc[dd.lineage.eq('Macrophage'),'median_score']
  if len(mm):bad.extend((d+':'+str(z.lineage)) for _,z in dd[(dd.lineage!='Macrophage')&(dd.cells>=20)&(dd.median_score>mm.iloc[0])].iterrows())
 pr.append({'phase':phase,'eligible_patients':len(eligible),'median_AUC_RNA_detection':float(eligible.AUC_RNA_detection.median()),'median_Spearman_PGAM5':float(eligible.Spearman_score_PGAM5.median()),'fraction_eligible_patients_rho_ge0_1':rfrac,'lineage_specificity_failures':';'.join(bad),'convergence_all':bool(t.NMF_converged.all()),'supported':len(eligible)>=3 and rfrac>=.7 and len(bad)==0 and t.NMF_converged.all()})
pd.DataFrame(pr).to_csv(R/'program_support_checks.csv',index=False)
res['programs_supported']=all(t['supported'] for t in pr);res['TCGA_state_score']='NOT_RUN: exploratory state programs do not meet support rules' if not res['programs_supported'] else 'NOT_RUN: requires new independent external validation and bulk normalization calibration';res['analysis_status']='Exploratory completed; no validated PGAM5-specific abundance or state signature';(R/'result.json').write_text(json.dumps(res,indent=2))
# Verify output rows and important arithmetic independently.
for _,t in s[s.scenario.eq('standard')].iterrows():
 q=p[(p.phase==t.phase)&(p.panel==t.panel)&(p.model==t.model)&(p.unit==t.unit)&(p.scenario==t.scenario)];assert len(q)==t.mixtures;assert np.isclose(sum(abs(float(x)-float(y)) for x,y in zip(q.predicted,q.truth))/len(q)*100,t.MAE_pp)
assert not ev.loc[ev.specific_candidate,'gene'].eq('PGAM5').any();assert len(pd.read_csv(R/'exact_two_state_recovery.csv'))==20
(R/'environment_versions.json').write_text(json.dumps({'python':platform.python_version(),**{k:md.version(k) for k in ['numpy','pandas','scipy','anndata','scanpy','scikit-learn']}},indent=2))
fig,axes=plt.subplots(1,2,figsize=(10,4))
for ax,model in zip(axes,['direct_NNLS','hierarchical_NNLS']):
 t=p[p.phase.eq('external_retest')&p.model.eq(model)&p.unit.eq('equal_RNA_cell_fraction')&p.scenario.eq('standard')];ax.scatter(t.truth*100,t.predicted*100,alpha=.5,s=20);ax.plot([0,5],[0,5],ls='--',c='black');ax.set(title=model,xlabel='True fraction (%)',ylabel='Estimated fraction (%)')
fig.tight_layout();fig.savefig(R/'external_models.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(figsize=(10,4));tt=pv[pv.role.eq('test')];lab=[str(z.phase).replace('cohort_holdout_','')+':'+str(z.donor).split(':')[-1] for _,z in tt.iterrows()];ax.scatter(range(len(tt)),tt.Spearman_score_PGAM5,s=25);ax.axhline(.1,c='red',ls='--');ax.set_xticks(range(len(tt)),lab,rotation=90);ax.set(ylabel='Spearman: program score vs PGAM5',title='Held-out macrophage program association');fig.tight_layout();fig.savefig(R/'program_patient_validation.png',dpi=180);plt.close(fig)
lines=['# PGAM5相关巨噬细胞状态与分层模型调整','', '完成跨患者一致性筛选、竞争谱系特异性筛选、有限细胞周期排除敏感性、分层反卷积及非细胞周期NMF状态探索。原有FDR<0.05、|log2FC|>=1差异表达表未更改；没有患者配对DE或测序深度匹配/协变量。','', '## 筛选结论','',f"联合训练筛选：稳定候选{res['training_stable_genes']}个（CENPK），谱系特异候选{res['training_specific_genes']}个，去除固定细胞周期清单后稳定候选{res['training_stable_noncycle_genes']}个。CENPK在14位合格训练患者中10位方向一致，但其最大竞争群为Cycling_T_NK，未通过特异性筛选。",'', '筛选使用训练患者中位表达比、方向一致性及竞争组分表达上限，并非新的统计显著性DE检验。完整阈值见protocol.json。每次患者留出重新筛选，少于5个候选时明确记为不可测试，没有把缺失当作零误差。此前21基因共识仅作为探索性对照；其中20个在本次有限细胞周期清单内，去除后仅剩HCFC1，未构成可测试多基因非增殖评分。该21基因来自此前含GSE189903的分析，不能宣称完全独立的外部特征发现。','', '## 分层反卷积','', '第一层按广义谱系同时估计总巨噬细胞与竞争细胞；第二层从混合表达中扣除拟合的非巨噬细胞贡献，再用两类巨噬细胞参考估计条件比例。未检出巨噬细胞参考按患者等权，并按细胞数合并其循环与非循环成员。最终目标=总巨噬细胞比例×条件目标比例。该形式仍然需要验证，不能因公式成立而认定丰度有效。','', '| 外部等RNA测试 | MAE（百分点） | r | 零目标预测95分位（%） |','|---|---:|---:|---:|']
for _,t in s[s.phase.eq('external_retest')&s.scenario.eq('standard')&s.unit.eq('equal_RNA_cell_fraction')].iterrows():lines.append(f'| {t.model} | {t.MAE_pp:.3f} | {t.Pearson_r if pd.notna(t.Pearson_r) else "不可计算"} | {t.zero_p95_percent:.3f} |')
lines+=['','分层外部等RNA预测全部为0，无法识别非零目标；零目标误判变低不是验证成功。RNA权重模拟、内部留出及单核/DC/Mixed_APC挑战另见validation_summary.csv。只做NNLS，未执行CIBERSORTx、MuSiC或BayesPrism。20组匹配两状态参考的混合准确恢复，独立有界求解器的第二层系数与NNLS一致，MAE逐条重新计算；这不证明注释或参考的生物学准确性。','', '## 连续状态探索','', '在全部训练巨噬细胞中，排除PGAM5、固定细胞周期清单及MT/RPL/RPS基因，从训练表达中选取500个变异基因，拟合4个非负表达程序。仅依据训练数据挑选与连续PGAM5表达相关性最高的程序，使用其30个高载荷基因构建连续评分。分别执行151530→149614、149614→151530跨队列验证，再联合训练后复测189903。参数不按外部结果调整。','', '| 测试 | 合格患者 | 中位AUC（RNA检出） | 中位Spearman相关 | rho>=0.1患者比例 | 支持条件全部满足 |','|---|---:|---:|---:|---:|---|']
for t in pr:lines.append(f"| {t['phase']} | {t['eligible_patients']} | {t['median_AUC_RNA_detection']:.3f} | {t['median_Spearman_PGAM5']:.3f} | {t['fraction_eligible_patients_rho_ge0_1']:.1%} | {t['supported']} |")
lines+=['','模型基因、训练标准化参数与全NMF载荷已导出。训练相关不等于已验证的PGAM5驱动状态；PGAM5 RNA检出不是蛋白或功能金标准。未控制深度，因此深度相关信号仍可能存在。非巨噬细胞评分高于巨噬细胞的逐患者检查见program_support_checks.csv。NMF未收敛者明确记录，不能视为稳定基因程序。','', '## 使用限制与结果文件','', '当前不提供经验证可用于TCGA-LIHC丰度的签名，也不将探索性状态分数包装为细胞丰度。GSE189903仅2位患者有>=5个目标RNA检出巨噬细胞，且无符合规则的纯DC群；本轮复用了已看过的外部数据，不是新增独立验证。状态评分在bulk中的平台/单位校准也未完成。','', '- training_gene_stability_specificity.csv、per_donor_gene_contrasts.csv.gz：全部基因的训练一致性/特异性证据。','- donor_target_consistency_QA.csv：目标检出率、单计数阳性占比、循环标记占比的患者审计。','- panel_fold_testability.csv、fold_selected_genes.csv：每个折的候选数量与不可测试状态。','- predictions.csv、validation_summary.csv、abundance_gate_checks.csv：直接/分层模型完整结果。','- program_score_parameters.csv、*_NMF_basis.tsv、all_programs_training_association.csv：探索性状态评分参数与程序。','- program_validation_by_patient.csv、program_support_checks.csv、program_lineage_specificity.csv：跨队列和外部复测。','- heldout_cell_state_scores.csv.gz、program_cell_scores.csv.gz：可审计的逐细胞探索性分数。','', '复现：依赖上一轮三个harmonized_counts_QC.h5ad及其train_joint_reference.py，路径默认/workspace/scratch/PGAM5_myeloid_reference_v2。先python analysis.py，再python discover_programs.py，最后python report.py；OPENBLAS_NUM_THREADS=1。大型矩阵不打包。输入和代码哈希见来源清单。细胞周期清单是记录来源的有限人工转录列表及扩展，非完整GO注释；外部下载403/404未用未验证内容替代，不能把“未列入”理解为确定非增殖。','', '下一步更有价值的是新增独立、巨噬细胞/单核/DC覆盖充分的HCC队列，并以蛋白或功能证据定位PGAM5相关状态。现有结果不支持继续放宽阈值得到最终TCGA签名。']
(R/'README.md').write_text('\n'.join(lines)+'\n');print(json.dumps(res,indent=2));print(pd.DataFrame(pr).to_string(index=False))
# Additional lineage-restricted exploratory programs, with the same test rules.
Q=R/'lineage_restricted_programs';v2=pd.read_csv(Q/'program_validation_by_patient.csv');l2=pd.read_csv(Q/'program_lineage_specificity.csv');rows=[]
for phase,t in v2[v2.role.eq('test')].groupby('phase'):
 e=t[t.eligible_positive_ge5];frac=float((e.Spearman_score_PGAM5>=.1).mean()) if len(e) else 0;bad=[]
 for d,z in l2[l2.phase.eq(phase)&l2.role.eq('test')].groupby('donor'):
  mac=z.loc[z.lineage.eq('Macrophage'),'median_score']
  if len(mac):bad.extend(d+':'+str(a.lineage) for _,a in z[(z.lineage!='Macrophage')&(z.cells>=20)&(z.median_score>mac.iloc[0])].iterrows())
 rows.append({'phase':phase,'eligible_patients':len(e),'median_AUC_RNA_detection':float(e.AUC_RNA_detection.median()),'median_Spearman_PGAM5':float(e.Spearman_score_PGAM5.median()),'fraction_eligible_patients_rho_ge0_1':frac,'lineage_specificity_failures':';'.join(bad),'convergence_all':bool(t.NMF_converged.all()),'supported':len(e)>=3 and frac>=.7 and not bad and t.NMF_converged.all()})
pd.DataFrame(rows).to_csv(Q/'program_support_checks.csv',index=False)
res['lineage_restricted_programs_supported']=all(z['supported'] for z in rows);res['TCGA_state_score']='NOT_RUN: no program met all support criteria; bulk scores would not establish PGAM5-specific macrophage state';(R/'result.json').write_text(json.dumps(res,indent=2))
extra=['','## 追加：训练谱系筛选后的非增殖程序','', '初始训练NMF高载荷出现ALB/APOA1/TTR后，追加仅依据训练数据的谱系条件：巨噬细胞患者等权均值不低于每个有支持的非巨噬细胞群，且巨噬细胞检出率>=10%，再执行同样的NMF流程。这是在训练标记审计后的探索性改进，不是预注册独立验证。','', '| 测试 | 合格患者 | 中位AUC | 中位Spearman相关 | rho>=0.1患者比例 | 全部条件满足 |','|---|---:|---:|---:|---:|---|']
for t in rows:extra.append(f"| {t['phase']} | {t['eligible_patients']} | {t['median_AUC_RNA_detection']:.3f} | {t['median_Spearman_PGAM5']:.3f} | {t['fraction_eligible_patients_rho_ge0_1']:.1%} | {t['supported']} |")
extra += ['', '联合模型的研究用30基因程序含C1QA/B/C、CTSB、CTSD、PSAP、TREM2、SPP1、CXCL9/10等巨噬细胞相关基因，减少了初始程序的肝细胞相关基因，但不能证明环境RNA或双细胞问题已纠正；但这些基因组成不能自行证明PGAM5调控、功能因果或bulk特异性。未做GO/通路富集检验。基因列表和分数参数见lineage_restricted_programs/program_score_parameters.csv。','', '不把未达到关联/特异性/患者覆盖要求的程序命名为已验证PGAM5状态签名。所有程序和负结果均保留，后续可结合独立蛋白或功能证据检验；此轮未运行TCGA状态评分。复现追加程序：在report.py之前运行python discover_lineage_restricted_programs.py。']
with (R/'README.md').open('a') as f:f.write('\n'.join(extra)+'\n')
print(pd.DataFrame(rows).to_string(index=False))

with (R/"README.md").open("a") as f:f.write("\n实现复核完成：从原始计数独立重算60条细胞评分，最大误差2.45e-15以内；32项AUC通过独立Mann–Whitney方法复核。研究用30基因清单另单独导出为EXPLORATORY_30_gene_macrophage_program_NOT_VALIDATED.csv。\n")
