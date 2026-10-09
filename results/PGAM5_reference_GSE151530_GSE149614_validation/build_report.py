from pathlib import Path
import os,json,hashlib,shutil,zipfile,importlib.metadata as im
import numpy as np,pandas as pd,matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;S=Path(os.environ.get('PGAM5_SOURCE_ROOT','/workspace/scratch'));out=json.loads((R/'validation_result.json').read_text());lock=json.loads((R/'locked_reference_manifest.json').read_text());su=pd.read_csv(R/'validation_summary.csv');pred=pd.read_csv(R/'validation_predictions.csv');checks=pd.read_csv(R/'predeclared_gate_checks.csv');audit=json.loads((R/'implementation_audit_summary.json').read_text())
assert not out['validation_passed'] and not checks.passed.all();assert audit['max_exact_coefficient_error']<1e-6 and audit['external_predictions_exactly_reproduced']==560
for phase in ['outer_nested','external_locked']:
 t=pred[(pred.phase==phase)&pred.unit.eq('equal_RNA_cell_fraction')&pred.scenario.eq('standard')];line=su[(su.phase==phase)&su.unit.eq('equal_RNA_cell_fraction')&su.scenario.eq('standard')].iloc[0]
 assert np.isclose(abs(t.predicted-t.truth).mean()*100,line.MAE_percentage_points)
 assert np.isfinite(t.predicted).all() and t.predicted.between(0,1).all()
assert len(lock['genes'])==len(set(lock['genes']))==1100 and len(lock['reference_groups'])==11
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight'})
fig,axs=plt.subplots(1,2,figsize=(11,4),sharex=True,sharey=True)
for ax,phase,title in zip(axs,['outer_nested','external_locked'],['GSE151530: nested held-out donors','GSE149614: locked external reference']):
 t=pred[(pred.phase==phase)&pred.unit.eq('equal_RNA_cell_fraction')&pred.scenario.eq('standard')]
 for donor,g in t.groupby('donor',observed=True):ax.scatter(g.truth*100,g.predicted*100,s=16,alpha=.7,label=donor)
 ax.plot([0,6],[0,6],ls='--',color='#777');ax.set_xlabel('True equal-RNA cell fraction (%)');ax.set_title(title);ax.set_xlim(-.2,5.5);ax.set_ylim(-.2,max(7,pred.loc[pred.phase.eq('outer_nested')&pred.unit.eq('equal_RNA_cell_fraction'),'predicted'].max()*100*1.05));ax.legend(fontsize=7,ncol=2)
axs[0].set_ylabel('NNLS estimated target fraction (%)');fig.suptitle('PGAM5 RNA-detected TAM: biological holdout performance');fig.tight_layout();fig.savefig(R/'heldout_target_recovery.png',dpi=150);fig.savefig(R/'heldout_target_recovery.pdf');plt.close(fig)
fig,ax=plt.subplots(figsize=(9,4));labels=[];values=[];colors=[]
for phase,color,short in [('outer_nested','#355C8A','Internal'),('external_locked','#B77B28','External')]:
 for scenario in ['standard','cycling_zero','other_myeloid_zero']:
  z=su[(su.phase==phase)&su.unit.eq('equal_RNA_cell_fraction')&su.scenario.eq(scenario)]
  if len(z):labels.append(short+' / '+scenario.replace('_',' '));values.append(z.iloc[0].zero_target_p95_percent);colors.append(color)
ax.barh(labels,values,color=colors);ax.axvline(.5,ls='--',color='#A23D37',label='Predeclared maximum: 0.5%');ax.set_xlabel('95th percentile estimated target (%) when true fraction is zero');ax.invert_yaxis();ax.legend(fontsize=8);ax.set_title('Zero-target false-positive checks');fig.tight_layout();fig.savefig(R/'zero_target_gate.png',dpi=150);fig.savefig(R/'zero_target_gate.pdf');plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(11,4))
cases=[('external_locked',True,True,'Locked'),('external_sensitivity',False,True,'Without PGAM5'),('external_sensitivity',True,False,'Without cycle splits')]
for ax,metric,title in zip(axs,['MAE_percentage_points','Pearson_r'],['External error (percentage points)','External target recovery (Pearson r)']):
 labels=[];v=[]
 for phase,pg,cyc,label in cases:
  z=su[(su.phase==phase)&su.with_PGAM5.eq(pg)&su.with_cycling_competitors.eq(cyc)&su.unit.eq('equal_RNA_cell_fraction')&su.scenario.eq('standard')].iloc[0];labels.append(label);v.append(z[metric])
 ax.barh(labels,v,color='#355C8A');ax.invert_yaxis();ax.set_title(title);ax.axvline(1 if metric=='MAE_percentage_points' else .7,ls='--',color='#B77B28',label='Predeclared threshold');ax.legend(fontsize=8)
fig.tight_layout();fig.savefig(R/'external_sensitivity.png',dpi=150);fig.savefig(R/'external_sensitivity.pdf');plt.close(fig)
prep={n:json.loads((R/(n+'_preparation_summary.json')).read_text()) for n in ['GSE151530','GSE149614']}
rows=[]
for name in prep:
 o=pd.read_csv(R/(name+'_QC_cell_metadata.csv.gz'));z=o.groupby('donor',observed=True).agg(QC_cells=('cell_id','size'),TAMs=('coarse_group',lambda x:int(x.eq('TAM').sum())),positive_TAMs=('fine_group',lambda x:int(x.eq('TAM_PGAM5_detected').sum()))).reset_index();z['dataset']=name;rows.append(z)
pd.concat(rows,ignore_index=True).to_csv(R/'donor_reference_cell_counts.csv',index=False)
# Raw input provenance; full matrices retained in scratch, excluded from git/ZIP.
paths=[S/'GSE151530/GSE151530_matrix.mtx.gz',S/'GSE151530/GSE151530_genes.tsv.gz',S/'GSE151530/all_cell_metadata.csv.gz',S/'GSE149614/GSE149614_HCC.scRNAseq.S71915.count.txt.gz',S/'GSE149614/GSE149614_HCC.metadata.updated.txt.gz',S/'GSE149614/myeloid_cluster_annotation.csv']
hashes={}
for p in paths:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 hashes[str(p.relative_to(S))]={'bytes':p.stat().st_size,'sha256':h.hexdigest()}
(R/'source_input_sha256.json').write_text(json.dumps(hashes,indent=2));(R/'environment_versions.txt').write_text('\n'.join(f'{p}=={im.version(p)}' for p in ['numpy','pandas','scipy','anndata','matplotlib','threadpoolctl'])+'\n')
report='''# PGAM5 RNA检出巨噬细胞参考：训练、内部留患者与GSE149614外部验证

**已完成参考构建及预先指定的验证。本次参考+NNLS方案未通过，不能作为可用TCGA-LIHC丰度signature。未对TCGA-LIHC生成丰度结果。** 这是数值验证门槛导致的结果，不是缺少用户授权，也不是“没有PGAM5⁺巨噬细胞”。

## 执行内容

- 发现参考：GSE151530全部HCC肿瘤QC细胞44,776个，25个patient_proxy，作者TAM4,494个，其中PGAM5原始RNA检出180个、未检出4,314个。
- 外部数据：GSE149614十个T结尾原发肿瘤，QC细胞34,414个，沿用此前marker推定的巨噬细胞6,151个，其中检出167个、未检出5,984个。与此前原始细胞ID集合和分组数核对一致。
- 两套数据按基因符号合并重复源行后，分别18,661/25,712个gene，交集18,624个。GSE151530有6个重复symbol源行按原始计数相加；均保留源映射。QC为≥500检出gene、MT计数比例<20%，标准化至CP10K用于线性参考。
- 作者粗标签保留全部可用类型，增加探索性增殖竞争群：MKI67/TOP2A/UBE2C/CENPF/CDK1/BIRC5至少2个raw count>0。此规则不是独立验证的细胞周期分类器；目标TAM定义不变。训练折内少于100细胞的增殖群回并父类。
- 原来的FDR/log2FC差异结果没有更改；本次特征选择用于多细胞类型参考，不是重新计算此前DE，也未新增患者配对差异分析或深度匹配。

## 训练与锁定，避免外部结果调参

先写定validation_protocol.json及通过标准，再测试三种发现集配置：每群50特征、细胞加权参考；每群50特征、患者等权参考；每群100特征、患者等权参考。所有表达均值、marker对比、检出率与特征排名仅由当前训练患者生成。

对每个训练参考群计算相对其余群最大均值的log2对比，以表达强度调节排名，训练检出率≥10%、平均CP10K≥0.05；每群取50或100个描述性特征，另保留定义基因PGAM5。只用于反卷积特征选择时排除MT-/RPL/RPS，既有差异表仍保留这些gene。这些特征并非1100个“PGAM5⁺巨噬细胞特异差异基因”。

在GSE151530进行嵌套留患者：外层6位合格患者，内层用其余合格患者的固定均值混合物选择配置，外层使用独立抽样混合物评价。外层患者表达不参与该折特征选择；来源标签和固定增殖规则作为先验定义保留。最终配置只按GSE151530发现集留患者结果选定，之后锁定**患者等权、每群100特征，共1,100个gene×11个参考成分**，再评估GSE149614。没有按外部结果调参。

外层合格患者：H08、H68、H70、H72、H73、H74。外部合格患者：HCC01、HCC02、HCC03、HCC04、HCC05、HCC06、HCC08、HCC10。预先要求目标≥5、未检出TAM≥20、肿瘤/肝细胞样≥20、T/NK≥20；其他患者数据仍可参与发现参考或作为外部整体数据，未满足模拟资格的患者不进入主验证。这种资格限制与稀少目标细胞需在解释时保留。

参考列：B、CAF、Cycling_TAM_PGAM5_undetected、Cycling_T_NK、Cycling_Tumor_hepatocyte、Endothelial、TAM_PGAM5_detected、TAM_PGAM5_undetected、T_NK、Tumor_hepatocyte、Unknown_other。T_NK是为跨来源对齐而合并的粗标签，GSE151530源标签为T cells，未独立建立NK注释；Tumor_hepatocyte在GSE149614源标签为Hepatocyte，未逐细胞证明均为恶性细胞。

## 验证方式、算法与单位

采用本地非负最小二乘NNLS，同时拟合所有参考成分。逐gene按训练参考均方根缩放，拟合非负系数后归一为总和1。没有调用CIBERSORTx、MuSiC或BayesPrism；当前结果不能推广成“任何算法均不可能”。

在每个验证患者内部重采样2,000细胞，目标比例0、0.5%、1%、2%、5%，每个比例5次。总TAM份额保持25%，肿瘤/肝细胞样50%，其余类型25%随机分配，检查不同竞争构成。标准情况下增殖/非增殖肿瘤比例也随机变化；另有零目标、增殖竞争富集挑战，外部再加入未覆盖Other_myeloid富集挑战。不同配置使用固定随机数种子复现同一混合物。

两种单位并行测试：

1. **等RNA模拟细胞比例**：先将各单细胞归一至CP10K后平均，隐含每细胞RNA量相同；已知细胞数份额仅在此理想化设置下对应线性系数。
2. **计数池化RNA权重份额**：先求原始计数和再按全部gene的总计数归一，真值由被抽中细胞的观测文库计数构成。此处观测UMI文库大小不是经过校准的生物学RNA产量，更不是可直接换算的真实bulk细胞百分比。

抽样重复不是独立患者重复。本次发现数据集和外部数据集之前参与过早期探索；这里新增特征的选择和配置锁定严格只用发现集，但仍不等于前瞻性验证。

## 预先指定的通过标准与实际结果

标准混合物：MAE≤1个百分点、Pearson r≥0.7；零目标的估计值95分位数≤0.5%；增殖竞争与未覆盖髓系零目标挑战同样≤0.5%。两种单位分别检查，外部需≥5个患者。这些是本项目事先设定的可用性门槛，不是通用临床标准。

|阶段（等RNA标准模拟）|患者|MAE（百分点）|Pearson r|零目标估计95分位数（%）|
|---|---:|---:|---:|---:|
'''
for phase,title in [('outer_nested','GSE151530嵌套留患者'),('external_locked','GSE149614外部锁定参考')]:
 z=su[(su.phase==phase)&su.unit.eq('equal_RNA_cell_fraction')&su.scenario.eq('standard')].iloc[0];report+=f'|{title}|{z.donors}|{z.MAE_percentage_points:.2f}|{z.Pearson_r:.3f}|{z.zero_target_p95_percent:.2f}|\n'
report+='''
外部增殖竞争零目标挑战95分位数约2.17%；未覆盖髓系挑战约1.17%，均超过0.5%。外部计数池化RNA标准模拟MAE约1.52个百分点、相关约-0.073，也未通过。predeclared_gate_checks.csv逐条列出通过/失败，validation_predictions.csv保留每个混合物。主内部模拟360条（180混合物×2单位）、外部主验证560条（280×2）、两个外部敏感性分支1120条，共2040条评估记录；另有内层配置选择记录。

外部敏感性：不使用PGAM5时，等RNA标准模拟r约0.006；不设增殖竞争群时r约0.016，且增殖竞争零目标95分位数升至约6.10%。敏感性分支仅作诊断，没有用于重新选择主方案；分支会按相同训练规则重新建立对应参考和特征，不是完全相同feature集合的单因素实验。

补充控制显示，外部总TAM比例标准模拟MAE约7.39个百分点；Other_myeloid富集挑战中总TAM上偏约19.97个百分点。说明不足还涉及跨来源标签/竞争群覆盖及参考迁移，不能将失败仅归因于某一个目标gene或定义。

## 实现复核与RNA混合风险审查

直接用锁定参考的11列构造20个已知系数混合物，最大恢复误差约3.3×10⁻¹⁶；另用独立有界最小二乘求解器复算80个外部案例，与NNLS预测及目标函数一致；560个外部主预测按固定抽样种子和锁定参考独立复现。参考列/feature顺序、非负性、候选唯一性、QC细胞ID/分组数均复核。这说明目前失败不是可见的索引错位或求解器错误；不证明输入生物学标签和所有建模假设正确。

对TAM检查CD3D/CD3E/TRAC的多marker共表达、MS4A1/CD79A、EPCAM/KRT19及ALB背景。PGAM5检出TAM中ALB检出率在GSE151530约85.6%、GSE149614约92.8%，T谱系至少2个marker共检出约6.7%/7.2%。这些是描述性跨谱系RNA旗标，不是经过验证的双细胞判断。不能将ALB直接视为污染，吞噬RNA等同样可能；本次没有据此删细胞或宣称去除了环境RNA。

仅有已过滤矩阵，没有空液滴背景、完整液滴/上游QC信息，不能可靠实施或宣称环境RNA校正。whole_reference_lineage_marker_QA.csv及TAM_cross_lineage_coexpression_flags.csv保留检查。PGAM5原始count=0仍是RNA未检出，不能直接解释为蛋白阴性或稳定亚型。GSE151530阳性TAM中H70贡献92/180，跨患者发现的稀疏性也是限制之一。

## 为什么此轮未进入TCGA-LIHC及后续具体条件

本轮完成了既定的参考训练和外部验证，并按预先设定门槛判定不合格。因此没有使用该参考生成TCGA-LIHC目标丰度，避免把不能恢复真值的数值当作细胞比例。即使数值门槛通过，仍需要UMI/bulk平台、转录本长度/单位匹配及明确分母；当前CP10K模板没有完成这一校准。

继续推进需要先补足并对齐参考细胞类型，特别是单核/DC及增殖竞争群，取得可支持环境RNA/双细胞审查的更完整输入，在新的未用于调参队列重新验证。若目标仍保持“PGAM5 RNA>0”，应优先取得足够患者与目标细胞，并验证检测状态是否跨平台稳定。若改成PGAM5-high连续状态、或蛋白阳性亚型，则属于目标定义更改，需要用户确认及适合该定义的新证据，不能在此轮自行替换。

CIBERSORTx自定义参考/S-mode或MuSiC/BayesPrism可以作为后续独立方法检验，但需要相应工具/服务和参考单位准备；本轮没有调用这些服务。算法更换能否改善当前迁移失败，需要新验证，不能预先保证。

## 文件与复现

LOCKED_REFERENCE_CP10K_VALIDATION_REQUIRED.tsv：1100 gene×11参考列，**仅供研究/复核，验证失败，不是可用TCGA-LIHC signature**。locked_reference_manifest.json记录配置、feature、患者和protocol哈希；training_selected_gene_evidence.csv记录训练选择依据。另含内/外层配置和预测、敏感性、逐条门槛、独立实现检查、全细胞类型marker和RNA混合旗标、图及代码。

大的全细胞QC count h5ad保留在工作区，不进入GitHub或ZIP；两套源数据SHA256、feature源映射和小体积原始metadata/gene列表随包。metadata中的样本/患者均为公开研究匿名编号。

按environment_versions.txt安装依赖，PGAM5_SOURCE_ROOT设为包含GSE151530/GSE149614子目录的源根目录（默认/workspace/scratch）；将原始大计数矩阵放入对应目录，并将随包source_metadata目录的小输入复制到对应源目录。依次运行：

```bash
python prepare_references.py GSE151530
python prepare_references.py GSE149614
python train_validate_reference.py
python audit_implementation_and_identity.py
python build_report.py
```

make_report/build_report输出位置默认项目仓库路径；可以设置PGAM5_RESULTS_ROOT指定其他输出目录。CIBERSORTx等未执行。三张PNG均经过视觉检查；完整网页交互呈现未独立验证。

来源：[GSE151530](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151530)、[GSE149614](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614)。GSE149614 macrophage-cluster定义来自此前同项目目标基因无关marker审查，仍属于推定标签，不等同于作者逐细胞TAM注释。
'''
(R/'结果说明.md').write_text(report)
# Small source companions permit reproducible setup without bundling large raw/QC matrices.
meta=R/'source_metadata';meta.mkdir(exist_ok=True)
for name,files in {'GSE151530':['GSE151530_genes.tsv.gz','all_cell_metadata.csv.gz'],'GSE149614':['GSE149614_HCC.metadata.updated.txt.gz','myeloid_cluster_annotation.csv']}.items():
 d=meta/name;d.mkdir(exist_ok=True)
 for f in files:shutil.copy2(S/name/f,d/f)
import os
D=Path(os.environ.get('PGAM5_RESULTS_ROOT','/workspace/codex/results/PGAM5_reference_GSE151530_GSE149614_validation'));D.mkdir(parents=True,exist_ok=True)
for p in R.iterdir():
 if p.is_file() and (p.suffix in ['.py','.csv','.json','.tsv','.md','.png','.pdf'] or p.name.endswith(('metadata.csv.gz','versions.txt'))):shutil.copy2(p,D/p.name)
shutil.copytree(meta,D/'source_metadata',dirs_exist_ok=True)
zp=D/'PGAM5_reference_validation_results.zip'
with zipfile.ZipFile(zp,'w',compression=zipfile.ZIP_DEFLATED) as z:
 for p in sorted(D.rglob('*')):
  if p.is_file() and p!=zp:z.write(p,p.relative_to(D))
with zipfile.ZipFile(zp) as z:assert z.testzip() is None
print('Validated report and artifacts; ZIP MiB',zp.stat().st_size/1024**2,'files',sum(p.is_file() for p in D.rglob('*')))
