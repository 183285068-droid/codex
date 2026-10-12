"""Build a source-backed Markdown results narrative from reviewed tables."""
import pandas as pd,json
from analyze import ROOT,PROGRAMS,COMPS,PANELS
ZH={'MHCII_presentation':'MHC-II抗原呈递','Inflammatory_myeloid':'炎症髓系RNA','CCL_CXCL8_chemokines':'CCL/CXCL8趋化因子','CXCR3_ligands':'CXCR3配体','Checkpoint_ligands':'检查点配体','Regulatory_mediators':'IL10/TGFB1/LGALS9调节因子','Costimulation_ligands':'共刺激配体','Scavenger_LPL_state':'清道夫/LPL状态','C1Q_complement':'C1Q补体（源门重叠）','Uptake_clearance_proxy':'摄取/清除代理（源门重叠）','SPP1_MMP12_remodeling':'SPP1/MMP12重塑','Hypoxia_response':'缺氧应答','Cycling':'增殖程序'}
s=pd.read_csv(ROOT/'program_statistics.csv');p=pd.read_csv(ROOT/'patient_equal_section_effects.csv');m=pd.read_csv(ROOT/'patient_direction_summary.csv');n=pd.read_csv(ROOT/'section_summary.csv');audit=json.loads((ROOT/'audit_results.json').read_text())
primary=s[s.comparison.eq(COMPS[0])];pp=p[p.comparison.eq(COMPS[0])].pivot(index='program',columns='patient',values='median_standardized_beta').reindex(PROGRAMS);pos=p[p.comparison.eq(COMPS[1])].pivot(index='program',columns='patient',values='median_standardized_beta').reindex(index=PROGRAMS,columns=pp.columns)
table='| RNA程序 | A009 | B013 | C004 | C008 | 主比较同向患者 | 仅检出对照同向患者 |\n| --- | ---: | ---: | ---: | ---: | --- | --- |\n'
def direction(v):
 v=v.dropna();return f'↑{int((v>0).sum())} / ↓{int((v<0).sum())}（{len(v)}例）'
for k in PROGRAMS:table+='| '+ZH[k]+' | '+' | '.join(f'{pp.loc[k,x]:+.3f}' if pd.notna(pp.loc[k,x]) else 'NA' for x in ['A009','B013','C004','C008'])+' | '+direction(pp.loc[k])+' | '+direction(pos.loc[k])+' |\n'
significant=primary[primary.primary_q_208.lt(.05)];full=s[s.joint_q_832.lt(.05)];consistent=[ZH[k] for k in PROGRAMS if ((pp.loc[k]>0).all() and (pos.loc[k]>0).all()) or ((pp.loc[k]<0).all() and (pos.loc[k]<0).all())]
summary=f'主比较208项中{int(primary.status.eq("tested").sum())}项可检验，{len(significant)}项q<0.05，最小主q={primary.primary_q_208.min():.4f}；全部832项中{int(s.status.eq("tested").sum())}项可检验，{len(full)}项联合q<0.05。'
summary+='4位患者的13个程序检验均未达到患者层面BH q<0.05。' if not m.patient_BH_q_13.lt(.05).any() else ''
if consistent:summary+='两种对照下4例方向均一致的描述性程序为：'+ '、'.join(consistent)+'；方向一致不等于显著或已证实功能改变。'
else:summary+='没有程序在两种对照下均具备4例可检验的一致方向证据，结果具有病例/对照依赖性。'
counttable=n.to_markdown(index=False) if False else '| 患者 | 切片 | 巨噬相关位点 | PGAM5-high | 其余检出 | 未检出 |\n| --- | ---: | ---: | ---: | ---: | ---: |\n'
for patient,z in n.groupby('patient'):counttable+=f'| {patient} | {len(z)} | {int(z.mac_spots.sum())} | {int(z.high.sum())} | {int(z.other_detected.sum())} | {int(z.undetected.sum())} |\n'
txt=f'''# GSE281759：PGAM5-high与其他巨噬细胞相关区域的功能RNA差异

**本轮比较的是Visium混合位点的功能相关RNA程序，无法确认单个巨噬细胞的功能改变，也不能推断PGAM5导致这些改变。** 采用13套预先固定的程序、55个唯一标志基因；PGAM5不进入任何被比较的功能分数。

## 统计结果与可解释范围

{summary}

本轮使用4患者/16肿瘤切片（15原发、1复发）、25,354 QC位点、4,091个巨噬RNA富集位点。382个PGAM5-high与3,709个其他巨噬RNA位点作主比较；后者包含1,123个其他检出、2,586个未检出源。功能表达比较无需有空间邻居，因此保留7个空间孤立源（1个其他检出、6个未检出）；前轮邻接分析分别为1,122个其他检出及2,580个未检出，分母差异来自是否排除孤立源。

{counttable}

### 患者效应方向

表中数值为逐片“调整后高表达组系数 / 当前比较总体程序分数SD”的患者等权中位数。正数表示high更高，负数表示更低；不是log2 fold change、Cohen d、细胞丰度、功能活性单位或因果效应。每位患者的多个切片不当作独立病例。

{table}

**最一致的描述性候选是清道夫/LPL相关RNA偏低**：主比较4例均负，检出对照3例可评估者均负；C004检出对照因源数/残余自由度不足未能检验，不补成第4例。主比较的共刺激配体4例均负，但检出对照中有正有负，不能称稳健降低。MHC-II主比较3例偏高，检出对照3例偏低，尤其依赖对照定义。不能据此宣称抗原呈递、免疫抑制或吞噬能力整体改变。

清道夫/LPL组合的方向并不代表每个基因都下降：主比较的LPL单基因患者效应为A009 +0.074、B013 −0.005、C004 +0.530、C008 −0.098；STAB1为−0.129、−0.041、−0.597、+0.074。MARCO在C004缺乏可估计变异。组合下降应作为待验证状态信号，**不能写成LPL普遍下调、脂质代谢下降或吞噬功能下降**。B013组合效应仅−0.003 SD，也不支持四患者有相同效应幅度。

患者只有4例，精确双侧符号翻转共16种分配，最小可达p为0.125；患者层面检验仅探索性，依赖患者效应近似符号对称性。不能把位点数用作总体患者样本量。各程序的患者效应与完整统计见 [patient_direction_summary.csv](patient_direction_summary.csv)、[program_statistics.csv](program_statistics.csv)。

## 可视化

![患者功能差异](figures/patient_adjusted_function_effects.png)

- [逐切片调整后功能差异](figures/section_adjusted_function_effects.png)
- [两种对照的患者效应](figures/patient_comparison_sensitivity.png)
- [未调整与调整结果对照](figures/unadjusted_adjusted_comparison.png)
- [排除T-RNA源敏感性](figures/nonT_function_sensitivity.png)
- [MT≤5%敏感性](figures/MT5_function_sensitivity.png)
- [55个标志基因的描述性差异](figures/marker_gene_effects.png)
- [测序质量及混合RNA代理](figures/quality_and_mixing_proxies.png)

8套图均同时提供PNG/PDF。热图固定对称色标，超过色标的值只在颜色上饱和，标注保留真实数值；NA不可当作0。点图的点是不同患者的效应，不是置信区间。

## 程序如何定义

这些是根据现有marker构建的定向多基因RNA程序，不是完整GO/KEGG通路、GSVA/ssGSEA/GSEA，也不是M1/M2二分类。每套分数为其基因log1p(CP10k)的算术均值，包含表达零值；CP10k分母来自原始全部36,601个feature的库大小。

'''
txt+='| 程序 | 固定基因 |\n| --- | --- |\n'
for k,gs in PANELS.items():txt+='| '+ZH[k]+' | '+', '.join(gs)+' |\n'
txt+=f'''
MHC-II分数反映转录层面抗原呈递相关信号，不能直接测量呈递能力；IL10、TGFB1、LGALS9及检查点配体具有情境依赖作用，不能把均值升高等同于免疫抑制。清道夫/LPL状态不是完整脂质代谢通路；SPP1/MMP12仅两基因代理。C1Q和摄取程序与最初巨噬门重叠，存在选择偏倚，不可据此独立证明补体/吞噬增强。缺氧和增殖标志也可能来自混合肿瘤或其他细胞。

**无法可靠评估线粒体自噬、氧化磷酸化、完整脂质代谢或全基因组富集**：现有AnnData只有179个选择标志，相关基因集不完整。PGAM5本身的已知功能不能替代本数据证据。未检出差异不等于真实功能相同。

## 比较、调整与稳健性

继续使用固定巨噬RNA门：C1QA/B/C、CD68、CSF1R、MERTK、MSR1均值≥切片Q75、≥3核心标志检出；Pan-myeloid≥切片中位数。PGAM5-high为门内PGAM5计数>0且归一化表达≥检出源Q75，包含阈值并列。未检出不称完全不表达。

主比较high vs所有其他巨噬RNA源，敏感性为①high vs其他PGAM5检出源；②排除同时满足T-RNA门的源后high vs其他；③限制到原作者MT≤5% QC子集。后两项保留原始源门及high标签，不重新求分位数，因此是固定标签的限制样本检验，不是重做全流程注释。另输出去除C004复发切片的主195项BH。

逐切片线性模型调整logUMI、MT%、巨噬核心、泛髓系、肝细胞、上皮、T、B、浆细胞、内皮、成纤维RNA代理及3×3粗空间区块。程序内基因从协变量的基因列表中删除，避免同一基因直接同时构成结果与调整项；原始源门保持不变。RNA代理不能当作反卷积细胞比例，回归调整不能使位点成为纯巨噬细胞。使用SVD去除冗余协变量，FWL获得调整后high系数。

每片、每比较4,999次高表达标签置换，种子281759103+GSM数字部分；只在3×3全片空间区块×当前源UMI二分位×巨噬核心二分位内交换标签，每次重新拟合调整后系数。双侧等尾p带+1校正，不假设条件零分布以0为中心。层内近似可交换性是前提；粗区块/协变量不能完全消除空间自相关、组织结构及混合成分的影响。

high≥5、对照≥10、残余自由度≥15且零分布非退化才报告检验p。主208项（16×13）和联合832项（16×13×4）分别BH；不可检验项p保持缺失，仅校正时计作1。缺失数据、不可检验和真实0分数明确区分。没有将置换零分布区间冒充效应置信区间。

55个标志基因额外提供调整系数、表达与检出率、分组原始计数求和/全库UMI求和的描述性伪bulk CP10k，以及加0.1 CP10k偏移后的log2比值；**没有基因级显著性检验，不能称为正式DEG列表**。见 [marker_gene_descriptive_contrasts.csv](marker_gene_descriptive_contrasts.csv)。所有变化均为关联；不存在虚拟敲除、蛋白或细胞功能实验。

## 来源与复现

官方数据：[GSE281759](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281759)，论文PMID 41950768。使用同仓库[前轮空间数据](../GSE281759_spatial_HCC_PGAM5/README.md)的16个AnnData，源提交 `43a6b17ed4da307f4dbf83250ce051a247ca2103`。输入SHA/字节见 [reference_manifest.json](reference_manifest.json)，原始GEO附件来源见前轮 [source_manifest.json](../GSE281759_spatial_HCC_PGAM5/source_manifest.json)。本轮NCBI FTP访问被环境策略拒绝，未取得全基因矩阵，未声称重新核验所有原始附件。

主QC为in_tissue、UMI≥500、基因≥200、MT≤25%，与论文MT≤5%不同；本轮增加MT5限制样本敏感性。A009N非肿瘤切片排除；B0013T1–T4与B013T5同患者B013，C004T1P/R是同患者原发/复发。没有病理ROI，不能把所有位点都称作恶性组织。

```bash
git clone https://github.com/183285068-droid/codex.git
cd codex/results/GSE281759_PGAM5_macrophage_function
python -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python analyze.py
OPENBLAS_NUM_THREADS=1 python plot.py
OPENBLAS_NUM_THREADS=1 python audit.py
python build_readme.py
python package_results.py
```

**审核{audit['checks']:,}项全部通过**：[audit_results.json](audit_results.json)。核查输入SHA、全库归一化、程序分数、原高表达定义、全部比较分母、独立完整OLS对FWL系数及前3次置换拟合、层内标签数量、等尾p、BH、患者符号翻转。审核复用了协变量构造代码，未声称模型选择经过外部验证。另人工检查主要图的数值、标注和排版。

完整结果包：[GSE281759_PGAM5_macrophage_function_results.zip](GSE281759_PGAM5_macrophage_function_results.zip)。原输入AnnData不重复打包，位于同仓库前轮目录。
'''
(ROOT/'README.md').write_text(txt)
print(summary)
