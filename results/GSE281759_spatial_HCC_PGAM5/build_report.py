from pathlib import Path
import json
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parent
q=pd.read_csv(ROOT/'sample_QC_summary.csv');s=pd.read_csv(ROOT/'PGAM5_spatial_autocorrelation.csv');lr=pd.read_csv(ROOT/'candidate_LR_codetection_summary.csv');regions=pd.read_csv(ROOT/'PGAM5_T_neighbor_source_summary.csv');audit=json.loads((ROOT/'audit_results.json').read_text())
def md(t):
    names=list(t.columns);rows=['| '+' | '.join(names)+' |','| '+' | '.join(['---']*len(names))+' |']
    for row in t.itertuples(index=False,name=None):rows.append('| '+' | '.join('—' if pd.isna(v) else str(v) for v in row)+' |')
    return '\n'.join(rows)
patient=q.groupby('patient',sort=False).agg(sections=('sample','size'),qc_spots=('qc_spots','sum'),PGAM5_detected=('PGAM5_detected','sum'),Mac_RNA_spots=('Mac_RNA_enriched','sum'),PGAM5_high_mac_spots=('PGAM5_high_mac_RNA_spots','sum'),Mac_Tex_neighbor_sources=('mac_sources_with_Tex_neighbors','sum'),PGAM5_high_Tex_neighbor_sources=('high_mac_sources_with_Tex_neighbors','sum')).reset_index()
patient.to_csv(ROOT/'patient_descriptive_summary.csv',index=False)
sections=q[['sample','tumor_status','qc_spots','expression_clusters','PGAM5_detection_pct','T_activation_supported','T_exhaustion_supported']].copy();sections['PGAM5_detection_pct']=sections.PGAM5_detection_pct.round(2)
stat=s[s.method.eq('UMI_MT_stratified')][['sample','Moran_I','permutation_p','BH_q_32_tests']].copy();stat.iloc[:,1:]=stat.iloc[:,1:].round(5)
d=lr[(lr.source_group=='PGAM5_high')&(lr.target_state=='Exhaustion_supported')].groupby('pair',sort=False)[['eligible_edges','supported_edges','supported_unique_sources']].sum().sort_values('supported_edges',ascending=False)
d.to_csv(ROOT/'PGAM5_high_Tex_candidate_summary.csv');examples=d.head(8).reset_index()
reg=regions[regions.region.eq('Mac_with_Tex_neighbor')][['sample','n_source_spots','PGAM5_detection_pct','PGAM5_mean_log1p_CP10k']].copy();reg.iloc[:,2:]=reg.iloc[:,2:].round(4)
text=f'''# GSE281759：HCC肿瘤切片空间RNA分析

已完成4位患者16张肿瘤切片（15张原发、1张复发）的质控、逐切片表达降维聚类、真实坐标组织图叠加、PGAM5空间检验，以及巨噬细胞RNA与活化/耗竭T细胞RNA邻近区域的描述分析。

**主要结果：**25,510个输入位点中保留25,354个；9,758个检出PGAM5 RNA（{100*q.PGAM5_detected.sum()/q.qc_spots.sum():.2f}%）。PGAM5全位点表达的Moran空间自相关在32项检验BH校正后全部不显著（最小q={s.BH_q_32_tests.min():.5f}）。这不能推出不存在局部生物学信号，也不代表进行了PGAM5巨噬细胞专属空间富集检验。

巨噬细胞RNA富集位点4,091个，PGAM5高表达子集382个；743个巨噬细胞RNA位点有耗竭支持的T-RNA直接邻居，其中74个属于PGAM5高表达子集。**这是混合位点RNA信号共定位，不是已证实的PGAM5⁺巨噬细胞身份、细胞浸润百分比或抑制T细胞的因果证据。**

## 图表入口

- [队列质控、PGAM5检出与空间检验](figures/cohort_QC_PGAM5_spatial_statistics.png)
- [A009患者空间总览](figures/A009_PGAM5_immune_spatial_overview.png)
- [B013患者空间总览](figures/B013_PGAM5_immune_spatial_overview.png)
- [C004原发/复发空间总览](figures/C004_PGAM5_immune_spatial_overview.png)
- [C008患者空间总览](figures/C008_PGAM5_immune_spatial_overview.png)
- [PGAM5-high→耗竭T-RNA候选配体受体共检出](figures/PGAM5_high_Tex_candidate_LR_codetection.png)
- 每张切片另有`{{sample}}_spatial_RNA_maps.png/pdf`（分群、PGAM5、巨噬、T、活化、耗竭六面板）及`{{sample}}_UMAP_marker_dotplot.png/pdf`，位于[figures](figures)。UMAP是表达嵌入，组织图使用真实Visium坐标。

## 患者及切片范围

{md(patient)}

{md(sections)}

排除A009N邻近非肿瘤切片。GEO中B0013T1–T4与B013T5命名不同，依据系列/论文统一患者编号为B013，保留原切片名称。C004T1P为原发，C004T1R为复发；不能把两张当作两位患者。系列共4位患者，多切片/位点不是独立病例，不进行位点充当病例的复发预后检验。

## 数据、质控与表达分群

全部144个官方附件已取得并记录原始字节数及SHA256：[source_manifest.json](source_manifest.json)。每张切片包含矩阵、barcode、36,601条feature、原始`tissue_positions_list.csv`、scalefactors、低/高分辨组织图及对齐图。条码一一匹配，10个重复symbol以唯一feature ID保留；PGAM5只有一个feature。

主质控采用前述空间分析规则：in_tissue=1、UMI≥500、检出基因≥200、线粒体计数比例≤25%。原论文使用基因≥200、线粒体≤5%，本轮不是其精确复现；5%会使A009T2仅保留131/1,565、C008T3仅保留74/1,494个位点。主规则尽量保持空间覆盖，同时[sample_QC_summary.csv](sample_QC_summary.csv)报告原论文5%规则（未额外加UMI500）与30%敏感性规则的保留数；这些是覆盖敏感性核查，**尚未重新计算所有5%阈值下的分群/空间统计**。较高MT切片的区域可能受应激或质量影响。

全36,601基因库总计数用于CP10k，再log1p。每切片独立选2,000个高变基因（≥10个位点检出，排除MT和典型核糖体基因），标准化、PCA40、30PC/15邻居、Leiden0.6、UMAP，随机种子281759。每片3–6个S组，**不同切片的S0等编号不对应同一群**。Wilcoxon给出每群top100标志；[expression_group_names.csv](expression_group_names.csv)以非MT/核糖体正向标志基因对命名为暂定混合RNA区域，非单细胞纯细胞类型注释，也不是肿瘤克隆鉴定。肝细胞/上皮、间质/内皮、髓系、T/B/浆细胞等程序均在[programs.json](programs.json)及各片cluster_profiles中报告。Dotplot点面积为该群原始计数>0的比例，颜色为包含零值的平均log1p(CP10k)。

## RNA程序、门控与空间邻接

19套固定基因程序每套使用其基因log1p(CP10k)算术均值，保留零值；所有请求基因都存在于feature表，仍可能全零或低检出。基因全集和覆盖表公开。

- 巨噬RNA门：Macrophage_core≥切片Q75、至少3个核心基因检出、Pan_myeloid≥切片中位数。核心为C1QA/B/C、CD68、CSF1R、MERTK、MSR1。
- T-RNA门：CD3D/E/G、TRAC、TRBC1/2均值≥切片Q75，至少2个谱系基因检出。
- 活化/耗竭/效应支持：T门内对应程序≥T门Q75、分数>0且至少2个程序基因检出；各状态可重叠。活化为CD69/IL2RA/TNFRSF9/ICOS/CD40LG/IL2；耗竭核心为PDCD1/LAG3/TIGIT/CTLA4/TOX/CXCL13。实际耗竭定义不是所有T细胞终末耗竭的充分证据。
- PGAM5-high：在巨噬RNA门内PGAM5原始计数>0的位点，归一化表达≥检出位点Q75；其余分Other_detected与Undetected。分位门逐切片计算，**不是可跨切片直接解释的绝对细胞阈值**。未检出不等于RNA完全不存在。

直接邻接来自原Visium六边形格：array row/col偏移(0,±2)、(±1,±1)，仅保留主QC位点；无自环，每位点最多6个邻居。按真实像素坐标乘原低分辨scalefactor叠加原组织图，未使用UMAP或人为生成位置。组织图中的肿瘤切片含混合组织，**未获得独立病理ROI，不能称每个位点均为恶性细胞或定义MVI区域**。

空间连续色标跨片共用；显示上限取各片99分位的最大值，超限仅绘图截断。低分辨图供复查，原始高分辨图可由下载脚本获取。便携AnnData保留原始计数、来自全基因库的归一化表达、真实坐标、图像与scale、邻接和UMAP；**仅为marker子集，不能用子集总计数重新归一化或复算全基因HVG/差异**。

## PGAM5空间自相关检验

{md(stat)}

Moran I使用所有主QC位点的PGAM5归一化表达与二元直接邻接。每片两类零假设：①全片随机置换；②仅在UMI五分位×MT三分位的最多15层内置换PGAM5表达，保留每层表达分布。各4,999次；双侧p比较距**各自精确置换期望值**的偏离，采用(p计数+1)/5,000。条件零模型的中心并不固定为-1/(n-1)；计算考虑同层无放回协方差。32项合并BH，所有q≥0.26453。条件零模型控制粗粒度测序质量分层，未完全消除组织结构、混合细胞比例和连续测序深度影响，也不是因果检验。

## 耗竭T-RNA邻近区域的PGAM5表达

下表每行源位点去重；包括PGAM5零值。无符合区域时表达统计缺失，不能记为零。

{md(reg)}

邻近源位点在A009共412个、B013共89个、C008共242个，C004两片本轮门控未识别符合耗竭支持的T位点。后者不证明组织中不存在耗竭T细胞；低检出和阈值可使其不可识别。

PGAM5-high源→耗竭支持T-RNA靶的直接边共105条（来自74个去重源位点）。21对来自固定LIANA1.10.0资源快照的候选中，下列RNA共检出较多：

{md(examples)}

配体与受体全部亚基均>0才计共检出；边分数为sqrt(源配体CP10k×靶受体亚基最小CP10k)。此处是固定候选的自定义RNA共检出描述，**没有运行LIANA/CellChat通讯推断，也没有互作显著性检验**。105条边重复共享源或靶，不能当作105次独立生物学重复。PVR/NECTIN2–TIGIT、SPP1–CD44等可用于后续候选检查，不能据此证明PGAM5改变通讯、T细胞活性或耗竭。PDCD1/TIGIT/CTLA4也用于T状态判定，候选支持并非独立耗竭验证。

## 核验和可复现

独立原始矩阵核验通过{audit['checks_passed']:,}项，见[audit_results.json](audit_results.json)：144文件长度/SHA、逐片QC条码、全库CP10k、19套程序、原坐标/scale/图像、独立KDTree格邻接、门控、PGAM5分层、置换完整性/分层不变/精确中心/p/BH、所有21候选受体亚基计分及汇总。图表已人工检查空间叠加、标注与dotplot。软件版本固定于[requirements.txt](requirements.txt)，Python3.12。

```bash
cd results/GSE281759_spatial_HCC_PGAM5
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
export GSE281759_WORK="$PWD/work"
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
export MPLCONFIGDIR="$GSE281759_WORK/matplotlib"
export XDG_CACHE_HOME="$GSE281759_WORK/cache"
mkdir -p "$GSE281759_WORK"
python download_inputs.py
python preprocess.py
python spatial_statistics.py
python plot_results.py
python audit.py
python build_report.py
python package_results.py
```

原始附件约454MiB、scratch全基因对象不加入Git；可按manifest重取。便携对象和分片完整表已发布。各患者包与公共总览包须解压到同一目录以重新绘图/复核，或直接克隆仓库。

## 来源

- [官方GEO系列GSE281759](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281759)及GSM8627601–GSM8627617；原始SOFT、文件清单在[provenance](provenance)。
- 原论文：[PMID41950768](https://pubmed.ncbi.nlm.nih.gov/41950768/)，[DOI10.1016/j.jhepr.2026.101790](https://doi.org/10.1016/j.jhepr.2026.101790)，PMC13087693；官方全文/摘要XML保存在provenance。论文重点为INTS6相关内皮/肿瘤区域与复发；本轮PGAM5分析独立于论文结论，未复现其蛋白、临床、CARD/CellTrek/CellChat流程。
'''
(ROOT/'README.md').write_text(text)
