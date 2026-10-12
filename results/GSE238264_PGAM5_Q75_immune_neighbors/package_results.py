from pathlib import Path
import json,hashlib,zipfile
import pandas as pd
R=Path(__file__).resolve().parent;s=pd.read_csv(R/'source_immune_neighborhoods.csv');t=pd.read_csv(R/'spatial_conditional_tests.csv');audit=pd.read_csv(R/'audit_checks.csv');assert audit.passed.all();he=json.loads((R/'HE_registration_checks.json').read_text());assert len(he)==7 and all(x['all_coordinates_inside_image'] for x in he)
h=s[s.source_group=='PGAM5_high'];edge=h.pivot(index='subtype',columns='sample',values='n_edges');types=['T','B','Plasma','NK','cDC1','cDC2','pDC','Monocyte','Neutrophil'];samples=['HCC1R','HCC2R','HCC3R','HCC4R','HCC5NR','HCC6NR','HCC7NR'];table=['| 目标RNA候选 | '+' | '.join(samples)+' |','|---|'+'---:|'*7]
for k in types:table.append('| '+k+' | '+' | '.join(str(int(edge.loc[k,x])) for x in samples)+' |')
body='''# GSE238264：Q75 PGAM5高表达巨噬RNA源与免疫亚群的空间邻近及候选通讯

**7张切片均完成Q75分析。PGAM5高表达源与T、B、浆细胞、NK、cDC1、cDC2、pDC、单核和中性粒细胞RNA候选区域存在不同程度的空间邻接，但本轮没有发现校正后显著的High相对其他PGAM5检出源的邻域富集。**候选配体–受体为RNA支持，不是已验证细胞通讯。

## 固定范围和观察单位

HCC1R–HCC4R为应答、HCC5NR–HCC7NR为未应答，来自新辅助卡博替尼＋纳武利尤单抗后的切除样本。本轮按用户要求分析各切片免疫邻域，治疗标签仅用于背景，未追加治疗应答组间检验。切片各对应1位患者，点和边不能作为独立患者重复。

源点完整沿用[原始基线](../GSE238264_PGAM5_spatial_response/README.md)：巨噬核心Q75、至少3核心标志、泛髓系≥中位数、真实空间有邻居；PGAM5原计数>0源内的log1p(CP10k)Q75定义High。7样本High数依次8/1/10/23/4/36/2，总84，与原门控逐点相同。本轮没有使用Q20/Q25/Q30/Q40等宽门控替代Q75。

所有程序从原始H5提取完整基因计数，以全库UMI归一化到CP10k后log1p；重复gene symbol合并计数。原质控不变，保留14,946点。每点混有多个细胞，PGAM5和配体不保证来自同一个巨噬细胞；“巨噬RNA源”是代理定义。

## 免疫候选门控

T门保留原CD3D/E/G、TRAC、TRBC1/2门。其余每个程序按各切片Q75、分数>0、至少2基因检出且有邻居定义，并要求谱系锚定标志。候选不互斥、不强制加和100%，同一混合空间点可能支持多个标签。

| 类别 | 程序基因 | 锚定条件 |
|---|---|---|
| B | MS4A1/CD79A/CD79B/CD19/CD22 | MS4A1或CD79A |
| Plasma | MZB1/JCHAIN/SDC1/XBP1/DERL3 | **MZB1且JCHAIN**，减少仅凭共享分泌程序命名 |
| NK | NKG7/GNLY/KLRD1/KLRF1/NCR1 | GNLY、KLRD1或NCR1；仍可能混合细胞毒T信号 |
| cDC1 | CLEC9A/XCR1/BATF3/IRF8/CADM1 | CLEC9A或XCR1 |
| cDC2 | CD1C/FCER1A/CLEC10A/CD1E | CD1C或FCER1A |
| pDC | CLEC4C/IL3RA/TCF4/GZMB | CLEC4C或IL3RA；不能独立确定纯pDC细胞 |
| Monocyte | FCN1/VCAN/S100A8/S100A9/LYZ | FCN1或VCAN |
| Neutrophil | FCGR3B/CSF3R/CXCR2/S100A8/S100A9 | FCGR3B、CSF3R或CXCR2 |

标志覆盖见marker_coverage.csv。单核与巨噬谱系连续、T/NK和树突状态重叠及环境RNA仍影响身份，不能据候选点数推断实际细胞数。没有匹配单细胞去卷积或蛋白来源验证；未将HE形态自行标注为肿瘤边界/TLS。

## High源→目标候选的真实第一圈有向边数

TABLE

无边不意味着细胞或通路不存在。cDC2虽在部分切片有候选点，Q75 High源没有合格直接邻接边；不能据此证明没有互作。边数可共享源/目标；免疫标签重叠会让同一边计入多类，不能跨行相加当作独立观察。

HCC4R具有较多High→T（83条）、单核（50条）及B（23条）边；HCC6NR以单核候选邻接为主（33条），High→T仅1条。两者差异仅为单病例描述，不构成应答/未应答总体机制。HCC2R只有1个High源，HCC7NR只有2个，保留定位与描述。

[源加权邻域对比图](figures/immune_neighborhood_coverage.png) · [全部邻域支持表](source_immune_neighborhoods.csv)

## 空间条件比较

每个源计算“该免疫候选邻居/全部保留邻居”，无此类邻居计0。主比较High vs Other detected，PGAM5零值作描述性参考。每源一个均值，避免直接对共享邻接边做独立t/Wilcoxon检验。

在切片3×3空间网格及格内候选源UMI中位数分层交换High/Other标签，保持每层High数；双侧统计围绕精确条件期望，避免将未调整全局差值直接当作零中心效应。配置≤4,999时全枚举，否则4,999次Monte Carlo；固定种子23826475起。该探索性方法依赖层内可交换性，粗分层不能保留全部局部相关，也不能解决混合细胞及临床混杂，不是患者总体因果估计。

7切片×9目标共63项固定BH家族，不可检验p仅校正时设1、结果保留NA。HCC2R/HCC7NR组数不足且没有可交换源标签，18项不检验；其余45项可计算（含10项所有结局0），其中35项有非零变异。最低p约0.05357，全部有效检验q=1。未发现Q75 High组特异富集，但检出能力低，不能认定不存在生物学关系。HCC3R仅12种条件配置、HCC5NR仅2种，尤其缺乏统计分辨率。

[p/q可视化](figures/spatial_immune_tests.png) · [完整p/q与状态表](spatial_conditional_tests.csv) · [配置支持度](permutation_support.csv)

## 14条候选配体–受体RNA支持

候选逐条核对LIANA 1.10.0 omni_resource快照：SPP1–CD44、LGALS9–HAVCR2、PVR–TIGIT、NECTIN2–TIGIT、CD274–PDCD1、CXCL9/10–CXCR3、CCL5–CCR5、CXCL12–CXCR4、CXCL13–CXCR5、TNFSF13B–TNFRSF13C、TNFSF13–TNFRSF17、TGFB1–TGFBR1/TGFBR2、IL10–IL10RA/IL10RB。未执行LIANA/CellChat算法，数据库来源只说明已知候选，不能证明本样本功能活性。

方向设为巨噬RNA源配体→相邻免疫候选受体。复合受体要求全部亚基检出，以最低亚基CP10k为代理；边RNA兼容得分=sqrt(配体CP10k×受体代理CP10k)，零表达边保留分母；每源先平均再组内平均。NaN/NA为无目标边，0为有边但没有共检出。候选边不等于独立样本或通讯概率，本轮LR仅描述，不另加稀疏亚群p值。

HCC4R示例：

- High→T的CXCL12–CXCR4为72/83共检出边，SPP1–CD44为58/83。
- High→B的CXCL12–CXCR4为22/23，SPP1–CD44为15/23；分别只是10个High源、21个B候选目标的邻接集合，不能把高百分比理解为强通讯概率。
- High→单核的CXCL12–CXCR4为44/50；HCC6NR同候选为6/33。不能用这两病例比例推断治疗应答因果。
- 浆细胞共检出需结合很小邻接边分母逐项查看，不据无边或低计数否定BAFF/APRIL通路。

RNA配体可来自源点内成纤维或其他细胞，尤其CXCL12不能自动归属巨噬细胞。共享标志和混合点使靶细胞归属不确定。

[候选LR共检出热图](figures/immune_LR_codetection.png) · [全部14候选各组/亚群表](LR_by_immune_subtype.csv) · selected_LR_pairs.csv

## HE融合图

HELINKS

每张图为九目标面板，蓝色半透明点为免疫RNA候选，橙圈为Q75 High源。用官方高分辨组织图像、原pixel坐标及各自tissue_hires_scalef，半径来自spot_diameter_fullres；未人工旋转、平移或推测病理ROI。所有源/目标均处于对应图像范围。PNG/PDF均可下载；HE图不替代单细胞身份验证。

## 核验和复现

AUDIT项核验通过：输入哈希、源分组冻结、完整库UMI、原始计数程序及锚定门控、独立KDTree几何邻接、邻域均值、全部LR逐边表达/复合亚基/聚合边数、条件p和63项BH。HE配准边界另见HE_registration_checks.json。

完整仓库中需保留相邻原基线、GSE245908固定数据库快照，以及7个已解压原包。设置GSE238264_INPUT为包含HCC1R至HCC7NR的目录；本结果包不重复大原始计数和图像。环境同requirements.txt，复现：

```bash
export GSE238264_INPUT=/你的路径/inputs
export MPLCONFIGDIR=/tmp/gse238264_mpl XDG_CACHE_HOME=/tmp/gse238264_cache
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
python analyze.py
python audit.py
python plot.py
python plot_statistics.py
python package_results.py
```

包内提供代码、逐点免疫门控、全部有向候选边/得分、邻域表、置换审计、统计及10组PNG/PDF，不声称已确认巨噬细胞对免疫亚群的功能调控。
'''
body=body.replace('TABLE','\n'.join(table)).replace('AUDIT',str(len(audit))).replace('HELINKS','\n'.join(f'- [{x}：HE免疫邻域图](figures/{x}_HE_immune_neighbors.png)' for x in samples));(R/'README.md').write_text(body)
z=R/'GSE238264_PGAM5_Q75_immune_neighbors_all_results.zip';files=[p for p in sorted(R.rglob('*')) if p.is_file() and p!=z and p.name!='package_manifest.json' and '__pycache__' not in p.parts];records=[dict(path=str(p.relative_to(R)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files];(R/'package_manifest.json').write_text(json.dumps(records,indent=2))
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as f:
 for p in files+[R/'package_manifest.json']:f.write(p,p.relative_to(R))
with zipfile.ZipFile(z) as f:
 assert f.testzip() is None
 for rec in records:assert hashlib.sha256(f.read(rec['path'])).hexdigest()==rec['sha256']
print('Verified',len(files),'files; ZIP',z.stat().st_size,'bytes')
