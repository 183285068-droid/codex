from pathlib import Path
import json,hashlib,zipfile
import pandas as pd
R=Path(__file__).resolve().parent;d=pd.read_csv(R/'patient_spatial_metrics.csv');a=pd.read_csv(R/'audit_checks.csv');assert a.passed.all();h=json.loads((R/'HE_registration_checks.json').read_text());assert len(h)==7 and all(x['all_coordinates_inside_image'] for x in h)
lines=['| 患者 | 应答 | QC空间点 | 巨噬RNA点 | PGAM5检出巨噬点 | 高表达点 | High源有T邻居 |','|---|---|---:|---:|---:|---:|---:|']
for row in d.itertuples():lines.append(f'| {row.sample} | {"应答" if row.response=="responder" else "未应答"} | {row.QC_spots} | {row.mac_spots} | {row.PGAM5_detected_mac_spots} | {row.PGAM5_high_spots} | {row.high_with_T_neighbors} |')
body='''# GSE238264：PGAM5高表达巨噬RNA空间特征与新辅助治疗应答

**7位患者全部纳入：HCC1R–HCC4R为应答，HCC5NR–HCC7NR为未应答。未发现PGAM5巨噬RNA检出比例与病理应答的显著关联。高表达源点稀少，部分患者的高表达邻域指标不能完成稳定的全7患者比较。**

## 病例和治疗定义

官方SOFT/XML一致确认7个不同患者对应7张Visium切片（GSM7661255–GSM7661261），均接受新辅助cabozantinib＋nivolumab，分析为治疗后切除标本。原试验摘要中的5/15应答是试验总体，不替代本空间子集的4/7应答。标签采用GEO的responder/non-responder，与复发或影像RECIST不混用。具体病理应答阈值和早期复发患者对应关系未在所上传元数据中提供，本轮不自行推定。

所有结果是治疗后关联，不能称为治疗前预测能力、疗效因果机制或已验证的临床生物标志。组织点混合多种细胞，PGAM5不能定位到其中某一个巨噬细胞；“巨噬/髓系RNA区域”是操作性代理标签。

## 完整输入与统一质控

四个分卷包按part00起顺序恢复并成功解包，其余3个原包直接解包。7样本均包含H5、MatrixMarket、条形码、特征表、真实空间位置、scale和HE图像。每个样本两种格式计数逐项相同，非负整数，原始矩阵为36,601个特征（不等于唯一gene symbol数）。重复symbol在程序计算时合并原计数；归一化分母始终为提供的完整矩阵总UMI，未用标志子集重新归一化。

QC为in_tissue=1、UMI≥500、检测特征≥200、MT≤25%。保留14,946个点。HCC2R仅保留1,207/2,766，测序/样本质量差异明显；没有按治疗组调整阈值。逐点total_counts、n_genes、pct_mt均公开。

TABLE

巨噬门控：C1QA/B/C、CD68、CSF1R、MERTK、MSR1平均log1p(CP10k)达到每张切片Q75，≥3标志检出；LST1/TYROBP/FCER1G/AIF1/LYZ/CTSS达到切片中位数；至少1个保留空间邻居。PGAM5原计数>0为检出，检出巨噬RNA源的归一化PGAM5 Q75为High，其他检出与零值分别保留。

切片内Q75是相对阈值，不是跨患者绝对高表达阈值。High/所有检出比例约25%是设计产生的，不能据此比较高表达亚群丰度；本轮主检出率的分母是全部巨噬RNA源，且报告原始数。表达稀疏与文库深度影响检出，巨噬门控也可能受混合RNA影响。

## 预设比较与统计

在读取七患者结果前固定4主指标：巨噬RNA点中PGAM5检出比例、High源邻居巨噬RNA比例、High源邻居T-RNA比例、High源邻居成纤维RNA平均程序。真实六角第一圈邻接，不用UMAP或条形码顺序推断坐标。每个源先汇总邻居，随后每患者平均；患者等权，不将重复边或数千组织点作为独立病例。

High邻域均值要求每患者至少3个High源。HCC2R只有1个，HCC7NR只有2个，这三项邻域主指标只能描述其他5位患者，不做选择性删除后显著性检验。图明确标注3/4及2/3可用病例，而不显示为完整组。表内NA是证据不足，不是表达0或邻居不存在。

PGAM5巨噬RNA检出率：应答组患者均值9.82%，未应答组24.90%；差值−15.08个百分点。全35种4/3标签分配的双侧绝对均值差检验p=0.22857，按4项预设主检验（不可检验p仅校正时设1）BH q=0.91429。未发现显著差异，也不证明两组相同。该无协变量精确检验依赖患者标签在零假设下可交换；本观察研究的临床/深度混杂、很小样本量及标签分辨率仍限制推断，不能据此建立治疗机制。

High占全QC点比例为次要描述/比较：应答0.422%、未应答0.600%，p=0.77143、6项预设次要家族q=1。其余5项High邻近T程序要求至少3个High源有T邻居，未应答病例均不足，不计算组间p。

可用5患者的High邻近T覆盖均值为应答29.17%、未应答6.48%；成纤维程序均值0.963 vs0.265；均为缺少两病例的描述，不能解读为显著应答机制或不受缺失影响的总体差异。未对这5例另做更有利的统计。没有建立肿瘤/间质病理ROI，成纤维RNA不能直接等同肿瘤边缘或治疗相关纤维化。

[患者层面对比图](figures/patient_response_primary_comparisons.png) · [全部指标](patient_spatial_metrics.csv) · [完整统计](response_association_tests.csv) · [35配置审计](exact_patient_permutation_audit.csv)

## HE融合可视化

每位患者提供4面板PNG/PDF：原始HE、PGAM5髓系源分组、T-RNA＋High源、成纤维RNA候选＋High源。读取对应tissue_hires_scalef，以x=pixel_col×scale、y=pixel_row×scale配准，点半径来自官方spot_diameter_fullres；无手动旋转/镜像/平移。7张图的坐标与原位置匹配，且均在原HE范围内。图像直观看到的区位差异没有独立病理ROI支持，不追加未经验证的病理命名。

HELINKS

橙圈=High源，蓝圈=其他PGAM5检出源，灰圈=巨噬RNA源中PGAM5未检出。T-RNA与成纤维候选用半透明蓝/紫区域显示，不能等同单细胞。T门控为CD3D/E/G、TRAC、TRBC1/2的Q75、≥2谱系标志和有邻居；耗竭/活化/效应/初始记忆/Treg程序用于探索性邻域描述。Treg支持另要求FOXP3检出。程序不是确定功能状态，不把细胞毒直接判定为CD8；本轮未进行配体–受体通讯推断。

## 核验与复现

AUDIT项输入与计算核验通过：输入哈希、原计数及两个矩阵格式、原始QC、全库归一化、门控、独立KDTree几何邻接、每源及每患者汇总、35配置精确p及预设BH。HE缩放/坐标边界校验另见HE_registration_checks.json。结果包保留每点注释及每条几何边，支持后续敏感性分析。

解压7个原包到同一目录，目录下应有HCC1R/…/HCC7NR/。原包较大，未重复放入ZIP；input_manifest.json记录全部使用文件的SHA256。新上传GEO SOFT和固定患者应答表随包附上。

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
export GSE238264_INPUT=/你的路径/inputs
export MPLCONFIGDIR=/tmp/gse238264_mpl XDG_CACHE_HOME=/tmp/gse238264_cache
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
python analyze.py
python audit.py
python plot.py
python build_report.py
```

此流程没有单细胞去卷积、恶性细胞CNV判定、病理复发标签或跨队列验证。本病例集可提供描述性空间RNA线索；当前证据不足以断言PGAM5高表达巨噬细胞促进/抑制治疗应答。
'''
body=body.replace('TABLE','\n'.join(lines)).replace('AUDIT',str(len(a))).replace('HELINKS','\n'.join(f'- [{sample}：HE空间特征](figures/{sample}_HE_spatial_features.png)' for sample in d['sample']));(R/'README.md').write_text(body)
z=R/'GSE238264_PGAM5_spatial_response_all_results.zip';files=[p for p in sorted(R.rglob('*')) if p.is_file() and p!=z and p.name!='package_manifest.json' and '__pycache__' not in p.parts];manifest=[dict(path=str(p.relative_to(R)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files];(R/'package_manifest.json').write_text(json.dumps(manifest,indent=2))
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as f:
 for p in files+[R/'package_manifest.json']:f.write(p,p.relative_to(R))
with zipfile.ZipFile(z) as f:
 assert f.testzip() is None
 for rec in manifest:assert hashlib.sha256(f.read(rec['path'])).hexdigest()==rec['sha256']
print('Verified ZIP:',len(files),'files;',z.stat().st_size,'bytes')
