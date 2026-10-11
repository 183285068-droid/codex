from pathlib import Path
import hashlib
import json
import zipfile
import pandas as pd

ROOT=Path(__file__).resolve().parent
q=pd.read_csv(ROOT/'sample_QC_summary.csv')
c=pd.read_csv(ROOT/'J259X_AVITI_cluster_profiles.csv')
n=pd.read_csv(ROOT/'expression_cluster_names.csv')
tech=pd.read_csv(ROOT/'paired_platform_summary.csv')
validation=json.loads((ROOT/'validation.json').read_text())
lines=['# GSE283406肿瘤病例：表达预分析与空间数据缺口',
'',
'已完成公开表达矩阵的质控、表达聚类、标志基因、PGAM5及免疫RNA程序可视化和同位点测序技术比较。**真正的组织坐标表达图、空间域、邻接巨噬–T细胞候选互作分析尚未完成：GEO公开文件及论文链接的作者代码仓库没有提供barcode对应坐标、比例因子或J259B逐位点病理标签。**UMAP不是组织位置，不据此推断空间相邻。',
'',
'## 病例和样本范围',
'',
'官方总体设计仅包含**一位匿名、HBV及饮酒背景的中分化HCC病例**，两张切片：J259X仅肿瘤；J259B包含肿瘤及邻近非肿瘤组织。每张切片均有NextSeq标准测序与AVITI深度测序。4个GSM是2张切片的技术比较，不是4位独立患者。',
'',
'主分析使用AVITI的J259X肿瘤切片；J259B作混合组织背景展示，不能将其全部点纳入“病理确认肿瘤点”，也未依据表达群替代病理标签。原论文两张切片共6,320个输入位点；本次AVITI质控后共6,294个。',
'',
'| GSM | 切片 | 测序 | 输入点 | 质控点 | PGAM5检出点 | 检出比例 |',
'|---|---|---|---:|---:|---:|---:|']
for r in q.itertuples():
    lines.append(f'| {r.gsm} | {r.section} | {r.platform} | {r.input_spots} | {r.qc_spots} | {r.PGAM5_detected} | {r.PGAM5_detection_pct:.2f}% |')
lines += ['',
'原研究：Yu et al., *Frontiers in Cell and Developmental Biology*, 2025，DOI [10.3389/fcell.2025.1600129](https://doi.org/10.3389/fcell.2025.1600129)，[PMID 40530332](https://pubmed.ncbi.nlm.nih.gov/40530332/)。已通过NCBI公开接口读取并保存正式全文XML；原论文报告16个整合空间表达群及CARD/CellChat分析，这些是作者结果，**不是本次重新计算或复现的结果**。',
'',
'## 可用文件与缺失文件',
'',
'已获取并SHA256核验全部20个GEO附件，总计约133 MiB；每个GSM只有barcodes.tsv.gz、features.tsv.gz、matrix.mtx.gz、组织hires和lowres PNG。官方归档[filelist](provenance/filelist.txt)同样仅列这20个附件。',
'',
'作者代码仓库：[Spatial_transcriptomics_deep_vs_standard](https://github.com/SilviaLiu12345/Spatial_transcriptomics_deep_vs_standard)，已核查commit `e9fcbdca77c9f0fad02618bbe54ca4e3591dbd5c`的全部跟踪文件。其Seurat.R从作者本地`../code/J259B/outs/`和`../code/J259X/outs/`调用Load10X_Spatial，这些outs目录、含图像坐标的Seurat对象以及逐点TumorLabel不在该公开仓库中。本次只读取作者文本定位方法和输入依赖，未执行外部脚本。',
'',
'要补齐本轮空间分析，需要J259X/J259B的：',
'',
'1. `spatial/tissue_positions.csv`或旧版`tissue_positions_list.csv`，含barcode、in_tissue、array行列及原图像素行列。',
'2. `spatial/scalefactors_json.json`，用于将原图坐标映射到已公开的组织图像。',
'3. J259B的逐位点肿瘤/非肿瘤病理标签，或可配准的病理ROI。',
'',
'也可提供带图像与坐标的Seurat RDS / 完整AnnData空间对象。不能从barcodes排序或UMAP恢复真实组织坐标，不能把没有标签的J259B表达群自动视为病理肿瘤边界。公开图像可供形态展示，但没有barcode-像素映射就无法叠加真实表达。',
'',
'## 质控与方法',
'',
'数据采用10x Visium FFPE探针检测，面板18,085个Ensembl特征。固定主质控为UMI≥500、检出特征≥200；宽松质控UMI≥250、基因≥100的保留数量也单独报告，未将宽松结果混入主分析。MT面板比例仅报告，未套用poly(A)单细胞的25%线粒体阈值。',
'',
'有3个重复gene symbol（TBCE、HSPA14、TMSB15B），对应不同Ensembl ID。保留原特征和计数，以唯一后缀区分，不删除或错误汇总；PGAM5只有1个唯一特征。所有归一化均以全18,085特征总UMI计算CP10k，再log1p；保留PGAM5零值。',
'',
'每张AVITI切片分别进行表达分析：至少10点检出的基因中排除MT及核糖体候选，选择2,000个HVG（Seurat flavor），scale截断10，PCA40维，15邻居/前30PC，Leiden分辨率0.6（igraph、2次迭代），UMAP固定种子283406。J259X为3群，J259B为4群。**这是本轮参数下的表达分群，未验证跨分辨率或重采样稳定性；不代表只有这些生物学群，也不对应作者原论文16群编号。**没有空间坐标参与聚类或地图绘制。',
'',
'每群前100标志使用Wilcoxon、tie correction、群内BH筛选排序。聚类和标志选择使用同一批位点，p/q用于探索性标志排序，不能解释为独立患者级验证；未经空间相关校正，也未将其作为肿瘤病例总体统计推断。',
'',
'## 肿瘤切片J259X的表达特征',
'',
'| 表达群 | 暂定混合RNA描述 | 点数 | PGAM5检出比例 | PGAM5平均表达 |',
'|---|---|---:|---:|---:|']
descriptions={'E0':'ALB/PPARGC1A：肝源代谢RNA，伴相对较高T谱系信号',
              'E1':'IGHG1/IGKC：免疫球蛋白与基质混合RNA',
              'E2':'REG3A/GLUL：肝源代谢与分泌相关RNA'}
for r in c.itertuples():
    lines.append(f'| {r.expression_cluster} | {descriptions[r.expression_cluster]} | {r.n_spots} | {r.PGAM5_detection_pct:.2f}% | {r.PGAM5_mean_log1p_CP10k:.4f} |')
lines += ['',
'平均表达为含零值的log1p(CP10k)均值。E0中PGAM5检出率及平均表达高于另两群，但本轮没有作组间空间统计检验；该观察不能证明PGAM5来自T细胞、巨噬细胞或肿瘤细胞，也不能证明促进或抑制免疫。',
'',
'J259B的4个表达群包括CYP2E1/HAMP代谢/基质混合、IGHG1/IGKC免疫球蛋白相关、ALB/HRG肝源分泌/急性期相关、COL1A1/COL3A1基质相关RNA；这些不是病理肿瘤/非肿瘤划分。每对命名标志均核实对应群平均表达高于其他群；程序相对信号及命名证据见[核查表](expression_cluster_naming_evidence.csv)。名称为暂定混合RNA描述，不是单细胞注释、独立克隆或功能验证。',
'',
'## 巨噬、T细胞、活化与耗竭RNA程序',
'',
'各面板见[programs.json](programs.json)，每个程序为**面板中实际存在基因**的log1p(CP10k)均值；未检出值仍为0，缺失于探针面板的基因单独记missing而不是0。覆盖表见[program_gene_coverage.csv](program_gene_coverage.csv)。',
'',
'重要覆盖限制：FCN1、KRT18、LDHA、BNIP3、SAA1、SAA2缺失；HLA-DRA/HLA-DRB1/HLA-DPA1/HLA-DPB1缺失，因此`Antigen_presentation`输出事实上仅为CD74单基因代理，**不能作为完整抗原呈递程序或MHC-II功能评价**。完整覆盖的巨噬核心、T谱系、早期活化、核心耗竭面板用于本轮主RNA门控。不要直接跨平台/跨研究比较这些程序的绝对值。',
'',
'巨噬RNA富集门控：C1QA/C1QB/C1QC/CD68/CSF1R/MERTK/MSR1均值达到切片75分位，≥3标志检出，泛髓系均值≥切片中位数。T谱系门控：CD3D/CD3E/CD3G/TRAC/TRBC1/TRBC2均值达到75分位，≥2标志检出。本次缺少坐标，因此不附加“有邻居”条件，不能据这些门控计算邻接关系。',
'',
'早期活化：CD69/IL2RA/TNFRSF9/ICOS/CD40LG/IL2。核心耗竭相关：PDCD1/LAG3/TIGIT/CTLA4/TOX/CXCL13。在T谱系门控内，程序≥75分位、>0且≥2程序标志检出才记“支持”；活化和耗竭可重叠。单次RNA程序不能确认终末耗竭。',
'',
'J259X AVITI：483个巨噬RNA富集点，其中83点检出PGAM5（17.18%），含零值PGAM5均值0.1536；575个T谱系富集点，36个活化支持点和37个耗竭相关支持点。点/程序可以重叠，**不是单细胞数或浸润细胞比例**。在巨噬RNA富集检出点中最高25%的PGAM5标为相对高表达，仅用于表达状态描述；尚无法定位其与耗竭T区域的空间关系。',
'',
'## NextSeq–AVITI同位点技术比较',
'',
'两种测序的原始barcode集合完全相同；解压后的hires/lowres图像均逐字节一致。技术表达按section+barcode配对，不合并为额外生物学重复。主表比较仅包含两种测序均通过QC的交集位点。',
'',
'| 切片 | 共同QC点 | PGAM5归一化Spearman rho | 两者均检出 | 仅AVITI检出 | 仅NextSeq检出 |',
'|---|---:|---:|---:|---:|---:|']
for r in tech.itertuples():
    lines.append(f'| {r.section} | {r.shared_QC_spots} | {r.PGAM5_normalized_spearman:.3f} | {r.PGAM5_both_detected} | {r.PGAM5_only_AVITI} | {r.PGAM5_only_NextSeq} |')
lines += ['',
'同位点一致性支持本次PGAM5测量的技术可重复性，不构成独立患者验证、细胞归属验证或通讯证据。未对空间缺失的点实施邻接检验、Moran’s I或配体–受体空间边分析，也没有运行CARD、CellChat或已验证细胞丰度反卷积。',
'',
'## 图表与核验',
'',
'- [表达群UMAP及名称](figures/AVITI_expression_clusters_UMAP.png)',
'- [PGAM5、巨噬、T与耗竭RNA表达UMAP](figures/PGAM5_immune_expression_UMAP.png)',
'- [标志基因dotplot](figures/AVITI_marker_dotplot.png)',
'- [各表达群PGAM5比较](figures/PGAM5_expression_group_comparison.png)',
'- [QC与同位点技术比较](figures/QC_and_paired_platform_comparison.png)',
'- [原组织图像（没有表达叠加）](figures/original_tissue_images_no_expression_overlay.png)',
'',
f'已通过{validation["passed"]}项原始输入独立检查：附件SHA256/大小、原始UMI与零值、全基因归一化、程序覆盖/分值、门控、portable子集计数、有限UMAP坐标、同位点配对和同图像确认、命名标志富集。见[audit_checks.csv](audit_checks.csv)。所有PNG实际检查；所有PDF已导出。',
'',
'## 复现和下载',
'',
'```bash',
'python3 -m venv .venv',
'source .venv/bin/activate',
'pip install -r requirements.txt',
'export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 NUMBA_NUM_THREADS=4',
'python download_inputs.py',
'python analyze.py',
'python plot.py',
'python audit.py',
'python build_report.py',
'```',
'',
'可设置`GSE283406_WORK`到现有可写数据目录（本次/workspace/scratch/GSE283406），否则默认本目录work/。先下载并核验公开输入；本结果包含图表PNG/PDF、表达和QC表、带原计数及全基因归一化的portable marker h5ad、方法、正式文章及GEO来源记录，**不包含133 MiB原始矩阵/图像附件或大型全基因缓存**。h5ad是表达标志子集，不可重新当全基因库归一化，也不含空间坐标。',
'',
'整个病例仅1位患者，空间研究的生物学重复有限。当前完成的范围是公开空间测序数据的表达预分析；等待坐标、比例因子和病理标签补齐后才能继续真实空间特征分析。',
]
(ROOT/'README.md').write_text('\n'.join(lines)+'\n')
zipname='GSE283406_HCC_expression_preanalysis_results.zip'
files=[p for p in ROOT.rglob('*') if p.is_file() and not any(x in p.parts for x in ['__pycache__','inputs','work']) and p.name not in [zipname,'package_manifest.csv']]
manifest=pd.DataFrame([{'file':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files])
manifest.to_csv(ROOT/'package_manifest.csv',index=False)
with zipfile.ZipFile(ROOT/zipname,'w',compression=zipfile.ZIP_DEFLATED) as z:
    for p in files+[ROOT/'package_manifest.csv']:z.write(p,p.relative_to(ROOT))
with zipfile.ZipFile(ROOT/zipname) as z:
    assert z.testzip() is None
    for row in manifest.itertuples():assert hashlib.sha256(z.read(row.file)).hexdigest()==row.sha256
print(f'Packaged {len(files)+1} files; {(ROOT/zipname).stat().st_size/2**20:.2f} MiB; CRC and SHA256 verified.')
