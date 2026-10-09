from pathlib import Path
import json,hashlib,shutil,zipfile,importlib.metadata as im
import pandas as pd
R=Path(__file__).resolve().parent
hashes={}
for p in sorted(R.glob('GSM*.tar.gz'))+[R/'GSE242889_family.soft.gz']:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 hashes[p.name]={'bytes':p.stat().st_size,'sha256':h.hexdigest()}
(R/'input_sha256.json').write_text(json.dumps(hashes,indent=2))
(R/'environment_versions.txt').write_text('\n'.join(f'{p}=={im.version(p)}' for p in ['numpy','pandas','scipy','anndata','scanpy','igraph','statsmodels','matplotlib','nbformat','nbclient','ipykernel','threadpoolctl'])+'\n')
s=json.loads((R/'analysis_summary.json').read_text())
u=pd.read_csv(R/'GSE242889_HCC_tumor_PGAM5_upregulated.csv',index_col=0)
coverage=u[['gene','log2FC','FDR']].copy()
for folder in sorted((R/'counts').iterdir()):
 genes=pd.read_csv(folder/'genes.tsv',sep='\t',header=None)
 assert genes[1].eq('PGAM5').sum()==1
 coverage['present_in_'+folder.name.split('_')[0]]=coverage.index.isin(set(genes[0]))
coverage.to_csv(R/'candidate_gene_presence_by_sample.csv')
assert coverage.filter(like='present_in_').all().all()

text='''# GSE242889：全部HCC肿瘤样本的PGAM5⁺巨噬细胞差异分析

纳入1T、2T、3T、4T、5T，共5位HCC患者的5个肿瘤样本。按前一轮仅肿瘤组织的范围，排除配对邻近非肿瘤NT。所有MVI状态均纳入，不按MVI筛选、分组或建模；GEO样本特征没有逐样本MVI标签，未自行猜测标签。

'''
text+=f"原始肿瘤细胞{s['source_tumor_cells']:,}个，基本质控后{s['QC_tumor_cells']:,}个。识别**{s['QC_inferred_macrophages']:,}个推定巨噬细胞**；PGAM5检出**{s['PGAM5_positive_cells']}个**，未检出**{s['PGAM5_undetected_cells']:,}个**，全部5个样本贡献两组细胞。\n\n"
text+=f"筛选得到**{s['upregulated_excluding_PGAM5']}个上调候选、{s['downregulated_excluding_PGAM5']}个下调候选**（排除PGAM5）。前10上调候选："+'、'.join(s['top20_up'][:10])+'。\n\n'
text+='''## 与之前数据集不同：细胞类型为探索性重新注释

GEO仅提供各样本计数、gene和barcode文件，没有作者细胞类型标签。先在全部肿瘤QC细胞中按10,000计数标准化+log1p，取2000个Seurat高变基因，并排除PGAM5，然后30个PC、15近邻、igraph Leiden resolution=1、随机种子0聚类。

依据不使用PGAM5的marker检查，选取C1QA/B/C检出率均≥75%、CSF1R≥60%、TYROBP≥90%、CD1C<40%、FCN1<50%的cluster：6、7、19、21、23、29。排除CD1C/FCER1A高表达的DC样cluster 5、14、20以及FCN1高表达的单核细胞样cluster 16。所有cluster的marker数值和纳入决定均保留。

这是**标志基因推定、以C1Q高表达为主的巨噬细胞群**，不等同于作者原始注释，也不保证覆盖C1Q低表达的全部巨噬细胞。阈值是目标基因无关marker审查后的探索性定义，不能称为预先验证的分类器。

## 差异分析方法

- 基本QC：检出基因数≥500、线粒体计数比例<20%。
- PGAM5原始UMI>0定义RNA检出组，=0为RNA未检出，不等同于蛋白阴性。
- 各样本通过Ensembl ID外连接，缺失feature行填0；不同样本feature名单不同，但此次82个上调候选均在全部5个输入gene名单中存在。具体检查见candidate_gene_presence_by_sample.csv。
- 所有选定巨噬细胞合并，标准化至每细胞10,000计数+log1p，Wilcoxon校正ties，全部58,336基因BH校正。
- FDR <0.05、|Scanpy近似log2FC| ≥1、至少一组检出率 ≥10%。
- 不进行患者配对、深度匹配、MVI/患者/样本/测序深度协变量校正。
- PGAM5保留在full_DE及DE_including_PGAM5表，从候选上/下调表排除。核糖体、线粒体基因保留。
- log2FC为平均log表达回变换后的近似比值，不是算术均值之比。

## 主要解释限制

这些结果只能作为探索性PGAM5 RNA检出相关候选。该数据集存在**广泛ALB RNA背景**，在巨噬细胞cluster中的检出率约89%–100%；尚未区分环境RNA、双细胞或吞噬来源，未进行环境RNA或双细胞校正。C1Q/CSF1R/TYROBP支持髓系巨噬细胞身份，但不能证明这些细胞或差异基因不受混合RNA影响。

另外，注释只覆盖推定C1Q较高巨噬细胞，患者/亚型/深度效应未调整，输入基因名单不同。不能据此认定已验证的PGAM5⁺巨噬细胞signature、患者间可重复性、细胞特异性或因果机制。

## 结果、复核与复现

包含完整差异表、含PGAM5的筛选表、上调/下调候选、top20、每样本分组数量、巨噬细胞元数据、全部cluster marker及纳入决定、候选gene输入覆盖检查、已执行notebook、三组PNG/PDF和代码。ZIP不含完整原始矩阵或QC计数h5ad。

复核检查样本范围、QC、原始PGAM5计数、marker-cluster归属、完整候选集合及BH；另从原始QC计数用SciPy重算前10候选的秩检验p值和log2FC。notebook已执行，3张PNG均已视觉检查，未检查完整notebook网页交互布局。

来源：[GEO GSE242889](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242889)、[官方family SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242889/soft/GSE242889_family.soft.gz)，PMID 37972953。下载的五个GSM tar.gz文件实际是未压缩的POSIX tar，不能强制按gzip解压。

按environment_versions.txt安装依赖，依次运行python download_extract.py、python prepare_annotation.py、python analyze_PGAM5.py、python build_notebook.py。脚本相对于自身目录定位输入。input_sha256.json保存原始文件SHA256，independent_calculation_checks.csv保存复核。download_extract.py对已有下载文件直接使用；如下载中断留下不完整文件，删除该文件后重新下载。下载脚本使用TLS并验证安全解压路径。
'''
(R/'结果说明.md').write_text(text)
D=Path('/workspace/codex/results/GSE242889_HCC_tumor_PGAM5');D.mkdir(parents=True,exist_ok=True)
names=['analysis_summary.json','all_tumor_cell_metadata_before_QC.csv.gz','all_tumor_QC_cluster_metadata.csv.gz','macrophage_analysis_metadata.csv.gz','global_cluster_marker_QA.csv','cluster_annotation_decisions.csv','TAM_marker_QA.csv','candidate_gene_presence_by_sample.csv','PGAM5_by_sample.csv','independent_calculation_checks.csv','input_sha256.json','environment_versions.txt','download_extract.py','prepare_annotation.py','analyze_PGAM5.py','build_notebook.py','package_results.py','GEO_series_record.txt','GSE242889_HCC_tumor_PGAM5_analysis.ipynb','结果说明.md']
for p in list(R.glob('GSE242889_HCC_tumor_PGAM5*.csv'))+list(R.glob('*.png'))+list(R.glob('*.pdf')):names.append(p.name)
for name in set(names):shutil.copy2(R/name,D/name)
zp=D/'GSE242889_HCC_tumor_PGAM5_results.zip'
with zipfile.ZipFile(zp,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(D.iterdir()):
  if p.is_file() and p!=zp:z.write(p,p.name)
with zipfile.ZipFile(zp) as z:assert z.testzip() is None
print(json.dumps({'files':len(list(D.iterdir())),'ZIP_MiB':zp.stat().st_size/1024**2},indent=2))
