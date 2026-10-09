from pathlib import Path
import json,hashlib,shutil,zipfile,importlib.metadata as im
import pandas as pd
R=Path(__file__).resolve().parent
hashes={}
for name in ['GSE146115_HCC1-2-5-9_count_with_ERCC.txt.gz','GSE146115_raw_to_prcessed_data_column.txt.gz','GSE146115_family.soft.gz']:
 p=R/name;h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 hashes[name]={'bytes':p.stat().st_size,'sha256':h.hexdigest()}
(R/'input_sha256.json').write_text(json.dumps(hashes,indent=2))
(R/'environment_versions.txt').write_text('\n'.join(f'{p}=={im.version(p)}' for p in ['numpy','pandas','scipy','anndata','scanpy','igraph','statsmodels','matplotlib','nbformat','nbclient','ipykernel','threadpoolctl'])+'\n')
s=json.loads((R/'analysis_summary.json').read_text())
prep=json.loads((R/'preparation_summary.json').read_text());e=json.loads((R/'exact_detection_summary.json').read_text())
u=pd.read_csv(R/'GSE146115_HCC_tumor_PGAM5_upregulated.csv',index_col=0)
assert len(u)==5 and not u.ambiguous_symbol.any()
text='''# GSE146115：四例HCC肿瘤的PGAM5检出相关巨噬细胞候选分析

全部HCC1、HCC2、HCC5、HCC9肿瘤样本纳入。原始3200个细胞来自4位患者，各800细胞；16个GEO样本记录是各200细胞的测序分部，不是16个独立患者。Fluidigm C1+HTSeq为read counts，不是UMI。

'''
text+=f"基本质控后**{prep['QC_cells']}个细胞**，重新推定**{s['QC_inferred_macrophages']}个巨噬细胞**。PGAM5检出**{s['PGAM5_positive_cells']}个**、未检出**{s['PGAM5_undetected_cells']}个**。检出组16个来自HCC2、1个来自HCC9；HCC1及HCC5没有检出组细胞。\n\n"
text+='''主分析得到**5个上调候选、0个下调候选**（排除PGAM5）：

**STK26、ZNF785、MS4A1、LINC01776、TOP3B。**

补充Fisher双侧精确基因检出率检验中，**无非PGAM5基因通过BH FDR<0.05，这5个候选也全部不通过，Fisher FDR均为1.0**。该检验比较检出率而不是全表达分布，不能直接替代Wilcoxon，但结果不支持将这5个候选当作稳健signature。

## 原始输入处理与质控差异

92条ERCC-数字外源spike-in行在QC、标准化和差异检验前移除。人类ERCC1、ERCC2等基因保留。184个缺失计数全部位于被移除的ERCC对照行（2个细胞各92项），人类计数没有缺失或负值。

源矩阵未提供完整线粒体基因集，只识别到旧名称RNR1、RNR2两条线粒体rRNA。本次使用**(RNR1+RNR2)/人类总计数**作为不完整的线粒体比例下限，执行<20%过滤，同时保留检出人类基因行数≥500。不能将其解释为完整线粒体比例，也不能假设未提供的基因计数为0。因而线粒体QC不与此前完整MT-基因统计完全等价。

源gene名中27条为日期格式损坏/歧义名称，如1-Mar、2-Mar，其中重复名称可能来自不同基因。所有人类原始行保留为唯一row_ID，不猜测修复、不合并不同基因；全表及候选表含ambiguous_symbol标记。这27条不用于聚类特征，全部5个上调候选的名称均无此歧义。

官方raw-to-processed映射表的COL/ROW编号有前导零。对3200条有效细胞映射去除编号前导零后，与计数矩阵列一一核对；不将该映射文件的其他说明/校验部分当作细胞记录。

## 巨噬细胞为探索性重新注释

没有作者细胞类型标签。全体QC细胞标准化至10000+log1p，2000个Seurat高变基因（排除PGAM5及歧义名称），30PC、15近邻、Leiden resolution1、seed0。

髓系富集global cluster4进一步以1500HVG、15PC、10近邻、Leiden0.6、seed0细分；依据C1QA≥90%、C1QB≥70%、C1QC≥50%、CD68≥70%、CSF1R≥25%、TYROBP≥60%、CD1C<10%、FCN1<75%的cluster检出率，纳入subcluster0、1、2、4，排除FCN1较高/混合subcluster3。

这是目标基因无关marker审查后的探索性定义，以C1Q较高巨噬细胞为主，不保证覆盖全部巨噬细胞或排除所有双细胞。所有global/subcluster marker数值及纳入决定保留，可审查和重算。PGAM5不用于聚类。

## 按请求的差异分析

PGAM5原始read count>0定义RNA检出，=0为未检出，不等同于蛋白阴性。全部选定巨噬细胞合并，10000计数标准化+log1p，Wilcoxon ties校正，全部27135人类源基因行作BH校正。

标准：FDR <0.05、|Scanpy近似log2FC| ≥1、至少一组检出率 ≥10%。PGAM5保留完整及含PGAM5筛选表，从上/下调候选中排除。核糖体等基因保留。按要求不作患者配对、深度匹配、患者/样本/深度协变量校正。

补充Fisher精确检验基于同一批原始计数中的基因检出/未检出2×2表，双侧p，对全部27135行BH。检验的对象是检出率，和Wilcoxon表达秩不同。

## 解释限制

检出组仅17细胞，16来自同一患者；本次无法证明患者间可重复性。前4候选均仅在3/17个检出组细胞中出现，在124个对照中未检出；近似log2FC约24–26由近零分母和变换伪计数驱动，不应解释为可靠的千万倍生物学差异。TOP3B在5/17与3/124细胞中检出，也未通过精确检出率FDR。

MS4A1通常为B细胞相关转录本，不能直接认定为巨噬细胞特异表达。推定巨噬细胞中ALB广泛检出，环境RNA、吞噬RNA和双细胞来源未区分或校正。重新注释、局部线粒体质控、日期格式损坏、患者偏重和小分组均限制signature结论。因此只保留为按请求产生的探索性候选，**不建立已验证signature**。

## 导出与复现

包括完整差异表、含PGAM5筛选表、5个上调及空下调候选、每样本计数、细胞元数据、global及myeloid marker审查、ERCC排除及歧义gene清单、所有基因Fisher结果及候选复核、已执行notebook、PNG/PDF、代码、输入SHA256及版本。

复核了范围、QC下限代理、原始PGAM5分组、cluster归属、BH校正和全部候选成员，并用SciPy从原始QC计数重算全部5候选的秩检验p和log2FC。Notebook已执行，3张PNG均已视觉检查，未检查完整notebook网页交互布局。

来源：[GEO GSE146115](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146115)、[官方补充文件](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146115/suppl/)、[family SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146115/soft/GSE146115_family.soft.gz)，PMID33531041。

将input_sha256.json列出的3个官方gz文件放在脚本目录，安装environment_versions.txt，依次运行prepare_annotation.py、refine_myeloid.py、analyze_PGAM5.py、exact_detection_check.py、build_notebook.py。脚本相对于自身定位输入。ZIP不包含完整原始或QC计数矩阵；重跑脚本会生成这些文件。source_cell_header_mapping.csv保存已核对的3200条有效源映射。
'''
(R/'结果说明.md').write_text(text)
D=Path('/workspace/codex/results/GSE146115_HCC_PGAM5');D.mkdir(parents=True,exist_ok=True)
names=['analysis_summary.json','preparation_summary.json','exact_detection_summary.json','all_cell_metadata_before_QC.csv.gz','all_tumor_QC_cluster_metadata.csv.gz','myeloid_cluster_metadata.csv.gz','macrophage_analysis_metadata.csv.gz','global_cluster_marker_QA.csv','myeloid_cluster_marker_QA.csv','cluster_annotation_decisions.csv','TAM_marker_QA.csv','PGAM5_by_sample.csv','ambiguous_source_gene_symbols.csv','excluded_spikeins.csv','source_cell_header_mapping.csv','independent_calculation_checks.csv','input_sha256.json','environment_versions.txt','prepare_annotation.py','refine_myeloid.py','analyze_PGAM5.py','exact_detection_check.py','build_notebook.py','package_results.py','GEO_series_record.txt','GSE146115_HCC_tumor_PGAM5_analysis.ipynb','结果说明.md']
for p in list(R.glob('GSE146115_HCC_tumor_PGAM5*.csv'))+list(R.glob('*exact_detection_check.csv'))+list(R.glob('exact_detection_Fisher_all_genes.csv'))+list(R.glob('*.png'))+list(R.glob('*.pdf')):names.append(p.name)
for name in set(names):shutil.copy2(R/name,D/name)
zp=D/'GSE146115_HCC_PGAM5_results.zip'
with zipfile.ZipFile(zp,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(D.iterdir()):
  if p.is_file() and p!=zp:z.write(p,p.name)
with zipfile.ZipFile(zp) as z:assert z.testzip() is None
print(json.dumps({'files':len(list(D.iterdir())),'ZIP_MiB':zp.stat().st_size/1024**2},indent=2))
