from pathlib import Path
import sys,json,hashlib,shutil,zipfile,importlib.metadata as im
import pandas as pd
R=Path(__file__).resolve().parent
repo_results=Path(sys.argv[1]) if len(sys.argv)>1 else Path('/workspace/codex/results')
D=repo_results/'GSE189903_HCC_tumor_PGAM5';D.mkdir(parents=True,exist_ok=True)
hashes={}
for name in ['GSE189903_Info.txt.gz','GSE189903_genes.tsv.gz','GSE189903_barcodes.tsv.gz','GSE189903_matrix.mtx.gz','GSE189903_family.soft.gz']:
 p=R/name;h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 hashes[name]={'bytes':p.stat().st_size,'sha256':h.hexdigest()}
(R/'input_sha256.json').write_text(json.dumps(hashes,indent=2))
(R/'environment_versions.txt').write_text('\n'.join(f'{p}=={im.version(p)}' for p in ['numpy','pandas','scipy','anndata','scanpy','statsmodels','matplotlib','nbformat','nbclient','ipykernel','threadpoolctl'])+'\n')
s=json.loads((R/'analysis_summary.json').read_text())
u=pd.read_csv(R/'GSE189903_HCC_tumor_PGAM5_upregulated.csv',index_col=0);new=set(u.gene)
sources={'GSE140228_Droplet':repo_results/'GSE140228_allTumor_PGAM5/Droplet/upregulated.csv',
 'GSE140228_Smartseq2':repo_results/'GSE140228_allTumor_PGAM5/Smartseq2/upregulated.csv',
 'GSE149614':repo_results/'GSE149614_PGAM5_primaryT/GSE149614_primaryT_PGAM5_upregulated.csv',
 'GSE151530':repo_results/'PGAM5_HCC/HCC_PGAM5_FDR005_log2FC1_upregulated.csv'}
overlap_summary={};status=u[['gene','log2FC','FDR','fraction_PGAM5_positive','fraction_PGAM5_undetected']].copy()
for key,p in sources.items():
 t=pd.read_csv(p,index_col=0);old=set(t.gene);common=new&old
 overlap_summary[key]={'count':len(common),'genes':sorted(common),'input_file':str(p.relative_to(repo_results)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
 status[key+'_upregulated']=status.gene.isin(old)
 # Symbol-level intersections; keep all source records if a symbol has multiple IDs.
 old_rows=t[t.gene.isin(common)].copy();old_rows.index.name='source_row_ID';old_rows=old_rows.reset_index()
 joined=u[u.gene.isin(common)].reset_index().merge(old_rows,on='gene',suffixes=('_GSE189903','_'+key),validate='one_to_many')
 joined.to_csv(R/('upregulated_shared_with_'+key+'.csv'),index=False)
status.to_csv(R/'upregulated_previous_dataset_membership.csv')
full=pd.read_csv(R/'GSE189903_HCC_tumor_PGAM5_full_DE.csv',index_col=0)
core=full[full.gene.isin(['OAF','SERINC2','SLC25A4','SQLE','TUBG1'])].copy()
assert len(core)==5
core['passes_requested_up_thresholds']=core.FDR.lt(.05)&core.log2FC.ge(1)&core[['fraction_PGAM5_positive','fraction_PGAM5_undetected']].max(axis=1).ge(.1)
core.to_csv(R/'previous_five_core_candidates_in_GSE189903.csv')
s['descriptive_previous_up_intersections']=overlap_summary
s['previous_five_candidates_passing_up_thresholds']=core[core.passes_requested_up_thresholds].gene.tolist()
(R/'analysis_summary.json').write_text(json.dumps(s,indent=2))
text='''# GSE189903 HCC肿瘤核心及肿瘤边缘：PGAM5⁺巨噬细胞差异分析

按用户明确选择，仅纳入官方诊断Hepatocellular carcinoma且组织为Tumor core或Tumor border的样本；排除邻近非肿瘤组织和ICC。GEO family SOFT的诊断和组织信息通过S_ID映射到细胞Info表，不仅靠样本名推测。使用原始作者TAM标签，不进行新的巨噬细胞聚类。

'''
text+=f"纳入**{s['HCC_tumor_samples']}个样本、{s['HCC_tumor_patients']}位患者**：12个肿瘤核心、4个肿瘤边缘；共有{s['HCC_tumor_cells']:,}个细胞，作者TAM注释{s['raw_TAMs']:,}个，质控后**{s['QC_TAMs']:,}个巨噬细胞**。其中**PGAM5检出{s['PGAM5_positive_cells']}个，未检出{s['PGAM5_undetected_cells']:,}个**，13个样本、4位患者贡献PGAM5检出细胞。\n\n"
text+=f"按既定筛选标准得到**{s['upregulated_excluding_PGAM5']}个上调候选、{s['downregulated_excluding_PGAM5']}个下调候选**（均排除PGAM5）。前10上调候选："+'、'.join(s['top20_up'][:10])+'。\n\n'
text+='''## 方法与筛选口径

- PGAM5原始UMI计数>0定义检出组，计数=0是RNA未检出，不能等同于蛋白阴性。
- QC：检出基因数≥500、线粒体计数比例<20%；每细胞计数标准化至10,000后log1p。
- 合并所有选定样本的QC巨噬细胞，Wilcoxon秩和检验、ties校正、全部25,885基因作BH校正。
- FDR <0.05、|Scanpy近似log2FC| ≥1，至少一组检出率≥10%。
- 按用户要求不进行患者配对、深度匹配、患者/样本/组织或测序深度协变量校正。
- PGAM5保留在full_DE.csv及DE_including_PGAM5.csv；从上、下调候选列表排除。核糖体与线粒体基因保留。
- log2FC为平均log表达回变换后的近似比值，不是算术均值之比；对照接近零时对伪计数敏感。

## 与此前结果的描述性交集

'''
for key,z in overlap_summary.items():text+=f"与{key}上调候选有**{z['count']}个**交集："+'、'.join(z['genes'])+'。\n\n'
text+='''此前的5个共同候选OAF、SERINC2、SLC25A4、SQLE、TUBG1在本次分析中均未达到完整上调筛选标准，具体log2FC、FDR及检出率见previous_five_core_candidates_in_GSE189903.csv。这不等同于证明没有生物学关联，也未证明此前候选错误；只是本次特定范围和筛选条件下不通过。交集按gene symbol匹配，不作别名转换，不把描述性交集当成独立验证。

## 文件与复现

主表：完整差异表、包含PGAM5的筛选表、上调/下调候选、top20上调；另有每样本/患者/组织的计数统计、QC前后细胞元数据、诊断与组织映射、marker检查、此前候选交集及5基因检查。

GSE189903_HCC_tumor_PGAM5_analysis.ipynb已执行，包含结果表与图形；提供PNG/PDF及分析/图形代码。input_sha256.json记录原始输入大小及SHA256，environment_versions.txt记录软件版本，independent_calculation_checks.csv记录独立复核。

复核已检查样本范围、原始barcode/细胞顺序、基因唯一标识、PGAM5原始计数、QC、完整候选集合及BH校正，另用SciPy从原始QC计数重算前10候选的秩检验p值和log2FC。3张PNG均已视觉检查；未进行完整notebook网页交互界面检查。

来源：[GEO GSE189903](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189903)、[官方补充文件](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189903/suppl/)、[family SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189903/soft/GSE189903_family.soft.gz)，PMID 36476645。

将input_sha256.json中的5个官方gz文件放在脚本目录，按environment_versions.txt安装依赖，运行python analyze_HCC_tumor.py，再运行python build_notebook.py。脚本相对于自身目录定位输入。ZIP包含结果表、图形、已执行notebook及代码，不含完整原始矩阵或QC计数h5ad；运行分析会重新生成这些计数输出。package_results.py为导出辅助脚本，需提供仓库已有结果目录以重算与此前数据集的交集。

## 解释与限制

89个上调基因是PGAM5 RNA检出相关候选，尚不是已验证signature。合并细胞FDR不代表患者间可重复性；本队列仅4位患者，多个样本来自同一患者。按要求未调整的患者、区域、巨噬细胞亚型与测序深度差异均可能贡献结果。

巨噬细胞标志C1QA/B/C、CD68、LST1、TYROBP支持原始TAM注释；部分细胞检出ALB等非髓系转录本，未厘清环境RNA、双细胞或吞噬来源，未进行双细胞或环境RNA校正。未进行通路富集、signature打分、患者间验证或独立功能验证。
'''
(R/'结果说明.md').write_text(text)
names=['analysis_summary.json','sample_diagnosis_tissue_mapping.csv','HCC_tumor_celltypes_by_sample.csv','HCC_tumor_TAM_metadata_before_QC.csv.gz','HCC_tumor_TAM_analysis_metadata.csv.gz','PGAM5_by_sample.csv','PGAM5_by_patient_tissue.csv','PGAM5_by_tissue.csv','TAM_marker_QA.csv','independent_calculation_checks.csv','input_sha256.json','environment_versions.txt','analyze_HCC_tumor.py','build_notebook.py','package_results.py','GEO_series_record.txt','结果说明.md','upregulated_previous_dataset_membership.csv','previous_five_core_candidates_in_GSE189903.csv','GSE189903_HCC_tumor_PGAM5_analysis.ipynb']
for f in list(R.glob('GSE189903_HCC_tumor_PGAM5*.csv'))+list(R.glob('*.png'))+list(R.glob('*.pdf'))+list(R.glob('upregulated_shared_with_*.csv')):names.append(f.name)
for name in set(names):shutil.copy2(R/name,D/name)
zp=D/'GSE189903_HCC_tumor_PGAM5_results.zip'
with zipfile.ZipFile(zp,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in sorted(D.iterdir()):
  if f.is_file() and f!=zp:z.write(f,f.name)
with zipfile.ZipFile(zp) as z:assert z.testzip() is None
print(json.dumps({'files':len(list(D.iterdir())),'ZIP_MiB':zp.stat().st_size/1024**2,'previous_intersection_counts':{k:z['count'] for k,z in overlap_summary.items()}},indent=2))
