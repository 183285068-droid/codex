from pathlib import Path
import json,hashlib,shutil,zipfile,importlib.metadata as im
import pandas as pd
R=Path(__file__).resolve().parent
hashes={}
for f in sorted(R.glob('GSE140228*.gz')):
 h=hashlib.sha256()
 with f.open('rb') as handle:
  for block in iter(lambda:handle.read(4*1024*1024),b''):h.update(block)
 hashes[f.name]={'bytes':f.stat().st_size,'sha256':h.hexdigest()}
(R/'input_sha256.json').write_text(json.dumps(hashes,indent=2))
(R/'environment_versions.txt').write_text('\n'.join(f'{pkg}=={im.version(pkg)}' for pkg in ['numpy','pandas','scipy','anndata','scanpy','statsmodels','matplotlib','nbformat','nbclient','ipykernel','threadpoolctl'])+'\n')
s=json.loads((R/'analysis_summary.json').read_text())
u={tech:pd.read_csv(R/tech/'upregulated.csv',index_col=0) for tech in ['Droplet','Smartseq2']}
ids=set(u['Droplet'].index)&set(u['Smartseq2'].index)
rows=[]
for i in sorted(ids):
 row={'Ensembl':i,'gene':u['Droplet'].loc[i,'gene']}
 for tech in u:
  for col in ['log2FC','FDR','fraction_PGAM5_positive','fraction_PGAM5_undetected']:row[f'{tech}_{col}']=u[tech].loc[i,col]
 rows.append(row)
pd.DataFrame(rows).to_csv(R/'shared_upregulated_between_platforms.csv',index=False)
s['descriptive_platform_overlap']={'shared_upregulated_genes':len(ids),'genes':[u['Droplet'].loc[i,'gene'] for i in sorted(ids)],'note':'Descriptive Ensembl overlap, not independent replication; both platforms include donor D20171109.'}
(R/'analysis_summary.json').write_text(json.dumps(s,indent=2,ensure_ascii=False))
text='''# GSE140228 所有 Tumor 样本：PGAM5⁺ 巨噬细胞差异分析

选择官方 cellinfo 的 Tissue == Tumor，包含全部原始 celltype_sub 以 Mφ- 开头的六类巨噬细胞；排除 Mono、DC 等标签。Droplet（UMI）与 Smart-seq2（read counts）分别合并各自所有 Tumor 样本检验。

**Droplet 的 8 个 tumor 样本中，7 个作者标注为 HCC、1 个标注为 CC。因此这是所有 Tumor 分析，不能称为纯 HCC 分析。Smart-seq2 的 10 个 tumor 样本均为 HCC，其中8个有QC巨噬细胞。样本数不等于患者数；D20171109 同时出现在两种技术中。**

| 技术 | Tumor样本 | Tumor患者 | QC巨噬细胞 | PGAM5检出 | PGAM5未检出 | 上调候选 | 下调候选 |
|---|---:|---:|---:|---:|---:|---:|---:|
'''
for tech in u:
 z=s[tech];text+=f"| {tech} | {z['tumor_samples']} | {z['tumor_donors']} | {z['QC_macrophages']} | {z['PGAM5_positive']} | {z['PGAM5_undetected']} | {z['upregulated_excluding_PGAM5']} | {z['downregulated_excluding_PGAM5']} |\n"
text+='''
PGAM5原始计数>0定义检出组，计数=0为RNA未检出，不能等同于蛋白阴性。QC：检出基因数≥500，线粒体计数比例<20%；每细胞标准化至10,000后log1p。合并细胞Wilcoxon检验，校正ties，对全部54,574基因作BH校正。按要求不进行患者配对、深度匹配、样本或测序深度协变量校正。

筛选：FDR <0.05、|Scanpy近似log2FC| ≥1，至少一组检出率≥10%。PGAM5保留在full_DE.csv及DE_including_PGAM5.csv，从候选upregulated.csv和downregulated.csv排除。核糖体和线粒体基因保留。log2FC是平均log表达回变换后的近似比值，不是算术均值之比；对照接近零时具有伪计数敏感性。

每个平台目录包含完整差异表、包含PGAM5的筛选表、上/下调候选、top20上调、每样本统计、QC前后细胞元数据及marker检查。根目录包含已执行notebook、PNG/PDF图形、分析代码、输入SHA256、环境版本和独立计算复核。

'''
text+=f'两平台上调交集共 **{len(ids)}** 个基因，见shared_upregulated_between_platforms.csv。部分患者跨技术重复，该交集不能视为独立队列验证。\n\n'
for tech in u:text+=tech+' 前10上调候选：'+'、'.join(s[tech]['top20_up'][:10])+'。\n\n'
text+='''## 解释与限制

这些是PGAM5 RNA检出相关候选，而非已验证signature。合并细胞FDR不代表患者间可重复性。按要求未调整的测序深度、肿瘤类型、巨噬细胞亚型或患者差异均可能影响结果；marker检查可见部分巨噬细胞检出ALB、CD3D、EPCAM等转录本，未厘清环境RNA、双细胞或吞噬来源，未进行双细胞去除或环境RNA校正。含CC的Droplet结果不能作为HCC特异结论。

复核验证细胞/基因唯一性、metadata与矩阵对齐、PGAM5原始计数、样本选择、QC、完整候选集合及BH校正；另用SciPy从原始QC计数重算各技术前5候选的秩检验p值和log2FC。notebook已执行，PNG图已检查；未进行完整notebook网页交互界面检查。

## 来源与复现

[GEO GSE140228](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140228)，[官方补充文件](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140228/suppl/)，PMID 31675496。

下载input_sha256.json列出的7个官方gz文件至脚本目录，按environment_versions.txt安装依赖，运行python analyze_tumor.py，再运行python build_notebook.py。代码以脚本目录定位输入。Smart-seq2对重复SYMBOL采用原矩阵的SYMBOL_ENSEMBL格式逐行核对，以Ensembl唯一标识基因。

为便于浏览器下载，ZIP包含结果表及复现代码，不含完整原始矩阵或QC计数h5ad。运行分析代码会生成QC计数h5ad及所有结果；已执行notebook保留复核输出及图形。下载ZIP并解压，可用Excel导入CSV查看。
'''
(R/'结果说明.md').write_text(text)
D=Path('/workspace/codex/results/GSE140228_allTumor_PGAM5');D.mkdir(parents=True,exist_ok=True)
for tech in u:
 (D/tech).mkdir(exist_ok=True)
 for f in (R/tech).iterdir():
  if f.is_file() and f.suffix!='.h5ad':shutil.copy2(f,D/tech/f.name)
for name in ['analysis_summary.json','input_sha256.json','environment_versions.txt','analyze_tumor.py','build_notebook.py','package_results.py','independent_calculation_checks.csv','shared_upregulated_between_platforms.csv','结果说明.md','GEO_series_record.txt','GSE140228_allTumor_PGAM5_analysis.ipynb']:
 shutil.copy2(R/name,D/name)
for f in list(R.glob('*.png'))+list(R.glob('*.pdf')):shutil.copy2(f,D/f.name)
zip_path=D/'GSE140228_allTumor_PGAM5_results.zip'
with zipfile.ZipFile(zip_path,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in sorted(D.rglob('*')):
  if f.is_file() and f!=zip_path:z.write(f,f.relative_to(D))
with zipfile.ZipFile(zip_path) as z:assert z.testzip() is None
print(json.dumps({'ZIP_MiB':zip_path.stat().st_size/1024**2,'platform_up_intersection':len(ids),'shared_genes':s['descriptive_platform_overlap']['genes']},ensure_ascii=False,indent=2))
