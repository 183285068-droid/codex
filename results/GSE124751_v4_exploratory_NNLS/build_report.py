"""Write the reviewed report and a self-contained package of completed work."""
from pathlib import Path
import hashlib, json, zipfile
import pandas as pd

ROOT = Path(__file__).resolve().parent


def main():
    summary = json.loads((ROOT/'summary.json').read_text())
    audit = json.loads((ROOT/'audit_results.json').read_text())
    assert audit['passed']
    sensitivity = pd.read_csv(ROOT/'transformation_sensitivity.csv')
    data = pd.read_csv(ROOT/'sample_coefficients_and_metadata.csv')
    counts = data.groupby('subtype',sort=True).agg(samples=('GSM','size'),
                nonzero=('TAM_PGAM5_detected',lambda x: (x>0).sum()))
    counts.to_csv(ROOT/'subtype_sample_counts.csv')
    n_checks = len(audit['checks'])
    text = f'''# GSE124751：v4探索性NNLS系数已完成，生存分析待取得临床随访

**本轮已完成206份HCC肿瘤表达样本的固定NNLS投影；尚未完成OS/RFS/PFS生存分析。**
GEO没有逐患者生存时间与事件状态。已找到同中心早期队列的公开cBioPortal临床文件，
但其实际Git LFS内容所在域名被当前云环境网络策略拦截。
这属于访问限制，不能据此判断该队列没有生存随访。当前没有可报告的KM曲线、logrank p值或Cox HR。

## 已完成结果

| 项目 | 结果 |
| --- | --- |
| 纳入表达样本 | 206份肿瘤：195 C-HCC，11 HCC-IS |
| 平台 | GPL11154，Illumina HiSeq 2000 RNA-seq |
| 表达矩阵 | 20,501基因 × 206样本，无缺失表达值 |
| 原参考 | 固定850基因 × 17组分v4联合参考及原行权重 |
| 本平台实际覆盖 | 763/850基因，包含PGAM5；不填补缺失87基因 |
| 目标组分 | TAM_PGAM5_detected |
| 非零目标系数 | 68份 |
| 目标系数为0 | 138份 |
| 固定中位数分组 | 中位数=0；High>0（68份），Low=0（138份） |
| 加权参考秩 / 条件数 | 17 / {summary['weighted_reference_condition']:.4f} |
| 加权相对拟合残差中位数 | {summary['primary_median_relative_fit_residual']:.4f} |

上述分组先在全部206份表达样本中固定，再考虑未来临床终点缺失。
因中位数为0，High/Low实际是非零/零系数的比较；**0不代表组织内没有PGAM5表达巨噬细胞**。

![NNLS分布与转换敏感性](NNLS_distribution_and_sensitivity.png)

### 表达值转换敏感性

GEO明确记录的是RSEM估计后固定上四分位归一化的log2表达值，未明确注明是否加1。
因此主分析预先假设log2(x+1)，以2^x−1转回线性尺度；不将这些数据称为原始读段计数。
对全20,501基因进行每样本和归一化至10,000，再按763个共有基因计算NNLS。
另外固定两种敏感性变换：零值保留为0、正值转为2^x；以及全部直接转为2^x。

| 变换 | 非零系数样本 | 与主分析相比分组变化 | 与主分析的Spearman r | 最大系数绝对差 |
| --- | --- | --- | --- | --- |
'''
    for row in sensitivity.itertuples():
        text += f'| {row.variant} | {row.nonzero} | {row.group_changes_vs_primary} | {row.Spearman_r_vs_primary:.6f} | {row.max_absolute_difference_vs_primary:.8f} |\n'
    text += f'''
变换的选择未利用任何生存结局，也没有改变参考、权重、求解器或筛选额外标志基因。
数值/分组对该转换假设稳定，不构成生物学丰度验证。

## 临床资料线索及身份核对

1. [GSE124751](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124751)有206份表达样本；未纳入另一系列GSE124752。
   195份为常规HCC，11份为免疫细胞丰富型HCC（HCC-IS）。
2. [GSE54504](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54504)为同中心早期可切除HCC队列。
   原完整肿瘤编号精确对应198份本轮表达样本：195 C-HCC和3 HCC-IS。
   例如本轮H060664T（GSM3553191）精确对应早期GSM1317251的Tumor (H060664T)。
   该早期GEO记录中的“early recurrence”仅是状态，不能替代时间和删失信息。
3. 其余8份HCC-IS在本轮表达矩阵中以L开头，但GEO标题为相同数字后缀的H编号。
   全部206列的顺序及这8例后缀相符，故已标注为元数据前缀别名，并保留两个原名。
   **此别名用于整理本轮GEO描述，不用于宣称临床患者身份已经验证。**
4. [原队列论文（PMID 24798001）](https://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=24798001&retmode=xml)
   和[后续论文（PMC5960090）](https://pmc.ncbi.nlm.nih.gov/articles/PMC5960090/)的数据可用性说明指向
   cBioPortal研究`lihc_amc_prv`。
   [官方公开datahub目录](https://github.com/cBioPortal/datahub/tree/{summary['source_datahub_commit']}/public/lihc_amc_prv)
   将其标为“Liver Hepatocellular Carcinoma (AMC, Hepatology 2014)”，231个肿瘤/配对正常患者队列。
   本轮仅需其中公开临床文件，不需原始受控测序数据或患者私密资料。

### 尚未取得的临床文件

| 文件 | 官方Git LFS对象SHA256 | 实际文件字节数 |
| --- | --- | --- |
| data_clinical_patient.txt | 629c2d2994250b55f069ddd4a7a5900f5ac0f3ab1d9da37b96686f2c8e8729c1 | 12050 |
| data_clinical_sample.txt | 5ecf933d8834175841f59efa2f8b9f251531877fd89a58a1da81b8042762ac49 | 20728 |

`inputs/*lfs_pointer.txt`是130字节的官方指针，**不是临床表**。两张真实表尚未读取，
因此尚不确定其具体终点字段、有效随访人数，以及与本轮样本的最终匹配人数。
198个GEO编号对应关系是临床链接机会，不等于198例已有可用随访。

尝试结果保存在`network_and_clinical_status.json`：GEO及GitHub普通文件可访问，
公开临床Git LFS实际内容下载所需的`media.githubusercontent.com`被当前代理以CONNECT HTTP403拦截。
cBioPortal API及其原生Git LFS下载主机也被当前网络策略拦截。
已保存允许`media.githubusercontent.com`的环境配置草稿，保留此前设置；
**保存草稿不等于修改当前网络，需在环境设置中保存并发布后再重试。**

取得真实临床表后的分析要求：核验字段定义和单位，以完整患者编号精确匹配，
并核对年龄等独立信息；审计重复、缺失和HCC-IS覆盖差异。
只有实际记录的生存时间/事件对才可计算KM、logrank和Cox。
原论文出现DSS/RFS等终点时按原定义报告，不改称OS/PFS。
该数据可能与此前TCGA系数的平台表现不同，不能移用TCGA绝对阈值。

## 方法与限制

参考及行权重与此前v4分析逐字节一致。只取实际覆盖的763个参考基因，
按原行尺度对参考和查询共同加权；使用SciPy固定NNLS，并以17组分系数之和归一化。
不能根据本队列结局筛选参考基因、改变权重或选择最佳cutpoint。

**该参考先前未通过目标缺失混合物等丰度验证，当前系数不能解释为已验证细胞比例。**
主分析相对加权拟合残差中位数约79.9%，提示该参考对这些bulk表达的整体解释有限。
PGAM5本身并非巨噬细胞特异基因，测到PGAM5也不建立bulk中的细胞来源特异性。
本轮变换稳定和数值审计通过仅说明计算可复现。

## 本轮核验

{n_checks}项检查通过。直接由原始发布的log2矩阵重新构建{audit['reconstructed_query_values']:,}个查询数值；
对{audit['independent_BVLS_samples']}份样本用独立BVLS有界最小二乘求解复核，
17组分系数最大差为{audit['BVLS_max_coefficient_error']:.3g}。
另外核查输入SHA256、行列顺序、基因覆盖、样本与别名数量、完整编号链接及全局中位数分组。
PNG/PDF已导出并检查。该审计不涵盖尚未取得的临床资料或生存分析。

## 文件及复现

- `sample_coefficients_and_metadata.csv`：206样本主分析目标系数、GSM、原矩阵名、GEO标题、前缀别名标记、年龄、亚型及预设分组。
- `coefficients_*.csv`：三种变换的全部17组分系数。
- `transformation_sensitivity.csv`、`fit_diagnostics_*.csv`：敏感性及拟合诊断。
- `earlier_GEO_exact_title_link_candidates.csv`：严格保持原矩阵编号的早期GEO对应候选，未填入未观察到的临床终点。
- `input_manifest.csv`、`analysis_protocol.json`：不可变输入来源、校验和及分析规则。
- `audit_results.json`：独立核验记录。
- `inputs/`：公共原表达、官方元数据、固定参考及临床LFS指针。论文正文没有重新分发。

在本目录用Python 3.12运行：

```bash
python -m pip install -r requirements.txt
export MPLCONFIGDIR=/tmp/pgam5-matplotlib
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
python analyze.py
python audit.py
python build_report.py
```

以上完整复现本轮已完成的NNLS与核验。
网络允许后，运行`python download_clinical.py`取得并核验两张公开临床真实表；
该脚本只下载源文件，不自动推测终点字段或生成尚未核验的KM分析。

结果包：`GSE124751_v4_exploratory_NNLS_results.zip`。打开GitHub文件后点击“Download raw file”。
'''
    if (ROOT/'survival_audit_results.json').exists():
        from write_survival_report import render
        text=render(summary,audit)
    (ROOT/'README.md').write_text(text)
    files = [p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts
             and p.suffix != '.zip' and p.name not in {'package_file_sha256.csv','zip_sha256.txt'}]
    inventory = pd.DataFrame([{'file':str(p.relative_to(ROOT)), 'bytes':p.stat().st_size,
                              'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(files)])
    inventory.to_csv(ROOT/'package_file_sha256.csv',index=False)
    files.append(ROOT/'package_file_sha256.csv')
    archive = ROOT/'GSE124751_v4_exploratory_NNLS_results.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(files):
            z.write(p,ROOT.name+'/'+str(p.relative_to(ROOT)),compress_type=
                    zipfile.ZIP_STORED if p.suffix=='.gz' else zipfile.ZIP_DEFLATED)
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for row in inventory.itertuples():
            assert hashlib.sha256(z.read(ROOT.name+'/'+row.file)).hexdigest()==row.sha256
    assert archive.stat().st_size < 30*1024*1024
    (ROOT/'zip_sha256.txt').write_text(hashlib.sha256(archive.read_bytes()).hexdigest()+'  '+archive.name+'\n')
    print(f'Package verified: {len(files)} files, {archive.stat().st_size / 1024**2:.2f} MiB')


if __name__ == '__main__':
    main()
