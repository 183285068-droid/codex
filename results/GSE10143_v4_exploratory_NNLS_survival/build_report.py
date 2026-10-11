"""Write the result report, pinned dependencies and checksum-verified download package."""
from pathlib import Path
import json, hashlib, zipfile
import pandas as pd

R = Path(__file__).resolve().parent


def main():
    preparation = json.loads((R / 'preparation_summary.json').read_text())
    audit = json.loads((R / 'audit_summary.json').read_text())
    stats = pd.read_csv(R / 'survival_statistics.csv').set_index('endpoint')
    cox = pd.read_csv(R / 'Cox_binary_and_continuous.csv')
    versions = json.loads((R / 'environment_versions.json').read_text())
    (R / 'requirements.txt').write_text('\n'.join(f'{name}=={version}' for name, version in versions.items()
                                                if name != 'python') + '\n')
    rows = []
    for endpoint in ['OS', 'TTR']:
        s = stats.loc[endpoint]
        rows.append(f'| {endpoint} | {int(s.patients)} | {int(s.high_patients)} / {int(s.low_patients)} | '
                    f'{int(s.events)} | {s.logrank_p:.6f} | {s.logrank_BH_q:.6f} | '
                    f'{s.HR_high_vs_low:.3f} ({s.HR_lower95:.3f}–{s.HR_upper95:.3f}) |')
    continuous_rows = []
    for s in cox[cox.model.eq('continuous_per_SD')].itertuples():
        continuous_rows.append(f'| {s.endpoint} | {s.HR:.3f} ({s.HR_lower95:.3f}–{s.HR_upper95:.3f}) | {s.p:.6f} | {s.PH_rank_p:.6f} |')
    report = f'''# GSE10143：固定v4探索性NNLS系数与生存分析

沿用此前确认的**850基因×17组分联合参考和固定加权NNLS**，以GSE10143中的**80份HCC肿瘤组织（GSE10141）**为主分析。在当前探索性投影下，高、低系数组之间的OS及复发时间差异均未达到统计学显著。该结果不能解释为经过验证的PGAM5⁺巨噬细胞真实浸润丰度。

## 生存分析结果

按全部80例肿瘤样本的系数中位数 **{preparation['median_cutoff']:.10f}** 分组：High > 中位数，Low ≤ 中位数，各40例。该分组在生存终点筛选前固定，未搜索最佳p值切点；复发分析沿用同一分组。

| 终点 | 有效患者 | High / Low | 事件数 | log-rank p | BH q（两个终点） | Cox HR High/Low（95% CI） |
|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

- OS：死亡事件High 11例、Low 21例。HR的置信区间跨1，不能据此认定高系数组有保护作用。
- TTR：复发事件High 24例、Low 25例。终点是GEO原字段“time to HCC recurrence”，**不将其改称PFS或DFS**。按原始复发事件状态计算KM；未进行竞争风险累积发生率分析。
- 二元Cox模型的PH秩检验：OS p={stats.loc['OS', 'PH_rank_p']:.6f}；TTR p={stats.loc['TTR', 'PH_rank_p']:.6f}，未检测到显著违反比例风险假设。
- 本轮未做临床变量校正，不能声称该系数是独立预后因素；官方随访资料没有年龄或分期，其他病理字段存在缺失。

以全部80例系数标准差为单位的连续变量Cox检验同样未显著：

| 终点 | 每增加1 SD的HR（95% CI） | Cox p | PH秩检验p |
|---|---:|---:|---:|
{chr(10).join(continuous_rows)}

![OS及TTR Kaplan–Meier曲线](KM_OS_TTR.png)

图中含95%置信区间、删失标记、风险人数表、log-rank p、BH q及二元Cox HR。

## 组织范围、随访来源与患者匹配

[GSE10143](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10143)是SuperSeries，并非387位独立HCC肿瘤患者：

| 子集 | 组织及样本数 | 本轮用途 |
|---|---|---|
| GSE10141 | 80份HCC肿瘤组织 | 唯一表达输入；官方OS随访 |
| GSE10140 | 82份肝炎/肝硬化背景肝组织 | 仅提取同一患者的复发随访，不混入表达输入 |
| GSE10142 | 225份非肿瘤肝组织 | 不纳入；逐患者生存资料仍在GEO标为未公开 |

总系列旧说明写有“Clinical data has been withheld …”。但GSE10141后来公开了[官方临床表](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10141/suppl/GSE10141_clinical_data_addition_081611.xls.gz)，且GSE10141与GSE10140的逐样本矩阵也含随访。核查子集后，本轮获得足够数据完成分析，并非用论文汇总数字推算个人随访。

OS以官方XLS逐GSM一对一匹配：80例时间均>0，事件仅0/1；全部时间和事件与肿瘤GEO矩阵一致。事件1为死亡，0为存活/删失；时间单位为**天**。

TTR按官方样本名称中的完整数字患者编号连接 `hcc_002` ↔ `bkg_002` 等。80位肿瘤患者中，**62位**能一对一匹配；每一对的OS时间、OS事件均完全一致，并独立核对XLS中的患者编号。其余18位无匹配复发随访，排除TTR分析且不插补。62位匹配患者的复发时间均>0、事件均0/1，因此全部62位有效。

患者027、083的复发删失时间晚于OS删失时间。两者的OS与复发记录均为删失，可能反映不同随访截止；按源数据保留并显式标记，没有将其更改为复发事件。背景肝子集中患者130的复发时间为0，但其不在本轮80位肿瘤患者的62位匹配对象中，未用于生存计算。逐患者链接与来源字段见 `clinical_link_audit.csv`，排除明细见 `endpoint_excluded_patients.csv`。

## 平台适配与方法

平台为**GPL5474，Illumina DASL约6,000基因面板**。官方处理为红/绿探针信号求和后按DAP进行cubic-spline归一化；矩阵是**线性强度**，本轮不进行2的幂变换。

1. 按官方Symbol精确匹配，排除空白/歧义符号；同一基因多个探针取线性强度中位数。6144个探针对应6100个基因。
2. 以全部映射基因的强度总和对每个样本做标量归一化（×10000），仅用于计算。这不是测序CP10K/TPM，也不校准细胞比例单位。
3. 原850基因中**416个可测（48.94%）**；缺测基因在查询和参考中同步删除，不补零、不重新筛选基因或训练参考。保留原17组分、原行权重。
4. 对每个样本求解 `min β≥0 || diag(1/row_scale)(Reference β − Query) ||²`，再以 `β/sum(β)` 作相对系数；目标列为 `TAM_PGAM5_detected`。
5. 80例中78例目标系数>0、2例=0；系数范围 {preparation['target_min']:.6f}–{preparation['target_max']:.6f}。中位数>0，主比较为40/40的中位数组，而非非零/零比较。
6. 无结局驱动特征选择或切点搜索；两侧log-rank检验，OS/TTR两项p进行BH校正，另报告未经校正的二元及连续Cox模型。

原参考和权重逐字节复制自此前已确认版本；来源与所有输入SHA256见 `source_provenance.json`。加权参考数值秩17，条件数 {preparation['weighted_reference_condition']:.4f}；中位相对拟合残差 **{preparation['median_relative_fit_residual']:.4f}**。数值可求解不代表生物学/平台匹配可靠。

## 解释边界

**GPL5474官方注释没有PGAM5符号、Entrez 192111或RefSeq NM_138575。** 这是芯片未覆盖PGAM5的证据，不表示这些组织内PGAM5没有表达。

本轮投影仅使用另外416个可测参考基因，缺失超过一半参考特征；原850×17参考也未通过独立混合物丰度验证。芯片强度与单细胞RNA参考存在平台差异，拟合残差较高，故这些系数只能作为探索性参考投影指标，**不能把0.05866换算解释为5.866%的PGAM5⁺巨噬细胞**，也不能凭此次阴性生存结果排除PGAM5相关状态的生物学作用。

## 验证与复现

`audit.py`已通过**{audit['checks_passed']}项核验**，失败0项，包括全部33,280个查询值从官方线性矩阵重建、80例OS官方XLS/矩阵核对、62例同患者复发链接核对、原参考/权重固定性、结局盲法分组、独立log-rank风险集、KM阶梯值及风险人数表。24例样本另外使用BVLS复核17个组分，最大系数差 {audit['independent_solver_max_difference']:.3g}。这些验证核查计算和来源，不验证真实细胞丰度。

Python {versions['python']}；依赖固定在 `requirements.txt`。结果包中已含所需公开输入，无需再次联网下载。在解压目录运行：

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python analyze.py
python audit.py
python build_report.py
```

## 文件

- `KM_OS_TTR.png` / `.pdf`：可分享及导出的KM图。
- `survival_statistics.csv`：log-rank、BH q及二元Cox结果。
- `Cox_binary_and_continuous.csv`：二元及连续Cox、PH检验。
- `patient_coefficients_clinical.csv`：所有80位患者系数、固定分组及原始临床字段。
- `all_component_coefficients.csv`：17组分系数；`fit_diagnostics.csv`：拟合残差。
- `OS_analysis_patients.csv` / `TTR_analysis_patients.csv`：实际用于每一终点的患者。
- `clinical_link_audit.csv` / `endpoint_excluded_patients.csv`：患者匹配与排除记录。
- `reference_gene_coverage.csv` / `used_probe_mapping.csv`：参考覆盖与官方探针映射。
- `KM_curve_values.csv` / `KM_numbers_at_risk.csv`：图形底层数据。
- `protocol.json` / `source_provenance.json` / `audit_summary.json`：固定方法、来源及核验。
- `inputs/`：公开原始矩阵、官方临床XLS、平台注释与固定参考。
- `SHA256SUMS`：包内每个文件的校验值。

[下载完整结果包](GSE10143_v4_exploratory_NNLS_survival_results.zip)：在GitHub打开后点击 **Download raw file**。
'''
    (R / 'README.md').write_text(report)
    excluded = {'SHA256SUMS', 'GSE10143_v4_exploratory_NNLS_survival_results.zip'}
    files = [f for f in sorted(R.rglob('*')) if f.is_file() and f.name not in excluded
             and '__pycache__' not in f.parts and '.venv' not in f.parts]
    manifest = '\n'.join(f'{hashlib.sha256(f.read_bytes()).hexdigest()}  {f.relative_to(R).as_posix()}' for f in files) + '\n'
    (R / 'SHA256SUMS').write_text(manifest)
    package = R / 'GSE10143_v4_exploratory_NNLS_survival_results.zip'
    with zipfile.ZipFile(package, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for f in files + [R / 'SHA256SUMS']:
            z.write(f, f.relative_to(R).as_posix())
    with zipfile.ZipFile(package) as z:
        assert z.testzip() is None
        for line in manifest.splitlines():
            expected, name = line.split('  ', 1)
            assert hashlib.sha256(z.read(name)).hexdigest() == expected
    print(json.dumps({'zip_MiB': package.stat().st_size / 1024**2, 'file_count': len(files) + 1,
                      'archive_CRC_and_all_SHA256_verified': True}, indent=2))


if __name__ == '__main__':
    main()
