"""Final report for the complete patient-linked GSE124751 analysis."""
from pathlib import Path
import json
import pandas as pd

ROOT=Path(__file__).resolve().parent


def render(projection,projection_audit):
    s=json.loads((ROOT/'survival_summary.json').read_text())
    audit=json.loads((ROOT/'survival_audit_results.json').read_text())
    assert audit['passed'] and s['survival_analysis_performed']
    stats=pd.read_csv(ROOT/'survival_statistics.csv')
    cox=pd.read_csv(ROOT/'Cox_binary_continuous_age_adjusted.csv')
    sens=pd.read_csv(ROOT/'survival_sensitivity.csv')
    text=f'''# GSE124751：固定v4探索性NNLS系数与OS / 源表DFS生存分析

**196例主分析中，非零系数与零系数组的OS和DFS均未发现统计学显著差异。**
沿用此前固定的850基因×17组分联合参考、原始行权重与NNLS；没有根据随访结果更改基因、系数或分组阈值。
临床表没有PFS，本报告使用原表的OS和DFS（其中事件标签为复发），不将DFS改称PFS。
该参考未通过既往丰度验证，结果是探索性RNA参考系数与随访的关联，不能解释为已验证的细胞丰度效应。

## 主结果

| 终点 | 患者数（High / Low） | 事件数（High / Low） | logrank p | BH q（两个主终点） | Cox HR，High / Low（95%CI） |
| --- | --- | --- | --- | --- | --- |
'''
    for r in stats.itertuples():
        text+=f'| {r.endpoint} | {r.patients}（{r.high_patients} / {r.low_patients}） | {r.events}（{r.high_events} / {r.low_events}） | {r.logrank_p:.6f} | {r.logrank_BH_q:.6f} | {r.HR_high_vs_low:.3f}（{r.HR_lower95:.3f}–{r.HR_upper95:.3f}） |\n'
    text+='''
![OS与DFS Kaplan–Meier曲线，含风险人数、置信区间及检验](KM_OS_DFS.png)

曲线显示95%置信区间、删失标记和各时间点开始时的风险人数。
主分析OS共29个死亡事件；置信区间较宽，结果不显著不能证明不存在关联。
OS两组的KM中位生存均未达到；DFS High组未达到中位数，Low组中位数46.1个月。

## 入组、患者链接与质量核验

原GSE124751包含206份肿瘤表达样本：195常规HCC（C-HCC）和11免疫细胞丰富型HCC（HCC-IS）。
矩阵有20,501个基因，包含PGAM5，覆盖原参考的763/850个基因。
按全部206份表达样本的预设中位数=0分组：High>0（68份），Low=0（138份）。
临床排除后仍保持原分组，主分析为High 61例、Low 135例。
**0是模型系数边界解，不代表组织内没有PGAM5表达巨噬细胞。**

两份用户上传RTF经严格ASCII正文解析后，转出的字节数和SHA256与官方Git LFS对象完全相同：

| 解码表 | 字节数 | 官方SHA256 |
| --- | --- | --- |
| data_clinical_patient.txt | 12050 | 629c2d2994250b55f069ddd4a7a5900f5ac0f3ab1d9da37b96686f2c8e8729c1 |
| data_clinical_sample.txt | 20728 | 5ecf933d8834175841f59efa2f8b9f251531877fd89a58a1da81b8042762ac49 |

临床源表各有231条一对一患者记录，无重复患者ID。实际链接过程为：

- 206份表达中198份的完整肿瘤编号去掉末尾T后，精确对应临床PATIENT_ID和SAMPLE_ID；
  同时在同中心早期GSE54504的`Tumor (完整编号T)`记录中逐一核对。
- 8份HCC-IS在矩阵中以L开头，在GEO标题中以H开头；元数据后缀和顺序相符，
  已保留两种原名，但没有借此推断临床身份，因缺乏完整编号的精确临床匹配而排除。
- 其余198份中196份年龄完全相符。H080196T的GEO年龄57、临床年龄49；
  H080204T的GEO年龄53、临床年龄57。按临床取得前已固定的年龄核验规则，主分析排除这2例。
- 主分析纳入193例C-HCC和3例HCC-IS。没有将8份缺失HCC-IS的随访或系数当作0补入分析。
- OS和DFS时间均为有限正值、状态标签有效，没有DFS时间晚于已知死亡时间的情况。

### 源表字段的限制

OS按原表`OS_MONTHS`（初始诊断后月数）、`0:LIVING / 1:DECEASED`报告。
DFS按原表`DFS_MONTHS`（初始治疗后月数）、`0:DiseaseFree / 1:Recurred`报告。
源表中确有死亡但DFS事件为0的患者，因此没有擅自把死亡补写为DFS事件，
也没有把本终点视作已重新验证的复发或死亡复合终点。所有结论限于源表定义。

H104572T的OS仍存活随访为26.4个月、DFS复发时间为28.3个月，存在终点随访时间不一致。
因为OS记录不是已知死亡时间，主分析保留源表两个终点并标记异常，另做剔除该例的DFS敏感性分析。
没有修改或替换原随访时间。

公开临床sample表的所有231行均写成`LIAD / Hepatocellular Adenoma`，与原HCC研究标题、GEO病理范围不一致。
本轮按原研究及GEO完整编号确定HCC入组，明确记录该字段冲突，没有使用该字段作入组过滤。
源表还存在18个`IIII`等级值和取1–4的`CIRRHOSIS`字段；均保留原值，未擅自转成分期或二元肝硬化状态。
年龄之外没有可靠的性别/TNM字段用于本次调整。

## 敏感性分析

| 分析 | 终点 | 人数 | logrank p | HR（95%CI） |
| --- | --- | --- | --- | --- |
'''
    labels={'all198_exact_IDs_including_age_conflicts':'纳入全部198个精确编号，含2例年龄冲突',
            'exclude_living_followup_discordance':'DFS剔除H104572T随访时间不一致例'}
    for r in sens[sens.analysis.isin(labels)].itertuples():
        text+=f'| {labels[r.analysis]} | {r.endpoint} | {r.patients} | {r.logrank_p:.6f} | {r.HR_high_vs_low:.3f}（{r.HR_lower95:.3f}–{r.HR_upper95:.3f}） |\n'
    text+='''
这些分析均未发现显著的高低分组差异。两种预设备选表达转换未改变分组，OS/DFS logrank结果与主分析相同。

### 补充Cox模型

连续系数按原全部206份表达样本的标准差缩放，未按终点或随访子集重新标准化。
另外提供只调整年龄的Cox模型。以下q对2终点×4模型的8个目标系数检验统一BH校正；
这一补充检验家族与两个主logrank终点的校正分别报告。

| 终点 | 补充模型 | 目标HR（95%CI） | 未校正p | BH q（8个目标检验） |
| --- | --- | --- | --- | --- |
'''
    labels={'continuous_per_SD':'连续系数每增加1 SD','binary_age_adjusted':'高低分组，调整年龄',
            'continuous_per_SD_age_adjusted':'连续系数每增加1 SD，调整年龄'}
    target=cox[cox.covariate.isin(['high','coefficient_per_SD']) & cox.model.isin(labels)]
    for r in target.itertuples():
        text+=f'| {r.endpoint} | {labels[r.model]} | {r.HR:.3f}（{r.HR_lower95:.3f}–{r.HR_upper95:.3f}） | {r.p:.6f} | {r.target_test_BH_q_across_8:.6f} |\n'
    text+=f'''
连续系数与DFS的未校正p约0.0177，调整年龄后约0.0219；但在8个目标检验统一BH校正后q约0.0875，
只能报告为探索性线索，不能用它替代不显著的预设高低组比较或宣称已确立保护作用。
OS年龄调整模型中年龄协变量的比例风险检验p约0.030，提示该补充模型的年龄效应假设需谨慎；
主二元系数的OS/DFS比例风险检验未低于0.05，详见模型表。没有根据结局继续选择其他模型。

## 固定NNLS方法与用途限制

参考和全部三种表达转换的系数在临床数据取得前已发布于Git提交
`5d1a42f2b4448213c74730d05c7489923baf8d4c`；本轮逐字节核验未改变这些文件。
原850×17参考和行权重也与更早的v4分析逐字节一致，缺失87个基因不作插补或重新筛选。
用共同763个基因计算加权NNLS，17组分系数以各自和归一化，目标是`TAM_PGAM5_detected`。

GEO明确发布的是上四分位归一化log2表达，但未明确是否加1。
主分析按预设log2(x+1)假设，以2^x−1反变换；对全20,501基因归一化至每样本总和10,000后投影。
敏感性变换为正值2^x/零保留0，以及全部直接2^x；两者都没有改变分组。
这些不是原始测序读段计数。

**该参考先前未通过丰度验证，且本队列加权相对拟合残差中位数为{projection['primary_median_relative_fit_residual']:.4f}。
不能将系数解释为真实PGAM5⁺巨噬细胞百分比。**
PGAM5不是巨噬细胞特异表达基因；该bulk关联不证明细胞来源特异性、抗肿瘤功能或因果作用。
其他病理因素、组间基线差异、临床误标及随访更新差异仍可能影响结果。

## 独立核验与复现

表达及投影共{len(projection_audit['checks'])}项检查通过：重建{projection_audit['reconstructed_query_values']:,}个查询数值，
28份样本独立BVLS系数最大差{projection_audit['BVLS_max_coefficient_error']:.3g}。
生存另有{len(audit['checks'])}项独立检查通过：官方临床SHA256、逐患者编号和年龄、原冻结系数、风险人数、
完整KM乘积估计、带并列事件的logrank、两个主终点BH，以及独立Efron部分似然的Cox HR、CI和Wald p。
PNG/PDF已导出并检查，ZIP通过CRC和逐文件SHA256核验。数值验证不等于丰度或临床终点重新验证。

在本目录运行：

```bash
python -m pip install -r requirements.txt
export MPLCONFIGDIR=/tmp/pgam5-matplotlib
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
python convert_clinical_rtf.py
python analyze.py
python survival.py
python audit.py
python audit_survival.py
python build_report.py
```

上传RTF和已验证的解码表都已在`inputs/`，所以本轮复现无需修改云环境网络配置。
`download_clinical.py`仅是可选公开源下载脚本；云环境原生下载限制仍存在，已通过提供文件解决本轮数据需求。

主要文件：

- `KM_OS_DFS.png / .pdf`：OS/DFS曲线与风险人数。
- `survival_statistics.csv`、`Cox_binary_continuous_age_adjusted.csv`、`survival_sensitivity.csv`：主比较、补充模型和敏感性结果。
- `patient_coefficients_clinical.csv`、`OS_analysis_patients.csv`、`DFS_analysis_patients.csv`：患者级表达系数及随访。
- `clinical_link_audit.csv`、`endpoint_excluded_patients.csv`：全部206样本的链接、冲突及排除原因。
- `KM_curve_values.csv`、`KM_numbers_at_risk.csv`、`independent_logrank_risk_sets.csv`：可核查的图表与检验数值。
- `clinical_RTF_conversion_audit.json`、`audit_results.json`、`survival_audit_results.json`：来源与数值核验。
- `preclinical_frozen_coefficients.csv`、`input_manifest.csv`、`survival_protocol.json`：冻结系数、输入来源及方法。
- `GSE124751_v4_exploratory_NNLS_results.zip`：完整自包含结果包，GitHub点击“Download raw file”。

官方来源：[GSE124751](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124751)、
[GSE54504](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54504)、
[cBioPortal官方研究目录](https://github.com/cBioPortal/datahub/tree/18763f8eee86d4886d7210068f745f1e509498ee/public/lihc_amc_prv)、
[后续论文数据可用性说明（PMC5960090）](https://pmc.ncbi.nlm.nih.gov/articles/PMC5960090/)。
'''
    return text
