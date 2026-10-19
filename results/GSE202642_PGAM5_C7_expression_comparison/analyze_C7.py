"""Sample-aware and exploratory C7 comparisons using archived GSE202642 values."""
from pathlib import Path
from itertools import product
import hashlib
import json
import os
import importlib.metadata as metadata

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
INPUT = Path(os.environ.get('GSE202642_SOURCE_DIR', ROOT.parent/'GSE202642_PGAM5_subpopulation_assessment'))
OUT = ROOT/'outputs'
OUT.mkdir(exist_ok=True)
ALL = list(range(15))
CORE = [0, 2, 3, 4, 6, 7, 9, 11, 12, 13]
UNCERTAIN = [1, 5, 8, 10, 14]
TARGET = 7
MIN_CELLS = 10
MIN_PATIENTS = 5
SOURCE_REF = '7a1ac2a0940b26df8fcfb5309eb47689dd20c145'
SOURCE_DIR = 'results/GSE202642_PGAM5_subpopulation_assessment'

def stars(q):
    if not np.isfinite(q): return 'NT'
    return '****' if q < .0001 else '***' if q < .001 else '**' if q < .01 else '*' if q < .05 else 'ns'

def bh(table):
    table['q_BH_14'] = stats.false_discovery_control(table.p_value.fillna(1).to_numpy(), method='bh')
    table.loc[table.p_value.isna(), 'q_BH_14'] = np.nan
    table['significance'] = table.q_BH_14.map(stars)
    return table

def exact_swap(delta):
    x = np.asarray(delta, float)
    null = np.asarray(list(product([-1., 1.], repeat=len(x)))) @ x / len(x)
    return float((np.abs(null) >= abs(x.mean()) - 1e-12).mean())

def paired(aggregated, endpoint='mean_log1p_CP10k', minimum_cells=MIN_CELLS):
    wide = aggregated[aggregated.cells >= minimum_cells].pivot(index='patient_id', columns='cluster', values=endpoint).reindex(columns=ALL)
    rows=[]
    for c in ALL:
        if c == TARGET: continue
        p = wide[[TARGET,c]].dropna()
        delta = p[TARGET] - p[c]
        n = len(p)
        rows.append(dict(target='C7', comparator=f'C{c}', comparator_number=c,
            comparator_scope='macrophage_related' if c in CORE else 'unresolved_mixed_RNA',
            endpoint=endpoint, min_cells_per_patient_cluster=minimum_cells,
            paired_patients=n,
            mean_difference_C7_minus_comparator=float(delta.mean()) if n else np.nan,
            median_difference_C7_minus_comparator=float(delta.median()) if n else np.nan,
            p_value=exact_swap(delta.to_numpy()) if n>=MIN_PATIENTS else np.nan,
            status='tested' if n>=MIN_PATIENTS else 'insufficient_paired_patients'))
    return bh(pd.DataFrame(rows))

def save(fig, stem):
    for ext in ['png','svg','pdf']:
        fig.savefig(OUT/f'{stem}.{ext}', dpi=190, facecolor='white')
    plt.close(fig)

def bracket_plot(summary, tests, groups, names, stem):
    all_groups = len(groups)==15
    fig, ax = plt.subplots(figsize=(17 if all_groups else 14, 10 if all_groups else 8.7))
    x = np.arange(len(groups)); target_x = groups.index(TARGET)
    color = ['#c28732' if c==TARGET else '#aeb3b8' if c in UNCERTAIN else '#386b8d' for c in groups]
    bars = ax.bar(x, summary.loc[groups,'mean_log1p_CP10k'], width=.66, color=color)
    for c, bar in zip(groups,bars):
        if c in UNCERTAIN: bar.set_hatch('//');bar.set_edgecolor('#747a80')
        value = summary.loc[c,'mean_log1p_CP10k']
        ax.text(bar.get_x()+bar.get_width()/2,value+.002,f'{value:.4f}',ha='center',fontsize=9.5)
    comps=[c for c in groups if c!=TARGET]
    comps.sort(key=lambda c:(abs(groups.index(c)-target_x),c))
    base = float(summary.loc[groups,'mean_log1p_CP10k'].max())+.016
    for level,c in enumerate(comps):
        r=tests.set_index('comparator_number').loc[c]
        left,right=sorted([target_x,groups.index(c)])
        y=base+level*.014
        ax.plot([left,left,right,right],[y,y+.0025,y+.0025,y],color='#293541',lw=.8)
        ax.text((left+right)/2,y+.0037,f'{stars(r.q_BH_14)}  q={r.q_BH_14:.2g}',
            ha='center',va='bottom',fontsize=9,bbox={'facecolor':'white','edgecolor':'none','pad':.3})
    labels=[]
    for c in groups:
        identifier=f'C{c}'+('?' if c in UNCERTAIN else '')
        if all_groups: labels.append(f'{identifier}\nn={summary.loc[c,"cells"]}')
        else:
            genes=names.loc[c,'marker_gene_1']+'/'+names.loc[c,'marker_gene_2']
            labels.append(f'{identifier}\n{genes}\nn={summary.loc[c,"cells"]}')
    ax.set_xticks(x,labels,fontsize=9)
    ax.set_ylim(0,base+(len(comps)-1)*.014+.023)
    ax.set_ylabel('Mean PGAM5 log1p(CP10k), including zero-expression cells',fontsize=11)
    ax.set_xlabel('Original RNA cluster | cell count',fontsize=11)
    ax.grid(axis='y',color='#e6e8eb');ax.set_axisbelow(True)
    title='All 15 original RNA clusters (unresolved groups flagged)' if all_groups else 'Macrophage-related RNA clusters'
    ax.set_title(f'GSE202642 | PGAM5 expression abundance\n{title} | C7 vs each other cluster',fontsize=16,pad=18)
    footer='Two-sided cell-level Mann-Whitney U; BH-adjusted q across all 14 planned C7 comparisons.\n'
    footer+='* q<0.05; ** q<0.01; *** q<0.001; **** q<0.0001; ns: q>=0.05.\n'
    if all_groups: footer+='? / hatched grey bars: unresolved mixed-RNA groups (C1, C5, C8, C10, C14).\n'
    else: footer+='C0 and C2 retain macrophage-related annotations with hepatic/tumor RNA caveats.\n'
    footer+='Exploratory cell-level significance; patient-level evidence is reported separately. Bar heights match the source UMAP metric.'
    fig.text(.08,.035,footer,fontsize=9.2,va='bottom',linespacing=1.45)
    fig.tight_layout(rect=[.015,.18,.995,.99])
    save(fig,stem)

def main():
    source=INPUT/'PGAM5_UMAP_expression_values.csv.gz'
    assert hashlib.sha256(source.read_bytes()).hexdigest() == '9b17deb293e7a1e220b0b9dfb9c8ca6e6729f70a1383d7789f8d79a402db14da'
    d=pd.read_csv(source)
    names=pd.read_csv(INPUT/'functional_subgroup_names.csv').set_index('primary_cluster')
    old=pd.read_csv(INPUT/'PGAM5_expression_by_cluster.csv').set_index('cluster').loc[ALL]
    assert len(d)==12135 and d.cell_id.is_unique and d.patient_id.nunique()==7
    assert set(d.primary_cluster)==set(ALL)
    assert not d[['cell_id','patient_id','primary_cluster','PGAM5_counts','total_counts']].isna().any().any()
    assert (d.total_counts>0).all() and (d.PGAM5_counts>=0).all() and (d.PGAM5_counts<=d.total_counts).all()
    d['CP10k']=d.PGAM5_counts/d.total_counts*10000
    d['log1p_CP10k']=np.log1p(d.CP10k)
    d['detected']=d.PGAM5_counts>0
    assert int(d.detected.sum())==868
    np.testing.assert_allclose(d.log1p_CP10k,d.PGAM5_log1p_CP10k,atol=1e-12,rtol=1e-12)
    np.testing.assert_allclose(d.CP10k,d.PGAM5_CP10k_recomputed,atol=1e-12,rtol=1e-12)
    summary=d.groupby('primary_cluster').agg(cells=('cell_id','size'),detected=('detected','sum'),
        patient_count=('patient_id','nunique'),mean_CP10k=('CP10k','mean'),mean_log1p_CP10k=('log1p_CP10k','mean'),
        median_CP10k=('CP10k','median')).reindex(ALL)
    summary.index.name='cluster_number'
    summary['cluster']=[f'C{c}' for c in ALL]
    summary['detection_percent']=summary.detected/summary.cells*100
    summary['functional_name_CN']=names.functional_name_CN
    summary['functional_name_EN']=names.functional_name_EN
    summary['mixed_RNA_flag']=names.mixed_RNA_flag
    summary['display_scope']=['unresolved_mixed_RNA' if c in UNCERTAIN else 'macrophage_related' for c in ALL]
    np.testing.assert_array_equal(summary.cells,old.cells)
    np.testing.assert_array_equal(summary.detected,old.PGAM5_detected_cells)
    np.testing.assert_allclose(summary.mean_CP10k,old.mean_PGAM5_CP10k,atol=1e-12,rtol=1e-12)
    np.testing.assert_allclose(summary.mean_log1p_CP10k,old.mean_PGAM5_log1p_CP10k,atol=1e-12,rtol=1e-12)
    np.testing.assert_allclose(d.primary_cluster.map(summary.mean_log1p_CP10k),d.cluster_mean_PGAM5_log1p_CP10k,atol=1e-12,rtol=1e-12)
    summary.to_csv(OUT/'cluster_expression_summary.csv')
    agg=d.groupby(['patient_id','primary_cluster']).agg(cells=('cell_id','size'),detected=('detected','sum'),
        mean_log1p_CP10k=('log1p_CP10k','mean'),mean_CP10k=('CP10k','mean'),
        PGAM5_counts=('PGAM5_counts','sum'),library_counts=('total_counts','sum')).reset_index().rename(columns={'primary_cluster':'cluster'})
    agg['detection_percent']=agg.detected/agg.cells*100
    agg['pseudobulk_CP10k']=agg.PGAM5_counts/agg.library_counts*10000
    agg['eligible_primary']=agg.cells>=MIN_CELLS
    # Patient/cell identifiers and patient-level detail tables are not republished.
    patient=paired(agg)
    patient.to_csv(OUT/'C7_vs_others_patient_primary.csv',index=False)
    pd.concat([paired(agg,endpoint,minimum) for endpoint in ['mean_log1p_CP10k','mean_CP10k','pseudobulk_CP10k']
        for minimum in [1,10,20]],ignore_index=True).to_csv(OUT/'C7_patient_sensitivity.csv',index=False)
    paired(agg,'detection_percent').to_csv(OUT/'C7_patient_detection_secondary.csv',index=False)
    rows=[];fisher_rows=[]
    target=d[d.primary_cluster.eq(TARGET)]
    for c in ALL:
        if c==TARGET: continue
        other=d[d.primary_cluster.eq(c)]
        result=stats.mannwhitneyu(target.log1p_CP10k,other.log1p_CP10k,alternative='two-sided',method='asymptotic',use_continuity=True)
        rows.append(dict(target='C7',comparator=f'C{c}',comparator_number=c,
            comparator_scope='macrophage_related' if c in CORE else 'unresolved_mixed_RNA',
            cells_C7=len(target),cells_comparator=len(other),
            mean_log1p_difference_C7_minus_comparator=float(target.log1p_CP10k.mean()-other.log1p_CP10k.mean()),
            mean_CP10k_difference_C7_minus_comparator=float(target.CP10k.mean()-other.CP10k.mean()),
            U=float(result.statistic),rank_biserial_C7_minus_comparator=float(2*result.statistic/(len(target)*len(other))-1),p_value=float(result.pvalue)))
        f=stats.fisher_exact([[target.detected.sum(),(~target.detected).sum()],[other.detected.sum(),(~other.detected).sum()]])
        fisher_rows.append(dict(target='C7',comparator=f'C{c}',odds_ratio=float(f.statistic),p_value=float(f.pvalue)))
    cell=bh(pd.DataFrame(rows));cell.to_csv(OUT/'C7_vs_others_cell_exploratory.csv',index=False)
    bh(pd.DataFrame(fisher_rows)).to_csv(OUT/'C7_detection_Fisher_exploratory.csv',index=False)
    # Independent checks of the custom exact test and plotted multiple-testing labels.
    assert exact_swap(np.zeros(5))==1 and np.isclose(exact_swap(np.ones(5)),.0625)
    w=agg[agg.cells>=MIN_CELLS].pivot(index='patient_id',columns='cluster',values='mean_log1p_CP10k')
    pair=w[[7,3]].dropna()
    scipy_perm=stats.permutation_test((pair[7].to_numpy(),pair[3].to_numpy()),lambda a,b:np.mean(a-b),
        permutation_type='samples',n_resamples=np.inf,alternative='two-sided')
    np.testing.assert_allclose(exact_swap((pair[7]-pair[3]).to_numpy()),scipy_perm.pvalue)
    order=np.argsort(cell.p_value.to_numpy());p=cell.p_value.to_numpy()[order]
    q=np.minimum.accumulate((p*14/np.arange(1,15))[::-1])[::-1].clip(0,1)
    np.testing.assert_allclose(q,cell.q_BH_14.to_numpy()[order])
    assert len(cell)==14 and set(cell.comparator_number)==set(ALL)-{TARGET}
    assert len(set(CORE)-{TARGET})==9
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
        'axes.spines.right':False,'axes.titleweight':'bold','svg.fonttype':'none','pdf.fonttype':42})
    bracket_plot(summary,cell,CORE,names,'01_C7_macrophage_barplot_significance')
    bracket_plot(summary,cell,ALL,names,'02_C7_all15_RNA_groups_barplot_significance')
    fig,axs=plt.subplots(1,2,figsize=(15,7.8),gridspec_kw={'width_ratios':[1.35,1]})
    rng=np.random.default_rng(202642)
    for i,c in enumerate(CORE):
        v=agg[agg.cluster.eq(c)];good=v[v.eligible_primary];low=v[~v.eligible_primary]
        color='#c28732' if c==7 else '#386b8d'
        if len(good):
            bp=axs[0].boxplot([good.mean_log1p_CP10k],positions=[i],widths=.5,manage_ticks=False,patch_artist=True,showfliers=False,
                medianprops={'color':'#17212b','lw':1.5})
            bp['boxes'][0].set(facecolor=color,alpha=.25)
            axs[0].scatter(i+rng.uniform(-.14,.14,len(good)),good.mean_log1p_CP10k,s=30,c=color,edgecolor='white',lw=.4,zorder=3)
        if len(low):axs[0].scatter(i+rng.uniform(-.14,.14,len(low)),low.mean_log1p_CP10k,s=28,facecolor='none',edgecolor='#a3a3a3',zorder=2)
    axs[0].set_xticks(range(len(CORE)),[f'C{c}' for c in CORE])
    axs[0].set_ylabel('Patient mean PGAM5 log1p(CP10k), zeros included')
    axs[0].set_xlabel('Macrophage-related RNA cluster');axs[0].set_ylim(bottom=0)
    axs[0].grid(axis='y',color='#e8e8e8');axs[0].set_axisbelow(True)
    axs[0].set_title('A  One point per patient-cluster',loc='left',fontsize=12)
    axs[1].axis('off');axs[1].set_title('B  Patient-paired C7 comparisons',loc='left',fontsize=12)
    tr=patient[patient.comparator_number.isin(CORE)].copy()
    tabledata=[]
    for r in tr.itertuples():
        tabledata.append([f'C7 vs {r.comparator}',str(r.paired_patients),f'{r.q_BH_14:.3g}' if np.isfinite(r.q_BH_14) else '-',r.significance])
    table=axs[1].table(cellText=tabledata,colLabels=['Comparison','Patients','BH q','Result'],cellLoc='center',bbox=[0,.13,1,.82])
    table.auto_set_font_size(False);table.set_fontsize(10)
    for (row,col),cell_obj in table.get_celld().items():
        cell_obj.set_edgecolor('#e1e4e8');cell_obj.set_facecolor('#eaf0f5' if row==0 else 'white')
    fig.suptitle('GSE202642 | patient-aware PGAM5 comparison',fontsize=16,y=.98)
    fig.text(.07,.045,'Filled points: >=10 cells per patient-cluster; hollow points: fewer than 10 cells (descriptive only).\n'
        'Two-sided exact paired label-swap test; >=5 paired patients required; BH over 14 planned C7 comparisons.\n'
        'NT: insufficient paired patients. Absent clusters are not assigned zero expression.',fontsize=9.5,linespacing=1.4)
    fig.tight_layout(rect=[.02,.18,.99,.94]);save(fig,'03_C7_patient_aware_comparison')
    tested=patient[patient.status.eq('tested')]
    coretests=cell[cell.comparator_number.isin(CORE)]
    verification=dict(cells=len(d),patients=d.patient_id.nunique(),clusters=15,PGAM5_detected_cells=int(d.detected.sum()),
        raw_counts_and_library_totals_recomputed=True,all_previous_cluster_statistics_match=True,source_UMAP_mean_metric_preserved=True,
        primary_macrophage_related_groups=[f'C{c}' for c in CORE],unresolved_mixed_groups=[f'C{c}' for c in UNCERTAIN],
        C7_cell_comparisons_planned=14,C7_primary_bar_brackets=9,cell_significant_macrophage_comparisons=int((coretests.q_BH_14<.05).sum()),
        cell_significant_all15_comparisons=int((cell.q_BH_14<.05).sum()),patient_testable_comparisons=len(tested),
        patient_significant_comparisons=int((tested.q_BH_14<.05).sum()),patient_minimum_q=float(tested.q_BH_14.min()) if len(tested) else None,
        exact_test_independently_checked=True,BH_independently_checked=True,
        packages={p:metadata.version(p) for p in ['numpy','pandas','scipy','matplotlib']})
    (OUT/'verification.json').write_text(json.dumps(verification,indent=2)+'\n')
    input_hashes={name:hashlib.sha256((INPUT/name).read_bytes()).hexdigest() for name in
        ['PGAM5_UMAP_expression_values.csv.gz','functional_subgroup_names.csv','PGAM5_expression_by_cluster.csv']}
    (ROOT/'source_provenance.json').write_text(json.dumps(dict(repository='183285068-droid/codex',source_commit=SOURCE_REF,
        source_directory=SOURCE_DIR,reference_figure='PGAM5_cluster_mean_expression_UMAP.pdf',
        source_cell_file_git_sha='96ddb713febde32ba59c359de854361fa161020e',
        input_sha256=input_hashes,patient_identity='Seven sample labels treated as seven independent patients, following the archived user-confirmed mapping.',
        analysis_date='2026-10-19',all_zero_values_retained=True,no_reclustering=True,no_imputation=True),ensure_ascii=False,indent=2)+'\n')
    rows='\n'.join(f'| C7 vs {r.comparator} | {r.p_value:.4g} | {r.q_BH_14:.4g} | {r.significance} |' for r in coretests.itertuples())
    readme=f'''# GSE202642：PGAM5表达丰度与C7定向比较

基于原始 `PGAM5_cluster_mean_expression_UMAP.pdf` 对应的12,135个细胞、7位患者、C0–C14原始分群；未重新聚类或更改细胞成员。

## 图和统计表

- [主柱状图：C7与其他巨噬细胞相关群](outputs/01_C7_macrophage_barplot_significance.png) · [PDF](outputs/01_C7_macrophage_barplot_significance.pdf) · [SVG](outputs/01_C7_macrophage_barplot_significance.svg)
- [对照柱状图：全部15个原始RNA群，混合群另行标记](outputs/02_C7_all15_RNA_groups_barplot_significance.png) · [PDF](outputs/02_C7_all15_RNA_groups_barplot_significance.pdf) · [SVG](outputs/02_C7_all15_RNA_groups_barplot_significance.svg)
- [患者层面图与比较结果](outputs/03_C7_patient_aware_comparison.png) · [PDF](outputs/03_C7_patient_aware_comparison.pdf) · [SVG](outputs/03_C7_patient_aware_comparison.svg)
- [各群表达汇总](outputs/cluster_expression_summary.csv)
- [C7定向细胞层面比较](outputs/C7_vs_others_cell_exploratory.csv) · [患者层面主分析](outputs/C7_vs_others_patient_primary.csv)
- [患者层面敏感性分析](outputs/C7_patient_sensitivity.csv) · [检出率次要分析](outputs/C7_patient_detection_secondary.csv) · [探索性Fisher检验](outputs/C7_detection_Fisher_exploratory.csv)

## 主要结果

C7共403个细胞、79个PGAM5检出细胞；平均CP10k={summary.loc[7,'mean_CP10k']:.6f}，平均log1p(CP10k)={summary.loc[7,'mean_log1p_CP10k']:.6f}，检出率={summary.loc[7,'detection_percent']:.2f}%。
在巨噬细胞相关群的9个定向比较中，细胞层面有{int((coretests.q_BH_14<.05).sum())}个q<0.05；C4与C7不显著。全部15群中有{int((cell.q_BH_14<.05).sum())}/14个定向比较显著，C4和C14不显著。
**这些是细胞层面的探索性显著性；同一患者内细胞不是独立生物学重复。**
患者层面的14个比较中，只有{len(tested)}个满足预设覆盖门槛；均未达到BH q<0.05。未通过门槛的比较标NT，不能解释成表达相同或不显著。

| 定向比较 | 原始P | BH q（14个比较） | 柱状图标注 |
|---|---:|---:|---|
{rows}

## 表达量与统计定义

1. CP10k = PGAM5原始计数 / 原始全转录组文库总计数 × 10,000。**柱高为每个细胞先取log1p(CP10k)后求群内算术平均，保留全部零值，与原UMAP均值一致。** 同时导出原始量纲CP10k均值与检出率；两种均值不混用。C7的平均CP10k为15群最高，但平均log1p(CP10k)低于混合RNA群C14，不声称C7在所有量纲均最高。
2. 用户指定C7作为参照，事先定义C7与其他14个原始群的双侧Mann–Whitney U检验（渐近法、并列值和连续性校正）。BH对整个14比较族校正。主图只展示其中9个巨噬细胞相关群对比，**没有重新缩小校正族或重算更有利的q值**。U检验比较表达分布，星号不应当作对算术均值的独立生物学重复检验。
3. 主图保留C0、C2、C3、C4、C6、C7、C9、C11、C12、C13。依据原注释，C1为低复杂度混合RNA，C5为内皮RNA混合，C8为T/NK RNA混合，C10为粒细胞RNA混合，C14为肝细胞RNA混合，均待身份复核；它们不作为已确认巨噬细胞亚群命名，保留在15群对照图并以灰色斜线和问号标明。C0/C2仍保留原巨噬细胞相关注释和肝源/肿瘤RNA限定。此为基于原注释的展示分层，不删除原始细胞。
4. 患者层面主分析：每患者-亚群至少10个细胞，每对至少5位共享合格患者；每位患者贡献一个含零值的平均log1p(CP10k)，等权配对比较。不存在的亚群不补零。用双侧精确配对标签交换检验，枚举2^n种符号翻转并比较绝对均值差；零假设要求配对标签可交换/差值符号对称。
5. 患者比较同样按14个预定比较做BH；不能检验的比较以P=1参与校正，输出P/q留空并标NT。次要检出率和敏感性分析分别组成各自14比较校正族。敏感性涵盖细胞数阈值1/10/20、平均log1p(CP10k)、平均CP10k、汇总计数pseudobulk CP10k，不据此挑选显著结论。
6. 7位患者的独立性沿用原分析中用户已确认的样本映射，详见源目录原有HCC_patient_sample_mapping.csv。这里没有重新识别或额外验证患者身份。10个细胞、5位患者为本次操作门槛，不是通用生物学标准；原分析候选群方向一致性使用的20细胞规则另在敏感性表中保留。

星号：`* q<0.05`、`** q<0.01`、`*** q<0.001`、`**** q<0.0001`；`ns`表示q>=0.05；`NT`表示配对患者不足。

## 解释限制和复现

许多群由单一患者主导；细胞层面极小P值不能取代跨患者重复。原数据尚未校正环境RNA/双细胞，未控制测序深度、细胞周期等混杂。本结果衡量RNA表达，不是蛋白丰度、功能活性或因果作用，也不自动验证TCGA反卷积用途。

```bash
python -m pip install -r requirements.txt
MPLCONFIGDIR=/tmp/gse202642-mpl python analyze_C7.py
```

本目录仅发布统计汇总、图和代码，不重复发布逐细胞数据、患者映射或含患者标识的明细。输入读取同仓库相邻的 `GSE202642_PGAM5_subpopulation_assessment` 原目录，source_provenance.json记录源提交和SHA-256；请克隆完整仓库后运行。也可通过环境变量 `GSE202642_SOURCE_DIR` 指定原数据目录。运行将重生成汇总outputs和此README。
原始计数/全库分母已重算，15群均值、检出数与原表及UMAP映射逐项一致；精确置换和BH实现经过独立数值核对。outputs/verification.json保存核查结果。

方法参考：[SciPy Mann–Whitney U](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html)；[BH校正](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html)；[单细胞伪重复问题](https://doi.org/10.1038/s41467-021-21038-1)。
'''
    (ROOT/'README.md').write_text(readme)
    print(json.dumps(verification,indent=2))
    print(coretests[['target','comparator','q_BH_14','significance']].to_string(index=False))

if __name__=='__main__':main()
