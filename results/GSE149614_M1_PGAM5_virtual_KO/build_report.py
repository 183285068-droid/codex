"""Publish the final model comparison, report and independently verified ZIPs."""
from pathlib import Path
import json
import hashlib
import zipfile
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

R = Path(__file__).resolve().parent
U = R / 'unsparsified'


def read(folder):
    return (json.loads((folder/'stability_summary.json').read_text()),
            pd.read_csv(folder/'run_metrics.csv'),
            json.loads((folder/'audit_summary.json').read_text()))


def package(folder, name, exclude_sensitivity=False):
    files = sorted(p for p in folder.rglob('*') if p.is_file()
        and '__pycache__' not in p.parts and p.suffix not in ['.zip', '.pyc']
        and p.name not in ['manifest.json', 'archive_verification.json']
        and (not exclude_sensitivity or 'unsparsified' not in p.relative_to(folder).parts))
    manifest = [{'file':str(p.relative_to(folder)), 'bytes':p.stat().st_size,
                 'SHA256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files]
    (folder/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    path = R/name
    with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in files+[folder/'manifest.json']:
            z.write(p,p.relative_to(folder))
    with zipfile.ZipFile(path) as z:
        assert z.testzip() is None and len(set(z.namelist()))==len(z.namelist())
        entries=len(z.namelist())
    assert path.stat().st_size < 100*1024**2
    report={'archive':name,'bytes':path.stat().st_size,'entries':entries,
            'SHA256':hashlib.sha256(path.read_bytes()).hexdigest(),'CRC_passed':True}
    (folder/'archive_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(name,round(path.stat().st_size/1024**2,2),'MiB; CRC verified')


def main():
    b,m,a=read(R)
    s,n,v=read(U)
    assert a['checks']==a['passed'] and v['checks']==v['passed']
    info=json.loads((R/'input_summary.json').read_text())
    summary=[]
    for name,folder,stability,metrics in [('Sparse q=0.95',R,b,m),('Unfiltered q=0',U,s,n)]:
        summary.append({'model':name,'significant_counts':metrics.significant_downstream_genes.tolist(),
            'recurrent_at_least4':stability['recurrent_ge4'], 'recurrent_all5':stability['recurrent_all5'],
            'min_pairwise_rho':stability['min_distance_rho'],
            'min_top50_Jaccard':stability['min_top50_Jaccard'],
            'passes_prespecified_technical_stability':stability['passes_prespecified_stability']})
    (R/'model_comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
    pd.DataFrame(summary).to_csv(R/'model_comparison.csv',index=False)
    left=pd.read_csv(R/'downstream_consensus.csv').set_index('gene')
    right=pd.read_csv(U/'downstream_consensus.csv').set_index('gene').reindex(left.index)
    overlap=set(left.nlargest(20,'median_distance').index)&set(right.nlargest(20,'median_distance').index)
    cross_rho=float(spearmanr(left.median_distance,right.median_distance).statistic)
    comparisons=left[['median_distance','significant_seed_count']].join(
        right[['median_distance','significant_seed_count']],lsuffix='_sparse',rsuffix='_unfiltered')
    comparisons.to_csv(R/'cross_parameter_gene_comparison.csv')
    recurrent=comparisons[(comparisons.significant_seed_count_sparse>=4)&
                          (comparisons.significant_seed_count_unfiltered>=4)]
    recurrent.to_csv(R/'cross_parameter_recurrent_candidates.csv')
    (R/'cross_parameter_summary.json').write_text(json.dumps({
        'median_displacement_rank_rho':cross_rho,'TOP20_intersection_count':len(overlap),
        'TOP20_intersection_genes':sorted(overlap),
        'comparison_valid_only_if_both_models_have_displacements_above_noise':
            b['meaningful_KO_seed_count']==5 and s['meaningful_KO_seed_count']==5},indent=2)+'\n')
    fig,axs=plt.subplots(1,2,figsize=(12,4.5))
    for i,(metric,color,label) in enumerate([(m,'#277ea7','Sparse q=0.95'),(n,'#e69e39','Unfiltered q=0')]):
        vals=metric.significant_downstream_genes.to_numpy()
        axs[0].scatter(i+np.linspace(-.08,.08,5),vals,s=45,color=color,label=label)
        axs[0].text(i,vals.max()+1,', '.join(map(str,vals)),ha='center',fontsize=9)
    axs[0].set_xticks([0,1],['Sparse q=0.95','Unfiltered q=0'])
    axs[0].set_xlim(-.5,1.5)
    axs[0].set_ylim(-1,max(m.significant_downstream_genes.max(),n.significant_downstream_genes.max())+5)
    axs[0].set_ylabel('Downstream genes with model q<0.05')
    axs[0].set_title('Five repeats per model; no-KO controls = 0')
    colors=np.where(left.significant_seed_count.ge(4)&right.significant_seed_count.ge(4),'#277ea7','#bbbbbb')
    axs[1].scatter(left.median_distance,right.median_distance,s=5,c=colors,rasterized=True)
    axs[1].set_xlabel('Sparse: median unsigned displacement')
    axs[1].set_ylabel('Unfiltered: median unsigned displacement')
    axs[1].set_title(f'Parameter comparison: rank rho={cross_rho:.3f}\nTOP20 overlap: {len(overlap)}/20')
    axs[1].legend(handles=[
        Line2D([0],[0],linestyle='',marker='o',markersize=5,color='#277ea7',
               label='Both models: q<0.05 in at least 4/5 repeats'),
        Line2D([0],[0],linestyle='',marker='o',markersize=5,color='#bbbbbb',label='Other modeled genes')],
        loc='upper left',fontsize=7,frameon=False)
    fig.suptitle('GSE149614 M1 | PGAM5 virtual KO | network-parameter sensitivity')
    fig.tight_layout()
    for ext in ['png','pdf']:fig.savefig(R/('06_model_configuration_comparison.'+ext),dpi=200)
    plt.close(fig)

    metric_counts=lambda frame:', '.join(map(str,frame.significant_downstream_genes))
    status=lambda obj:'通过' if obj['passes_prespecified_stability'] else '未通过'
    number=lambda x:'不可计算（无有效位移）' if x is None else f'{x:.3f}'
    top=pd.read_csv(R/'TOP20_downstream_genes.csv')
    rows=['| 排名 | 下游基因 | 中位无方向网络位移 | 模型q<0.05次数 |',
          '|---:|---|---:|---:|']
    for x in top.itertuples():rows.append(f'| {x.rank} | {x.gene} | {x.median_distance:.6g} | {x.significant_seed_count}/5 |')
    report=f'''# GSE149614 M1：PGAM5虚拟敲除与可视化

固定GSE149614的10个原发HCC肿瘤T样本、最近9群版本中的**M1：MKI67/TOP2A增殖型巨噬细胞**。共208个细胞，19个检出PGAM5，189个零值全部保留；覆盖10个样本/患者编号，最大样本贡献{info['largest_sample_fraction']:.2%}。不重新聚类、不按PGAM5检出筛除细胞，也不将原始观测UMAP称为敲除后UMAP。

## 模型结果

完成5次稀疏主模型（网络保留阈值q=0.95）、5次关闭全局稀疏化的敏感性模型（网络阈值q=0），以及10次相同WT网络的无敲除数值对照。网络阈值q与统计学BH q是不同参数。

| 模型 | 各次模型BH q<0.05的下游基因数 | 至少4/5次显著 | 5/5次显著 | 预设技术稳定性 |
|---|---|---:|---:|---|
| 稀疏主模型 | {metric_counts(m)} | {b['recurrent_ge4']} | {b['recurrent_all5']} | {status(b)} |
| 关闭稀疏化 | {metric_counts(n)} | {s['recurrent_ge4']} | {s['recurrent_all5']} | {status(s)} |

主模型有效干预{b['meaningful_KO_seed_count']}/5次；最低/中位两两位移排序rho={number(b['min_distance_rho'])}/{number(b['median_distance_rho'])}，前50 Jaccard最低/中位={number(b['min_top50_Jaccard'])}/{number(b['median_top50_Jaccard'])}。敏感性模型有效干预{s['meaningful_KO_seed_count']}/5次；最低/中位rho={number(s['min_distance_rho'])}/{number(s['median_distance_rho'])}，Jaccard最低/中位={number(s['min_top50_Jaccard'])}/{number(s['median_top50_Jaccard'])}。

两方案中位位移的排序rho={cross_rho:.3f}，TOP20交集{len(overlap)}/20；[两方案均至少4/5次显著的共同候选](cross_parameter_recurrent_candidates.csv)共{len(recurrent)}个。技术种子重复和网络参数敏感性不能替代独立患者或真实敲除实验验证。通过本轮技术稳定标准也不证明真实因果调控、PGAM5蛋白作用或抗肿瘤功能变化。

## 图形与统计

- [M1定位及敲除前实际PGAM5表达UMAP](01_M1_target_and_observed_PGAM5_UMAP.png)、[样本构成](02_M1_sample_representation.png)。
- **[主模型TOP20条形图](03_PGAM5_virtual_KO_TOP20.png)**、[主模型TOP20重复热图](05_TOP20_cross_seed_heatmap.png)、[主模型稳定性和无敲除对照](04_repetition_stability_and_controls.png)。
- [关闭稀疏化TOP20](unsparsified/03_PGAM5_virtual_KO_TOP20.png)、[重复热图](unsparsified/05_TOP20_cross_seed_heatmap.png)、[稳定性及对照](unsparsified/04_repetition_stability_and_controls.png)。
- [网络参数敏感性对比](06_model_configuration_comparison.png)。所有PNG均有同名PDF。
- [主模型完整基因汇总](all_genes_perturbation_summary.csv)、[TOP20精确数值](TOP20_downstream_genes.csv)、[重复候选](recurrent_candidates.csv)、[跨参数逐基因对比](cross_parameter_gene_comparison.csv)。

TOP20条形表示5次模型的中位**无方向网络位移**，误差线为Q1–Q3；蓝色表示5/5次模型BH q<0.05，橙色4/5，灰色≤3/5。Q1–Q3不是患者置信区间。位移与模型FC均不是实际RNA表达fold change，不能据此将基因判为上调或下调。若有效信号少于20个，只展示通过数值噪声过滤的基因。

{chr(10).join(rows)}

PGAM5始终保留在输入、完整统计和[被强制干预靶基因表](PGAM5_forced_target_summary.csv)中；由于其被直接干预，不计入独立下游TOP20或稳定候选。

## 方法与预设标准

scTenifoldpy 0.5.1：PC回归构建30个有放回重采样网络（每个208个细胞、3成分），CP张量秩5、最大1,000次、容差1e-5，网络归一化后保留6位小数并转置、去自环，再将PGAM5对应调控行设为0；其他行与WT完全不变。30维流形对齐计算无方向位移，模型卡方p及BH q描述模型内相对扰动，不是实际敲除实验差异表达检验。

原始矩阵包含25,712基因。按检出率≥5%排除MT-/RPL/RPS，取log1p(CP10k)的Seurat标准化离散度前1,500基因，再确保加入PGAM5，共1,501基因。网络使用CPM，分母始终为完整25,712基因原始文库总计数。没有额外回归患者、测序深度或细胞周期；M1仅19个PGAM5检出细胞，存在稀疏表达和网络推断限制。

两模型使用相同细胞、基因、种子149614/42/101/202/303及其余参数。预设标准：全部张量收敛、全部WT PGAM5行非零且位移高于数值噪声、全部无敲除对照无显著、所有两两rho≥0.8、所有前50 Jaccard≥0.5。模型BH q<0.05在至少4/5次出现作为重复候选；重复频率及中位q不是合并p值。

独立审计：主模型{a['passed']}/{a['checks']}项、敏感性模型{v['passed']}/{v['checks']}项通过。包括全部208个细胞和原始计数、完整文库分母、19检出/189零值、实际KO矩阵SHA256、仅PGAM5行置零、对齐位移、独立卡方p和BH q重算、数值噪声过滤、10次无敲除对照及TOP20排序。审计证明实现正确，不自动证明生物学结论。

## 下载与复现

- [主模型完整包](GSE149614_M1_PGAM5_virtual_KO_results.zip)：输入、5次网络及数值结果、图、代码、运行依赖源码快照。
- [关闭稀疏化敏感性包](GSE149614_M1_PGAM5_virtual_KO_unsparsified_results.zip)：独立输入、5次网络、图、统计和代码。两包独立解压，无须二进制合并。

GitHub打开ZIP点击Download raw file。原始计数不改变；输入M1_raw.h5ad为完整25,712基因，network_raw_counts.csv.gz仅为选中1,501基因，不可用其行和重新归一化。原分析[完整巨噬细胞矩阵](../GSE149614_primaryT_hierarchical_clustering/macrophages_annotated_full.h5ad)的哈希记录在input_summary.json。

Python3.12，建议4线程及≥8 GiB内存。依赖版本见requirements.txt。vendor/保存本轮实际执行的scTenifold源码及许可证、SHA256，run.py优先加载该快照，避免未来同名版本变化。

```bash
pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
export MPLCONFIGDIR="$PWD/.matplotlib"
python visualize.py
python audit.py
```

从头重算请在新目录复制inputs/、vendor/、脚本、sample_summary.csv和requirements.txt后运行prepare_inputs.py、run.py、visualize.py、audit.py。run.py跳过已有metrics.json的完成种子；重新运行全部模型需使用无models/的新目录。manifest.json与archive_verification.json提供完整性记录。
'''
    (R/'README.md').write_text(report)
    (U/'README.md').write_text('# GSE149614 M1 PGAM5：关闭稀疏化敏感性模型\n\n'
        f'固定208个细胞、1,501基因、5个种子，网络q=0。各次模型BH q<0.05下游基因数：{metric_counts(n)}。预设技术稳定性：{status(s)}。\n\n'
        '[TOP20](03_PGAM5_virtual_KO_TOP20.png) · [热图](05_TOP20_cross_seed_heatmap.png) · [稳定性与对照](04_repetition_stability_and_controls.png)。\n\n'
        '位移无方向，不是表达上调/下调。误差线为技术重复Q1–Q3。完整输入及vendor/依赖源码支持独立重算；安装requirements.txt，运行python visualize.py及python audit.py。\n\n'
        '[完整方法与主模型比较](https://github.com/183285068-droid/codex/blob/main/results/GSE149614_M1_PGAM5_virtual_KO/README.md)。\n')
    package(U,'GSE149614_M1_PGAM5_virtual_KO_unsparsified_results.zip')
    package(R,'GSE149614_M1_PGAM5_virtual_KO_results.zip',exclude_sensitivity=True)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
