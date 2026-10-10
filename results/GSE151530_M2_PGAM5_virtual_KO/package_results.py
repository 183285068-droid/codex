from pathlib import Path
import json,shutil,hashlib,zipfile,sys
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;U=R/'unsparsified'
b=json.loads((R/'stability_summary.json').read_text());s=json.loads((U/'stability_summary.json').read_text());a=json.loads((R/'audit_summary.json').read_text());v=json.loads((U/'audit_summary.json').read_text());m=pd.read_csv(R/'run_metrics.csv');n=pd.read_csv(U/'run_metrics.csv');pilot=json.loads((R/'high_rank_pilot/151530/metrics.json').read_text());top=pd.read_csv(U/'TOP20_downstream_genes.csv');rec=pd.read_csv(U/'recurrent_candidates.csv')
assert a['checks']==a['passed'] and v['checks']==v['passed']
counts=lambda t:', '.join(map(str,t.significant_downstream_genes))
table='\n'.join(f'|{r.rank}|{r.gene}|{r.median_distance*1e7:.3f}|{r.significant_seed_count}/5|' for r in top.itertuples())
fig,ax=plt.subplots(figsize=(10,4));groups=[m.significant_downstream_genes.to_numpy(),n.significant_downstream_genes.to_numpy(),np.array([pilot['significant_downstream_genes']])];labels=['Sparse q=0.95, rank=5\n5 seeds','Unfiltered q=0, rank=5\n5 seeds','Sparse high-rank pilot\n1 seed'];colors=['#888888','#277ea7','#ca983c']
for i,(g,col) in enumerate(zip(groups,colors)):
 ax.scatter(i+(np.linspace(-.08,.08,len(g)) if len(g)>1 else np.array([0])),g,c=col,s=45);ax.text(i,float(g.max())+1.2,', '.join(map(str,g)),ha='center',fontsize=10)
ax.set_xlim(-.5,2.5);ax.set_xticks(range(3),labels);ax.set_ylabel('Downstream genes with model q<0.05');ax.set_ylim(-2,float(n.significant_downstream_genes.max())+5);ax.set_title('GSE151530 M2 | sensitivity to network sparsification\nNumerical controls: zero significant genes in all 11 runs');fig.tight_layout();fig.savefig(R/'07_model_configuration_comparison.png',dpi=180);fig.savefig(R/'07_model_configuration_comparison.pdf');plt.close(fig)
readme=f'''# GSE151530 M2：PGAM5虚拟敲除与参数敏感性

固定前一轮GSE151530 HCC巨噬细胞16群版本中的**M2：MKI67/TOP2A增殖型，171个细胞**。来自21个样本、17个患者编号代理；PGAM5检出22个细胞，149个零值全部保留。H70贡献88/171个细胞（51.46%）。本轮不按PGAM5检出进一步选择细胞，不把样本数当作独立患者数。

## 结果与可解释范围

完成**5次稀疏主模型、5次关闭稀疏化的敏感性模型、1次高阶试算**，以及对应11次相同WT网络无敲除数值对照。

|方案|每次模型q<0.05的非靶基因数|≥4/5次显著|5次均显著|整体稳定标准|
|---|---|---:|---:|---|
|主模型：q=0.95，秩5|{counts(m)}|0|0|未通过：无可分辨敲除信号|
|敏感性：q=0，秩5|{counts(n)}|{s['recurrent_ge4']}|{s['recurrent_all5']}|{'通过' if s['passes_prespecified_stability'] else '未通过'}|
|高阶试算：q=0.95，秩10，1次|{pilot['significant_downstream_genes']}|不适用|不适用|未产生有效敲除信号|

稀疏主模型5次均没有超出数值噪声阈值的下游位移，不能将它理解为PGAM5没有生物学作用，也不能对约10⁻¹⁵的数值误差做TOP20生物学排名。4次WT PGAM5行完全为零；种子101张量拟合产生弱出边，但最大下游位移5.11×10⁻¹⁰仍低于约1.05×10⁻⁸的数值分辨阈值。种子151530的30个原始重采样网络中，PGAM5出边均被q=0.95阈值过滤，见[网络过滤诊断](06_PGAM5_edge_filtering_diagnostic.png)及[逐网络表](PGAM5_network_sparsification_diagnostic.csv)。独立重建稀疏网络确认每个保留PGAM5出边数均为0。

随后固定关闭全局稀疏化（q=0，对全部边一致处理），保留原始回归权重，其他主模型参数、细胞、基因和种子不变，完成5次敏感性重复。该方案每次得到26–28个模型候选，18个至少4/5次显著，15个5次均显著；**这是条件于该网络参数的探索性候选**。最低/中位两两排序rho={s['min_distance_rho']:.3f}/{s['median_distance_rho']:.3f}，前50 Jaccard最低/中位={s['min_top50_Jaccard']:.3f}/{s['median_top50_Jaccard']:.3f}，仍未达到预设最低rho≥0.8及Jaccard≥0.5。稀疏与不稀疏方案结果差异明显，不能声称候选对算法选择稳健。

所有11次分解均达到各自残差变化收敛标准，所有11次无敲除对照均无显著基因。技术重复不等于独立患者验证。高阶试算仅完成种子151530，秩10、容差1e-7的分解在1,577次迭代后收敛，但没有PGAM5出边；其结果保留于high_rank_pilot，不与5次主模型合并候选。

## 图形

- [M2定位及敲除前PGAM5表达UMAP](01_M2_target_and_observed_PGAM5_UMAP.png)、[样本构成](02_M2_sample_representation.png)。UMAP沿用上一轮观测坐标，右侧只显示M2实际PGAM5 RNA表达，不是敲除后细胞坐标。
- **[敏感性模型TOP20条形图](unsparsified/03_PGAM5_virtual_KO_TOP20.png)**：条形为5次中位无方向网络位移，误差线为Q1–Q3；蓝色5/5、橙色4/5、灰色≤3/5次模型q<0.05。误差线不是患者层面置信区间。
- [TOP20跨种子热图](unsparsified/05_TOP20_cross_seed_heatmap.png)、[排序相关及对照](unsparsified/04_repetition_stability_and_controls.png)。
- [两种方案与高阶试算的对比](07_model_configuration_comparison.png)。
- [稀疏模型的有效干预诊断](04_repetition_stability_and_controls.png)、[无可解释TOP20的说明图](03_PGAM5_virtual_KO_TOP20.png)。

全部提供同名PDF。TOP20按中位位移排序，不等同于18个≥4/5次显著的候选名单。

|排名|敏感性模型非靶基因|中位无方向网络位移（×10⁻⁷）|模型q<0.05次数|
|---:|---|---:|---:|
{table}

[TOP20精确统计](unsparsified/TOP20_downstream_genes.csv) · [18个重复候选](unsparsified/recurrent_candidates.csv)。PGAM5始终包含在输入和全部模型统计中；其被强制干预，单独保留在各方案的PGAM5_forced_target_summary.csv，不计入独立下游证据。

候选中出现补体/巨噬标志C1QC、C3，抗原呈递相关HLA-DPA1，以及GBP1、IFI44L、ISG15等干扰素相关RNA；同时出现CHGA、PGA5、TFF3、PIGR等上皮/分泌相关RNA，提示须谨慎解释环境RNA或混合信号。这些是名单内容的描述，未进行新的通路富集检验，不能据此判定敲除后巨噬细胞抗肿瘤能力增强或减弱。

## 方法及核查

使用scTenifoldpy 0.5.1的网络构建、CP张量分解、PGAM5网络行置零、30维流形对齐和差异调控函数，复现scTenifoldKnk方法步骤。输入来自上一轮[完整巨噬矩阵](https://github.com/183285068-droid/codex/blob/main/results/GSE151530_HCC_hierarchical_clustering/macrophages_annotated_full.h5ad)，哈希见input_summary.json。输入筛选取原始检出率≥5%，排除MT-/RPL/RPS，按log1p(CP10k)的Seurat标准化离散度取前1,500基因，再加入PGAM5，最终1,501基因。网络用CPM，分母是原始18,661基因完整文库总计数，不是子集总计数。

每次建立30个网络，每个有放回抽取171个细胞，PC回归3成分。主模型/敏感性模型CP秩5、最多1,000次迭代、容差1e-5；高阶试算秩10、最多2,000次、容差1e-7。张量归一化后保留6位小数并转置、去除对角自环。PGAM5行置零，其他行保持完全一致。预设技术稳定标准记录于protocol.json：所有两两rho≥0.8、前50 Jaccard≥0.5、全部分解收敛、全部WT PGAM5行非零、全部无敲除对照无显著。没有在看敏感性候选后修改标准；不对低于数值分辨率的位移计算生物学排序一致性。

**Distance无方向；模型FC是位移平方/背景期望，不是RNA表达fold change。**模型卡方p和BH q不是实际敲除实验的差异表达检验。重复频率与median_model_q不是合并p值。未额外回归测序深度、患者或细胞周期；稀疏PGAM5检出、样本不均衡、可能的环境/混合RNA及网络参数依赖限制因果解释。

独立核查：稀疏主模型及高阶试算**{a['checks']}/{a['passed']}项通过**；敏感性模型**{v['checks']}/{v['passed']}项通过**。包含全部171个细胞及原始计数、完整文库分母、22个检出/149零值、网络行置零与其他行完全不变、实际KO矩阵哈希、对齐位移、独立卡方p及BH q重算、11次无敲除对照及TOP20噪声过滤/排序。审计证明计算实现正确，不改变整体稳定性未通过的结论。

网络文件为WT、KO_delta和实际KO SHA256；为保留浮点正负零的字节一致性，KO由WT复制后将PGAM5行设为正零，再核对KO_delta和哈希。该存储优化没有修改模型数值或统计结果。完整分解残差、对齐坐标和单次结果均保存。

## 下载与复现

为使每份小于100 MiB，结果分成两个**独立ZIP包，无须二进制合并**：

- [稀疏主模型及高阶试算包](GSE151530_M2_PGAM5_virtual_KO_results.zip)：输入、诊断图、5次主模型、1次高阶试算及代码。
- [关闭稀疏化敏感性模型包](GSE151530_M2_PGAM5_virtual_KO_unsparsified_results.zip)：5次敏感性模型、TOP20图、统计及独立输入/代码。解压后进入unsparsified目录。

GitHub打开ZIP后点击 **Download raw file**。每个包包含自己的manifest.json并通过ZIP CRC检查。Python3.12，建议4线程、≥8 GiB内存；敏感性模型峰值约3.2 GiB。各方案目录运行：

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
export MPLCONFIGDIR="$PWD/.matplotlib"
python visualize.py
python audit.py
```

从头复现，在新目录复制该方案inputs/、prepare_inputs.py、run.py、visualize.py、audit.py、requirements.txt及sample_summary.csv，运行python prepare_inputs.py、python run.py、python visualize.py、python audit.py。run.py仅跳过已有metrics.json的完成种子。稀疏模型的高阶试算可另运行python run_high_rank_pilot.py 151530；audit.py在存在高阶结果时额外检查该结果。network_diagnostic.py复现种子151530的30个原始与过滤网络诊断。

全部输入原始RNA保持不变。M2_raw.h5ad含完整18,661基因；network_raw_counts.csv.gz为1,501基因子集，勿用子集求和重新归一化。不同聚类版本编号不可直接套用M2标签。PNG/PDF和CSV可直接下载查看，所有候选均按探索性网络预测解释。
'''
(R/'README.md').write_text(readme)
(U/'README.md').write_text('# 关闭全局稀疏化的PGAM5虚拟敲除敏感性模型\n\n使用同一GSE151530 M2的171个细胞、1,501基因，q=0、CP秩5；5次模型显著非靶基因数为'+counts(n)+'。18个至少4/5次显著，15个5次均显著；整体预设稳定标准未通过，且结果与q=0.95主模型不同。\n\n[TOP20条形图](03_PGAM5_virtual_KO_TOP20.png) · [跨种子热图](05_TOP20_cross_seed_heatmap.png) · [稳定性及对照](04_repetition_stability_and_controls.png)。位移无方向，不是RNA上调/下调；误差线为5次Q1–Q3，不是患者置信区间。PGAM5在完整结果及forced_target表保留。\n\n本包内可独立运行python visualize.py和python audit.py；安装requirements中的依赖，并设置OPENBLAS_NUM_THREADS=4和MPLCONFIGDIR为可写目录。新目录复制inputs及脚本后用python prepare_inputs.py、python run.py、python visualize.py、python audit.py从头复现。所有原始完整文库计数分母均保存。网络文件以WT及KO_delta保存，KO由WT复制后PGAM5行置零并核对差分和SHA256。\n\n完整方法、参数依赖及输入来源见[主报告](https://github.com/183285068-droid/codex/blob/main/results/GSE151530_M2_PGAM5_virtual_KO/README.md)。\n')
D=Path(sys.argv[1]).expanduser().resolve() if len(sys.argv)>1 else R.parent/(R.name+'_packaged');assert D!=R and R not in D.parents;D.mkdir(parents=True,exist_ok=True)
for item in R.iterdir():
 if item.name in ['__pycache__','.matplotlib','package.log'] or item.suffix=='.zip' or item.name=='manifest.json':continue
 if item.is_dir():shutil.copytree(item,D/item.name,dirs_exist_ok=True,ignore=shutil.ignore_patterns('__pycache__','.matplotlib'))
 else:shutil.copy2(item,D/item.name)
for log in D.rglob('*.log'):log.write_text('\n'.join(line.rstrip() for line in log.read_text().splitlines())+'\n')
for sub in [D,D/'unsparsified']:
 files=sorted(f for f in sub.rglob('*') if f.is_file() and f.suffix!='.zip' and f.name!='manifest.json' and (sub!=D or 'unsparsified' not in f.relative_to(D).parts))
 (sub/'manifest.json').write_text(json.dumps([dict(file=str(f.relative_to(sub)),bytes=f.stat().st_size,SHA256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files],indent=2));name='GSE151530_M2_PGAM5_virtual_KO_results.zip' if sub==D else 'GSE151530_M2_PGAM5_virtual_KO_unsparsified_results.zip';z=D/name
 with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as out:
  for f in files+[sub/'manifest.json']:out.write(f,str(f.relative_to(D)))
 with zipfile.ZipFile(z) as out:assert out.testzip() is None
 assert z.stat().st_size<100*1024**2
 print(name,round(z.stat().st_size/1024**2,2),'MiB, CRC verified',flush=True)
print('PACKAGED',D,flush=True)
