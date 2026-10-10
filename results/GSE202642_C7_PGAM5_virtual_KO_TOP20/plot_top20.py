from pathlib import Path
import pandas as pd,numpy as np,matplotlib,json,hashlib,zipfile,shutil
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
R=Path(__file__).resolve().parent;src=R/'inputs';src.mkdir(exist_ok=True)
for seed in [202642,42,101,202,303]:
 f=src/f'seed{seed}.csv'
 if not f.exists():shutil.copy2(R.parent/'GSE202642_C7_PGAM5_virtual_KO_v2/ensemble'/str(seed)/'differential_regulation.csv',f)
data={s:pd.read_csv(src/f'seed{s}.csv').set_index('Gene') for s in [202642,42,101,202,303]};distance=pd.DataFrame({s:v.Distance for s,v in data.items()});q=pd.DataFrame({s:v['adjusted p-value'] for s,v in data.items()});t=pd.DataFrame({'gene':distance.index,'median_distance':distance.median(axis=1),'Q1_distance':distance.quantile(.25,axis=1),'Q3_distance':distance.quantile(.75,axis=1),'significant_seeds':q.lt(.05).sum(axis=1),'min_model_q':q.min(axis=1),'max_model_q':q.max(axis=1)}).reset_index(drop=True);t['is_forced_KO_target']=t.gene.eq('PGAM5');t.to_csv(R/'all_genes_perturbation_summary.csv',index=False);target=t[t.is_forced_KO_target];target.to_csv(R/'PGAM5_forced_target_summary.csv',index=False);top=t[~t.is_forced_KO_target].sort_values(['median_distance','gene'],ascending=[False,True]).head(20).copy();top.insert(0,'rank',range(1,21));top.to_csv(R/'TOP20_downstream_genes.csv',index=False)
v=top.iloc[::-1];scale=1e7;fig,ax=plt.subplots(figsize=(10,10));colors=['#2878a5' if n==5 else '#e99b3d' if n==4 else '#aaaaaa' for n in v.significant_seeds];ax.barh(v.gene,v.median_distance*scale,color=colors);ax.errorbar(v.median_distance*scale,np.arange(20),xerr=np.vstack(((v.median_distance-v.Q1_distance)*scale,(v.Q3_distance-v.median_distance)*scale)),fmt='none',ecolor='#333333',elinewidth=1.1,capsize=3)
limit=float(v.Q3_distance.max()*scale);ax.set_xlim(0,limit*1.27)
for i,r in enumerate(v.itertuples()):ax.text(max(r.Q3_distance,r.median_distance)*scale+limit*.015,i,f'{r.significant_seeds}/5',va='center',fontsize=9)
ax.set_xlabel('Median unsigned network displacement (×10⁻⁷)\nBars: median across 5 seeds; error bars: interquartile range');ax.set_title('PGAM5 virtual knockout | C7 MKI67/TOP2A macrophages\nTOP20 downstream genes by median perturbation',pad=18);ax.legend(handles=[Patch(color='#2878a5',label='Model FDR < 0.05 in 5/5 seeds'),Patch(color='#e99b3d',label='4/5 seeds'),Patch(color='#aaaaaa',label='≤3/5 seeds')],loc='lower right',fontsize=8);ax.spines[['top','right']].set_visible(False);fig.text(.5,.015,'Unsigned model prediction; no expression up/down direction. Overall model stability criteria not met.\nPGAM5 is the forced KO target and is reported separately.',ha='center',fontsize=9);fig.tight_layout(rect=(0,.06,1,1));fig.savefig(R/'PGAM5_KO_TOP20_barplot.png',dpi=180);fig.savefig(R/'PGAM5_KO_TOP20_barplot.pdf');plt.close(fig)
lines=['|排名|基因|中位网络位移（×10⁻⁷）|模型FDR<0.05次数|','|---:|---|---:|---:|']
for r in top.itertuples():lines.append(f'|{r.rank}|{r.gene}|{r.median_distance*scale:.3f}|{r.significant_seeds}/5|')
(R/'README.md').write_text('# C7 PGAM5虚拟敲除：TOP20下游网络扰动\n\n使用v2中增加重采样方案的5次运行（种子202642、42、101、202、303），按非靶基因中位流形位移Distance降序选TOP20。条形长度为中位值，误差线为5次运行的25%–75%分位范围；蓝色5/5次模型FDR<0.05、橙色4/5、灰色≤3/5。频率不是合并p值，误差线不是患者层面的置信区间。\n\n'+ '\n'.join(lines)+ '\n\nPGAM5本身被强制干预，单独保存在PGAM5_forced_target_summary.csv，并在all_genes_perturbation_summary.csv保留。主图比较TOP20下游基因，避免把强制目标自身的位移当成下游证据。网络位移无方向，不能解释为RNA上调或下调；模型FC不是表达fold change。本轮整体模型未通过预设稳定性标准。这里TOP20按中位位移排序，与按显著频率优先排序的21个重复候选名单不同。\n\n候选主要涉及热休克/蛋白稳态（HSP、DNAJB1、PPP1R15A）、即时应答（FOS/JUN、DUSP1、ATF3、NFKBIA），另有SOD2、PLIN2等。此处为描述性汇总，不新增通路检验或推断抗肿瘤功能方向。\n\n输入5张完整模型统计表已随包保存。复现：安装numpy、pandas、matplotlib后运行python plot_top20.py。PNG及PDF供查看和导出，CSV供复核；SHA256SUMS记录文件完整性。\n');assert len(top)==20 and top.gene.nunique()==20;assert top.median_distance.is_monotonic_decreasing;(R/'validation.json').write_text(json.dumps({'genes':20,'seeds':5,'ranking':'median unsigned Distance','uncertainty':'Q1 to Q3 across seeds, not confidence interval','PGAM5_retained_separately':True},indent=2))
files=sorted(p for p in R.rglob('*') if p.is_file() and p.name!='SHA256SUMS' and p.suffix!='.zip');(R/'SHA256SUMS').write_text('\n'.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(R)) for p in files)+'\n');z=R/'GSE202642_C7_PGAM5_virtual_KO_TOP20_results.zip'
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as h:
 for p in sorted(R.rglob('*')):
  if p.is_file() and p!=z:h.write(p,str(p.relative_to(R)))
with zipfile.ZipFile(z) as h:assert h.testzip() is None
print(top[['rank','gene','significant_seeds']].to_string(index=False))
