# 关闭全局稀疏化的PGAM5虚拟敲除敏感性模型

使用同一GSE151530 M2的171个细胞、1,501基因，q=0、CP秩5；5次模型显著非靶基因数为26, 26, 28, 26, 28。18个至少4/5次显著，15个5次均显著；整体预设稳定标准未通过，且结果与q=0.95主模型不同。

[TOP20条形图](03_PGAM5_virtual_KO_TOP20.png) · [跨种子热图](05_TOP20_cross_seed_heatmap.png) · [稳定性及对照](04_repetition_stability_and_controls.png)。位移无方向，不是RNA上调/下调；误差线为5次Q1–Q3，不是患者置信区间。PGAM5在完整结果及forced_target表保留。

本包内可独立运行python visualize.py和python audit.py；安装requirements中的依赖，并设置OPENBLAS_NUM_THREADS=4和MPLCONFIGDIR为可写目录。新目录复制inputs及脚本后用python prepare_inputs.py、python run.py、python visualize.py、python audit.py从头复现。所有原始完整文库计数分母均保存。网络文件以WT及KO_delta保存，KO由WT复制后PGAM5行置零并核对差分和SHA256。

完整方法、参数依赖及输入来源见[主报告](https://github.com/183285068-droid/codex/blob/main/results/GSE151530_M2_PGAM5_virtual_KO/README.md)。
