# GSE149614 M1 PGAM5：关闭稀疏化敏感性模型

固定208个细胞、1,501基因、5个种子，网络q=0。各次模型BH q<0.05下游基因数：26, 27, 27, 26, 23。预设技术稳定性：通过。

[TOP20](03_PGAM5_virtual_KO_TOP20.png) · [热图](05_TOP20_cross_seed_heatmap.png) · [稳定性与对照](04_repetition_stability_and_controls.png)。

位移无方向，不是表达上调/下调。误差线为技术重复Q1–Q3。完整输入及vendor/依赖源码支持独立重算；安装requirements.txt，运行python visualize.py及python audit.py。

[完整方法与主模型比较](https://github.com/183285068-droid/codex/blob/main/results/GSE149614_M1_PGAM5_virtual_KO/README.md)。
