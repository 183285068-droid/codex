# 恢复原始全基因表达矩阵

179-marker 子集无法逆推出删去的表达列。恢复指重新获取 GEO 已发布的完整计数附件，不必先从 FASTQ 重跑 Space Ranger。当前环境的 NCBI 网络访问受限；以下在能够访问 NCBI 的机器执行。

```bash
git clone https://github.com/183285068-droid/codex.git
cd codex
python -m pip install requests anndata pandas scipy
export GSE281759_WORK="$PWD/data/GSE281759"
python results/GSE281759_spatial_HCC_PGAM5/download_inputs.py
```

下载脚本获取 `source_manifest.json` 列出的附件，按字节数及 SHA256 核验，包括16张肿瘤切片。每张的三个表达文件是 `GSM编号_切片_matrix.mtx.gz`、`_features.tsv.gz`、`_barcodes.tsv.gz`。完整空间分析另需 `_tissue_positions_list.csv.gz`、`_scalefactors_json.json.gz` 及组织图像，下载脚本也会获取。

以 A009T1 为例读取完整计数并保存独立 h5ad：

```python
from pathlib import Path
import os, gzip
import numpy as np
import pandas as pd
import scipy.io
import anndata as ad

work = Path(os.environ['GSE281759_WORK'])
prefix = work / 'inputs' / 'GSM8627602_A009T1'
features = pd.read_csv(str(prefix) + '_features.tsv.gz', sep='\t', header=None)
barcodes = pd.read_csv(str(prefix) + '_barcodes.tsv.gz', header=None)[0]
with gzip.open(str(prefix) + '_matrix.mtx.gz', 'rb') as handle:
    counts = scipy.io.mmread(handle).T.tocsr()
assert counts.shape == (len(barcodes), len(features))
assert barcodes.is_unique and features[0].is_unique
a = ad.AnnData(
    X=counts,
    obs=pd.DataFrame(index=pd.Index(barcodes.astype(str), name='barcode')),
    var=pd.DataFrame({'gene_id': features[0].values,
                      'original_gene_symbol': features[1].values},
                     index=pd.Index(features[1].values, name='gene')),
)
a.var_names_make_unique()
a.obs['total_counts'] = np.asarray(counts.sum(axis=1)).ravel()
a.write_h5ad(work / 'A009T1_full_gene_raw_counts.h5ad')
```

此对象 `.X` 为完整原始计数，尚未 QC、log归一化或加入空间坐标，不能直接输入本轮脚本。对其余切片依 `sample_scope.json` 的 GSM/切片对应关系重复；不要按名字排序拼接位点。将已有基线对象的 QC 后条形码、空间信息及固定源门按条形码对齐恢复，然后用全基因总 UMI 做 CP10k/log1p，提取新增 T 标志。验证共享基因计数与已有对象一致。

现有 `preprocess.py` 的最终便携输出仍是 marker 子集，直接重跑不会自动补齐 CD4/CD8/Treg/记忆标志。应修改为另存完整对象，或扩展标志列表后重新输出到新目录，以保留既有分析。重复基因符号需保留基因 ID 并检查，不可随意把不同 ID 合并。
