# OEP000321：来源已核查，RNA表达及临床输入尚待提供

已确认该项目对应Gao等2019年Cell论文 **Integrated Proteogenomic Characterization of HBV-Related Hepatocellular Carcinoma**，PMID31585088、DOI10.1016/j.cell.2019.08.052。

**本轮尚未计算OEP000321的NNLS系数或生存统计量。** 当前云环境代理拦截官方NODE项目页、原论文补充资料下载主机和一个候选公开Zenodo档案。由于用户无法保存/发布环境，本轮不再依赖修改环境设置；可通过提供公开源文件继续。

已核查NCBI原始论文记录及公开PMC再分析文献、GitHub公开仓库。可访问来源主要是研究说明、蛋白矩阵或派生统计表，尚未找到能够核验原患者身份的完整基因级RNA表达矩阵及配套临床随访。不会用蛋白数据替代RNA参考，也不会从论文汇总HR或KM图推断逐患者随访。

## 所需文件

请从[OEP000321官方项目](https://www.biosino.org/node/project/detail/OEP000321)及[原论文](https://doi.org/10.1016/j.cell.2019.08.052)提供：

1. **基因级RNA表达矩阵**：行有原基因标识，列有患者/样本编号。优先原发布的处理后矩阵，注明是counts、FPKM、TPM或log2值；保留原样本名称和组织类型。无需下载原始FASTQ或蛋白矩阵。
2. **临床随访表**：有相同患者编号、生存时间、时间单位及死亡/删失或复发事件。原论文Table S1是临床资料线索，但该表实际内容和终点尚未在本环境读取，不能先宣称其具有某个终点或有效患者数。
3. 若表达矩阵列名本身无法区分tumor/adjacent，附上原样本对应表。

可上传CSV/TSV/TXT、Excel及其gzip/zip文件；保留原表完整字段。使用原公开患者编号，不需要患者姓名。只有实际输入核验后，才能确定本次RNA矩阵覆盖数及各终点的可用人数。

## 已固定的下一步方法

沿用此前850基因×17组分v4联合参考和原行权重，固定NNLS；只取本平台实际共有基因，不依据生存结局更改参考。输入后首先核验基因ID、表达单位/变换、组织选择、重复及完整编号链接。以全部合格肿瘤表达样本的中位数分High/Low，再处理终点缺失；使用原记录的OS/RFS/DFS/PFS定义，生成KM风险人数、logrank、Cox及预定校正。不会将RFS/DFS改称PFS。

参考文件SHA256及稳定GitHub来源保存在`frozen_reference_manifest.csv`；来源核查和缺失输入状态在`source_status.json`。该参考未通过此前丰度验证，即使取得数据，输出也仍是探索性系数，不是已验证的真实细胞比例。

## 可核查来源

- [官方NODE项目](https://www.biosino.org/node/project/detail/OEP000321)：当前云环境CONNECT403，未取得数据目录或RNA文件。
- [NCBI原论文记录](https://pubmed.ncbi.nlm.nih.gov/31585088/)：元数据已通过允许访问的NCBI ESummary接口读取，保存于`original_study_pubmed_metadata.json`。
- [PMC8226228](https://pmc.ncbi.nlm.nih.gov/articles/PMC8226228/)：再分析文献明确引用OEP000321及原研究，不是本次RNA/随访源表。
- [PMC12524807](https://pmc.ncbi.nlm.nih.gov/articles/PMC12524807/)：TAM再分析报告引用OEP000321；未从该文取得本次需要的完整原始队列表。
- GitHub上一个二次研究的来源说明指向Cell2019 Table S1和Zenodo14553766；说明本身不能替代数据核验。Zenodo也被本云环境网络策略拦截，故没有采用其未读取的数据或统计结果。

当前阻碍是**表达矩阵及逐患者随访未取得**，不是已经证明该项目没有相关数据。
