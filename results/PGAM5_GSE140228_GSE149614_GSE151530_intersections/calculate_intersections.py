from pathlib import Path
import json,hashlib,zipfile,itertools
import pandas as pd
D=Path(__file__).resolve().parent
root=D.parent
paths={
 'GSE140228_Droplet':root/'GSE140228_allTumor_PGAM5/Droplet/upregulated.csv',
 'GSE140228_Smartseq2':root/'GSE140228_allTumor_PGAM5/Smartseq2/upregulated.csv',
 'GSE149614':root/'GSE149614_PGAM5_primaryT/GSE149614_primaryT_PGAM5_upregulated.csv',
 'GSE151530':root/'PGAM5_HCC/HCC_PGAM5_FDR005_log2FC1_upregulated.csv'
}
# Existing repository filenames are checked rather than silently replacing source versions.
for k,p in paths.items():
 if not p.exists():
  options=list(p.parent.glob('*upregulated.csv'))
  assert len(options)==1,(k,options)
  paths[k]=options[0]
tables={k:pd.read_csv(p,index_col=0) for k,p in paths.items()}
sets={}
for k,t in tables.items():
 assert t.FDR.lt(.05).all() and t.log2FC.ge(1).all()
 assert t[['fraction_PGAM5_positive','fraction_PGAM5_undetected']].max(axis=1).ge(.1).all()
 assert not t.gene.eq('PGAM5').any() and t.gene.notna().all()
 sets[k]=set(t.gene)
common_previous=sets['GSE149614']&sets['GSE151530']
droplet=common_previous&sets['GSE140228_Droplet']
smart=common_previous&sets['GSE140228_Smartseq2']
both=droplet&smart;either=droplet|smart
# Preserve every source gene row. This avoids choosing a row for duplicated symbols.
cols=['log2FC','FDR','fraction_PGAM5_positive','fraction_PGAM5_undetected']
def export(filename,genes,keys):
 rows=[]
 for gene in sorted(genes):
  row={'gene':gene}
  for key in keys:
   source=tables[key].loc[tables[key].gene.eq(gene)]
   assert len(source)>=1
   row[key+'_row_IDs']=';'.join(map(str,source.index))
   for c in cols:
    row[key+'_'+c]=float(source[c].iloc[0]) if len(source)==1 else ';'.join(map(str,source[c]))
  rows.append(row)
 pd.DataFrame(rows).to_csv(D/filename,index=False)
export('triple_intersection_GSE140228_Droplet.csv',droplet,['GSE140228_Droplet','GSE149614','GSE151530'])
export('triple_intersection_GSE140228_Smartseq2.csv',smart,['GSE140228_Smartseq2','GSE149614','GSE151530'])
export('strict_shared_in_all_four_analysis_lists.csv',both,list(tables))
presence=pd.DataFrame({'gene':sorted(either)})
for k in sets:presence[k+'_upregulated']=presence.gene.isin(sets[k])
presence.to_csv(D/'triple_intersection_GSE140228_either_platform.csv',index=False)
pairs=[]
for a,b in itertools.combinations(sets,2):
 pair=sets[a]&sets[b];pairs.append({'list_A':a,'list_B':b,'shared_upregulated_symbols':len(pair)})
 export('pair_'+a+'_'+b+'.csv',pair,[a,b])
pd.DataFrame(pairs).to_csv(D/'pairwise_intersection_counts.csv',index=False)
summary={'criteria':'FDR<0.05, log2FC>=1, >=10% detection in at least one group; PGAM5 excluded; symbol-level intersection',
 'sources':{k:{'path':str(p.relative_to(root)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'rows':len(tables[k]),'unique_symbols':len(sets[k]),'duplicate_symbols':sorted(tables[k].loc[tables[k].gene.duplicated(False),'gene'].unique().tolist())} for k,p in paths.items()},
 'triple_Droplet':{'count':len(droplet),'genes':sorted(droplet)},
 'triple_Smartseq2':{'count':len(smart),'genes':sorted(smart)},
 'triple_either_GSE140228_platform':{'count':len(either),'genes':sorted(either)},
 'strict_all_four_lists':{'count':len(both),'genes':sorted(both)},
 'notes':['GSE140228 Droplet includes 7 HCC and 1 CC Tumor samples.','GSE140228 technologies are analyzed separately and share donor D20171109.','Smartseq2 has 1637 gene rows but 1636 symbols because HULC has two Ensembl entries; symbol-level intersections deduplicate this symbol.','Descriptive intersections; no reanalysis, enrichment significance, cross-patient reproducibility or validated signature is claimed.']}
(D/'intersection_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
text='''# GSE140228、GSE149614、GSE151530：PGAM5巨噬细胞上调基因交集

沿用此前上调结果：FDR <0.05、log2FC ≥1、至少一组检出率 ≥10%，候选表均排除定义分组的PGAM5。不重新进行差异检验，按官方gene symbol匹配；没有进行别名转换。

'''
text+=f'GSE140228 Droplet∩GSE149614∩GSE151530：**{len(droplet)}个**。\n\n'+ '、'.join(sorted(droplet))+'。\n\n'
text+=f'GSE140228 Smart-seq2∩GSE149614∩GSE151530：**{len(smart)}个**。\n\n'+ '、'.join(sorted(smart))+'。\n\n'
text+=f'三个数据集且GSE140228两种技术均上调：**{len(both)}个**。\n\n'+ '、'.join(sorted(both))+'。\n\n'
text+=f'若GSE140228任一技术达到上调标准即可，与另外两个数据集的共同上调基因有**{len(either)}个**。该口径是两技术候选集合并集；未混合两技术的原始计数进行检验。\n\n'
text+='''GSE140228 Droplet有671个上调基因；Smart-seq2有1637条上调基因记录，对应1636个不同symbol（HULC为两条不同Ensembl记录）。GSE149614有512个、GSE151530有740个上调基因。交集统计按不同symbol计数。导出表保留各来源原始行ID、log2FC、FDR及两组检出率；如一个symbol对应多条记录，使用分号保存全部值，不擅自挑选其中一条。

## 文件

- triple_intersection_GSE140228_Droplet.csv：14个交集及三份分析中的效应值和FDR。
- triple_intersection_GSE140228_Smartseq2.csv：16个交集及三份分析中的效应值和FDR。
- strict_shared_in_all_four_analysis_lists.csv：5个基因及四份分析中的效应值和FDR。
- triple_intersection_GSE140228_either_platform.csv：25个基因及各分析上调状态。
- pair_*.csv：所有两两交集的基因及对应数值。
- pairwise_intersection_counts.csv：两两交集计数。
- intersection_summary.json：口径、全部交集基因、输入相对路径及SHA256。
- calculate_intersections.py：从仓库已有结果表重新计算本目录全部文件。

## 解释

5个共同上调基因可作为后续验证的候选核心集合，尚不能称为已验证的PGAM5⁺巨噬细胞signature。未计算交集的富集显著性，也未评估这些基因的巨噬细胞特异性、患者间可重复性或功能因果关系。

GSE140228的Droplet结果包含7个HCC及1个CC样本，不能称为HCC特异结果。GSE140228两种技术共享患者D20171109，技术交集不能视为额外独立患者队列验证。本交集继承原分析不进行患者配对、深度匹配或协变量校正的条件。
'''
(D/'结果说明.md').write_text(text)
zp=D/'PGAM5_three_dataset_upregulated_intersections.zip'
with zipfile.ZipFile(zp,'w',compression=zipfile.ZIP_DEFLATED) as z:
 for f in sorted(D.iterdir()):
  if f.is_file() and f!=zp:z.write(f,f.name)
with zipfile.ZipFile(zp) as z:assert z.testzip() is None
assert len(droplet)==14 and len(smart)==16 and len(both)==5 and len(either)==25
print(json.dumps({k:summary[k] for k in ['triple_Droplet','triple_Smartseq2','strict_all_four_lists','triple_either_GSE140228_platform']},indent=2,ensure_ascii=False))
