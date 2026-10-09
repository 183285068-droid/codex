from pathlib import Path
import pandas as pd,json,hashlib,itertools,zipfile
D=Path(__file__).resolve().parent;R=D.parent
paths={'GSE242889':R/'GSE242889_HCC_tumor_PGAM5/GSE242889_HCC_tumor_PGAM5_upregulated.csv','GSE149614':R/'GSE149614_PGAM5_primaryT/GSE149614_primaryT_PGAM5_upregulated.csv','GSE151530':R/'PGAM5_HCC/HCC_PGAM5_FDR005_log2FC1_upregulated.csv'}
tables={k:pd.read_csv(p,index_col=0) for k,p in paths.items()}
sets={}
for k,t in tables.items():
 assert t.gene.is_unique and t.gene.notna().all() and not t.gene.eq('PGAM5').any()
 assert t.FDR.lt(.05).all() and t.log2FC.ge(1).all()
 assert t[['fraction_PGAM5_positive','fraction_PGAM5_undetected']].max(axis=1).ge(.1).all()
 sets[k]=set(t.gene)
def export(filename,genes,keys):
 rows=[]
 for gene in sorted(genes):
  row={'gene':gene}
  for key in keys:
   t=tables[key].set_index('gene').loc[gene]
   for c in ['log2FC','FDR','fraction_PGAM5_positive','fraction_PGAM5_undetected']:row[key+'_'+c]=float(t[c])
  rows.append(row)
 pd.DataFrame(rows).to_csv(D/filename,index=False)
common=set.intersection(*sets.values())
export('three_dataset_shared_upregulated.csv',common,list(tables))
long=[]
for gene in sorted(common):
 for key,t in tables.items():
  row=t[t.gene.eq(gene)].iloc[0]
  long.append({'gene':gene,'dataset':key,**{c:float(row[c]) for c in ['log2FC','FDR','fraction_PGAM5_positive','fraction_PGAM5_undetected']}})
pd.DataFrame(long).to_csv(D/'three_dataset_shared_upregulated_long.csv',index=False)
pairs=[]
for a,b in itertools.combinations(sets,2):
 pair=sets[a]&sets[b]
 export('pair_'+a+'_'+b+'.csv',pair,[a,b])
 pairs.append({'dataset_A':a,'dataset_B':b,'count':len(pair),'genes':sorted(pair)})
summary={'criteria':'FDR<0.05, log2FC>=1, detection>=10% in either group; PGAM5 excluded. Official symbol intersection, no alias conversion.','source_inputs':{k:{'path':str(p.relative_to(R)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'upregulated_genes':len(sets[k])} for k,p in paths.items()},'triple_count':len(common),'triple_genes':sorted(common),'pairwise':pairs,'scope':{'GSE242889':'All five HCC tumor samples regardless of MVI; marker-inferred C1Q-rich macrophages','GSE149614':'Primary HCC tumor samples ending T','GSE151530':'All selected HCC samples'},'interpretation':'Descriptive shared association. Shared genes cannot establish a validated or macrophage-specific signature; no enrichment test or independent patient validation performed.'}
(D/'intersection_summary.json').write_text(json.dumps(summary,indent=2))
text='''# GSE242889、GSE149614、GSE151530上调基因交集

沿用既定标准：FDR <0.05、log2FC ≥1、至少一组检出率 ≥10%，排除定义分组的PGAM5。按gene symbol匹配，不重新拟合差异分析、不做别名转换。

GSE242889为全部5个HCC肿瘤样本（不按MVI筛选），采用marker推定C1Q较高巨噬细胞；GSE149614为原发HCC肿瘤T样本；GSE151530为此前选择的全部HCC样本。分别有82、512、740个上调候选。

'''
text+=f'三者共有**{len(common)}个共同上调基因**：'+ '、'.join(sorted(common))+'。\n\n'
text+='| 基因 | 数据集 | log2FC | FDR | PGAM5检出组中该基因检出率 | PGAM5未检出组中该基因检出率 |\n|---|---|---:|---:|---:|---:|\n'
for row in long:text+=f"| {row['gene']} | {row['dataset']} | {row['log2FC']:.4f} | {row['FDR']:.4g} | {row['fraction_PGAM5_positive']*100:.2f}% | {row['fraction_PGAM5_undetected']*100:.2f}% |\n"
text+='\n两两交集：\n\n'
for p in pairs:text+=f"- {p['dataset_A']} ∩ {p['dataset_B']}：{p['count']}个。\n"
text+='''
共同上调是描述性一致性，可以作为进一步验证的候选；这些共同基因不能直接构成已验证的PGAM5⁺巨噬细胞signature，也未证明巨噬细胞特异性、患者间可重复性或因果机制。继承此前合并细胞分析、不做患者配对及深度/样本协变量校正的条件。未计算交集富集显著性或独立患者验证。GSE242889未提供作者细胞类型标签，本次使用探索性marker推定巨噬细胞，且存在广泛ALB RNA背景，交集继承这些注释及混合RNA的不确定性。

CSV保留各来源的log2FC、FDR及两组检出率；JSON保留输入路径和SHA256。calculate_intersections.py从仓库已有上调表重新计算。three_dataset_shared_upregulated_long.csv适合在Excel查看，每行对应一个数据集。
'''
(D/'结果说明.md').write_text(text)
zp=D/'PGAM5_three_dataset_upregulated_intersections.zip'
with zipfile.ZipFile(zp,'w',compression=zipfile.ZIP_DEFLATED) as z:
 for f in sorted(D.iterdir()):
  if f.is_file() and f!=zp:z.write(f,f.name)
with zipfile.ZipFile(zp) as z:assert z.testzip() is None
print(json.dumps({'triple':sorted(common),'pairwise':[(p['dataset_A'],p['dataset_B'],p['count']) for p in pairs]},ensure_ascii=False))
