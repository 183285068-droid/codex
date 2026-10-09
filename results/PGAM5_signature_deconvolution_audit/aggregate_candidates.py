from pathlib import Path
import pandas as pd,json,hashlib,itertools,os
R=Path(__file__).resolve().parent;B=Path(os.environ.get('PGAM5_REPOSITORY_RESULTS_ROOT','/workspace/codex/results'))
files={
'GSE151530':B/'PGAM5_HCC/HCC_PGAM5_FDR005_log2FC1_upregulated.csv',
'GSE149614':B/'GSE149614_PGAM5_primaryT/GSE149614_primaryT_PGAM5_upregulated.csv',
'GSE189903':B/'GSE189903_HCC_tumor_PGAM5/GSE189903_HCC_tumor_PGAM5_upregulated.csv',
'GSE242889':B/'GSE242889_HCC_tumor_PGAM5/GSE242889_HCC_tumor_PGAM5_upregulated.csv',
'GSE202642':B/'GSE202642_HCC_tumor_PGAM5/GSE202642_HCC_tumor_PGAM5_upregulated.csv',
'GSE140228_Droplet':B/'GSE140228_allTumor_PGAM5/Droplet/upregulated.csv',
'GSE140228_Smartseq2':B/'GSE140228_allTumor_PGAM5/Smartseq2/upregulated.csv',
'GSE125449':B/'GSE125449_HCC_PGAM5/GSE125449_HCC_PGAM5_upregulated.csv',
'GSE146115':B/'GSE146115_HCC_PGAM5/GSE146115_HCC_tumor_PGAM5_upregulated.csv'}
sets={};tables=[];hashes={}
for name,p in files.items():
 d=pd.read_csv(p,index_col=0);assert d.gene.notna().all() and d.FDR.lt(.05).all() and d.log2FC.ge(1).all() and ~d.gene.eq('PGAM5').any()
 sets[name]=set(d.gene);e=d[['gene','log2FC','FDR']].copy();e['dataset']=name;tables.append(e.reset_index(drop=True));hashes[name]={'path':str(p.relative_to(B)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'rows':len(d),'unique_symbols':d.gene.nunique()}
long=pd.concat(tables,ignore_index=True);long.to_csv(R/'candidate_evidence_long.csv',index=False)
primary=['GSE151530','GSE149614','GSE189903','GSE242889','GSE202642'];genes=sorted(set.union(*(sets[n] for n in primary)))
a=pd.DataFrame({'gene':genes})
for name in sets:a[name]=a.gene.isin(sets[name])
a['primary_support_count']=a[primary].sum(axis=1);a['GSE140228_any_technology']=a[['GSE140228_Droplet','GSE140228_Smartseq2']].any(axis=1)
a=a.sort_values(['primary_support_count','gene'],ascending=[False,True]);a.to_csv(R/'candidate_recurrence.csv',index=False)
a[a.primary_support_count.ge(3)].to_csv(R/'primary_consensus_at_least3.csv',index=False)
pairs=[{'dataset1':x,'dataset2':y,'intersection':len(sets[x]&sets[y]),'genes':';'.join(sorted(sets[x]&sets[y]))} for x,y in itertools.combinations(primary,2)];pd.DataFrame(pairs).to_csv(R/'primary_pairwise_intersections.csv',index=False)
summary={'primary_datasets':primary,'sources':hashes,'all5_intersection':sorted(set.intersection(*(sets[n] for n in primary))),'primary_support_histogram':a.primary_support_count.value_counts().sort_index().to_dict(),'at_least3_genes':a.loc[a.primary_support_count.ge(3),'gene'].tolist(),'excluded_from_primary':'GSE125449 overlaps GSE151530 and has 4 positive TAMs; GSE146115 has 17 positive inferred macrophages and zero non-PGAM5 exact-detection candidates; GSE140228 includes one CC tumor and two partly overlapping technologies. Supplementary evidence shown, not independent votes; GSE154906 cancelled.'}
(R/'consensus_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
