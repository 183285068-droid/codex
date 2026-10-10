from pathlib import Path
import numpy as np,pandas as pd,json,itertools
from scipy.stats import spearmanr,fisher_exact
from statsmodels.stats.multitest import multipletests
from scTenifold.core._networks import manifold_alignment,d_regulation
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;protocol=json.loads((R/'protocol.json').read_text());seeds=protocol['seeds'];all_pairs=[];summaries=[];controlrows=[]
import sys
for scenario in (sys.argv[1:] or ['precision_only','ensemble']):
 results={seed:pd.read_csv(R/scenario/str(seed)/'differential_regulation.csv').set_index('Gene').drop(index='PGAM5') for seed in seeds};genes=results[seeds[0]].index
 distance=pd.DataFrame({str(k):v.Distance.reindex(genes) for k,v in results.items()});q=pd.DataFrame({str(k):v['adjusted p-value'].reindex(genes) for k,v in results.items()});frequency=q.lt(.05).sum(1);consensus=pd.DataFrame({'gene':genes,'significant_seed_count':frequency,'median_distance':distance.median(1),'median_model_q':q.median(1),'worst_model_q':q.max(1)});consensus['median_rank']=distance.rank(ascending=False).median(1);consensus=consensus.sort_values(['significant_seed_count','median_rank'],ascending=[False,True]);consensus.to_csv(R/(scenario+'_consensus.csv'),index=False)
 pairs=[]
 for s1,s2 in itertools.combinations(seeds,2):
  top1=set(results[s1].nlargest(50,'Distance').index);top2=set(results[s2].nlargest(50,'Distance').index);sig1=set(results[s1].index[results[s1]['adjusted p-value']<.05]);sig2=set(results[s2].index[results[s2]['adjusted p-value']<.05]);pairs.append(dict(scenario=scenario,seed1=s1,seed2=s2,distance_rho=float(spearmanr(distance[str(s1)],distance[str(s2)]).statistic),top50_Jaccard=len(top1&top2)/len(top1|top2),significant_Jaccard=len(sig1&sig2)/len(sig1|sig2) if sig1|sig2 else np.nan))
 all_pairs+=pairs;pair=pd.DataFrame(pairs)
 # Same WT network alignment must not create signal; all five models checked.
 for seed in seeds:
  f=R/scenario/str(seed);cache=f/'no_KO_control.csv'
  if cache.exists():null=pd.read_csv(cache)
  else:
   d=np.load(f/'WT_network.npz');wt=pd.DataFrame(d['weights'],index=d['genes'],columns=d['genes']);np.random.seed(seed)
   with threadpool_limits(limits=4):null=d_regulation(manifold_alignment(wt,wt,d=30),ko_genes=['PGAM5'])
   null.to_csv(cache,index=False)
  controlrows.append(dict(scenario=scenario,seed=seed,no_KO_significant=int(null['adjusted p-value'].lt(.05).sum()),max_null_distance=float(null.Distance.max())))
 nullpass=all(r['no_KO_significant']==0 for r in controlrows if r['scenario']==scenario);metrics=pd.DataFrame([json.loads((R/scenario/str(seed)/'metrics.json').read_text()) for seed in seeds]);converged=metrics.loc[metrics.scenario.eq(scenario),'tensor_converged'].all();passed=bool(pair.distance_rho.min()>=.8 and pair.top50_Jaccard.min()>=.5 and nullpass and converged)
 summaries.append(dict(scenario=scenario,min_rho=float(pair.distance_rho.min()),median_rho=float(pair.distance_rho.median()),min_top50_Jaccard=float(pair.top50_Jaccard.min()),median_top50_Jaccard=float(pair.top50_Jaccard.median()),genes_significant_ge4_seeds=int(frequency.ge(4).sum()),genes_significant_all5=int(frequency.eq(5).sum()),no_KO_controls_pass=nullpass,all_tensor_converged=bool(converged),passes_prespecified_stability=passed))
 # Program exploration based on recurring hits; preserve all original PGAM5 memberships.
 sets=json.loads((R/'inputs/program_definitions.json').read_text());bg=set(genes);hit=set(frequency.index[frequency.ge(4)]);rows=[]
 for term,gs in sets.items():
  g=set(gs)&bg
  if not g:continue
  h=g&hit;_,p=fisher_exact([[len(h),len(g-h)],[len(hit-g),len(bg-g-hit)]],alternative='greater');rows.append(dict(program=term,tested_non_target_genes=len(g),recurrent_hits=len(h),p_value=p,genes=';'.join(sorted(h))))
 t=pd.DataFrame(rows);t['BH_q']=multipletests(t.p_value,method='fdr_bh')[1];t.sort_values('BH_q').to_csv(R/(scenario+'_recurrent_program_enrichment.csv'),index=False)
 fig,axs=plt.subplots(1,2,figsize=(12,6));c=consensus.head(20).iloc[::-1];axs[0].barh(c.gene,c.significant_seed_count,color='#397bab');axs[0].set_xlim(0,5.5);axs[0].set_xlabel('Seeds with model FDR < 0.05 (out of 5)');axs[0].set_title(scenario+' | recurrent-ranked candidates');matrix=np.eye(5)
 for row in pairs:i=seeds.index(row['seed1']);j=seeds.index(row['seed2']);matrix[i,j]=matrix[j,i]=row['distance_rho']
 im=axs[1].imshow(matrix,vmin=0,vmax=1,cmap='viridis');axs[1].set_xticks(range(5),seeds,rotation=45);axs[1].set_yticks(range(5),seeds);axs[1].set_title('Downstream distance rank correlations')
 for i in range(5):
  for j in range(5):axs[1].text(j,i,f'{matrix[i,j]:.2f}',ha='center',va='center',color='white' if matrix[i,j]<.65 else 'black')
 fig.colorbar(im,ax=axs[1],shrink=.7);fig.tight_layout();fig.savefig(R/(scenario+'_stability.png'),dpi=180);fig.savefig(R/(scenario+'_stability.pdf'));plt.close(fig)
pd.DataFrame(all_pairs).to_csv(R/'pairwise_stability.csv',index=False);pd.DataFrame(controlrows).to_csv(R/'no_KO_control_summary.csv',index=False);pd.DataFrame(summaries).to_csv(R/'stability_summary.csv',index=False);(R/'summary.json').write_text(json.dumps(summaries,indent=2));print(json.dumps(summaries,indent=2),flush=True)
