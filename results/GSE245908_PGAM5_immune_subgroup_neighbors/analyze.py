"""Fixed RNA gates on mixed Visium spots, not resolved immune cell identities."""
from pathlib import Path
import os,json,gzip,itertools,math,hashlib
os.environ.setdefault('MPLCONFIGDIR','/tmp/GSE245908_immune_mpl')
import numpy as np,pandas as pd,scanpy as sc,scipy.sparse as sp
from scipy.stats import false_discovery_control
from PIL import Image
from io import BytesIO
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent/'GSE245908_spatial_PGAM5'
SAMPLES={'CHC20':'GSM7850822','CHC23':'GSM7850823'};SEED=24590881;B=4999
PANELS={
 'Macrophage_core':['C1QA','C1QB','C1QC','CD68','CSF1R','MERTK','MSR1'],
 'Pan_myeloid':['LST1','TYROBP','FCER1G','AIF1','LYZ','CTSS'],
 'All_T':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2'],
 'T_exhaustion':['PDCD1','LAG3','TIGIT','CTLA4','TOX','CXCL13'],
 'T_activation':['CD69','IL2RA','TNFRSF9','ICOS','CD40LG','IL2'],
 'CD8_cytotoxic_T':['IFNG','PRF1','GZMB','GZMA','CCL5','TBX21'],
 'Naive_memory_T':['CCR7','TCF7','LEF1','IL7R','SELL'],
 'Treg_related':['FOXP3','IL2RA','IKZF2','CTLA4','TIGIT'],
 'NK_related':['NKG7','GNLY','KLRD1','NCR1','KLRF1','FCGR3A'],
 'B_related':['MS4A1','CD79A','CD79B','CD19','CD22'],
 'Plasma_related':['MZB1','JCHAIN','SDC1','IGKC','DERL3'],
 'cDC2_related':['CD1C','FCER1A','CLEC10A','CD1E'],
 'cDC1_related':['CLEC9A','XCR1','BATF3','CADM1','IRF8'],
 'pDC_related':['CLEC4C','IL3RA','TCF4','GZMB','BST2'],
 'Inflammatory_monocyte':['FCN1','VCAN','S100A8','S100A9','LYZ'],
 'Neutrophil_related':['FCGR3B','CSF3R','CXCR2','S100A8','S100A9'],
 'Mast_related':['TPSAB1','TPSB2','CPA3','KIT','MS4A2']}
TARGETS=[k for k in PANELS if k not in ['Macrophage_core','Pan_myeloid']]
ZH={'All_T':'全部T谱系（基线）','T_exhaustion':'耗竭支持T相关','T_activation':'早期活化支持T相关','CD8_cytotoxic_T':'CD8/细胞毒T相关','Naive_memory_T':'CCR7/TCF7初始记忆T相关','Treg_related':'FOXP3/IKZF2调节T相关','NK_related':'KLRD1/GNLY NK相关','B_related':'MS4A1/CD79A B相关','Plasma_related':'MZB1/JCHAIN浆细胞相关','cDC2_related':'CD1C/FCER1A cDC2相关','cDC1_related':'CLEC9A/XCR1 cDC1相关','pDC_related':'CLEC4C/IL3RA pDC相关','Inflammatory_monocyte':'FCN1/VCAN炎症单核相关','Neutrophil_related':'FCGR3B/CSF3R粒细胞相关','Mast_related':'TPSB2/CPA3肥大细胞相关'}
PAIRS=[('SPP1','CD44'),('TGFB1','TGFBR1_TGFBR2'),('IL10','IL10RA_IL10RB'),('PVR','TIGIT'),('LGALS9','HAVCR2'),('CD274','PDCD1'),('CD80','CD28'),('CXCL9','CXCR3'),('CCL2','CCR2'),('CSF1','CSF1R'),('TNF','TNFRSF1A'),('TNF','TNFRSF1B'),('TNFSF13B','TNFRSF13C'),('TNFSF13','TNFRSF17'),('CXCL12','CXCR4')]

def dense(x):return np.asarray(x.toarray() if sp.issparse(x) else x)
def load(s,gsm):
    o=pd.read_csv(BASE/f'{s}_spot_annotations.csv',index_col=0)
    a=sc.read_10x_h5(BASE/'inputs'/s/f'{gsm}_{s}_filtered_feature_bc_matrix.h5');a.var_names_make_unique();a=a[o.index].copy();a.obs=o
    a.layers['counts']=a.X.copy();sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
    a.obsm['spatial']=o[['pixel_col','pixel_row']].to_numpy()
    pref=BASE/'inputs'/s/f'{gsm}_{s}'
    sf=json.loads(gzip.decompress(Path(str(pref)+'_scalefactors_json.json.gz').read_bytes()))
    im=np.asarray(Image.open(BytesIO(gzip.decompress(Path(str(pref)+'_tissue_lowres_image.png.gz').read_bytes()))).convert('RGB'))
    a.uns['spatial']={s:{'images':{'lowres':im},'scalefactors':sf}}
    positions=o[['array_row','array_col']].to_numpy();index={tuple(v):i for i,v in enumerate(positions)};edges=[]
    for i,(r,c) in enumerate(positions):
        for dr,dc in [(0,2),(0,-2),(1,1),(1,-1),(-1,1),(-1,-1)]:
            j=index.get((r+dr,c+dc))
            if j is not None:edges.append((i,j))
    ij=np.asarray(edges);a.obsp['spatial_connectivities']=sp.csr_matrix((np.ones(len(edges)),(ij[:,0],ij[:,1])),shape=(len(a),len(a)))
    return a

def annotate(a):
    o=a.obs;W=a.obsp['spatial_connectivities'];deg=np.asarray(W.sum(1)).ravel()
    coverage=[];cuts={}
    for k,gs in PANELS.items():
        found=[g for g in gs if g in a.var_names];assert len(found)==len(gs),(k,found)
        o['score_'+k]=dense(a[:,gs].X).mean(1);o['detected_'+k]=(dense(a[:,gs].layers['counts'])>0).sum(1)
        coverage.append(dict(sample=o['sample'].iloc[0],program=k,requested=len(gs),present=len(found),genes=';'.join(found)))
    m=(o.score_Macrophage_core>=o.score_Macrophage_core.quantile(.75))&(o.detected_Macrophage_core>=3)&(o.score_Pan_myeloid>=o.score_Pan_myeloid.median())&(deg>0)
    t=(o.score_All_T>=o.score_All_T.quantile(.75))&(o.detected_All_T>=2)&(deg>0)
    pos=m&o.PGAM5_raw.gt(0);cut=float(o.loc[pos,'PGAM5_log1p_CP10k'].quantile(.75));high=pos&o.PGAM5_log1p_CP10k.ge(cut)
    o['source_group']=np.where(m,np.where(high,'PGAM5_high',np.where(pos,'Other_detected','Undetected')),'Outside_mac_gate')
    o['Mac_RNA_enriched']=m;o['gate_All_T']=t;o['spatial_degree']=deg.astype(int)
    def raw(g):return dense(a[:,[g]].layers['counts']).ravel()
    def anchor(gs):return np.any(np.vstack([raw(g)>0 for g in gs]),axis=0)
    for k in TARGETS[1:]:
        lineage=t if k in ['T_exhaustion','T_activation','CD8_cytotoxic_T','Naive_memory_T','Treg_related'] else (deg>0)
        minmarker=3 if k in ['Plasma_related','pDC_related'] else 2
        q=float(o.loc[t,'score_'+k].quantile(.75)) if k.startswith('T_') or k in ['CD8_cytotoxic_T','Naive_memory_T','Treg_related'] else float(o['score_'+k].quantile(.75))
        mask=lineage&(o['score_'+k]>=q)&(o['score_'+k]>0)&(o['detected_'+k]>=minmarker)
        if k=='CD8_cytotoxic_T':mask&=anchor(['CD8A','CD8B'])
        if k=='Treg_related':mask&=raw('FOXP3')>0
        if k=='NK_related':mask&=anchor(['KLRD1','NCR1','KLRF1'])&~t
        if k=='B_related':mask&=anchor(['MS4A1','CD19'])
        if k=='Plasma_related':mask&=anchor(['MZB1','JCHAIN'])
        if k=='cDC2_related':mask&=anchor(['CD1C','CLEC10A'])
        if k=='cDC1_related':mask&=anchor(['CLEC9A','XCR1'])
        if k=='pDC_related':mask&=raw('CLEC4C')>0
        if k=='Inflammatory_monocyte':mask&=anchor(['FCN1','VCAN'])&(o.score_Pan_myeloid>=o.score_Pan_myeloid.median())
        if k=='Neutrophil_related':mask&=anchor(['FCGR3B','CSF3R'])
        if k=='Mast_related':mask&=anchor(['TPSB2','CPA3'])
        o['gate_'+k]=np.asarray(mask,bool);cuts[k]=q
    o['immune_gate_count']=o[['gate_'+k for k in TARGETS[1:]]].sum(1)
    a.uns['PGAM5_high_threshold']=cut;a.uns['target_score_thresholds']=cuts
    return coverage

def stratify(o,ix):
    xy=o[['array_col','array_row']].to_numpy();tile=np.minimum(((xy[ix]-xy.min(0))/(np.ptp(xy,axis=0)+1)*3).astype(int),2);block=tile[:,0]+3*tile[:,1];strata=[];groups=np.zeros(len(ix),int)
    for b in sorted(set(block)):
        members=np.flatnonzero(block==b);depth=o.total_counts.to_numpy()[ix[members]];level=(depth>np.median(depth)).astype(int)
        for lv in [0,1]:
            part=members[level==lv]
            if len(part):strata.append(part);groups[part]=b*2+lv
    return strata,groups

def evaluate(o,eligible,Y,labels,kind,rng,s,targetcounts):
    ix=np.flatnonzero(eligible);h=o.source_group.eq('PGAM5_high').to_numpy()[ix].astype(float);y=Y[ix];n1=h.sum();n0=len(h)-n1
    groups,code=stratify(o,ix);space=math.prod(math.comb(len(g),int(h[g].sum())) for g in groups)
    info=dict(sample=s,kind=kind,n_eligible=len(ix),n_high=int(n1),n_other=int(n0),n_strata=len(groups),configuration_space=str(space),movable_sources=sum(len(g) for g in groups if 0<h[g].sum()<len(g)))
    if n1<3 or n0<3 or space<=1:
        return [dict(**info,target=k,statistic=np.nan,permutation_p=np.nan,status='not_testable',reason='Insufficient source groups / no exchangeable labels',n_target_spots=targetcounts[k]) for k in labels]
    if space<=B:
        choices=[list(itertools.combinations(g,int(h[g].sum()))) for g in groups];mat=np.zeros((space,len(h)))
        for row,parts in enumerate(itertools.product(*choices)):
            for part in parts:mat[row,list(part)]=1
        mode='exact_enumeration'
    else:
        mat=np.tile(h,(B,1))
        for g in groups:mat[:,g]=rng.permuted(mat[:,g],axis=1)
        mode='Monte_Carlo'
    eh=np.zeros(len(h))
    for g in groups:eh[g]=h[g].mean()
    observed=(h@y)/n1-((1-h)@y)/n0;center=(eh@y)/n1-((1-eh)@y)/n0;null=(mat@y)/n1-((1-mat)@y)/n0
    audit_name=kind+('_'+labels[0] if len(labels)==1 else '')
    np.savez_compressed(ROOT/f'{s}_{audit_name}_permutation_audit.npz',source_barcodes=o.index[ix].to_numpy(str),high=h,outcomes=y,labels=np.array(labels),groups=code,assignments=mat.astype(np.uint8),observed=observed,conditional_expected=center,null_deltas=null,mode=np.array(mode))
    rows=[]
    for j,k in enumerate(labels):
        d=dict(**info,target=k,statistic=float(observed[j]),conditional_expected=float(center[j]),excess_over_conditional_expectation=float(observed[j]-center[j]),mean_high=float(y[h>0,j].mean()),mean_other=float(y[h==0,j].mean()),n_target_spots=targetcounts[k],permutations=len(mat),mode=mode,permutation_p=np.nan,status='not_testable',reason='')
        if targetcounts[k]<3:d['reason']='Fewer than 3 target spots'
        elif np.ptp(y[:,j])==0:d['reason']='Constant source-neighbor outcome'
        elif np.ptp(null[:,j])<1e-12:d['reason']='Degenerate conditional null'
        else:
            exceed=(np.abs(null[:,j]-center[j])>=abs(observed[j]-center[j])-1e-12).sum();p=exceed/len(mat) if mode=='exact_enumeration' else (1+exceed)/(len(mat)+1)
            d.update(permutation_p=float(p),status='tested')
        rows.append(d)
    return rows

def run():
    resource=pd.read_csv(ROOT/'source_inputs/omni_resource.csv');lrpairs=[]
    for l,r in PAIRS:
        rr=resource[(resource.source_genesymbol==l)&(resource.target_genesymbol==r)];assert len(rr)
        lrpairs.append(dict(ligand=l,receptor_complex=r,pair=l+'--'+r,database_resources=';'.join(sorted(set(rr.resource))),snapshot='LIANA wheel 1.10.0 omni_resource.csv'))
    pd.DataFrame(lrpairs).to_csv(ROOT/'selected_candidate_LR_pairs.csv',index=False)
    (ROOT/'programs.json').write_text(json.dumps(PANELS,indent=2)+'\n');(ROOT/'subgroup_names.json').write_text(json.dumps(ZH,ensure_ascii=False,indent=2)+'\n')
    config=dict(seed=SEED,permutations=B,primary='30 neighbor-fraction tests (2 sections x 15 RNA gates, includes All-T baseline)',joint='90 planned tests (fraction, binary adjacency, source-not-target sensitivity); missing p filled with 1 only for BH',minimum_source_group=3,minimum_target_spots=3,SPP1_not_required=True,source='Unchanged prior Mac-RNA gate / PGAM5-positive Q75 high vs other detected; zeros descriptive',strata='3x3 full-section coordinate tiles x within-tile candidate-source UMI median; preserve local high counts',null='Exact label enumeration when configuration space <=4999; otherwise 4999 permutations',p='Two-sided deviations about exact conditional expectation, not zero',causal_claim=False)
    (ROOT/'analysis_config.json').write_text(json.dumps(config,indent=2)+'\n')
    alltests=[];counts=[];summaries=[];lrs=[];coverage=[];sources=[]
    for s,gsm in SAMPLES.items():
        a=load(s,gsm);coverage+=annotate(a);o=a.obs;W=a.obsp['spatial_connectivities'];deg=o.spatial_degree.to_numpy();groups=o.source_group.to_numpy();high=groups=='PGAM5_high';positive=np.isin(groups,['PGAM5_high','Other_detected']);targets={k:o['gate_'+k].to_numpy(bool) for k in TARGETS};targetcounts={k:int(v.sum()) for k,v in targets.items()};rng=np.random.default_rng(SEED+(s=='CHC23'))
        neighborcount=np.column_stack([np.asarray(W@targets[k].astype(int)).ravel() for k in TARGETS]);fraction=neighborcount/np.where(deg[:,None]>0,deg[:,None],np.nan)
        alltests+=evaluate(o,positive,fraction,TARGETS,'neighbor_fraction',rng,s,targetcounts)
        alltests+=evaluate(o,positive,(neighborcount>0).astype(float),TARGETS,'neighbor_present',rng,s,targetcounts)
        for j,k in enumerate(TARGETS):alltests+=evaluate(o,positive&~targets[k],fraction[:,[j]],[k],'source_not_target_fraction',rng,s,targetcounts)
        src=o.loc[o.Mac_RNA_enriched,['source_group','PGAM5_raw','PGAM5_log1p_CP10k','spatial_degree','total_counts','pct_mt','array_row','array_col']].copy()
        for j,k in enumerate(TARGETS):
            src[k+'_neighbor_count']=neighborcount[o.Mac_RNA_enriched,j];src[k+'_neighbor_fraction']=fraction[o.Mac_RNA_enriched,j]
            for group in ['PGAM5_high','Other_detected','Undetected']:
                keep=groups==group;n=int(keep.sum());nn=neighborcount[keep,j]
                counts.append(dict(sample=s,target=k,target_name_zh=ZH[k],source_group=group,n_sources=n,n_target_spots=targetcounts[k],sources_with_neighbor=int((nn>0).sum()),source_adjacency_pct=100*(nn>0).mean() if n else np.nan,directed_edges=int(nn.sum()),unique_neighbor_targets=int(np.sum((np.asarray(W[keep].sum(0)).ravel()>0)&targets[k])),mean_neighbor_fraction=float(fraction[keep,j].mean()) if n else np.nan))
        src.to_csv(ROOT/f'{s}_source_neighbor_subgroups.csv',index_label='barcode',float_format='%.17g')
        for k in TARGETS:
            mask=targets[k];summaries.append(dict(sample=s,target=k,target_name_zh=ZH[k],n_spots=int(mask.sum()),multi_gate_spots=int((mask&o.immune_gate_count.gt(1)).sum()),score_mean=float(o.loc[mask,'score_'+k].mean()) if mask.any() else np.nan))
        g=W.tocoo();keep=o.Mac_RNA_enriched.to_numpy()[g.row]&np.any(np.column_stack(list(targets.values())),axis=1)[g.col];ri,ci=g.row[keep],g.col[keep]
        e=pd.DataFrame(dict(source_barcode=o.index[ri],target_barcode=o.index[ci],source_group=groups[ri]))
        for k in TARGETS:e['target_'+k]=targets[k][ci]
        for l,r in PAIRS:
            ligand=dense(a[:,[l]].layers['counts']).ravel()/o.total_counts.to_numpy()*1e4;receptor=(dense(a[:,r.split('_')].layers['counts'])/o.total_counts.to_numpy()[:,None]*1e4).min(1)
            score=np.sqrt(ligand[ri]*receptor[ci]);name=l+'--'+r;e[name]=score
            for k in TARGETS:
                for group in ['PGAM5_high','Other_detected','Undetected']:
                    keep=(groups[ri]==group)&targets[k][ci];supported=keep&(score>0)
                    lrs.append(dict(sample=s,target=k,source_group=group,pair=name,eligible_edges=int(keep.sum()),unique_sources=len(set(ri[keep])),unique_targets=len(set(ci[keep])),codetected_edges=int(supported.sum()),supported_sources=len(set(ri[supported])),supported_targets=len(set(ci[supported])),codetection_pct=float(100*supported.sum()/keep.sum()) if keep.any() else np.nan,mean_geometric_CP10k=float(score[keep].mean()) if keep.any() else np.nan))
        e.to_csv(ROOT/f'{s}_directed_immune_LR_edges.csv.gz',index=False)
        o.to_csv(ROOT/f'{s}_spot_annotations.csv.gz',index_label='barcode',float_format='%.17g')
        panel=list(dict.fromkeys(['PGAM5','MYO19','CD8A','CD8B']+[g for v in PANELS.values() for g in v]+[g for l,r in PAIRS for g in [l]+r.split('_')]))
        p=a[:,panel].copy();p.write_h5ad(ROOT/f'{s}_immune_marker_subset.h5ad',compression='gzip')
        sources.append(dict(sample=s,qc_spots=len(a),Mac_RNA_spots=int(o.Mac_RNA_enriched.sum()),PGAM5_high=int(high.sum()),other_detected=int((groups=='Other_detected').sum()),undetected=int((groups=='Undetected').sum()),PGAM5_high_threshold=a.uns['PGAM5_high_threshold'],immune_union_spots=int(np.any(np.column_stack(list(targets.values())),axis=1).sum())))
        print('Analyzed',s,sources[-1],targetcounts,flush=True)
    tests=pd.DataFrame(alltests);assert len(tests)==90
    tests['BH_q_90_planned_tests']=false_discovery_control(tests.permutation_p.fillna(1).to_numpy());tests['primary_BH_q_30']=np.nan
    prim=tests.kind.eq('neighbor_fraction');tests.loc[prim,'primary_BH_q_30']=false_discovery_control(tests.loc[prim,'permutation_p'].fillna(1).to_numpy())
    tests.to_csv(ROOT/'spatial_conditional_tests.csv',index=False);pd.DataFrame(counts).to_csv(ROOT/'subgroup_neighborhood_counts.csv',index=False);pd.DataFrame(summaries).to_csv(ROOT/'target_gate_summary.csv',index=False);pd.DataFrame(lrs).to_csv(ROOT/'candidate_LR_by_immune_subgroup.csv',index=False);pd.DataFrame(sources).to_csv(ROOT/'source_group_summary.csv',index=False);pd.DataFrame(coverage).to_csv(ROOT/'program_coverage.csv',index=False)
    manifest=[]
    for p in [BASE/record['path'] for record in json.loads((ROOT/'source_inputs/base_requirements.json').read_text())]+[ROOT/'source_inputs/omni_resource.csv']:
        manifest.append(dict(path=str(p.relative_to(ROOT.parent)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    (ROOT/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

if __name__=='__main__':run()
