from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,subprocess
R=Path(__file__).resolve().parent; base=Path('/workspace/codex/results'); old=base/'GSE202642_PGAM5_subpopulation_assessment'; src=base/'GSE202642_HCC_tumor_PGAM5'
panels={'IFN_gamma_response':'STAT1 IRF1 GBP1 GBP2 GBP5 SOCS1 IDO1','T_cell_recruitment':'CXCL9 CXCL10 CXCL11 CCL5','Type1_effector_support':'IL12A IL12B IL18 CD40 IL15','MHC_II_presentation':'HLA-DRA HLA-DRB1 HLA-DPA1 HLA-DPB1 HLA-DQA1 HLA-DQB1 CD74 CIITA CTSS','MHC_I_processing':'HLA-A HLA-B HLA-C B2M TAP1 TAP2 PSMB8 PSMB9 NLRC5','Inhibitory_signaling_RNA':'IL10 TGFB1 CD274 PDCD1LG2 IDO1 LGALS9 VSIR','Lipid_associated_TAM':'SPP1 TREM2 APOE APOC1 LPL LGALS3 GPNMB CD9 CTSB','Matrix_remodeling':'SPP1 MMP9 MMP12 CTSB CTSL PLAUR FN1 TIMP1','Inflammatory_NFkB':'IL1B TNF NFKBIA NFKBIZ TNFAIP3 PTGS2 CCL3 CCL4','Neutrophil_recruitment':'CXCL1 CXCL2 CXCL3 CXCL5 CXCL8','Oxidative_stress_response':'HMOX1 NQO1 TXNRD1 GCLM SOD2 SRXN1','Proliferation':'MKI67 TOP2A CENPF CDK1 TYMS RRM2 PCLAF BIRC5'}
panels={k:['PGAM5']+v.split() for k,v in panels.items()};(R/'program_definitions.json').write_text(json.dumps(panels,indent=2))
a=pd.read_csv(old/'cell_annotations.csv.gz',index_col=0); names=pd.read_csv(old/'functional_subgroup_names.csv'); a=a.merge(names[['primary_cluster','functional_name_EN','functional_name_CN','mixed_RNA_flag']],on='primary_cluster',how='left').set_axis(a.index)
g=pd.read_csv(src/'GSE202642_features.tsv.gz',sep='\t',header=None)[1];b=pd.read_csv(src/'GSE202642_barcodes.tsv.gz',header=None)[0]
selected=sorted(set(sum(panels.values(),[])));assert set(selected)<=set(g)
# Duplicate source gene symbols are summed, as in the original harmonization.
pd.DataFrame([(i+1,selected.index(v)+1) for i,v in enumerate(g) if v in selected]).to_csv(R/'gene_map.tsv',sep='\t',index=False,header=False)
pd.DataFrame([(i+1,a.index.get_loc(v)+1) for i,v in enumerate(b) if v in a.index]).to_csv(R/'cell_map.tsv',sep='\t',index=False,header=False)
script='BEGIN{FS=" "} FILENAME==ARGV[1]{gm[$1]=$2;next} FILENAME==ARGV[2]{cm[$1]=$2;next} /^%/{next} !header{header=1;next} ($2 in cm){tot[cm[$2]]+=$3;if($1 in gm)print cm[$2],gm[$1],$3} END{for(c in tot)print c,0,tot[c]}'
with (R/'selected_counts.tsv').open('w') as out:
 p=subprocess.Popen(['gzip','-cd',str(R/'GSE202642_matrix.mtx.gz')],stdout=subprocess.PIPE)
 subprocess.run(['awk',script,str(R/'gene_map.tsv'),str(R/'cell_map.tsv'),'-'],stdin=p.stdout,stdout=out,check=True);p.stdout.close();assert p.wait()==0
v=pd.read_csv(R/'selected_counts.tsv',sep=' ',header=None).to_numpy();x=np.zeros((len(a),len(selected)),dtype=np.int32);sub=v[v[:,1]>0];np.add.at(x,(sub[:,0]-1,sub[:,1]-1),sub[:,2]);t=v[v[:,1]==0];tot=np.zeros(len(a));tot[t[:,0]-1]=t[:,2];assert np.array_equal(tot,a.total_counts);assert np.array_equal(x[:,selected.index('PGAM5')],a.PGAM5_counts)
np.savez_compressed(R/'program_raw_counts.npz',counts=x);pd.Series(selected,name='gene').to_csv(R/'program_genes.csv',index=False);a.to_csv(R/'cell_metadata.csv.gz');names.to_csv(R/'functional_subgroup_names.csv',index=False)
h=hashlib.sha256();
with (R/'GSE202642_matrix.mtx.gz').open('rb') as f:
 for chunk in iter(lambda:f.read(8388608),b''):h.update(chunk)
assert h.hexdigest()=='855789bdd50f129a4131d84e285d0b52a1b481644accc747acc6f99dcc28e2e4'
(R/'input_validation.json').write_text(json.dumps({'matrix_sha256':h.hexdigest(),'cells':len(a),'PGAM5_detected':int((a.PGAM5_counts>0).sum()),'total_counts_exact_match':True,'PGAM5_exact_match':True,'all_program_genes_present':True,'PGAM5_in_every_program':True},indent=2));print('Prepared',x.shape,flush=True)
