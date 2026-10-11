"""Post hoc descriptive marker-specificity checks; no added inferential tests."""
import numpy as np,pandas as pd,anndata as ad
from analyze import ROOT,SAMPLES
rows=[];sensitivity=[]
for s in SAMPLES:
    a=ad.read_h5ad(ROOT/f'{s}_immune_marker_subset.h5ad');o=a.obs;G=a.obsp['spatial_connectivities'];high=o.source_group.eq('PGAM5_high').to_numpy()
    def raw(g):return a[:,[g]].layers['counts'].toarray().ravel()
    for g in ['FOXP3','IKZF2','IL2RA','CTLA4','CD8A','CD8B','FCN1','VCAN','FCGR3B','CSF3R','MZB1','JCHAIN']:
        z=raw(g);t=o.gate_All_T.to_numpy(bool)
        rows.append(dict(sample=s,gene=g,n_spots=len(a),detected_spots=int((z>0).sum()),detected_pct=100*(z>0).mean(),T_gate_spots=int(t.sum()),T_gate_detected=int((t&(z>0)).sum()),median_raw=float(np.median(z)),maximum_raw=int(z.max())))
    gates={
      'Treg_add_IKZF2_and_IL2RA_or_CTLA4':o.gate_Treg_related.to_numpy(bool)&(raw('IKZF2')>0)&((raw('IL2RA')>0)|(raw('CTLA4')>0)),
      'CD8_cytotoxic_require_CD8B':o.gate_CD8_cytotoxic_T.to_numpy(bool)&(raw('CD8B')>0),
      'Monocyte_require_FCN1':o.gate_Inflammatory_monocyte.to_numpy(bool)&(raw('FCN1')>0),
      'Neutrophil_require_FCGR3B':o.gate_Neutrophil_related.to_numpy(bool)&(raw('FCGR3B')>0)}
    for name,m in gates.items():
        n=np.asarray(G@m.astype(int)).ravel();edges=G.tocoo();keep=high[edges.row]&m[edges.col]
        sensitivity.append(dict(sample=s,rule=name,target_spots=int(m.sum()),PGAM5_high_sources_with_neighbor=int((high&(n>0)).sum()),directed_edges=int(keep.sum()),unique_targets=len(set(edges.col[keep])),post_hoc_descriptive_only=True))
pd.DataFrame(rows).to_csv(ROOT/'marker_detection_specificity.csv',index=False)
pd.DataFrame(sensitivity).to_csv(ROOT/'strict_marker_descriptive_sensitivity.csv',index=False)
print(pd.DataFrame(sensitivity).to_string(index=False))
