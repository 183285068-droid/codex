"""Descriptive marker directions across the seven user-confirmed patients."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse

R=Path(__file__).resolve().parent

def main():
    obs=pd.read_csv(R/'cell_annotations.csv.gz')
    genes=pd.read_csv(R/'portable_feature_genes.csv').gene.to_numpy()
    raw=sparse.load_npz(R/'portable_raw_feature_counts.npz')
    norm=raw.multiply((1e4/obs.total_counts.to_numpy())[:,None]).tocsr()
    markers=pd.read_csv(R/'cluster_descriptive_markers.csv',dtype={'cluster':str})
    markers=markers[markers.geometry.eq('primary') & markers['rank'].le(10)]
    rows=[]
    for _,marker in markers.iterrows():
        j=int(np.flatnonzero(genes==marker.gene)[0])
        for patient in sorted(obs.patient_id.unique()):
            ii=obs.patient_id.eq(patient).to_numpy() & obs.primary_cluster.astype(str).eq(marker.cluster).to_numpy()
            oo=obs.patient_id.eq(patient).to_numpy() & ~obs.primary_cluster.astype(str).eq(marker.cluster).to_numpy()
            n1,n0=int(ii.sum()),int(oo.sum())
            mean1=float(norm[ii,j].mean()) if n1 else np.nan
            mean0=float(norm[oo,j].mean()) if n0 else np.nan
            rows.append({'cluster':marker.cluster,'gene':marker.gene,'rank':marker['rank'],'patient_id':patient,
                         'cluster_cells':n1,'other_cells':n0,'eligible':n1>=20 and n0>=20,
                         'cluster_mean_CP10k':mean1,'other_mean_CP10k':mean0,
                         'descriptive_log2_ratio_pseudocount0p05':np.log2((mean1+.05)/(mean0+.05)),
                         'cluster_detection_fraction':float((raw[ii,j]>0).mean()) if n1 else np.nan,
                         'other_detection_fraction':float((raw[oo,j]>0).mean()) if n0 else np.nan})
    table=pd.DataFrame(rows)
    table.to_csv(R/'cluster_marker_directions_by_patient.csv',index=False)
    summaries=[]
    for (cluster,gene),t in table.groupby(['cluster','gene'],observed=True):
        eligible=t[t.eligible]
        summaries.append({'cluster':cluster,'gene':gene,'rank':int(t['rank'].iloc[0]),
            'evaluable_patients':len(eligible),'patients_mean_expression_enriched':int(eligible.descriptive_log2_ratio_pseudocount0p05.gt(0).sum()),
            'enriched_direction_fraction':float(eligible.descriptive_log2_ratio_pseudocount0p05.gt(0).mean()),
            'median_patient_descriptive_log2_ratio':float(eligible.descriptive_log2_ratio_pseudocount0p05.median()),
            'recurs_in_at_least3_patients_and_ge70pct_direction':bool(len(eligible)>=3 and eligible.descriptive_log2_ratio_pseudocount0p05.gt(0).mean()>=.7)})
    pd.DataFrame(summaries).to_csv(R/'cluster_marker_patient_consistency.csv',index=False)
    print('Saved descriptive marker directions for',len(summaries),'cluster/gene pairs across seven confirmed patients.')

if __name__=='__main__':main()
