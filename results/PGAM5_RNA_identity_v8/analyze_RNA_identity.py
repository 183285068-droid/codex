"""RNA-only annotation evidence and training-only competitive gene screening.

Existing lineage/RNA memberships are preserved. Corroboration is descriptive.
No protein observations, patient-paired DE, or depth adjustments are used.
"""
from pathlib import Path
import hashlib
import json
import math
import numpy as np
import pandas as pd
from scipy import sparse

R = Path(__file__).resolve().parent
P = json.loads((R/'protocol.json').read_text())
S = Path(P['source_root'])
COHORTS = P['primary_cohorts']
PANELS = {
    'macrophage_core': 'CD68 CD163 CSF1R C1QA C1QB C1QC MERTK MSR1'.split(),
    'broad_myeloid': 'LST1 TYROBP FCER1G AIF1 SPI1 CTSS'.split(),
    'monocyte': 'FCN1 S100A8 S100A9 VCAN'.split(),
    'DC': 'CD1C FCER1A CLEC10A CLEC9A XCR1'.split(),
    'pDC': 'CLEC4C IL3RA TCF4 GZMB JCHAIN'.split(),
    'epithelial_hepatocyte': 'EPCAM KRT8 KRT18 KRT19 ALB APOA1 TTR'.split(),
    'lymphoid': 'CD3D CD3E TRAC CD79A MS4A1'.split(),
    'cycle_diagnostic': 'MKI67 TOP2A UBE2C CENPF CDK1 BIRC5'.split(),
    'stress_diagnostic': 'FOS FOSB JUN JUNB ATF3 EGR1 IER2 IER3 PPP1R15A HSPA1A HSPA1B'.split(),
}
STRESS = set(PANELS['stress_diagnostic'])
CYCLE = set((S/'cell_cycle_genes.txt').read_text().split())
G = pd.Index(pd.read_csv(S/'common_genes.csv').gene)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for part in iter(lambda:f.read(1024*1024),b''):
            h.update(part)
    return h.hexdigest()


def annotate(obs, raw, genes):
    assert len(obs) == raw.shape[0] and obs.cell_id.is_unique
    pg = raw[:,genes.get_loc('PGAM5')].toarray().ravel()
    assert np.array_equal(pg,obs.PGAM5_counts)
    out = obs.copy()
    macro = out.harmonized_group.eq('Macrophage').to_numpy()
    out['macrophage_identity'] = np.where(macro,'Existing_macrophage_annotation','Other_existing_lineage')
    out['PGAM5_RNA_status'] = np.where(pg>0,'RNA_detected','RNA_undetected')
    out['PGAM5_macrophage_annotation'] = np.where(macro,np.where(pg>0,
        'Macrophage_PGAM5_RNA_detected','Macrophage_PGAM5_RNA_undetected'),'Not_macrophage')
    for panel, symbols in PANELS.items():
        measured = [x for x in symbols if x in genes]
        x = raw[:,genes.get_indexer(measured)].toarray().astype(float)
        out[panel+'_measured_genes'] = len(measured)
        out[panel+'_detected_genes'] = (x>0).sum(1)
        out[panel+'_mean_log1p_CP10k'] = np.log1p(x*10000/obs.total_counts.to_numpy()[:,None]).mean(1) if measured else np.nan
    core = out.macrophage_core_detected_genes.ge(2).to_numpy()
    broad = out.broad_myeloid_detected_genes.ge(1).to_numpy()
    out['macrophage_RNA_marker_corroboration'] = core & broad
    out['other_lineage_RNA_review_flag'] = out.epithelial_hepatocyte_detected_genes.ge(2) | out.lymphoid_detected_genes.ge(2)
    out['RNA_lineage_evidence_level'] = np.where(macro,np.where(core&broad,
        'Existing_macrophage_with_marker_corroboration','Existing_macrophage_low_marker_detection'),
        'Other_lineage_marker_evidence_only')
    out['PGAM5_detected_macrophage_with_marker_corroboration'] = macro & (pg>0) & core & broad
    out['PGAM5_related_state'] = np.where(macro,'Unresolved_no_validated_RNA_state','Not_applicable')
    out['RNA_state_validation_passed'] = False
    out['usable_for_validated_TCGA_PGAM5_cell_abundance'] = False
    return out


def main():
    provenance, coverage, marker_rows = [],[],[]
    for c in COHORTS+P['supplementary']:
        primary = c in COHORTS
        paths = [S/f'{c}_metadata.csv.gz',S/f'{c}_raw_common.npz' if primary else S/f'{c}_raw_measured.npz']
        genes = G if primary else pd.Index(pd.read_csv(S/f'{c}_measured_genes.csv').gene)
        obs = pd.read_csv(paths[0])
        raw = sparse.load_npz(paths[1]).tocsr()
        out = annotate(obs,raw,genes)
        out.to_csv(R/f'{c}_RNA_annotation.csv.gz',index=False)
        macro = out.harmonized_group.eq('Macrophage')
        det = macro & out.PGAM5_counts.gt(0)
        coverage.append({'dataset':c,'scope':'primary_all_QC_HCC' if primary else 'supplementary_macrophages_only',
            'cells':len(out),'macrophages':int(macro.sum()),'PGAM5_RNA_detected_macrophages':int(det.sum()),
            'marker_corroborated_macrophages':int((macro & out.macrophage_RNA_marker_corroboration).sum()),
            'marker_corroborated_PGAM5_RNA_detected_macrophages':int((det & out.macrophage_RNA_marker_corroboration).sum()),
            'macrophages_other_lineage_review_flag':int((macro & out.other_lineage_RNA_review_flag).sum()),
            'PGAM5_detected_nonmacrophages':int((~macro & out.PGAM5_counts.gt(0)).sum()),
            'validated_TCGA_abundance':False})
        symbols = sorted(set(sum(PANELS.values(),[])) | {'PGAM5'})
        for group, part in obs.groupby('harmonized_group',sort=True):
            ix = part.index.to_numpy()
            for gene in symbols:
                measured = gene in genes
                x = raw[ix,genes.get_loc(gene)].toarray().ravel() if measured else None
                marker_rows.append({'dataset':c,'group':group,'gene':gene,'cells':len(ix),'measured':measured,
                    'detection_fraction':np.mean(x>0) if measured else np.nan,
                    'mean_CP10k':np.mean(x*10000/part.total_counts.to_numpy()) if measured else np.nan})
        provenance += [{'path':str(p),'bytes':p.stat().st_size,'sha256':digest(p)} for p in paths]
        if not primary:
            p=S/f'{c}_measured_genes.csv'
            provenance.append({'path':str(p),'bytes':p.stat().st_size,'sha256':digest(p)})
        print('Annotation',coverage[-1],flush=True)
        if c=='GSE149614':
            # Actual cells spanning RNA status and competing lineage; portable import/audit.
            pick=[]
            for mask in [macro & det,macro & ~det,~macro]: pick.extend(out.index[mask][:20])
            pick=np.array(pick,dtype=int)
            testgenes=genes[genes.isin(symbols)]
            pd.DataFrame(raw[pick][:,genes.get_indexer(testgenes)].toarray(),columns=testgenes).to_csv(R/'actual_cell_marker_counts.csv',index=False)
            out.iloc[pick].to_csv(R/'actual_cell_annotations.csv',index=False)
            pd.DataFrame({'gene':testgenes}).to_csv(R/'actual_cell_measured_genes.csv',index=False)
    pd.DataFrame(coverage).to_csv(R/'annotation_coverage.csv',index=False)
    pd.DataFrame(marker_rows).to_csv(R/'marker_expression_by_lineage.csv',index=False)
    roles=[]
    for panel,genes in PANELS.items():
        roles += [{'gene':g,'role':panel,'PGAM5_specific_evidence':False} for g in genes]
    roles.append({'gene':'PGAM5','role':'RNA_detection_definition_anchor','PGAM5_specific_evidence':False})
    pd.DataFrame(roles).to_csv(R/'gene_roles.csv',index=False)

    de={c:pd.read_csv(S/f'{c}_tie_corrected_DE.csv.gz') for c in COHORTS}
    profiles, donor_contrasts = {},{}
    for c in COHORTS:
        assert np.array_equal(de[c].gene,G)
        index=pd.read_csv(S/f'{c}_group_stats_index.csv')
        stats=np.load(S/f'{c}_group_statistics.npz')['values']
        index['broad_group']=index.group.str.replace(r'^Cycling_','',regex=True)
        index.loc[index.broad_group.str.startswith('TAM_'),'broad_group']='Macrophage'
        profiles[c]={}
        for group,part in index.groupby('broad_group',sort=True):
            if group=='Macrophage': continue
            donorprofiles=[]
            for donor,d in part.groupby('donor',sort=True):
                donorprofiles.append(stats[d.index,0].astype(float).sum(0)/d.cells.sum())
            profiles[c][group]={'profile':np.stack(donorprofiles).mean(0),'cells':int(part.cells.sum()),'donors':len(donorprofiles)}
        target=index[index.group.str.contains('PGAM5_detected',regex=False)]
        donorprofiles=[]
        for donor,d in target.groupby('donor',sort=True):
            donorprofiles.append(stats[d.index,0].astype(float).sum(0)/d.cells.sum())
        profiles[c]['Target_RNA_detected_macrophage']={'profile':np.stack(donorprofiles).mean(0),'cells':int(target.cells.sum()),'donors':len(donorprofiles)}
        ind=pd.read_csv(S/f'{c}_donor_contrasts_index.csv')
        values=np.load(S/f'{c}_donor_means.npz')['values']
        contrasts=[]
        for donor,d in ind.groupby('donor',sort=True):
            pos=d[d.status.eq('detected')].index[0]; neg=d[d.status.eq('undetected')].index[0]
            if ind.loc[pos,'cells']>=5 and ind.loc[neg,'cells']>=20:
                contrasts.append(np.log2((values[pos].astype(float)+.05)/(values[neg].astype(float)+.05)))
        donor_contrasts[c]=np.stack(contrasts)
        for suffix in ['tie_corrected_DE.csv.gz','group_stats_index.csv','group_statistics.npz','donor_contrasts_index.csv','donor_means.npz']:
            path=S/f'{c}_{suffix}'
            provenance.append({'path':str(path),'bytes':path.stat().st_size,'sha256':digest(path)})

    summaries=[]; support_rows=[]; compound=[]
    for held in COHORTS+['ALL_TRAINING']:
        training=[c for c in COHORTS if c!=held]
        votes=np.column_stack([de[c].strict_upregulated.to_numpy() for c in training])
        support=votes.sum(1)
        required=max(2,math.ceil(.6*len(training)))
        contrasts=np.concatenate([donor_contrasts[c] for c in training])
        direction=(contrasts>0).mean(0)
        base=np.array([g!='PGAM5' and not g.startswith(('MT-','RPL','RPS')) for g in G])
        recurrent=base & (support>=required) & (direction>=.7)
        cycle=np.array([g in CYCLE for g in G]); stress=np.array([g in STRESS for g in G])
        target=np.stack([profiles[c]['Target_RNA_detected_macrophage']['profile'] for c in training]).mean(0)
        groups=sorted(set().union(*[set(profiles[c]) for c in training])-{'Target_RNA_detected_macrophage'})
        competitors=[]; names=[]
        for group in groups:
            available=[profiles[c][group] for c in training if group in profiles[c]]
            if sum(x['cells'] for x in available)>=50 and sum(x['donors'] for x in available)>=2:
                competitors.append(np.stack([x['profile'] for x in available]).mean(0)); names.append(group)
        competitors=np.stack(competitors)
        win=competitors.argmax(0); max_other=competitors.max(0)
        specificity=np.log2((target+.05)/(max_other+.05))
        eligible=recurrent & ~cycle & ~stress & (specificity>=1)
        table=pd.DataFrame({'gene':G,'strict_cohort_votes':support,'required_votes':required,
            'eligible_donor_labels':len(contrasts),'donor_positive_direction_fraction':direction,
            'recurrent_candidate':recurrent,'cycle_list_member':cycle,'stress_list_member':stress,
            'training_target_mean_CP10k':target,'training_max_competitor_mean_CP10k':max_other,
            'strongest_training_competitor':[names[x] for x in win],
            'training_competitive_log2_ratio':specificity,'RNA_state_feature_eligible':eligible})
        for c in training:
            table[c+'_strict_up']=de[c].strict_upregulated.to_numpy()
        folder=R/held; folder.mkdir(exist_ok=True)
        table.to_csv(folder/'training_gene_evidence.csv.gz',index=False)
        selected=table[eligible].copy()
        selected.to_csv(folder/'eligible_RNA_state_genes.csv',index=False)
        sufficient=int(eligible.sum())>=5
        definition={'training_cohorts':training,'held_cohort':held,'RNA_state_features':G[eligible].tolist(),
            'independent_RNA_state_feature_count':int(eligible.sum()),'minimum_required_features':5,
            'feature_count_gate_passed':sufficient,'protein_validation_required':False,
            'interpretation':'Feature sufficiency only; no biological state or cellular abundance validation implied'}
        (folder/'definition.json').write_text(json.dumps(definition,indent=2)+'\n')
        summaries.append({'held_cohort':held,'training_cohorts':','.join(training),'required_cohort_votes':required,
            'recurrent_candidates':int(recurrent.sum()),'recurrent_cycle_candidates':int((recurrent&cycle).sum()),
            'recurrent_stress_candidates':int((recurrent&stress).sum()),'noncycle_nonstress_candidates':int((recurrent&~cycle&~stress).sum()),
            'competitive_specific_candidates':int(eligible.sum()),'feature_count_gate_passed':sufficient,
            'eligible_genes':';'.join(G[eligible]),'recurrent_genes':';'.join(G[recurrent])})
        report=table[recurrent | table.gene.eq('PGAM5')].copy(); report['held_cohort']=held
        compound.append(report)
        for c in COHORTS:
            for group,info in profiles[c].items():
                for i in np.flatnonzero(recurrent | np.asarray(G=='PGAM5')):
                    support_rows.append({'fold':held,'dataset':c,'role':'heldout' if c==held else 'training',
                        'group':group,'gene':G[i],'mean_CP10k':info['profile'][i],
                        'cells':info['cells'],'donor_labels':info['donors']})
        print('Screen',summaries[-1],flush=True)
    pd.DataFrame(summaries).to_csv(R/'RNA_state_screen_by_fold.csv',index=False)
    pd.concat(compound,ignore_index=True).to_csv(R/'recurrent_gene_competitor_evidence.csv',index=False)
    pd.DataFrame(support_rows).to_csv(R/'candidate_expression_by_competitor.csv',index=False)
    prior7=Path(P['v7_root'])
    previous=pd.read_csv(prior7/'all_validation_summary.csv')
    previous.to_csv(R/'prior_v7_validation_summary.csv',index=False)
    for path in [prior7/'all_validation_summary.csv',prior7/'all_validation_gates.csv',Path(P['v6_root'])/'validation_result.json',S/'common_genes.csv',S/'cell_cycle_genes.txt']:
        provenance.append({'path':str(path),'bytes':path.stat().st_size,'sha256':digest(path)})
    (R/'source_provenance.json').write_text(json.dumps({'inputs':provenance,'protocol_sha256':digest(R/'protocol.json'),
        'history':'Retrospective training-only holdout; not fresh blinded biological validation'},indent=2)+'\n')
    assert not any(x['feature_count_gate_passed'] for x in summaries), 'Sufficient feature panel requires implementing additional locked state testing, not automatic acceptance'
    (R/'validation_result.json').write_text(json.dumps({'protein_validation_required':False,
        'scope':'RNA-only target','RNA_detection_annotation_reproducible':True,
        'RNA_lineage_marker_evidence_added':True,'stable_RNA_state_feature_gate_passed':False,
        'previous_heldout_mixture_validation_passed':False,'independent_untouched_RNA_validation_available':False,
        'real_known_mixture_validation_available':False,'bulk_platform_cell_RNA_calibrated':False,
        'usable_for_validated_TCGA_PGAM5_cell_abundance':False,
        'RNA_zero_cells_reclassified_as_positive':0,
        'decision':'RNA detection annotation available; independent multigene state and calibrated abundance unsupported'},indent=2)+'\n')


if __name__=='__main__': main()
