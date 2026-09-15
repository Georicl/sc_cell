"""独立复核270个来源混合、供者隔离、CAF配对及全部梯度指标。"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy.io import mmread

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data/processed/Wu2021_tcell_validation_v1'
OUT=ROOT/'results/Wu2021_tcell_validation_v1'
DONORS=['CID3941','CID3948','CID4461']
checks={};maxdiff=0.
for manifest in ['frozen_inputs.json','simulation_hashes.json']:
    for path,expected in json.loads((OUT/manifest).read_text()).items():
        digest=hashlib.sha256()
        with (ROOT/path).open('rb') as f:
            for chunk in iter(lambda:f.read(1024*1024),b''):digest.update(chunk)
        checks['hash_'+path]=digest.hexdigest()==expected
metrics=pd.read_csv(OUT/'gradient_metrics.csv')
pure_metrics=pd.read_csv(OUT/'pure_control_summary.csv')
all_pairs=pd.read_csv(OUT/'paired_caf_effects.csv')
for donor in DONORS:
    out=OUT/donor
    a=ad.read_h5ad(DATA/f'{donor}_counts_before_qc.h5ad')
    s=ad.read_h5ad(DATA/donor/'simulated_spots.h5ad')
    d=pd.read_csv(out/'simulation_design.csv',index_col=0)
    m=pd.read_csv(out/'mixture_membership.csv')
    truths={'rna':pd.read_csv(out/'truth_rna_fraction.csv',index_col=0),'cells':pd.read_csv(out/'truth_cell_fraction.csv',index_col=0)}
    labels=a.obs.celltype_major.astype(str).where(~a.obs.celltype_major.eq('T-cells'),a.obs.celltype_minor.astype(str))
    checks[donor+'_author_umi']=np.array_equal(np.asarray(a.X.sum(axis=1)).ravel(),a.obs.nCount_RNA)
    checks[donor+'_author_genes']=np.array_equal(np.asarray((a.X>0).sum(axis=1)).ravel(),a.obs.nFeature_RNA)
    checks[donor+'_sample_identity']=set(a.obs['orig.ident'])=={donor} and d.donor_id.eq(donor).all()
    checks[donor+'_membership_labels']=np.array_equal(labels.loc[m.source_barcode],m.reference_label)
    checks[donor+'_unique_within_spot']=not m.duplicated(['spot_id','source_barcode']).any()
    checks[donor+'_gene_order']=a.var_names.equals(s.var_names)
    checks[donor+'_mtx']=(mmread(DATA/donor/'spatial_counts.mtx').T.tocsr()!=s.X).nnz==0
    checks[donor+'_spot_order']=s.obs_names.equals(d.index)
    checks[donor+'_umi_totals']=np.array_equal(np.asarray(s.X.sum(axis=1)).ravel(),d.nUMI) and np.array_equal(m.groupby('spot_id').component_umi.sum().reindex(d.index),d.nUMI)
    for name,values,agg in [('rna','component_umi','sum'),('cells','source_barcode','count')]:
        target=truths[name]
        table=m.pivot_table(index='spot_id',columns='reference_label',values=values,aggfunc=agg,fill_value=0).reindex(index=target.index,columns=target.columns,fill_value=0)
        checks[donor+'_'+name+'_truth']=np.allclose(table.div(table.sum(axis=1),axis=0),target)
    recipes=pd.read_csv(out/'frozen_recipes.csv',index_col=0)
    ct=pd.crosstab(m.spot_id,m.reference_label).reindex(index=d.index,columns=truths['cells'].columns,fill_value=0)
    checks[donor+'_recipe_counts']=np.array_equal(ct,recipes.loc[d.base_mixture,truths['cells'].columns].to_numpy())
    exact=paired=bounded=True
    for _,pair in d.groupby('base_mixture'):
        high=pair.index[pair.depth_fraction.eq(1)][0];low=pair.index[pair.depth_fraction.eq(.25)][0]
        barcodes=m.loc[m.spot_id.eq(high),'source_barcode']
        expected=np.asarray(a[barcodes].X.sum(axis=0)).ravel()
        actual=s[high].X.toarray().ravel()
        exact &= np.array_equal(expected,actual)
        paired &= set(barcodes)==set(m.loc[m.spot_id.eq(low),'source_barcode'])
        bounded &= bool((s[low].X.toarray().ravel()<=actual).all())
    checks[donor+'_90_original_mixtures']=exact and len(recipes)==90
    checks[donor+'_paired_depths']=paired and bounded
    cafequal=True
    for _,group in d[d.scenario.eq('gradient')].groupby(['depth_fraction','pair_id']):
        rows={int(r.n_caf):i for i,r in group.iterrows()}
        t0=m[m.spot_id.eq(rows[0])].set_index('source_barcode').component_umi.sort_index()
        for caf in [1,3]:
            t=m[m.spot_id.eq(rows[caf])&m.reference_label.isin(['T cells CD4+','T cells CD8+'])].set_index('source_barcode').component_umi.sort_index()
            cafequal &= t0.equals(t)
        cafequal &= bool((s[rows[0]].X.toarray()<=s[rows[1]].X.toarray()).all() and (s[rows[1]].X.toarray()<=s[rows[3]].X.toarray()).all())
    checks[donor+'_paired_caf_backgrounds']=cafequal
    for variant in ['baseline','donor_balanced','epithelial_split','balanced_split']:
        ref=ROOT/'data/processed/Wu2021_known_composition_v1' if variant=='baseline' else ROOT/'data/processed/Wu2021_reference_comparison_v1'/variant
        meta=pd.read_csv(ref/'reference_metadata.csv',index_col=0)
        checks[donor+'_'+variant+'_separation']=set(meta.index).isdisjoint(a.obs_names) and meta.donor_id.isin(['CID4471','CID4535']).all()
        for method in ['RCTD','NNLS']:
            folder=out/variant
            pred=pd.read_csv(folder/f'{method.lower()}_parent_weights.csv',index_col=0)
            raw=pd.read_csv(folder/f'{method.lower()}_relative_weights.csv',index_col=0)
            if variant!='baseline':
                mapping=pd.read_csv(ROOT/'results/Wu2021_reference_comparison_v1'/variant/'label_mapping.csv',index_col=0).parent_label
                raw=raw.T.groupby(mapping).sum().T
            tag=donor+'_'+variant+'_'+method
            checks[tag+'_weights']=set(pred.index)==set(d.index) and len(pred)==180 and np.isfinite(pred.values).all() and (pred.values>=0).all() and np.allclose(pred.sum(axis=1),1) and np.allclose(raw.loc[pred.index,pred.columns],pred)
            group=metrics[metrics.donor_id.eq(donor)&metrics.variant.eq(variant)&metrics.method.eq(method)]
            for row in group.itertuples():
                ids=d.index[d.scenario.eq('gradient')&d.depth_fraction.eq(row.depth_fraction)&d.n_caf.eq(row.n_caf)]
                actual=truths[row.truth].loc[ids,row.cell_type].to_numpy();estimated=pred.loc[ids,row.cell_type].to_numpy()
                error=(estimated-actual)*100;absent=truths['cells'].loc[ids,row.cell_type].to_numpy()==0
                maxdiff=max(maxdiff,abs(np.abs(error).mean()-row.mae_pp),abs(np.sqrt((error**2).mean())-row.rmse_pp),abs(error.mean()-row.bias_pp))
                if absent.any():
                    maxdiff=max(maxdiff,abs(estimated[absent].mean()*100-row.absent_mean_prediction_pct),abs((estimated[absent]>=.01).mean()-row.absent_ge1pct_fraction),abs((estimated[absent]>=.05).mean()-row.absent_ge5pct_fraction))
            pure_ok=True
            for row in pure_metrics[pure_metrics.donor_id.eq(donor)&pure_metrics.variant.eq(variant)&pure_metrics.method.eq(method)].itertuples():
                ids=d.index[d.scenario.eq('pure')&d.depth_fraction.eq(row.depth_fraction)&truths['cells'][row.cell_type].eq(1)]
                pure_ok &= int(pred.loc[ids].idxmax(axis=1).eq(row.cell_type).sum())==row.dominant_matches
                maxdiff=max(maxdiff,abs(pred.loc[ids,row.cell_type].mean()-row.mean_weight_on_true_type),abs(pred.loc[ids,'T cells CD8+'].mean()-row.mean_cd8_weight))
            checks[tag+'_pure_matches']=pure_ok
            for row in all_pairs[all_pairs.donor_id.eq(donor)&all_pairs.variant.eq(variant)&all_pairs.method.eq(method)].itertuples():
                e0=(pred.loc[row.base_spot,row.cell_type]-truths['rna'].loc[row.base_spot,row.cell_type])*100
                e1=(pred.loc[row.caf_spot,row.cell_type]-truths['rna'].loc[row.caf_spot,row.cell_type])*100
                maxdiff=max(maxdiff,abs((e1-e0)-row.error_change_pp),abs((abs(e1)-abs(e0))-row.absolute_error_change_pp))
checks['all_metrics_and_pair_differences']=maxdiff<1e-10
result={'checks':{k:bool(v) for k,v in checks.items()},'max_recalculation_difference':float(maxdiff)}
(OUT/'validation_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print('Checks:',sum(checks.values()),'/',len(checks),'max recalculation difference:',maxdiff)
assert all(checks.values())
