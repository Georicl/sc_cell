"""复核粗分类映射、参考特征、两种预测及全部指标。"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import anndata as ad

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/Wu2021_coarse_T_v1';DATA=ROOT/'data/processed/Wu2021_coarse_T_v1'
protocol=json.loads((OUT/'protocol.json').read_text());mapping=pd.read_csv(OUT/'label_mapping.csv',index_col=0).reference_label
checks={};maxdiff=0.
for path,expected in json.loads((OUT/'input_hashes.json').read_text()).items():
    h=hashlib.sha256()
    with (ROOT/path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    checks['hash_'+path]=h.hexdigest()==expected
meta=pd.read_csv(DATA/'reference_metadata.csv',index_col=0)
old=pd.read_csv(ROOT/'data/processed/Wu2021_known_composition_v1/reference_metadata.csv',index_col=0)
checks['same_reference_cells_and_counts']=meta.index.equals(old.index) and np.array_equal(meta.nUMI,old.nUMI)
checks['T_mapping']=set(mapping.index[mapping.eq('T cells')])==set(protocol['merge_to_T_cells']) and mapping.loc['NK cells']=='NK cells'
checks['metadata_mapping']=np.array_equal(old.reference_label.map(mapping),meta.reference_label)
ref=ad.read_h5ad(ROOT/'data/processed/Wu2021_known_composition_v1/reference_counts.h5ad')
profiles=pd.read_csv(DATA/'nnls_reference_profiles.csv',index_col=0)
norm=ref[:,profiles.index].X.astype(float).multiply((1/np.asarray(ref.X.sum(axis=1)).ravel())[:,None]).tocsr()
expected=np.column_stack([np.asarray(norm[meta.reference_label.eq(t).to_numpy()].mean(axis=0)).ravel() for t in profiles.columns])
checks['coarse_nnls_profiles']=np.allclose(expected,profiles.values,rtol=1e-10,atol=1e-14)
metrics=pd.read_csv(OUT/'comparison_metrics.csv');pure=pd.read_csv(OUT/'pure_control_summary.csv')
for entry in protocol['datasets']:
    donor=entry['donor_id'];res=ROOT/entry['result_dir'];out=OUT/donor
    design=pd.read_csv(res/'simulation_design.csv',index_col=0)
    checks[donor+'_not_in_reference']=donor not in set(meta.donor_id)
    fine_truths={n:pd.read_csv(res/f'truth_{f}_fraction.csv',index_col=0) for n,f in [('rna','rna'),('cells','cell')]}
    truths={}
    for name,truth in fine_truths.items():
        aggregate=truth.T.groupby(mapping).sum().T
        saved=pd.read_csv(out/f'truth_{name}_coarse.csv',index_col=0)
        checks[donor+'_'+name+'_truth']=np.allclose(saved,aggregate.loc[saved.index,saved.columns])
        truths[name]=saved
    for method in ['RCTD','NNLS']:
        fine=pd.read_csv(ROOT/entry['fine_prediction_dir']/f'{method.lower()}_relative_weights.csv',index_col=0)
        aggregate=fine.T.groupby(mapping).sum().T
        for approach in ['fine_fit_then_sum','coarse_reference_refit']:
            pred=pd.read_csv(out/f'{method.lower()}_{approach}_weights.csv',index_col=0)
            source=aggregate if approach=='fine_fit_then_sum' else pd.read_csv(out/f'{method.lower()}_relative_weights.csv',index_col=0)
            tag=donor+'_'+method+'_'+approach
            checks[tag+'_values']=set(pred.index)==set(design.index) and pred.index.is_unique and np.isfinite(pred.values).all() and (pred.values>=0).all() and np.allclose(pred.sum(axis=1),1) and np.allclose(pred,source.loc[pred.index,pred.columns])
            for row in metrics[metrics.donor_id.eq(donor)&metrics.method.eq(method)&metrics.approach.eq(approach)].itertuples():
                ids=design.index[design.scenario.eq(row.scenario)&design.depth_fraction.eq(row.depth_fraction)]
                estimate=pred.loc[ids,row.cell_type].to_numpy();actual=truths[row.truth].loc[ids,row.cell_type].to_numpy();e=(estimate-actual)*100
                maxdiff=max(maxdiff,abs(np.abs(e).mean()-row.mae_pp),abs(np.sqrt((e**2).mean())-row.rmse_pp),abs(e.mean()-row.bias_pp))
                absent=truths['cells'].loc[ids,row.cell_type].eq(0).to_numpy()
                if absent.any():maxdiff=max(maxdiff,abs(estimate[absent].mean()*100-row.absent_mean_prediction_pct),abs((estimate[absent]>=.01).mean()-row.absent_ge1pct_fraction))
            pure_ok=True
            for row in pure[pure.donor_id.eq(donor)&pure.method.eq(method)&pure.approach.eq(approach)].itertuples():
                ids=design.index[design.scenario.eq('pure')&design.depth_fraction.eq(row.depth_fraction)&fine_truths['cells'][row.original_type].eq(1)]
                pure_ok &= int(pred.loc[ids].idxmax(axis=1).eq(row.coarse_type).sum())==row.dominant_matches
                maxdiff=max(maxdiff,abs(pred.loc[ids,'T cells'].mean()-row.mean_T_weight))
            checks[tag+'_pure_controls']=pure_ok
checks['all_metrics_recomputed']=maxdiff<1e-10
sp=OUT/'spatial_CID4535'
raw=pd.read_csv(ROOT/'data/processed/CID4535_rctd_v1/rctd_weights_spot_by_celltype.csv',index_col=0)
saved=pd.read_csv(sp/'coarse_relative_weights.csv',index_col=0)
expected=raw.div(raw.sum(axis=1),axis=0).T.groupby(mapping).sum().T
checks['spatial_aggregation']=np.allclose(expected.loc[saved.index,saved.columns],saved) and len(saved)==1102
stats=pd.read_csv(sp/'spot_stability.csv')
checks['spatial_frequency_bounds']=stats.top10_frequency.between(-1e-12,1+1e-12).all()
result=dict(checks={k:bool(v) for k,v in checks.items()},max_metric_difference=float(maxdiff))
(OUT/'validation_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print('Checks:',sum(checks.values()),'/',len(checks),'max metric difference:',maxdiff)
assert all(checks.values())
