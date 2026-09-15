"""核对参考抽样、标签映射、固定spot和评价计算。"""
from pathlib import Path
import hashlib
import json
import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data/processed/Wu2021_reference_comparison_v1'
OUT=ROOT/'results/Wu2021_reference_comparison_v1'
BASE=ROOT/'results/Wu2021_known_composition_v1'
ref=ad.read_h5ad(ROOT/'data/processed/Wu2021_known_composition_v1/reference_counts.h5ad')
checks={}
metas={}
for variant in ['donor_balanced','epithelial_split','balanced_split']:
    meta=pd.read_csv(DATA/variant/'reference_metadata.csv',index_col=0)
    metas[variant]=meta
    x=mmread(DATA/variant/'reference_counts.mtx').T.tocsr()
    source=ref[meta.index].X
    checks[variant+'_exact_source_counts']=(x!=source).nnz==0
    checks[variant+'_unique_reference_cells']=meta.index.is_unique and meta.donor_id.isin(['CID4471','CID4535']).all()
    checks[variant+'_parent_labels']=np.array_equal(meta.parent_label,ref.obs.loc[meta.index,'reference_label'].astype(str))
    checks[variant+'_minimum_25']=meta.loc[meta.nUMI>=100,'reference_label'].value_counts().min()>=25
    if variant!='donor_balanced':
        actual=meta.loc[meta.parent_label.eq('Cancer Epithelial'),'reference_label']
        expected=ref.obs.loc[actual.index,'celltype_minor'].astype(str).map({'Cancer LumA SC':'Cancer LumA','Cancer LumB SC':'Cancer LumB'}).fillna('Cancer Other')
        checks[variant+'_author_subtype_mapping']=actual.equals(expected.rename('reference_label'))
    profile=pd.read_csv(DATA/variant/'nnls_reference_profiles.csv',index_col=0)
    idx=ref.var_names.get_indexer(profile.index)
    normalized=source[:,idx].multiply(1/np.asarray(source.sum(axis=1))).tocsr()
    expected=np.column_stack([np.asarray(normalized[meta.reference_label.eq(t).to_numpy()].mean(axis=0)).ravel() for t in profile.columns])
    checks[variant+'_nnls_profiles']=np.allclose(expected,profile.values,rtol=1e-10,atol=1e-14)
checks['split_uses_all_reference_cells']=metas['epithelial_split'].index.equals(ref.obs_names)
checks['balanced_variants_share_cells']=metas['donor_balanced'].index.equals(metas['balanced_split'].index)
counts=metas['donor_balanced'].groupby(['parent_label','donor_id']).size().unstack()
checks['donor_counts_balanced_except_cycling']=counts.drop(index='Cycling T-cells').nunique(axis=1).eq(1).all()
checks['cycling_pool_retained']=counts.loc['Cycling T-cells'].to_dict()=={'CID4471':8,'CID4535':19}
for name in ['spatial_counts.mtx','spatial_metadata.csv','genes.tsv']:
    old_hashes=json.loads((BASE/'input_hashes.json').read_text())
    path=ROOT/'data/processed/Wu2021_known_composition_v1'/name
    checks['original_spot_input_'+name]=hashlib.sha256(path.read_bytes()).hexdigest()==old_hashes[name]
for path,digest in json.loads((OUT/'input_hashes.json').read_text()).items():
    h=hashlib.sha256()
    with (ROOT/path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    checks['hash_'+path]=h.hexdigest()==digest
metrics=pd.read_csv(OUT/'comparison_metrics.csv')
design=pd.read_csv(BASE/'simulation_design.csv',index_col=0)
truths={'rna':pd.read_csv(BASE/'truth_rna_fraction.csv',index_col=0),'cells':pd.read_csv(BASE/'truth_cell_fraction.csv',index_col=0)}
max_error=0
for (variant,method),group in metrics.groupby(['variant','method']):
    folder=BASE if variant=='baseline' else OUT/variant
    suffix='relative' if variant=='baseline' else 'parent'
    pred=pd.read_csv(folder/f'{method.lower()}_{suffix}_weights.csv',index_col=0)
    if variant!='baseline':
        raw=pd.read_csv(folder/f'{method.lower()}_relative_weights.csv',index_col=0)
        mapping=pd.read_csv(folder/'label_mapping.csv',index_col=0).parent_label
        aggregate=raw.T.groupby(mapping).sum().T.loc[pred.index,pred.columns]
        checks[variant+'_'+method+'_aggregation']=np.allclose(aggregate,pred)
    checks[variant+'_'+method+'_valid_weights']=len(pred)==488 and set(pred.index)==set(design.index) and np.isfinite(pred.values).all() and (pred.values>=0).all() and np.allclose(pred.sum(axis=1),1)
    for row in group.itertuples():
        ids=design.index[design.scenario.eq(row.scenario)&design.depth_fraction.eq(row.depth_fraction)]
        cols=['Cancer Epithelial','Normal Epithelial'] if row.cell_type=='Epithelial total' else [row.cell_type]
        e=(pred.loc[ids,cols].sum(axis=1).to_numpy()-truths[row.truth].loc[ids,cols].sum(axis=1).to_numpy())*100
        max_error=max(max_error,abs(np.abs(e).mean()-row.mae_pp),abs(np.sqrt(np.square(e).mean())-row.rmse_pp),abs(e.mean()-row.bias_pp))
checks['all_mae_rmse_bias_recomputed']=max_error<1e-10
old=pd.read_csv(BASE/'benchmark_metrics.csv')
new=metrics[metrics.variant.eq('baseline')&~metrics.cell_type.eq('Epithelial total')]
keys=['method','truth','scenario','depth_fraction','cell_type']
compare=old.merge(new,on=keys,suffixes=('_old','_new'),validate='one_to_one')
checks['baseline_matches_previous_run']=len(compare)==len(old) and all(np.allclose(compare[f'{col}_old'],compare[f'{col}_new'],equal_nan=True) for col in ['mae_pp','rmse_pp','bias_pp','present_mae_pp','absent_mean_prediction_pct','absent_ge1pct_fraction'])
result={'checks':{key:bool(value) for key,value in checks.items()},'max_metric_difference_pp':float(max_error)}
(OUT/'validation_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print('Checks:',len(checks),'passed:',sum(checks.values()),'max metric difference:',max_error)
assert all(checks.values())
