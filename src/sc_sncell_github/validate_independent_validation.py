"""核对独立供者来源、模拟组成及全部评价指标。"""
from pathlib import Path
import hashlib
import json
import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data/processed/Wu2021_independent_validation_v1'
OUT=ROOT/'results/Wu2021_independent_validation_v1'
test=ad.read_h5ad(DATA/'CID4463_counts_before_qc.h5ad')
spots=ad.read_h5ad(DATA/'simulated_spots.h5ad')
design=pd.read_csv(OUT/'simulation_design.csv',index_col=0)
members=pd.read_csv(OUT/'mixture_membership.csv')
truths={'rna':pd.read_csv(OUT/'truth_rna_fraction.csv',index_col=0),'cells':pd.read_csv(OUT/'truth_cell_fraction.csv',index_col=0)}
checks={}
checks['author_umi_counts']=np.array_equal(np.asarray(test.X.sum(axis=1)).ravel(),test.obs.nCount_RNA)
checks['author_detected_genes']=np.array_equal(np.asarray((test.X>0).sum(axis=1)).ravel(),test.obs.nFeature_RNA)
checks['gene_order']=spots.var_names.equals(test.var_names)
checks['matrix_market_counts']=(mmread(DATA/'spatial_counts.mtx').T.tocsr()!=spots.X).nnz==0
checks['spot_order']=spots.obs_names.equals(design.index)
checks['donor_is_CID4463']=set(test.obs['orig.ident'])=={'CID4463'} and design.donor_id.eq('CID4463').all()
checks['no_duplicate_cells_within_spot']=not members.duplicated(['spot_id','source_barcode']).any()
checks['source_barcode_membership']=members.source_barcode.isin(test.obs_names).all()
author_labels=test.obs.celltype_major.astype(str).where(~test.obs.celltype_major.eq('T-cells'),test.obs.celltype_minor.astype(str))
checks['membership_labels_match_author']=np.array_equal(members.reference_label,author_labels.loc[members.source_barcode])
checks['component_umi_totals']=np.array_equal(members.groupby('spot_id').component_umi.sum().reindex(design.index),design.nUMI)
checks['matrix_umi_totals']=np.array_equal(np.asarray(spots.X.sum(axis=1)).ravel(),design.nUMI)
for name,values,agg in [('rna','component_umi','sum'),('cells','source_barcode','count')]:
    target=truths[name]
    table=members.pivot_table(index='spot_id',columns='reference_label',values=values,aggfunc=agg,fill_value=0).reindex(index=target.index,columns=target.columns,fill_value=0)
    checks[name+'_truth']=np.allclose(table.div(table.sum(axis=1),axis=0),target)
exact=paired=thinned=True
for _,pair in design.groupby('base_mixture'):
    original=pair.index[pair.depth_fraction.eq(1)][0];low=pair.index[pair.depth_fraction.eq(.25)][0]
    barcodes=members.loc[members.spot_id.eq(original),'source_barcode']
    actual=spots[original].X.toarray().ravel()
    exact &= np.array_equal(actual,np.asarray(test[barcodes].X.sum(axis=0)).ravel())
    paired &= set(barcodes)==set(members.loc[members.spot_id.eq(low),'source_barcode'])
    thinned &= bool((spots[low].X.toarray().ravel()<=actual).all())
checks.update(all_237_original_counts_reconstructed=exact,paired_depth_sources=paired,thinned_counts_bounded=thinned)
recipes=pd.read_csv(OUT/'frozen_recipes.csv',index_col=0)
checks['recipe_count_and_spot_count']=len(recipes)==237 and len(design)==474
count_table=pd.crosstab(members.spot_id,members.reference_label).reindex(index=design.index,columns=truths['cells'].columns,fill_value=0)
expected=recipes.loc[design.base_mixture,truths['cells'].columns].to_numpy()
checks['frozen_recipe_counts']=np.array_equal(count_table.to_numpy(),expected)
for manifest in ['frozen_reference_hashes.json','simulation_hashes.json']:
    for path,digest in json.loads((OUT/manifest).read_text()).items():
        h=hashlib.sha256()
        with (ROOT/path).open('rb') as f:
            for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
        checks['hash_'+path]=h.hexdigest()==digest
metrics=pd.read_csv(OUT/'validation_metrics.csv')
maxdiff=0
for (variant,method),group in metrics.groupby(['variant','method']):
    ref=ROOT/'data/processed/Wu2021_known_composition_v1' if variant=='baseline' else ROOT/'data/processed/Wu2021_reference_comparison_v1'/variant
    meta=pd.read_csv(ref/'reference_metadata.csv',index_col=0)
    checks[variant+'_reference_test_disjoint']=set(meta.index).isdisjoint(test.obs_names) and meta.donor_id.isin(['CID4471','CID4535']).all()
    checks[variant+'_reference_gene_order']=pd.read_csv(ref/'genes.tsv',header=None)[0].tolist()==test.var_names.tolist()
    folder=OUT/variant
    pred=pd.read_csv(folder/f'{method.lower()}_parent_weights.csv',index_col=0)
    native=pd.read_csv(folder/f'{method.lower()}_relative_weights.csv',index_col=0)
    if variant!='baseline':
        mapping=pd.read_csv(ROOT/'results/Wu2021_reference_comparison_v1'/variant/'label_mapping.csv',index_col=0).parent_label
        native=native.T.groupby(mapping).sum().T
    checks[variant+'_'+method+'_parent_aggregation']=np.allclose(pred,native.loc[pred.index,pred.columns])
    checks[variant+'_'+method+'_valid_weights']=set(pred.index)==set(design.index) and len(pred)==474 and np.isfinite(pred.values).all() and (pred.values>=0).all() and np.allclose(pred.sum(axis=1),1)
    for row in group.itertuples():
        ids=design.index[design.scenario.eq(row.scenario)&design.depth_fraction.eq(row.depth_fraction)]
        cols=['Cancer Epithelial','Normal Epithelial'] if row.cell_type=='Epithelial total' else [row.cell_type]
        error=(pred.loc[ids,cols].sum(axis=1).to_numpy()-truths[row.truth].loc[ids,cols].sum(axis=1).to_numpy())*100
        maxdiff=max(maxdiff,abs(np.abs(error).mean()-row.mae_pp),abs(np.sqrt((error**2).mean())-row.rmse_pp),abs(error.mean()-row.bias_pp))
checks['all_mae_rmse_bias_recomputed']=maxdiff<1e-10
result={'checks':{k:bool(v) for k,v in checks.items()},'max_metric_difference_pp':float(maxdiff)}
(OUT/'validation_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print('Checks:',sum(checks.values()),'/',len(checks),'max metric difference:',maxdiff)
assert all(checks.values())
