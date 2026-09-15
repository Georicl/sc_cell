"""从原始对象复核逐细胞标签与全部查看基因。"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import anndata as ad

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/Wu2021_tcell_label_review_v1'
review=pd.read_csv(OUT/'tcell_label_expression_review.csv',index_col=0)
provenance=json.loads((OUT/'provenance.json').read_text());markers=provenance['markers']
checks={}
for relative in provenance['source_hashes']:
    a=ad.read_h5ad(ROOT/relative)
    donor=str(a.obs['orig.ident'].iloc[0]);r=review[review['orig.ident'].eq(donor)]
    original=a[r.index]
    checks[donor+'_raw_marker_counts']=np.array_equal(original[:,markers].X.toarray(),r[markers].to_numpy())
    for col in ['celltype_major','celltype_minor','celltype_subset']:
        checks[donor+'_'+col]=np.array_equal(original.obs[col].astype(str),r[col])
    checks[donor+'_full_umi']=np.array_equal(np.asarray(original.X.sum(axis=1)).ravel(),r.nCount_RNA)
    checks[donor+'_T_marker_detection']=np.array_equal(original[:,['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2']].X.sum(axis=1).A1>0,r.T_markers_detected)
    checks[donor+'_CD8_detection']=np.array_equal(original[:,['CD8A','CD8B']].X.sum(axis=1).A1>0,r.CD8_detected)
checks['unique_barcodes']=review.index.is_unique
features=pd.read_csv(OUT/'marker_feature_coverage.csv')
checks['12_runs_CD8AB_retained']=bool(features[features.gene.isin(['CD8A','CD8B'])].in_RCTD.all())
checks['12_runs_CD4_FOXP3_absent']=bool((~features[features.gene.isin(['CD4','FOXP3'])].in_RCTD).all())
base=pd.read_csv(OUT/'baseline_rctd_profiles.csv',index_col=0)
checks['CD4_FOXP3_below_regression_threshold']=bool((base.loc[['CD4','FOXP3']].max(axis=1)<.0002).all())
for donor in ['CID3941','CID3948','CID4461']:
    s=pd.read_csv(OUT/'zero_CD8_mixture_source_review.csv')
    s=s[s.donor_id.eq(donor)&s.variant.eq('baseline')]
    parent=ROOT/'results/Wu2021_tcell_validation_v1'/donor
    members=pd.read_csv(parent/'mixture_membership.csv')
    selected=members[members.spot_id.isin(s.spot_id)]
    checks[donor+'_zero_CD8_sources_all_CD4']=selected.reference_label.eq('T cells CD4+').all()
    pred=pd.read_csv(parent/'baseline/rctd_relative_weights.csv',index_col=0)
    checks[donor+'_spot_prediction_link']=np.allclose(s.predicted_CD8_pct.to_numpy(),pred.loc[s.spot_id,'T cells CD8+'].to_numpy()*100)
checks={k:bool(v) for k,v in checks.items()}
(OUT/'validation_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
print('Checks:',sum(checks.values()),'/',len(checks))
assert all(checks.values())
