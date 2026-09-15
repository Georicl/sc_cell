"""复核每张切片的矩阵对齐、汇总权重、表达面板及新旧参考比较。"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy.io import mmread
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data/processed/Wu2021_multislice_CAF_T_v1';OUT=ROOT/'results/Wu2021_multislice_CAF_T_v1'

def main(samples,hash_inputs=True):
    protocol=json.loads((OUT/'protocol.json').read_text());checks=json.loads((OUT/'input_checks.json').read_text())
    if hash_inputs:
        for path,expected in json.loads((OUT/'input_hashes.json').read_text()).items():
            h=hashlib.sha256()
            with (ROOT/path).open('rb') as f:
                for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
            checks['hash_'+path]=h.hexdigest()==expected
    reference=ad.read_h5ad(ROOT/'data/processed/Wu2021_known_composition_v1/reference_counts.h5ad')
    genes=pd.read_csv(DATA/'genes.tsv',header=None)[0].astype(str)
    meta=pd.read_csv(DATA/'reference_metadata.csv',index_col=0)
    ref_mtx=mmread(DATA/'reference_counts.mtx').T.tocsr()
    checks['reference_exact_source_counts']=(ref_mtx!=reference[meta.index,genes].X).nnz==0
    checks['reference_labels_preserved']=np.array_equal(meta.reference_label,reference.obs.loc[meta.index,'reference_label'].astype(str))
    mapping=pd.read_csv(ROOT/'results/Wu2021_coarse_T_v1/label_mapping.csv',index_col=0).reference_label
    for sample in samples:
        out=OUT/sample;inp=DATA/sample
        a=ad.read_h5ad(inp/'spatial_counts_full.h5ad');sp_meta=pd.read_csv(inp/'spatial_metadata.csv',index_col=0)
        sample_genes=pd.read_csv(inp/'genes.tsv',header=None)[0].astype(str)
        x=mmread(inp/'spatial_counts.mtx').T.tocsr()
        checks[sample+'_exact_model_matrix']=(x!=a[:,sample_genes].X).nnz==0
        checks[sample+'_spatial_order']=a.obs_names.equals(sp_meta.index)
        checks[sample+'_coordinates']=np.array_equal(sp_meta[['x','y']],a.obs[['pxl_col_in_fullres','pxl_row_in_fullres']])
        qc=pd.read_csv(out/'all_spot_qc.csv',index_col=0)
        expected=qc.index[~(qc.in_tissue.eq(0)|qc.Classification.eq('Artefact').fillna(False))]
        checks[sample+'_QC_union']=a.obs_names.equals(expected)
        raw=pd.read_csv(out/'rctd_raw_weights.csv',index_col=0);fine=pd.read_csv(out/'rctd_relative_weights.csv',index_col=0);coarse=pd.read_csv(out/'coarse_relative_weights.csv',index_col=0)
        checks[sample+'_valid_weights']=np.isfinite(raw.values).all() and (raw.values>=0).all() and (raw.sum(axis=1)>0).all() and np.allclose(raw.div(raw.sum(axis=1),axis=0),fine)
        aggregate=fine.T.groupby(mapping).sum().T
        checks[sample+'_coarse_sum']=np.allclose(aggregate.loc[coarse.index,coarse.columns],coarse) and np.allclose(coarse.sum(axis=1),1)
        sensitivity=pd.read_csv(out/'T_cycling_sensitivity_by_spot.csv',index_col=0)
        checks[sample+'_cycling_sensitivity_sum']=np.allclose(sensitivity.T_without_cycling,fine[['T cells CD4+','T cells CD8+','NKT cells']].sum(axis=1)) and np.allclose(sensitivity.T_all,sensitivity.T_without_cycling+sensitivity.Cycling_T)
        removed=pd.read_csv(out/'model_removed_spots.csv').barcode.astype(str)
        checks[sample+'_removed_accounted']=set(a.obs_names)-set(raw.index)==set(removed) and set(raw.index).issubset(a.obs_names)
        scores=pd.read_csv(out/'marker_scores.csv',index_col=0);correlations=pd.read_csv(out/'marker_concordance.csv')
        for name,panel in protocol['marker_checks'].items():
            counts=a[coarse.index,panel].X.toarray();total=np.asarray(a[coarse.index].X.sum(axis=1)).ravel()
            calculated=np.log1p(counts/total[:,None]*10000).mean(axis=1)
            checks[sample+'_'+name+'_score']=np.allclose(scores[name+'_score'],calculated)
            checks[sample+'_'+name+'_detection']=np.array_equal(scores[name+'_genes_detected'],(counts>0).sum(axis=1))
            checks[sample+'_'+name+'_panel_UMI']=np.array_equal(scores[name+'_panel_umi'],counts.sum(axis=1))
        sensitivity_summary=pd.read_csv(out/'T_cycling_sensitivity_summary.csv').set_index('definition')
        for definition,col in [('main_T_all','T_all'),('diagnostic_without_cycling','T_without_cycling'),('Cycling_T_component','Cycling_T')]:
            checks[sample+'_'+definition+'_diagnostic_rho']=bool(np.isclose(spearmanr(sensitivity[col],scores.T_score).statistic,sensitivity_summary.loc[definition,'spearman'],equal_nan=True))
        for row in correlations.itertuples():
            rho=spearmanr(coarse[row.cell_type],scores[row.marker_panel+'_score']).statistic
            checks[sample+'_'+row.cell_type+'_'+row.marker_panel+'_rho']=bool(np.isclose(rho,row.spearman,equal_nan=True))
        review=pd.read_csv(out/'spot_review.csv',index_col=0)
        high=coarse['T cells'].sort_values(ascending=False,kind='stable').head(int(np.ceil(len(coarse)*.1))).index
        checks[sample+'_high_T_flag']=np.array_equal(review.high_T,review.index.isin(high))
        checks[sample+'_expression_flag']=np.array_equal(review.high_T_low_T_expression,review.high_T&scores.T_score.le(scores.T_score.quantile(.25)))
        if sample=='CID4535':
            oldraw=pd.read_csv(ROOT/'data/processed/CID4535_rctd_v1/rctd_weights_spot_by_celltype.csv',index_col=0).loc[coarse.index]
            old=oldraw.div(oldraw.sum(axis=1),axis=0).T.groupby(mapping).sum().T
            comparisons=pd.read_csv(out/'reference_comparison_metrics.csv')
            for row in comparisons.itertuples():
                e=(coarse[row.cell_type]-old[row.cell_type])*100
                checks['CID4535_'+row.cell_type+'_change']=bool(np.isclose(e.abs().mean(),row.mean_abs_change_pp) and np.isclose(e.mean(),row.mean_change_pp) and np.isclose(spearmanr(old[row.cell_type],coarse[row.cell_type]).statistic,row.spearman))
    checks={k:bool(v) for k,v in checks.items()}
    result=dict(samples=samples,checks=checks,n_passed=sum(checks.values()),n_checks=len(checks))
    filename='validation_all.json' if len(samples)==6 else 'validation_'+ '_'.join(samples)+'.json'
    (OUT/filename).write_text(json.dumps(result,indent=2)+'\n')
    print('Checks:',sum(checks.values()),'/',len(checks))
    assert all(checks.values())

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('samples',nargs='*',default=['CID4535']);args=parser.parse_args();main(args.samples)
