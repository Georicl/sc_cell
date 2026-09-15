"""构建T细胞粗分类参考，复用五位供者的固定模拟输入。"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import anndata as ad
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data/processed/Wu2021_coarse_T_v1'
OUT=ROOT/'results/Wu2021_coarse_T_v1'
BASE=ROOT/'data/processed/Wu2021_known_composition_v1'
DATA.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
TYPES=['T cells CD4+','T cells CD8+','NKT cells','Cycling T-cells']
datasets=[]
for donor,folder,child in [('CID4067','Wu2021_known_composition_v1',''),('CID4463','Wu2021_independent_validation_v1',''),('CID3941','Wu2021_tcell_validation_v1','CID3941'),('CID3948','Wu2021_tcell_validation_v1','CID3948'),('CID4461','Wu2021_tcell_validation_v1','CID4461')]:
    inp=ROOT/'data/processed'/folder/child;res=ROOT/'results'/folder/child
    fine=res if donor=='CID4067' else res/'baseline'
    datasets.append(dict(donor_id=donor,input_dir=str(inp.relative_to(ROOT)),result_dir=str(res.relative_to(ROOT)),fine_prediction_dir=str(fine.relative_to(ROOT))))
protocol=dict(frozen_at_utc=datetime.now(timezone.utc).isoformat(),reference_donors=['CID4471','CID4535'],merge_to_T_cells=TYPES,NK_policy='Separate NK cells; author major T-cells includes NK, so major labels alone are not used.',n_coarse_types=10,datasets=datasets,
    comparisons=['fine_fit_then_sum','coarse_reference_refit'],reference_cells='All original 12570 reference cells, no filtering or donor reweighting',
    rctd=dict(seed=9,mode='full',max_cores=2,ref_UMI_min=100,ref_n_cells_min=25),
    nnls=dict(features=2000,selection='variance/mean across coarse reference means; reference only'),
    truths=['captured RNA fraction','cell-number fraction'],evaluation='Per donor, scenario and depth; T-zero controls and NK pure controls; all coarse classes retained.',
    status='Reference development comparison using previously inspected donors, not a new independent validation.')
(OUT/'protocol.json').write_text(json.dumps(protocol,ensure_ascii=False,indent=2)+'\n')
meta=pd.read_csv(BASE/'reference_metadata.csv',index_col=0)
meta['original_label']=meta.reference_label
meta['reference_label']=meta.reference_label.where(~meta.reference_label.isin(TYPES),'T cells')
assert meta.reference_label.nunique()==10
meta.to_csv(DATA/'reference_metadata.csv')
mapping=meta[['original_label','reference_label']].drop_duplicates().set_index('original_label')
mapping.to_csv(OUT/'label_mapping.csv')
meta.groupby(['donor_id','reference_label']).size().rename('n_cells').to_csv(OUT/'reference_counts_by_type.csv')
ref=ad.read_h5ad(BASE/'reference_counts.h5ad')
assert ref.obs_names.equals(meta.index)
normalized=ref.X.astype(float).multiply((1/meta.nUMI.to_numpy())[:,None]).tocsr()
types=sorted(meta.reference_label.unique())
profiles=np.column_stack([np.asarray(normalized[meta.reference_label.eq(t).to_numpy()].mean(axis=0)).ravel() for t in types])
score=profiles.var(axis=1)/(profiles.mean(axis=1)+1e-12)
order=np.lexsort((ref.var_names.to_numpy(),-score));selected=order[score[order]>0][:2000]
pd.DataFrame(profiles[selected],index=ref.var_names[selected],columns=types).to_csv(DATA/'nnls_reference_profiles.csv')
hashes={}
paths=[BASE/'reference_counts.mtx',BASE/'genes.tsv',DATA/'reference_metadata.csv',DATA/'nnls_reference_profiles.csv',OUT/'protocol.json',OUT/'label_mapping.csv']
for entry in datasets:
    inp=ROOT/entry['input_dir'];res=ROOT/entry['result_dir']
    paths += [inp/'spatial_counts.mtx',inp/'spatial_metadata.csv',inp/'genes.tsv',res/'simulation_design.csv',res/'truth_rna_fraction.csv',res/'truth_cell_fraction.csv']
for path in paths:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    hashes[str(path.relative_to(ROOT))]=h.hexdigest()
(OUT/'input_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
print(meta.reference_label.value_counts())
