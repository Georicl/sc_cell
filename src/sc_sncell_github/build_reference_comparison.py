"""固定三种参考调整，复用上一轮488个模拟spot。"""
from pathlib import Path
import hashlib
import json
import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmwrite

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'data/processed/Wu2021_known_composition_v1'
DATA = ROOT / 'data/processed/Wu2021_reference_comparison_v1'
OUT = ROOT / 'results/Wu2021_reference_comparison_v1'
DATA.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)
protocol = dict(seed=20260916, rctd_seed=9,
    reference_donors=['CID4471', 'CID4535'], development_donor='CID4067',
    variants=['donor_balanced', 'epithelial_split', 'balanced_split'],
    balancing='Within each original type, sample equally from each donor without replacement. If balanced total <25, retain original pool. Single fixed seed.',
    epithelial_split='Author Cancer LumA SC and Cancer LumB SC become separate model labels; remaining cancer labels form Cancer Other. Aggregate all three back to Cancer Epithelial for 13-type evaluation.',
    fixed_spots='Wu2021_known_composition_v1, all 488, paired depths',
    rctd=dict(mode='full', max_cores=2, ref_UMI_min=100, ref_n_cells_min=25),
    nnls_features=2000, metrics=['MAE', 'RMSE', 'bias'],
    comparison='Exploratory reference development after inspecting CID4067 baseline; no independent validation donor in this comparison.')
(OUT / 'protocol.json').write_text(json.dumps(protocol, ensure_ascii=False, indent=2)+'\n')
ref = ad.read_h5ad(BASE / 'reference_counts.h5ad')
assert set(ref.obs.donor_id) == set(protocol['reference_donors'])
rng = np.random.default_rng(protocol['seed'])
selected = []
for _, group in ref.obs.groupby('reference_label', observed=True, sort=True):
    pools = [g.index.to_numpy() for _, g in group.groupby('donor_id', observed=True, sort=True)]
    n = min(map(len, pools))
    for pool in pools:
        selected.extend(rng.choice(pool, n, replace=False) if n*len(pools) >= 25 else pool)
assert len(selected) == len(set(selected))
inventory = []
for name in protocol['variants']:
    inp, out = DATA/name, OUT/name
    inp.mkdir(exist_ok=True); out.mkdir(exist_ok=True)
    a = ref[selected].copy() if name != 'epithelial_split' else ref.copy()
    a.obs['parent_label'] = a.obs.reference_label.astype(str)
    a.obs['reference_label'] = a.obs.reference_label.astype(str)
    if name != 'donor_balanced':
        cancer = a.obs.parent_label.eq('Cancer Epithelial')
        a.obs.loc[cancer, 'reference_label'] = a.obs.loc[cancer, 'celltype_minor'].astype(str).map(
            {'Cancer LumA SC': 'Cancer LumA', 'Cancer LumB SC': 'Cancer LumB'}).fillna('Cancer Other')
    meta = a.obs[['donor_id', 'parent_label', 'reference_label', 'celltype_minor']].copy()
    meta['nUMI'] = np.asarray(a.X.sum(axis=1)).ravel()
    assert meta.loc[meta.nUMI >= 100, 'reference_label'].value_counts().min() >= 25
    meta.to_csv(inp/'reference_metadata.csv', index_label='barcode')
    mmwrite(inp/'reference_counts.mtx', a.X.T)
    pd.Series(a.var_names).to_csv(inp/'genes.tsv', header=False, index=False, sep='\t')
    mapping = meta[['reference_label', 'parent_label']].drop_duplicates().set_index('reference_label')
    assert mapping.index.is_unique
    mapping.to_csv(out/'label_mapping.csv')
    inventory.append(meta.groupby(['donor_id','parent_label','reference_label'], observed=True).size().rename('n_cells').reset_index().assign(variant=name))
    normalized = a.X.astype(float).multiply((1/meta.nUMI.to_numpy())[:, None]).tocsr()
    types = sorted(meta.reference_label.unique())
    profiles = np.column_stack([np.asarray(normalized[meta.reference_label.eq(t).to_numpy()].mean(axis=0)).ravel() for t in types])
    score = profiles.var(axis=1)/(profiles.mean(axis=1)+1e-12)
    order = np.lexsort((a.var_names.to_numpy(), -score))
    keep = order[score[order] > 0][:2000]
    pd.DataFrame(profiles[keep], index=a.var_names[keep], columns=types).to_csv(inp/'nnls_reference_profiles.csv')
    print(name, a.shape, meta.reference_label.value_counts().to_dict(), flush=True)
pd.concat(inventory).to_csv(OUT/'reference_inventory.csv', index=False)
paths = [BASE/'spatial_counts.mtx', BASE/'spatial_metadata.csv'] + list(DATA.glob('*/*.csv')) + list(DATA.glob('*/*.mtx')) + list(DATA.glob('*/*.tsv'))
hashes = {}
for path in paths:
    digest = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''): digest.update(chunk)
    hashes[str(path.relative_to(ROOT))] = digest.hexdigest()
(OUT/'input_hashes.json').write_text(json.dumps(hashes, indent=2)+'\n')
