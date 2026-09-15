"""从来源细胞和保存的预测独立核对模拟测试。"""
from pathlib import Path
import hashlib
import json
import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
data = ROOT / 'data/processed/Wu2021_known_composition_v1'
out = ROOT / 'results/Wu2021_known_composition_v1'
spots = ad.read_h5ad(data / 'simulated_spots.h5ad')
ref = ad.read_h5ad(data / 'reference_counts.h5ad')
test = ad.read_h5ad(ROOT / 'data/processed/Wu2021_ER_donor_qc_v1/CID4067_counts_before_qc.h5ad')
design = pd.read_csv(out / 'simulation_design.csv', index_col=0)
members = pd.read_csv(out / 'mixture_membership.csv')
truth = pd.read_csv(out / 'truth_rna_fraction.csv', index_col=0)
cells = pd.read_csv(out / 'truth_cell_fraction.csv', index_col=0)
checks = {}
checks['gene_order'] = ref.var_names.equals(spots.var_names) and test.var_names.equals(spots.var_names)
checks['no_train_test_overlap'] = not bool(set(ref.obs_names) & set(test.obs_names))
checks['no_duplicate_cells_within_spot'] = not members.duplicated(['spot_id', 'source_barcode']).any()
checks['test_donor_only'] = members.donor_id.eq('CID4067').all()
for label, values, agg, target in [('rna', 'component_umi', 'sum', truth), ('cells', 'source_barcode', 'count', cells)]:
    table = members.pivot_table(index='spot_id', columns='reference_label', values=values, aggfunc=agg, fill_value=0)
    table = table.reindex(index=target.index, columns=target.columns, fill_value=0)
    checks[label + '_truth'] = np.allclose(table.div(table.sum(axis=1), axis=0), target)
checks['component_umi_totals'] = np.array_equal(members.groupby('spot_id').component_umi.sum().reindex(design.index), design.nUMI)
checks['matrix_umi_totals'] = np.array_equal(np.asarray(spots.X.sum(axis=1)).ravel(), design.nUMI)
exact = paired = thinning = True
for _, pair in design.groupby('base_mixture'):
    original = pair.index[pair.depth_fraction.eq(1.0)][0]
    low = pair.index[pair.depth_fraction.eq(0.25)][0]
    barcodes = members.loc[members.spot_id.eq(original), 'source_barcode']
    expected = np.asarray(test[barcodes].X.sum(axis=0)).ravel()
    actual = spots[original].X.toarray().ravel()
    exact &= np.array_equal(expected, actual)
    paired &= set(barcodes) == set(members.loc[members.spot_id.eq(low), 'source_barcode'])
    thinning &= bool((spots[low].X.toarray().ravel() <= actual).all())
checks.update(original_244_count_reconstruction=exact, paired_source_cells=paired, thinned_counts_bounded=thinning)
for name, expected in json.loads((out / 'input_hashes.json').read_text()).items():
    digest = hashlib.sha256()
    with (data / name).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    checks['input_hash_' + name] = digest.hexdigest() == expected
metrics = pd.read_csv(out / 'benchmark_metrics.csv')
for method in ['RCTD', 'NNLS']:
    pred = pd.read_csv(out / f'{method.lower()}_relative_weights.csv', index_col=0).loc[truth.index, truth.columns]
    checks[method + '_valid_weights'] = np.isfinite(pred.values).all() and (pred.values >= 0).all() and np.allclose(pred.sum(axis=1), 1)
    errors = []
    for row in metrics[metrics.method.eq(method)].itertuples():
        ids = design.index[design.scenario.eq(row.scenario) & design.depth_fraction.eq(row.depth_fraction)]
        target = truth if row.truth == 'rna' else cells
        e = (pred.loc[ids, row.cell_type] - target.loc[ids, row.cell_type]).to_numpy() * 100
        errors.extend([abs(np.abs(e).mean() - row.mae_pp), abs(np.sqrt((e ** 2).mean()) - row.rmse_pp), abs(e.mean() - row.bias_pp)])
    checks[method + '_metric_recalculation'] = max(errors) < 1e-10
checks = {key: bool(value) for key, value in checks.items()}
(out / 'validation_checks.json').write_text(json.dumps(checks, indent=2) + '\n')
print(json.dumps(checks, indent=2))
assert all(checks.values())
