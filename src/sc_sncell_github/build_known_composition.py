"""固定供者划分，生成CID4067的已知组成测试输入。"""
from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.io import mmwrite

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data/processed/Wu2021_known_composition_v1"
RESULT = ROOT / "results/Wu2021_known_composition_v1"
DATA.mkdir(parents=True, exist_ok=True)
RESULT.mkdir(parents=True, exist_ok=True)
TRAIN = ["CID4471", "CID4535"]
TEST = "CID4067"
SEED = 20260915

# 1. 先保存本轮配方与参数，再生成预测输入
protocol = {
    "reference_donors": TRAIN, "test_donor": TEST,
    "sampling_seed": SEED, "rctd_seed": 9,
    "factorial_CAF_cells": [0, 1, 3], "factorial_CD8_cells": [0, 1, 3],
    "factorial_cells_per_spot": 10, "fixed_PVL_cells": 1,
    "background_types": ["T cells CD4+", "Myeloid", "Endothelial", "B-cells", "NK cells", "NKT cells", "Cycling T-cells"],
    "background_cells": 1, "repeats_per_recipe": 3,
    "remaining_cells": "Cancer Epithelial",
    "pure_control_cells_per_spot": 5, "pure_repeats_per_present_type": 5,
    "depth_fractions": [1.0, 0.25],
    "sampling": "without replacement within each spot; cells can recur across spots",
    "thinning": "independent binomial thinning of each component cell count",
    "primary_target": "fraction of captured RNA UMIs contributed by each type after thinning",
    "secondary_target": "known cell-number fraction",
    "primary_metrics": ["MAE", "RMSE", "bias"],
    "diagnostic_positive_weight_threshold": 0.01,
    "nnls_feature_count": 2000,
    "nnls_gene_selection": "between-reference-type variance / mean of full-library-normalized profiles, reference only",
    "rctd": {"rctd_mode": "full", "max_cores": 2, "ref_UMI_min": 100, "ref_n_cells_min": 25},
    "reference_weighting": "pool cells within each reference label, matching initial RCTD reference construction",
    "gene_universe": "shared measured gene list from original Wu count matrix; zero rows retained at input",
    "independent_test_donors": 1,
}
(RESULT / "protocol.json").write_text(json.dumps(protocol, ensure_ascii=False, indent=2), encoding="utf-8")

# 2. 使用QC保留条码读取原始计数，各样本保留原测量基因全集
qc = pd.read_csv(ROOT / "results/Wu2021_ER_donor_qc_v1/cell_qc_records.csv", index_col=0)
objects = {}
for donor in TRAIN + [TEST]:
    p = ROOT / f"data/processed/Wu2021_ER_donor_qc_v1/{donor}_counts_before_qc.h5ad"
    a = ad.read_h5ad(p)
    use = qc.loc[a.obs_names, "qc_keep"].to_numpy()
    a = a[use].copy()
    a.obs["reference_label"] = qc.loc[a.obs_names, "reference_label"].astype(str)
    a.obs["donor_id"] = donor
    objects[donor] = a
genes = objects[TRAIN[0]].var_names
assert all(a.var_names.equals(genes) for a in objects.values())
reference = ad.concat([objects[d] for d in TRAIN], join="inner", merge="same")
test = objects[TEST]
assert set(reference.obs_names).isdisjoint(test.obs_names)
types = sorted(reference.obs.reference_label.unique())
assert set(test.obs.reference_label).issubset(types)
assert reference.obs.reference_label.value_counts().min() >= 25
split = pd.concat([
    reference.obs[["donor_id", "reference_label"]].assign(role="reference"),
    test.obs[["donor_id", "reference_label"]].assign(role="heldout_test"),
])
split.to_csv(RESULT / "cell_split.csv", index_label="barcode")
split.groupby(["role", "donor_id", "reference_label"]).size().rename("n_cells").to_csv(RESULT / "split_counts.csv")

# 3. 建立配方，CAF和CD8组合之外的细胞填入癌上皮及其他背景类型
recipes = []
for n_caf in protocol["factorial_CAF_cells"]:
    for n_cd8 in protocol["factorial_CD8_cells"]:
        for background in protocol["background_types"]:
            composition = {"CAFs": n_caf, "T cells CD8+": n_cd8, "PVL": 1,
                           background: 1, "Cancer Epithelial": 8 - n_caf - n_cd8}
            for repeat in range(3):
                recipes.append(("factorial", f"caf{n_caf}_cd8{n_cd8}_{background}", repeat, composition))
for cell_type in sorted(test.obs.reference_label.unique()):
    for repeat in range(5):
        recipes.append(("pure", cell_type, repeat, {cell_type: 5}))
rng = np.random.default_rng(SEED)
pools = {t: np.flatnonzero(test.obs.reference_label.eq(t).to_numpy()) for t in types}
type_index = {t: i for i, t in enumerate(types)}
counts_rows, truth_cells, truth_rna, spot_rows, memberships = [], [], [], [], []
for number, (scenario, recipe_name, repeat, composition) in enumerate(recipes):
    selected = []
    for cell_type, n in composition.items():
        if n:
            selected.extend(rng.choice(pools[cell_type], size=n, replace=False).tolist())
    assert len(selected) == len(set(selected))
    labels = test.obs.reference_label.iloc[selected].to_numpy()
    cell_fraction = np.array([(labels == t).sum() / len(selected) for t in types])
    source_counts = test.X[selected].tocsr()
    for depth in protocol["depth_fractions"]:
        spot_id = f"mix_{number:04d}_depth{int(depth * 100)}"
        component_counts = source_counts.copy()
        if depth < 1:
            component_counts.data = rng.binomial(component_counts.data, depth)
            component_counts.eliminate_zeros()
        umi_by_cell = np.asarray(component_counts.sum(axis=1)).ravel()
        n_umi = int(umi_by_cell.sum())
        assert n_umi > 0
        rna_by_type = np.zeros(len(types))
        for k, (cell_i, label, umi) in enumerate(zip(selected, labels, umi_by_cell)):
            rna_by_type[type_index[label]] += umi
            memberships.append({"spot_id": spot_id, "base_mixture": number,
                                "source_barcode": test.obs_names[cell_i], "donor_id": TEST,
                                "reference_label": label, "component_umi": int(umi)})
        counts_rows.append(sparse.csr_matrix(component_counts.sum(axis=0)))
        truth_cells.append(cell_fraction)
        truth_rna.append(rna_by_type / n_umi)
        spot_rows.append({"spot_id": spot_id, "base_mixture": number, "scenario": scenario,
                          "recipe": recipe_name, "repeat": repeat, "depth_fraction": depth,
                          "n_cells": len(selected), "nUMI": n_umi,
                          "x": len(spot_rows), "y": 0, "coordinate_type": "synthetic_index"})
spots = pd.DataFrame(spot_rows).set_index("spot_id")
counts_matrix = sparse.vstack(counts_rows, format="csr")
cell_truth = pd.DataFrame(truth_cells, index=spots.index, columns=types)
rna_truth = pd.DataFrame(truth_rna, index=spots.index, columns=types)
assert len(spots) == 488 and spots.base_mixture.nunique() == 244
assert np.allclose(cell_truth.sum(axis=1), 1) and np.allclose(rna_truth.sum(axis=1), 1)
assert np.array_equal(np.asarray(counts_matrix.sum(axis=1)).ravel(), spots.nUMI)
pd.DataFrame(memberships).to_csv(RESULT / "mixture_membership.csv", index=False)
cell_truth.to_csv(RESULT / "truth_cell_fraction.csv")
rna_truth.to_csv(RESULT / "truth_rna_fraction.csv")
spots.to_csv(RESULT / "simulation_design.csv")

# 4. 保存两个方法使用的输入，真值文件与模型输入分开
reference_meta = reference.obs[["donor_id", "reference_label"]].copy()
reference_meta["nUMI"] = np.asarray(reference.X.sum(axis=1)).ravel()
reference_meta.to_csv(DATA / "reference_metadata.csv", index_label="barcode")
spots[["x", "y", "nUMI"]].to_csv(DATA / "spatial_metadata.csv")
pd.Series(genes).to_csv(DATA / "genes.tsv", sep="\t", header=False, index=False)
mmwrite(DATA / "reference_counts.mtx", reference.X.T)
mmwrite(DATA / "spatial_counts.mtx", counts_matrix.T)
reference.write_h5ad(DATA / "reference_counts.h5ad", compression="gzip")
ad.AnnData(X=counts_matrix, obs=spots.copy(), var=pd.DataFrame(index=genes)).write_h5ad(
    DATA / "simulated_spots.h5ad", compression="gzip")

# 5. NNLS参考特征由参考侧计算，使用全库UMI归一化
normalized = reference.X.astype(float).multiply((1 / reference_meta.nUMI.to_numpy())[:, None]).tocsr()
profiles = np.column_stack([
    np.asarray(normalized[reference.obs.reference_label.eq(t).to_numpy()].mean(axis=0)).ravel()
    for t in types
])
score = profiles.var(axis=1) / (profiles.mean(axis=1) + 1e-12)
order = np.lexsort((genes.to_numpy(), -score))
selected_genes = order[score[order] > 0][:protocol["nnls_feature_count"]]
pd.DataFrame(profiles[selected_genes], index=genes[selected_genes], columns=types).to_csv(DATA / "nnls_reference_profiles.csv")
input_paths = [DATA / n for n in ["reference_counts.mtx", "spatial_counts.mtx", "genes.tsv", "reference_metadata.csv", "spatial_metadata.csv"]]
(RESULT / "input_hashes.json").write_text(json.dumps({p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in input_paths}, indent=2))
print("Reference:", reference.shape, "test donor:", test.shape)
print("Simulated spots:", len(spots), "base mixtures:", spots.base_mixture.nunique())
print(spots.groupby(["scenario", "depth_fraction"]).nUMI.agg(["size", "min", "median", "max"]))
