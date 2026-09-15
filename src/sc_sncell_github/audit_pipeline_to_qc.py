"""复算当前阶段的关键输入、QC与重抽样指标，仅输出审查记录。"""
from pathlib import Path
import json
import ast
import hashlib
import math
import numpy as np
import pandas as pd
import anndata as ad
import nbformat
from scipy.io import mmread
from scipy.stats import spearmanr, pearsonr

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/pipeline_audit_20260915"
OUT.mkdir(parents=True, exist_ok=True)
checks = []


def check(name, condition, detail=""):
    assert bool(condition), name
    checks.append({"check": name, "status": "pass", "detail": str(detail)})


def same_matrix(a, b):
    return a.shape == b.shape and (a != b).nnz == 0


# 1. 文档结构、语法和保存的运行状态
for path in sorted((ROOT / "src/sc_sncell_github").glob("*.ipynb")):
    n = nbformat.read(path, as_version=4)
    nbformat.validate(n)
    code = [c for c in n.cells if c.cell_type == "code"]
    for c in code:
        ast.parse(c.source)
    check(path.name + ": saved execution", all(c.execution_count is not None for c in code)
          and not any(o.output_type == "error" for c in code for o in c.outputs), len(code))

# 2. 单细胞计数、作者标签和标准化对象
source = ROOT / "data/raw/wu2021_scrna/Wu_etal_2021_BRCA_scRNASeq"
metadata = pd.read_csv(source / "metadata.csv", index_col=0)
barcodes = pd.read_csv(source / "count_matrix_barcodes.tsv", header=None, sep="\t")[0]
check("full cohort barcode coverage", metadata.index.is_unique and barcodes.is_unique
      and set(barcodes) == set(metadata.index), len(metadata))
inventory = pd.read_csv(ROOT / "results/Wu2021_donor_reference_inventory/donor_by_reference_celltype.csv", index_col=0)
label = metadata.celltype_major.astype("string").copy()
label.loc[metadata.celltype_major.eq("T-cells")] = metadata.loc[metadata.celltype_major.eq("T-cells"), "celltype_minor"].astype("string")
count_table = pd.crosstab(metadata["orig.ident"], label).reindex(index=inventory.index, columns=inventory.columns[2:], fill_value=0)
check("donor inventory counts", np.array_equal(count_table, inventory.iloc[:, 2:]), "26 IDs, 13 labels")
check("donor subtype mapping", metadata.groupby("orig.ident").subtype.nunique().eq(1).all()
      and metadata.groupby("orig.ident").subtype.first().reindex(inventory.index).equals(inventory.subtype))
sc = ad.read_h5ad(ROOT / "data/processed/CID4535_baseline_v1/single_cell_analysis.h5ad")
raw = ad.read_h5ad(ROOT / "data/processed/CID4535_baseline_v1/qc_counts.h5ad")
check("single cell counts layer", same_matrix(sc.layers["counts"], raw.X))
for col in ["celltype_major", "celltype_minor", "celltype_subset"]:
    check("single cell author " + col, sc.obs[col].astype(str).equals(metadata.loc[sc.obs_names, col].astype(str)))
expected = sc.layers["counts"].astype(float).copy()
total = np.asarray(expected.sum(axis=1)).ravel()
expected = expected.multiply((10000 / total)[:, None]).tocsr()
expected.data = np.log1p(expected.data)
check("single cell log normalization", np.max(np.abs((expected - sc.X).data), initial=0) < 1e-5)
check("single cell PCA UMAP finite", np.isfinite(sc.obsm["X_pca"]).all() and np.isfinite(sc.obsm["X_umap"]).all())
check("single cell HVG count", sc.var.highly_variable.sum() == 3000)

# 3. 空间矩阵从原始文件到过滤对象
sp = ad.read_h5ad(ROOT / "data/processed/CID4535_spatial/spatial_qc_filtering.h5ad")
sp_root = ROOT / "data/raw/wu2021_visium"
sp_meta = pd.read_csv(sp_root / "metadata/CID4535_metadata.csv", index_col=0)
pos = pd.read_csv(sp_root / "spatial/CID4535_spatial/tissue_positions_list.csv", header=None, index_col=0)
raw_dir = sp_root / "filtered_count_matrices/CID4535_filtered_count_matrix"
sp_barcodes = pd.read_csv(raw_dir / "barcodes.tsv.gz", header=None, compression=None)[0]
with (raw_dir / "matrix.mtx.gz").open("rb") as handle:
    sp_matrix = mmread(handle, spmatrix=True).T.tocsr()
retained = (~sp_meta.loc[sp_barcodes, "Classification"].eq("Artefact").fillna(False)
            & pos.loc[sp_barcodes, 1].ne(0)).to_numpy()
check("spatial retained barcodes", list(sp_barcodes[retained]) == list(sp.obs_names), "1127 -> 1102")
check("spatial retained raw counts", same_matrix(sp_matrix[retained], sp.X))
check("spatial coordinate order", np.array_equal(sp.obsm["spatial"], pos.loc[sp.obs_names, [5, 4]].to_numpy()))
check("spatial pathology annotations", sp.obs.Classification.astype("string").fillna("Missing").equals(
    sp_meta.loc[sp.obs_names, "Classification"].astype("string").fillna("Missing")))
check("spatial retained missing pathology", sp.obs.Classification.isna().sum() == 1)

# 4. 共同基因对象和R导出文件
combined = ROOT / "data/processed/CID4535_combine_v1"
ref = ad.read_h5ad(combined / "single_cell_reference_counts.h5ad")
spi = ad.read_h5ad(combined / "spatial_common_gene_counts.h5ad")
check("shared genes aligned", ref.var_names.equals(spi.var_names) and ref.var_names.is_unique, ref.n_vars)
check("reference raw source", same_matrix(ref.X, sc[:, ref.var_names].layers["counts"]))
check("spatial shared gene source", same_matrix(spi.X, sp[:, spi.var_names].X))
export = combined / "rctd_input"
export_genes = pd.read_csv(export / "genes.tsv", header=None)[0]
for name, obj in [("reference", ref), ("spatial", spi)]:
    exported = mmread(export / f"{name}_counts.mtx", spmatrix=True).T.tocsr()
    em = pd.read_csv(export / f"{name}_metadata.csv", index_col=0)
    check(name + " export counts", same_matrix(exported, obj.X))
    check(name + " export IDs", list(em.index) == list(obj.obs_names) and list(export_genes) == list(obj.var_names))
    check(name + " count validity", np.isfinite(obj.X.data).all() and (obj.X.data >= 0).all() and (obj.X.data == np.floor(obj.X.data)).all())

# 5. 独立重算全部重抽样数值
base_dir = ROOT / "data/processed/CID4535_rctd_v1"
boot_dir = ROOT / "data/processed/CID4535_reference_bootstrap_v1"
metrics_dir = ROOT / "results/CID4535_reference_bootstrap_v1"
b = pd.read_csv(base_dir / "rctd_weights_spot_by_celltype.csv", index_col=0).loc[sp.obs_names]
b = b.div(b.sum(axis=1), axis=0)
stored_metrics = pd.read_csv(metrics_dir / "metrics_by_run.csv").set_index(["seed", "cell_type"])
stored_overlap = pd.read_csv(metrics_dir / "high_spot_overlap.csv").set_index(["seed", "cell_type", "top_fraction"])
arrays = []
freq = np.zeros(b.shape)
max_metric_delta = 0.0
for seed in range(1001, 1011):
    w = pd.read_csv(boot_dir / f"seed_{seed}/weights.csv", index_col=0).loc[b.index, b.columns]
    check(f"bootstrap {seed} counts and values", w.shape == (1102, 13) and np.isfinite(w.to_numpy()).all()
          and (w.to_numpy() >= 0).all() and (w.sum(axis=1) > 0).all())
    sampled = pd.read_csv(boot_dir / f"seed_{seed}/sampled_cells.csv", index_col=0)
    expected_types = ref.obs.reference_label.value_counts().sort_index()
    got_types = sampled.reference_label.value_counts().sort_index()
    check(f"bootstrap {seed} sample strata", expected_types.equals(got_types)
          and sampled.index.is_unique and sampled.source_barcode.isin(ref.obs_names).all())
    check(f"bootstrap {seed} source labels", all(
        ref.obs.reference_label.astype(str).loc[sampled.source_barcode].to_numpy() == sampled.reference_label.to_numpy()))
    w = w.div(w.sum(axis=1), axis=0)
    arrays.append(w.to_numpy())
    for j, t in enumerate(b.columns):
        stored = stored_metrics.loc[(seed, t)]
        computed = [spearmanr(b[t], w[t]).statistic, pearsonr(b[t], w[t]).statistic,
                    np.abs(b[t] - w[t]).mean() * 100, (w[t] - b[t]).mean() * 100]
        max_metric_delta = max(max_metric_delta, np.max(np.abs(np.array(computed) - stored[["spearman", "pearson", "mae_pp", "bias_pp"]].to_numpy())))
        for fraction in [0.05, 0.1, 0.2]:
            k = math.ceil(len(b) * fraction)
            ib = np.lexsort((b.index.to_numpy(), -b[t].to_numpy()))[:k]
            iw = np.lexsort((b.index.to_numpy(), -w[t].to_numpy()))[:k]
            overlap = len(set(ib) & set(iw))
            saved = stored_overlap.loc[(seed, t, fraction)]
            assert overlap == saved.n_overlap
            assert abs(overlap / k - saved.retention) < 1e-12
            if fraction == 0.1:
                freq[iw, j] += 0.1
check("bootstrap metrics recomputation", max_metric_delta < 1e-10, max_metric_delta)
check("all 390 hotspot comparisons", True)
cube = np.stack(arrays)
spots = pd.read_csv(metrics_dir / "spot_stability.csv").set_index(["barcode", "cell_type"])
for j, t in enumerate(b.columns):
    saved = spots.xs(t, level="cell_type").loc[b.index]
    np.testing.assert_allclose(saved.bootstrap_mean, cube[:, :, j].mean(axis=0), atol=1e-12)
    np.testing.assert_allclose(saved.bootstrap_sd, cube[:, :, j].std(axis=0, ddof=1), atol=1e-12)
    np.testing.assert_allclose(saved.top10_frequency, freq[:, j], atol=1e-12)
check("all 14326 spot stability records", True)

# 6. 三供者QC与复核使用表
qc = pd.read_csv(ROOT / "results/Wu2021_ER_donor_qc_v1/cell_qc_records.csv", index_col=0)
for donor in ["CID4471", "CID4067", "CID4535"]:
    obj = ad.read_h5ad(ROOT / f"data/processed/Wu2021_ER_donor_qc_v1/{donor}_counts_after_qc.h5ad")
    before = ad.read_h5ad(ROOT / f"data/processed/Wu2021_ER_donor_qc_v1/{donor}_counts_before_qc.h5ad")
    check(donor + " before/after cell and counts", before.obs_names.equals(obj.obs_names)
          and same_matrix(before[:, obj.var_names].X, obj.X))
    check(donor + " gene filter", ((obj.X > 0).sum(axis=0) >= 1).all())
    total = np.asarray(obj.X.sum(axis=1)).ravel()
    n_genes = np.asarray((obj.X > 0).sum(axis=1)).ravel()
    mt = np.asarray(obj.X[:, obj.var.gene_symbol.str.startswith("MT-").to_numpy()].sum(axis=1)).ravel() / total * 100
    records = qc.loc[obj.obs_names]
    check(donor + " QC recomputation", np.array_equal(total, records.total_counts) and np.array_equal(n_genes, records.n_genes_by_counts)
          and np.allclose(mt, records.pct_counts_mt, atol=1e-10))
    check(donor + " author labels", obj.obs.celltype_subset.astype(str).equals(metadata.loc[obj.obs_names, "celltype_subset"].astype(str)))
    panels = {"epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
              "T_lineage": ["CD3D", "CD3E", "TRAC"],
              "CAF": ["COL1A1", "COL1A2", "DCN", "LUM"],
              "PVL": ["RGS5", "CSPG4", "MCAM"]}
    panel_hits = {k: (obj[:, genes].X.toarray() >= 2).sum(axis=1) >= 2
                  for k, genes in panels.items()}
    check(donor + " panel coexpression flags", np.array_equal(
        panel_hits["epithelial"] & panel_hits["T_lineage"], records.review_epithelial_T)
        and np.array_equal(panel_hits["CAF"] & panel_hits["PVL"], records.review_CAF_PVL))
for flag, metric, high in [("review_low_umi", "total_counts", False),
                          ("review_high_umi", "total_counts", True),
                          ("review_low_genes", "n_genes_by_counts", False),
                          ("review_high_genes", "n_genes_by_counts", True),
                          ("review_high_mt", "pct_counts_mt", True)]:
    expected_flags = pd.Series(False, index=qc.index)
    for _, group in qc.groupby(["donor_id", "reference_label"]):
        values = group[metric].to_numpy(dtype=float)
        if metric != "pct_counts_mt":
            values = np.log1p(values)
        center = np.median(values)
        mad = np.median(np.abs(values - center))
        if len(group) >= 20 and mad > 0:
            expected_flags.loc[group.index] = values > center + 5 * mad if high else values < center - 5 * mad
    check("recomputed " + flag, expected_flags.equals(qc[flag]))
use = pd.read_csv(ROOT / "results/Wu2021_ER_target_review_v1/target_cell_use_table.csv", index_col=0)
check("target cell use records", len(use) == 3989 and use.index.is_unique and use.review_any.sum() == 303
      and use.exclude_after_review.sum() == 0)

pd.DataFrame(checks).to_csv(OUT / "calculation_checks.csv", index=False)
(OUT / "check_summary.json").write_text(json.dumps({"passed_checks": len(checks), "status": "calculation_checks_passed"}, indent=2))
print("通过核对项：", len(checks))
print("重抽样指标最大差异：", max_metric_delta)
