"""从Wu全队列矩阵分块提取三例，保留作者细胞条码与注释。"""
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

PROJECT = Path(__file__).resolve().parents[2]
SOURCE = PROJECT / "data/raw/wu2021_scrna/Wu_etal_2021_BRCA_scRNASeq"
OUTPUT = PROJECT / "data/processed/Wu2021_ER_donor_qc_v1"
DONORS = ["CID4471", "CID4067", "CID4535"]


def extract_donors():
    # 1. 根据条码确定需要提取的细胞列
    metadata = pd.read_csv(SOURCE / "metadata.csv", index_col=0)
    barcodes = pd.read_csv(SOURCE / "count_matrix_barcodes.tsv", sep="\t", header=None)[0].astype(str)
    genes = pd.read_csv(SOURCE / "count_matrix_genes.tsv", sep="\t", header=None)[0].astype(str)
    assert metadata.index.is_unique and barcodes.is_unique and genes.is_unique
    metadata = metadata.loc[barcodes].copy()
    selected = np.flatnonzero(metadata["orig.ident"].isin(DONORS).to_numpy())
    lookup = np.full(len(barcodes) + 1, -1, dtype=np.int32)
    lookup[selected + 1] = np.arange(len(selected), dtype=np.int32)
    selected_metadata = metadata.iloc[selected].copy()

    # 2. 分块读取坐标格式矩阵，每次只保留这三例的条目
    matrix_path = SOURCE / "count_matrix_sparse.mtx"
    header_lines = 0
    with matrix_path.open() as stream:
        for line in stream:
            header_lines += 1
            if not line.startswith("%"):
                n_genes, n_cells, n_entries = map(int, line.split())
                break
    assert n_genes == len(genes) and n_cells == len(barcodes)
    rows, cols, values = [], [], []
    observed_entries = 0
    for number, chunk in enumerate(pd.read_csv(
        matrix_path, sep=r"\s+", header=None, names=["gene", "cell", "count"],
        skiprows=header_lines, dtype=np.int32, chunksize=2_000_000
    ), 1):
        new_rows = lookup[chunk["cell"].to_numpy()]
        keep = new_rows >= 0
        rows.append(new_rows[keep])
        cols.append(chunk["gene"].to_numpy()[keep] - 1)
        values.append(chunk["count"].to_numpy()[keep])
        observed_entries += len(chunk)
        if number % 10 == 0:
            print(f"已读取 {observed_entries:,}/{n_entries:,} 条计数", flush=True)
    assert observed_entries == n_entries
    matrix = sparse.coo_matrix(
        (np.concatenate(values), (np.concatenate(rows), np.concatenate(cols))),
        shape=(len(selected), n_genes)
    ).tocsr()
    assert matrix.data.min() >= 0

    # 3. 保存各供者原始counts对象，逐细胞核对作者QC计数
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for donor in DONORS:
        keep = selected_metadata["orig.ident"].eq(donor).to_numpy()
        a = ad.AnnData(
            X=matrix[keep].copy(), obs=selected_metadata.loc[keep].copy(),
            var=pd.DataFrame({"gene_symbol": genes.to_numpy()}, index=pd.Index(genes, name="gene"))
        )
        a.obs.index.name = "barcode"
        np.testing.assert_array_equal(np.asarray(a.X.sum(axis=1)).ravel(), a.obs["nCount_RNA"])
        np.testing.assert_array_equal(np.asarray((a.X > 0).sum(axis=1)).ravel(), a.obs["nFeature_RNA"])
        a.uns["source_matrix"] = str(matrix_path.relative_to(PROJECT))
        a.uns["source_metadata"] = str((SOURCE / "metadata.csv").relative_to(PROJECT))
        a.write_h5ad(OUTPUT / f"{donor}_counts_before_qc.h5ad", compression="gzip")
        print(donor, a.shape, "计数与作者metadata一致", flush=True)


if __name__ == "__main__":
    extract_donors()
