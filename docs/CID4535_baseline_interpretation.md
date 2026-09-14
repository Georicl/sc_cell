# CID4535 单细胞基线判读

当前采用简单的 Scanpy 单样本流程，用于展示主要细胞组成并为后续空间分析准备。初始结果在这个用途下合理，尚不支持把每个聚类视为独立细胞类型。

## 固定参数

| 环节 | 设置 |
|---|---|
| 输入 | 作者处理后的 CID4535 原始 UMI 矩阵及 metadata |
| 细胞过滤 | 不额外过滤，保留 3,961 个细胞 |
| 基因过滤 | 保留至少在 1 个细胞中表达的基因，22,061 个 |
| 归一化 | 每细胞 10,000，随后 log1p |
| 高变基因 | Scanpy seurat 离散度方法，3,000 个 |
| PCA | 中心化，不做单位方差缩放；计算 50 个，使用前 30 个 |
| 邻居图 | 15 个邻居 |
| Leiden | igraph，resolution=0.5，n_iterations=2 |
| UMAP | min_dist=0.5；仅用于表达空间展示 |
| 随机种子 | 0 |
| Marker | 按作者 celltype_major 分组；全部保留基因的 log 表达，Wilcoxon、BH 校正，并输出表达比例 |

## 初始结果是否合理

- 重新计算的 UMI 数和检测基因数与作者 metadata 逐细胞一致，未发现矩阵和注释错配。
- 16 个聚类中，内皮、B 细胞、浆母细胞、髓系及上皮等主要谱系有 marker 与作者标签共同支持。
- 原 cluster 13 的增殖相关基因突出，适合解释为增殖状态，不必单独当作新细胞类型。
- 原 cluster 5 有上皮/PVL 特征冲突；cluster 8 主要是低计数的混合群；cluster 15 混有作者标注的正常和癌上皮。保留为待复核，不自动删除或强制定名。
- CAF 与 PVL、CD4/CD8 T 与 NK/NKT 尚未充分分开。这不妨碍当前的大类展示，但不能直接据此定义 CAF/CD8 空间分析所需的细分参考群。

作者标签下，癌上皮为 2,223 个，PVL 为 592 个，T-cells 为 396 个，髓系为 255 个，内皮为 219 个，CAF 为 102 个，浆母细胞为 96 个，B 细胞为 56 个，正常上皮为 22 个。图中明确写明这些类型来自作者；细胞捕获比例不等于组织面积比例。

## 本轮调整

主 notebook 从 resolution=0.8 恢复为 0.5；0.8 对照未改善主要混合问题。其余核心分析参数沿用初始基线，不追加回归、批次整合、严格 QC 或自动双细胞删除。

直接采用作者 metadata 中的 celltype_major、celltype_minor、celltype_subset，逐细胞核对并保存来源，不再创建空白人工注释。Marker 检验和点图按作者大类分组，用于表达合理性检查；基因组合的分组标题仅描述检查用途，不代表新亚群。主要 UMAP 对照 Leiden 聚类与作者标签，另保存作者标签下的细胞数量与比例。

## 原文方法的关系

Wu 2021 不只是调整聚类分辨率。作者公开代码包含标准化、部分步骤的协变量回归、跨患者谱系整合及根据低计数/异谱系信号等证据排除问题群后的重分析。本项目当前是单患者的简化展示，不能称为论文全部方法的严格复现。

- [原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC9044823/)
- [作者免疫重分析代码，固定 commit](https://github.com/Swarbricklab-code/BrCa_cell_atlas/blob/4ef33fc58e5cae97f9bf51e58d45bb49688274f6/reclustering_analysis/immune_cells/01B_reclustering_with_filtering_cells.R)
- [作者问题群过滤与重分析](https://github.com/Swarbricklab-code/BrCa_cell_atlas/blob/4ef33fc58e5cae97f9bf51e58d45bb49688274f6/reclustering_analysis/immune_cells/02B_filtering_doublet_clusters_and_reprocessing.R)

## 运行与后续

打开 src/sc_sncell_github/read_data.ipynb，选择项目 .venv 内核，按顺序运行。当前输出写入 results/CID4535_baseline_review_v1 和 data/processed/CID4535_baseline_review_v1；原 baseline 和 0.8 对照保留。

下一步读取 CID4535 的空间计数、坐标与组织图像，检查 spot 对齐和质量。细胞类型参考优先保留来源明确的作者标签；在正式开展 CAF–CD8 关联前，再明确参考样本范围和目标细胞定义。
