# 分析文件使用

当前保留 notebook 作为单细胞与空间分析入口，R 脚本负责 RCTD 拟合，Python 脚本负责论文成图和临床信息整理。

## 1. 环境与输入

在项目根目录运行 `uv sync`，然后选择 `.venv` 的 Python 内核。原始数据下载入口见 [data/README.md](../data/README.md)，下载和派生对象目录保留在本地。

R 分析使用本项目记录的 spacexr 1.2.0 接口，并加载 Matrix、SummarizedExperiment 和 SpatialExperiment；版本记录见 [sessionInfo.txt](../results/CID4535_rctd_v1/sessionInfo.txt)。Python 的 `uv.lock` 不管理 R 包。

原 notebook 和部分 R 脚本保留了当前作者本机的绝对 `project_dir`，在其他机器运行时应改为实际项目根目录。使用 `../../` 的 notebook 应从 `src/sc_sncell_github` 目录运行。下面列出的文件均为当前实际存在的入口。

## 2. 按顺序运行

| 文件 | 主要输出 |
|---|---|
| `read_data.ipynb` | 单细胞 QC、降维、作者注释与原始 counts 层 |
| `read_spatial.ipynb` | 空间坐标和图像对齐、过滤记录与空间对象 |
| `prepare_reference.ipynb` | 共同基因参考、13类标签与 RCTD 输入 |
| `run_rctd.R` | 原生 RCTD 结果及逐 spot 权重 |
| `insect_rctd.ipynb` | 反卷积权重的空间展示（保留现有文件名） |
| `rctd_reference_bootstrap.R` | 10次参考重抽样结果与抽样条码 |
| `inspect_reference_bootstrap.ipynb` | 稳定性指标、图形及结果解读 |

在项目根目录运行 R 脚本：

```bash
Rscript src/sc_sncell_github/run_rctd.R
Rscript src/sc_sncell_github/rctd_reference_bootstrap.R
```

完成上述分析后重绘论文图：

```bash
.venv/bin/python src/sc_sncell_github/make_cid4535_paper_figures.py
```

图组输出到 `results/CID4535_paper_v1`，包括 PDF、PNG 和成图 CSV。该命令读取本地已生成的单细胞及空间对象，不会重新运行 RCTD。

## 3. 后续队列规划文件

`donor_reference_inventory.ipynb` 统计 Wu 全队列作者标签；`build_er_donor_plan.py` 对照 GEO 和临床补充表生成 ER+ 候选纳入表。后者另需 `openpyxl`，原始补充表和 SOFT 来源 URL 记录在 `results/Wu2021_ER_donor_plan/provenance.json`，按其中相对路径放置后运行。

## 4. 公开结果与复现记录

仓库保留正式图组、支撑结果表、单细胞参数、RCTD权重与R版本信息，以及10轮的实际抽样条码和参数记录。原始矩阵、`.h5ad`、原生模型 `.rds` 和本地环境不随 Git 提交；它们由上述流程在本地生成。
