# 已保留的原始数据

`downloads/` 是原始下载包，`raw/` 是解压后的来源文件，`provenance/` 保留数据库元数据与下载记录。矩阵、作者注释、坐标和图像均以来源文件为准。

原始文件保留在本地；大型下载包和解压数据不会被 Git 默认跟踪。

## 数据清单

| 来源 | 当前本地内容 | 位置 | 计划用途 |
|---|---|---|---|
| Wu 2021，GSE176078 | 全队列原始计数与作者细胞元数据 | `downloads/wu2021_scrna/`；`raw/wu2021_scrna/Wu_etal_2021_BRCA_scRNASeq/` | 参考图谱、单细胞分析 |
| Wu 单样本，GSM5354538 | CID4535 下载包与解压文件，属于上面同一研究 | `downloads/wu2021_scrna_smoke/`；`raw/wu2021_scrna_smoke/` | 小规模读入练习，不是独立队列 |
| Wu Visium，Zenodo 4739739 | 6 张切片的计数、作者病理标签、坐标、H&E 和缩放信息 | `downloads/wu2021_visium/`；`raw/wu2021_visium/` | 开发与复现 |
| Honda 2024，GSE243022 | A1/B1/C1/D1 四个新样本的 H5/MTX 计数和空间资产 | `downloads/honda2024_visium/`；`raw/honda2024_visium/` | 独立空间案例候选 |
| Wang 2024 TNBC，Zenodo 14204217 | 作者 README、Clinical.tar 与解压内容、数据库文件清单 | `downloads/wang2024_tnbc/`；`raw/wang2024_tnbc/Clinical/`；`provenance/wang2024_tnbc.zenodo.json` | 外部筛选；表达与逐 spot 病理数据尚未下载 |
| Pal 2021，GSE161529 | GEO 系列/样本元数据 | `provenance/pal2021_scrna.family.soft.gz` | 独立参考候选；表达矩阵尚未下载 |
| Li 2025 CTA，Zenodo 15211538 | 数据库说明与文件清单 | `provenance/li2025_cta.zenodo.json` | 病理辅助验证候选；数据包尚未下载 |

## Wu 空间数据入口

```text
raw/wu2021_visium/
├── filtered_count_matrices/   各切片的 matrix、features、barcodes
├── metadata/                 作者元数据与病理分类
└── spatial/                  图像、spot 坐标与缩放信息
```

部分来源文件虽然以 `.gz` 结尾，实际内容可能是明文。读入时应检查格式，避免仅凭后缀判断压缩方式。

Wu 的 `T-cells` 大类含更细的 T/NK 等标签；CD8 分析应核对 minor/subset 字段。部分病理分类是癌、间质和淋巴细胞的混合标签，需要保留其原意。

## 独立数据使用前的检查

- **Honda：** GEO 原始包提供计数和标准空间资产，没有发现可直接用于本项目的独立病理分区表；需要核对附录或另行标注。四个新样本与论文复用的 Wu 六张切片应分开处理。
- **Wang：** 已读取的 Clinical.RDS 有 94 行，阵列映射有 281 行；论文的分析患者数还涉及纳入/排除，需要依据作者说明核对。README 区分 `Robjects/counts`（批次校正后）与 `Robjects/countsNonCorrected`（未作该批次校正），输入不能混用。
- **Pal：** 包含正常、肿瘤、淋巴结、全细胞和上皮富集等来源，需要逐样本筛选，并核对治疗与标签信息。

## 来源入口

- [Wu 单细胞：GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078)
- [Wu 空间：Zenodo 4739739](https://zenodo.org/records/4739739)
- [Honda 新空间样本：GSE243022](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243022)
- [Wang TNBC：Zenodo 14204217](https://zenodo.org/records/14204217)
- [Pal 单细胞：GSE161529](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161529)
- [Li CTA：Zenodo 15211538](https://zenodo.org/records/15211538)

来源记录也保存在 `provenance/`。引用、许可和后续下载以相应数据库及原论文的说明为准。
