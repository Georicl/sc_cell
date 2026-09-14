# Wu ER+供者纳入与用途表（草案v1）

## 1. 本次核对范围

以作者细胞元数据中11个ER+样本为范围，对照GEO逐样本记录与原论文Supplementary Table 1。纳入盘点不等于已纳入正式模型。只将CID4535记录为已用于开发，其他供者的用途均为建议，未冻结测试集。

## 2. 来源、组织学、治疗和目标细胞覆盖

| donor_id | GEO样本 | 组织学 | 治疗状态 | 总细胞 | CAF | CD8 T |
| --- | --- | --- | --- | --- | --- | --- |
| CID3941 | GSM5354516 | IDC | Naïve | 631 | 8 | 126 |
| CID3948 | GSM5354518 | IDC | Naïve | 2327 | 15 | 498 |
| CID4040 | GSM5354520 | IDC | Naïve | 2531 | 129 | 528 |
| CID4067 | GSM5354522 | IDC | Naïve | 3764 | 135 | 243 |
| CID4290A | GSM5354523 | IDC（候选映射） | Naïve（候选映射） | 5789 | 280 | 115 |
| CID4398 | GSM5354524 | IDC | Treated | 4451 | 178 | 926 |
| CID4461 | GSM5354526 | IDC | Naïve | 631 | 41 | 11 |
| CID4463 | GSM5354527 | IDC | Naïve | 1138 | 25 | 75 |
| CID4471 | GSM5354529 | ILC | Naïve | 8609 | 1292 | 159 |
| CID4530N | GSM5354537 | IDC（候选映射） | Naïve（候选映射） | 4409 | 368 | 108 |
| CID4535 | GSM5354538 | ILC | Naïve | 3961 | 102 | 98 |

IDC为浸润性导管癌，ILC为浸润性小叶癌；Naïve表示补充表标记未治疗。完整CSV保留Treatment status和Details of treatment原文、工作表行号及逐样本GEO链接。

所有11个GEO条目均写Primary Breast Tumor，均为整细胞10x Chromium数据。GEO统一描述3′/5′试剂，不能据此确定每个样本的具体化学版本或富集方式，相关字段保留待核。CID4040和CID4398缺少作者标记的癌上皮细胞，但这本身不能证明做过免疫富集。

## 3. 用途建议

| donor_id | 当前实际用途 | 用途建议 | 依据 |
| --- | --- | --- | --- |
| CID3941 | 未分配 | 低覆盖补充候选 | 未治疗IDC，总细胞631，CAF仅8个 |
| CID3948 | 未分配 | 免疫补充候选 | 未治疗IDC，CD8有498个，CAF仅15个 |
| CID4040 | 未分配 | 免疫/基质补充候选 | 未治疗IDC，CD8丰富；作者细胞标签无正常及癌上皮 |
| CID4067 | 未分配 | 保留测试候选，尚未冻结 | 未治疗IDC，CAF/CD8及癌上皮均有覆盖；IDC与当前ILC不同 |
| CID4290A | 未分配 | 患者映射核对后再分配 | 临床4290为候选对应；空间队列有CID4290，后缀A关系需记录 |
| CID4398 | 未分配 | 治疗相关敏感性候选 | 临床表为Treated，接受新辅助FEC-D；作者细胞标签无癌上皮 |
| CID4461 | 未分配 | 低覆盖补充候选 | 未治疗IDC，总细胞631，CD8仅11个 |
| CID4463 | 未分配 | 补充参考/测试候选 | 未治疗IDC，CAF25、CD8为75，多类较少 |
| CID4471 | 未分配 | 优先QC；同ILC的参考候选 | CAF和PVL丰富，与CID4535同为未治疗ILC；癌上皮仅212个 |
| CID4530N | 未分配 | 患者映射核对后再分配 | 临床4530为候选对应；不根据N后缀推断正常组织 |
| CID4535 | 已用于开发 | 已用于开发；不能作为未见测试供者 | 基线拟合与重抽样已使用；未治疗ILC |

## 4. 需要保留的区别

- CID4471和CID4535都是未治疗ILC，可优先检查同组织型参考覆盖。CID4471尚未进行本项目逐供者QC。
- CID4067是未治疗IDC，可保留为测试候选；与ILC基线存在组织学差异，不能把这种差异藏在同ER+标签下。
- CID4398在临床表中为Treated，治疗为Neoadjuvant FEC-D。先列为治疗敏感性候选，不直接混入未治疗起步方案。
- CID4290A与临床4290、CID4530N与临床4530仅做显式候选映射。后缀含义尚未获得直接说明，字段标记为条件性。本地空间队列存在CID4290；在确认是否同患者前，不能把两者分在参考和独立测试两侧。
- 发现CID3963在细胞metadata中为TNBC，GEO和临床补充表为ER+；临床表另记录既往治疗和可能复发。冲突单独保存，不自动加入这11个ER+样本，也不修改作者细胞标签。

## 5. 下一步操作入口

优先核对CID4471、CID4067的逐供者QC与类型覆盖，再确认参考/保留测试角色。临床同号后缀、样本富集信息及亚型冲突保留为明确待核项。没有运行新的表达分析或自动分配训练测试集。

## 6. 来源与复现

- [原论文患者补充表](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-021-00911-1/MediaObjects/41588_2021_911_MOESM4_ESM.xlsx)，Supplementary Table 1第4行为表头，逐样本Excel行号在CSV中。
- [GEO family SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/soft/GSE176078_family.soft.gz)，用Sample_title与orig.ident精确连接。
- 本地数量表：results/Wu2021_donor_reference_inventory/donor_by_reference_celltype.csv。
- 本地原始来源副本和SHA256记录保留于data/provenance及本目录provenance.json。
- 生成脚本：src/sc_sncell_github/build_er_donor_plan.py（需要pandas和openpyxl）。
