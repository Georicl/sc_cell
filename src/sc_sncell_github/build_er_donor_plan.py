"""建立ER+候选供者的来源核对与用途草案，不执行QC或划分训练集。"""

# 1. 读取之前的数量表，以及本次核对的GEO和临床补充表
from pathlib import Path
from datetime import datetime, timezone
import gzip
import hashlib
import json

import openpyxl
import pandas as pd

PROJECT = Path(__file__).resolve().parents[2]
OUT = PROJECT / "results/Wu2021_ER_donor_plan"
OUT.mkdir(parents=True, exist_ok=True)
COUNTS = PROJECT / "results/Wu2021_donor_reference_inventory/donor_by_reference_celltype.csv"
SOFT = PROJECT / "data/provenance/wu2021_scrna.family.soft.gz"
SUPP = PROJECT / "data/provenance/wu2021_supplementary_tables.xlsx"
SUPP_URL = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-021-00911-1/MediaObjects/41588_2021_911_MOESM4_ESM.xlsx"
SOFT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/soft/GSE176078_family.soft.gz"

counts = pd.read_csv(COUNTS).set_index("donor_id")
text = gzip.open(SOFT, "rt").read()
geo_rows = []
for block in text.split("^SAMPLE = ")[1:]:
    values = {}
    for line in block.splitlines():
        if line.startswith("!Sample_") and " = " in line:
            key, value = line.split(" = ", 1)
            values.setdefault(key, []).append(value)
    characteristics = dict(
        value.split(": ", 1)
        for value in values.get("!Sample_characteristics_ch1", [])
        if ": " in value
    )
    geo_rows.append({
        "donor_id": values["!Sample_title"][0],
        "gsm": block.splitlines()[0].strip(),
        "geo_subtype": characteristics.get("clinical_subtype", "未提供"),
        "geo_tissue": characteristics.get("tissue", "未提供"),
        "geo_source": values.get("!Sample_source_name_ch1", ["未提供"])[0],
        "geo_protocol": " | ".join(values.get("!Sample_extract_protocol_ch1", [])),
        "geo_treatment": " | ".join(values.get("!Sample_treatment_protocol_ch1", [])) or "未逐样本提供",
    })
geo = pd.DataFrame(geo_rows).set_index("donor_id")
assert geo.index.is_unique and counts.index.is_unique
assert set(counts.index) == set(geo.index)

workbook = openpyxl.load_workbook(SUPP, read_only=True, data_only=True)
sheet = workbook["Supplementary Table 1"]
rows = list(sheet.iter_rows(values_only=True))
clinical = pd.DataFrame(rows[4:], columns=rows[3]).dropna(subset=["Case ID"])
clinical["excel_row"] = clinical.index + 5
clinical["case_id_text"] = clinical["Case ID"].astype(str)
clinical = clinical.set_index("case_id_text")
assert clinical.index.is_unique
clinical.reset_index().to_csv(OUT / "clinical_table1_extracted.csv", index=False)
geo.to_csv(OUT / "geo_sample_fields.csv")

# 2. 显式记录用途建议，不把候选角色当成已冻结划分
# 对有后缀的样本仅记录候选临床映射；不通用删字母、不合并样本。
alias_candidates = {"CID4290A": "4290", "CID4530N": "4530"}
roles = {
    "CID4535": ("已用于开发", "已用于开发；不能作为未见测试供者", "基线拟合与重抽样已使用；未治疗ILC", "后续逐供者QC；开发用途已确定"),
    "CID4471": ("未分配", "优先QC；同ILC的参考候选", "CAF和PVL丰富，与CID4535同为未治疗ILC；癌上皮仅212个", "核对QC及正常/癌上皮标签覆盖，再决定参考用途"),
    "CID4067": ("未分配", "保留测试候选，尚未冻结", "未治疗IDC，CAF/CD8及癌上皮均有覆盖；IDC与当前ILC不同", "QC后评估测试组成覆盖；不可称为同组织型验证"),
    "CID4398": ("未分配", "治疗相关敏感性候选", "临床表为Treated，接受新辅助FEC-D；作者细胞标签无癌上皮", "不直接混入未治疗起步方案；富集或缺失原因待核"),
    "CID4040": ("未分配", "免疫/基质补充候选", "未治疗IDC，CD8丰富；作者细胞标签无正常及癌上皮", "核对细胞覆盖原因；不能构造该供者的癌上皮混合真值"),
    "CID3941": ("未分配", "低覆盖补充候选", "未治疗IDC，总细胞631，CAF仅8个", "QC后重新查看CAF覆盖，不直接按数量删除供者"),
    "CID3948": ("未分配", "免疫补充候选", "未治疗IDC，CD8有498个，CAF仅15个", "CAF测试覆盖有限，需明确目标类型"),
    "CID4461": ("未分配", "低覆盖补充候选", "未治疗IDC，总细胞631，CD8仅11个", "QC后重新查看CD8覆盖"),
    "CID4463": ("未分配", "补充参考/测试候选", "未治疗IDC，CAF25、CD8为75，多类较少", "QC后检查所有目标类别是否有足够覆盖"),
    "CID4290A": ("未分配", "患者映射核对后再分配", "临床4290为候选对应；空间队列有CID4290，后缀A关系需记录", "核实4290A与4290对应；若同患者，必须放在同一划分侧"),
    "CID4530N": ("未分配", "患者映射核对后再分配", "临床4530为候选对应；不根据N后缀推断正常组织", "核实后缀含义和患者对应，之后再定用途"),
}
er_ids = counts.index[counts["subtype"].eq("ER+")].tolist()
assert set(er_ids) == set(roles)
plan_rows = []
for donor in er_ids:
    raw_case = donor.removeprefix("CID")
    is_alias = donor in alias_candidates
    clinical_id = alias_candidates.get(donor, raw_case)
    clinical_row = clinical.loc[clinical_id]
    actual, proposed, reason, todo = roles[donor]
    entry = counts.loc[donor].to_dict()
    entry.update(geo.loc[donor].to_dict())
    entry.update({
        "donor_id": donor,
        "clinical_case_id": clinical_id,
        "clinical_mapping_status": "候选同号映射；后缀含义待核" if is_alias else "CID前缀标准化后与Case ID精确匹配",
        "patient_group_for_split": "待核实" if is_alias else "CID" + clinical_id,
        "clinical_subtype": clinical_row["Subtype by IHC"],
        "histology": clinical_row["Cancer Type"],
        "grade": clinical_row["Grade"],
        "treatment_status": clinical_row["Treatment status"],
        "treatment_details": clinical_row["Details of treatment"],
        "clinical_evidence_status": "临床字段依赖候选患者映射" if is_alias else "补充表1已核对",
        "sampling_method": "整细胞scRNA-seq；10x Chromium",
        "per_sample_enrichment": "GEO未明确逐样本富集/分选情况，待方法记录核对",
        "chemistry_per_sample": "GEO统一写v2 3′及5′，未明确逐样本归属",
        "current_role": actual,
        "proposed_role": proposed,
        "proposal_reason": reason,
        "inclusion_status": "盘点纳入，未冻结正式分析纳入",
        "qc_status": "单样本基线已做；全队列统一QC未做" if donor == "CID4535" else "未执行本项目逐供者QC",
        "next_check": todo,
        "clinical_sheet": "Supplementary Table 1",
        "clinical_excel_row": int(clinical_row["excel_row"]),
        "geo_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + geo.loc[donor, "gsm"],
        "clinical_source_url": SUPP_URL,
    })
    plan_rows.append(entry)
plan = pd.DataFrame(plan_rows).set_index("donor_id")
assert len(plan) == 11 and plan["subtype"].eq("ER+").all()
assert plan["geo_subtype"].eq("ER+").all() and plan["clinical_subtype"].eq("ER+").all()
assert plan["current_role"].eq("已用于开发").sum() == 1
plan.to_csv(OUT / "er_donor_inclusion_and_role.csv")

# 3. 单独记录亚型冲突，不擅自改写前面的数量表
conflicts = []
for donor in counts.index:
    author_subtype = counts.loc[donor, "subtype"]
    geo_subtype = geo.loc[donor, "geo_subtype"]
    # HER2+/ER+与图谱合并为HER2+不是ER+/TNBC冲突
    geo_broad = "HER2+" if "HER2+" in geo_subtype else geo_subtype
    if author_subtype != geo_broad:
        raw_id = donor.removeprefix("CID")
        c = clinical.loc[raw_id] if raw_id in clinical.index else None
        conflicts.append({
            "donor_id": donor, "metadata_subtype": author_subtype,
            "geo_subtype": geo_subtype,
            "clinical_subtype": c["Subtype by IHC"] if c is not None else "待核",
            "treatment_status": c["Treatment status"] if c is not None else "待核",
            "treatment_details": c["Details of treatment"] if c is not None else "待核",
            "pathology_note": c["Notable Pathological features"] if c is not None else "待核",
            "action": "保留原metadata标签；暂不并入当前ER+方案，需确定分型定义",
            "geo_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + geo.loc[donor, "gsm"],
            "clinical_excel_row": int(c["excel_row"]) if c is not None else None,
            "clinical_source_url": SUPP_URL,
        })
pd.DataFrame(conflicts).to_csv(OUT / "subtype_conflicts.csv", index=False)

# 4. 生成方便阅读的表，明确事实、用途建议和待核项
def md_table(frame):
    frame = frame.reset_index()
    rows = ["| " + " | ".join(map(str, frame.columns)) + " |",
            "| " + " | ".join(["---"] * len(frame.columns)) + " |"]
    rows += ["| " + " | ".join(str(x).replace("|", "/") for x in row) + " |"
             for row in frame.itertuples(index=False, name=None)]
    return "\n".join(rows)

clinical_view = plan[["gsm", "histology", "treatment_status", "total_cells", "CAFs", "T cells CD8+"]].copy()
clinical_view.loc[["CID4290A", "CID4530N"], "treatment_status"] += "（候选映射）"
clinical_view.loc[["CID4290A", "CID4530N"], "histology"] += "（候选映射）"
clinical_view.columns = ["GEO样本", "组织学", "治疗状态", "总细胞", "CAF", "CD8 T"]
role_view = plan[["current_role", "proposed_role", "proposal_reason"]].copy()
role_view.columns = ["当前实际用途", "用途建议", "依据"]
report = (
    "# Wu ER+供者纳入与用途表（草案v1）\n\n"
    "## 1. 本次核对范围\n\n"
    "以作者细胞元数据中11个ER+样本为范围，对照GEO逐样本记录与原论文Supplementary Table 1。"
    "纳入盘点不等于已纳入正式模型。只将CID4535记录为已用于开发，其他供者的用途均为建议，未冻结测试集。\n\n"
    "## 2. 来源、组织学、治疗和目标细胞覆盖\n\n" + md_table(clinical_view) +
    "\n\nIDC为浸润性导管癌，ILC为浸润性小叶癌；Naïve表示补充表标记未治疗。"
    "完整CSV保留Treatment status和Details of treatment原文、工作表行号及逐样本GEO链接。\n\n"
    "所有11个GEO条目均写Primary Breast Tumor，均为整细胞10x Chromium数据。"
    "GEO统一描述3′/5′试剂，不能据此确定每个样本的具体化学版本或富集方式，相关字段保留待核。"
    "CID4040和CID4398缺少作者标记的癌上皮细胞，但这本身不能证明做过免疫富集。\n\n"
    "## 3. 用途建议\n\n" + md_table(role_view) +
    "\n\n## 4. 需要保留的区别\n\n"
    "- CID4471和CID4535都是未治疗ILC，可优先检查同组织型参考覆盖。CID4471尚未进行本项目逐供者QC。\n"
    "- CID4067是未治疗IDC，可保留为测试候选；与ILC基线存在组织学差异，不能把这种差异藏在同ER+标签下。\n"
    "- CID4398在临床表中为Treated，治疗为Neoadjuvant FEC-D。先列为治疗敏感性候选，不直接混入未治疗起步方案。\n"
    "- CID4290A与临床4290、CID4530N与临床4530仅做显式候选映射。后缀含义尚未获得直接说明，字段标记为条件性。"
    "本地空间队列存在CID4290；在确认是否同患者前，不能把两者分在参考和独立测试两侧。\n"
    "- 发现CID3963在细胞metadata中为TNBC，GEO和临床补充表为ER+；临床表另记录既往治疗和可能复发。"
    "冲突单独保存，不自动加入这11个ER+样本，也不修改作者细胞标签。\n\n"
    "## 5. 下一步操作入口\n\n"
    "优先核对CID4471、CID4067的逐供者QC与类型覆盖，再确认参考/保留测试角色。"
    "临床同号后缀、样本富集信息及亚型冲突保留为明确待核项。没有运行新的表达分析或自动分配训练测试集。\n\n"
    "## 6. 来源与复现\n\n"
    f"- [原论文患者补充表]({SUPP_URL})，Supplementary Table 1第4行为表头，逐样本Excel行号在CSV中。\n"
    f"- [GEO family SOFT]({SOFT_URL})，用Sample_title与orig.ident精确连接。\n"
    "- 本地数量表：results/Wu2021_donor_reference_inventory/donor_by_reference_celltype.csv。\n"
    "- 本地原始来源副本和SHA256记录保留于data/provenance及本目录provenance.json。\n"
    "- 生成脚本：src/sc_sncell_github/build_er_donor_plan.py（需要pandas和openpyxl）。\n"
)
(OUT / "er_donor_plan.md").write_text(report, encoding="utf-8")
provenance = {
    "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
    "sources": [{"path": str(path.relative_to(PROJECT)),
                 "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "url": url}
                for path, url in [(SOFT, SOFT_URL), (SUPP, SUPP_URL), (COUNTS, "local inventory")]],
    "candidate_aliases_not_confirmed": alias_candidates,
    "n_er_samples": len(plan), "no_split_frozen": True,
}
(OUT / "provenance.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False), encoding="utf-8")
print(clinical_view.to_string())
print("亚型冲突：", conflicts)
print("已保存：", OUT)
