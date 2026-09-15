"""汇总六张切片的已完成拟合、表达核对和审核记录。"""
from pathlib import Path
import json
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/Wu2021_multislice_CAF_T_v1'
inventory=pd.read_csv(OUT/'slice_inventory.csv').set_index('sample_id')
summary=pd.read_csv(OUT/'completed_slice_summary.csv').set_index('sample_id')
markers=pd.read_csv(OUT/'all_marker_concordance.csv')
comparison=pd.read_csv(OUT/'CID4535/reference_comparison_metrics.csv').set_index('cell_type')
concordance=pd.read_csv(OUT/'CID4535/reference_expression_concordance.csv').set_index('cell_type')
assert set(summary.index)==set(inventory.index)
validation=json.loads((OUT/'validation_all.json').read_text())
previous=json.loads((OUT/'previous_result_review/review_status.json').read_text())
assert validation['n_passed']==validation['n_checks'] and all(r['exit_code']==0 for r in previous)
lines=['# 双供者参考与多切片CAF/T分析','',
    '## 1. 本轮完成的工作','',
    '固定CID4471和CID4535作为参考，保留原13类拟合，再汇总CD4、CD8、NKT和Cycling T为T细胞大类；NK单独保留。没有修改原作者细胞标签，也没有重新使用CD8特异权重进行空间解释。',
    '先在CID4535上核对与旧模型完全相同的1102个spot和15770个输入基因，完成双供者拟合、新旧分布对照与表达核验，再扩展到其余5张本地Wu切片。',
    f'六张切片共有{int(inventory.n_before.sum())}个原始spot，按作者Artefact与组织外标记的并集排除{int(inventory.n_excluded.sum())}个，保留{int(inventory.n_retained.sum())}个。RCTD最终拟合{int(summary.n_fitted.sum())}个，模型另外排除{int(summary.n_model_removed.sum())}个。',
    '', '## 2. CID4535的新旧参考对照','',
    '| 类型 | 空间排序相关 | 平均绝对变化（百分点） | 最高10%重叠 | 旧/新表达相关 |',
    '|---|---:|---:|---:|---:|']
for t in ['CAFs','T cells']:
    r=comparison.loc[t];e=concordance.loc[t]
    lines.append(f'| {t} | {r.spearman:.3f} | {r.mean_abs_change_pp:.2f} | {int(r.top_overlap)}/{int(r.top_k)}（{r.top_retention:.1%}） | {e.old_spearman:.3f}/{e.new_spearman:.3f} |')
lines+=['',
    '新旧参考保留了主要空间格局，但权重大小及部分高权重spot发生变化。双供者参考没有提高这两类的表达面板相关性；本轮作用是统一参考、记录敏感区域，不把换参考本身称为准确性提升。',
    'reference_high_spot_agreement.png展示两套参考共同较高及仅一套较高的区域；reference_comparison_by_spot.csv保留全部变化。CID4535本身也参与参考构建，这不是未见供者的空间验证。',
    '逐spot对照表还关联了既有单供者参考的10次重抽样高权重频率和标准差，供共同阅读。本轮没有重跑双供者参考的bootstrap，不把旧参考重抽样当作新参考的不确定性估计。',
    '', '## 3. 六张切片的表达核对','',
    'T面板为CD3D、CD3E、CD3G、TRAC、TRBC1、TRBC2；CAF面板为COL1A1、COL1A2、DCN、LUM；NK相关面板为NKG7、GNLY、KLRD1。分数为按完整空间文库归一化后的平均log1p(CP10k)。',
    '', '| 切片 | 作者亚型 | QC保留/拟合 | T表达相关 | CAF表达相关 | 高T但低T表达 |',
    '|---|---|---:|---:|---:|---:|']
for sample in inventory.index:
    s=summary.loc[sample];m=markers[markers.sample_id.eq(sample)]
    t=m[m.cell_type.eq('T cells')&m.marker_panel.eq('T')].iloc[0].spearman
    c=m[m.cell_type.eq('CAFs')&m.marker_panel.eq('CAF')].iloc[0].spearman
    lines.append(f'| {sample} | {inventory.loc[sample,"subtype"]} | {int(s.n_qc_retained)}/{int(s.n_fitted)} | {t:.3f} | {c:.3f} | {int(s.T_high_low_marker)}/{int(s.T_high_count)} |')
lines+=['',
    '“高T但低T表达”定义为预测T权重最高10%，同时T表达分数不高于该片25%分位阈值；零值并列时低表达集合可超过25%。该标记用于逐点复核，未据此删除spot或修改拟合。每张切片的spot_review.csv保存所有标记和原病理标签。',
    'CID4290需要重点复核：其242个T高权重spot中166个未检出这6个T标志，平均预测T权重4.73%，主要来自Cycling T和CD4分量。该片CAF与表达面板一致性较好，但不能把整张T预测图直接作为已确认的浸润分布。',
    '表达核对使用同一RNA数据，属于表达一致性，不是独立蛋白或细胞身份验证；NK相关表达也可来自细胞毒性T细胞。没有用单一基因或权重负相关判断CAF/T空间排斥。',
    '', '## 4. Cycling T敏感性诊断','',
    '看到分歧后，另行计算CD4/CD8/NKT汇总、不包含Cycling T的T表达相关。此为事后诊断，主结果仍使用原先固定的T细胞定义。',
    '', '| 切片 | 主定义T相关 | 不含Cycling T的相关 |',
    '|---|---:|---:|']
for sample in inventory.index:
    s=pd.read_csv(OUT/sample/'T_cycling_sensitivity_summary.csv').set_index('definition')
    lines.append(f'| {sample} | {s.loc["main_T_all","spearman"]:.3f} | {s.loc["diagnostic_without_cycling","spearman"]:.3f} |')
lines+=['',
    '去掉Cycling T后，CID4465、CID44971、1160920F的相关提高，但CID4535、CID4290和1142243F降低。因此不能全局删除Cycling T并称为修复；它在部分切片中伴随T表达，在另一些切片中可能混入其他来源的信号。',
    '', '## 5. 切片间如何比较','',
    'CID4535与CID4290为ER；CID4465、CID44971、1142243F和1160920F为TNBC。后四张使用ER+参考进行固定流程迁移检查，单列报告，不与ER合并成一个亚型差异检验。',
    '前四张使用15770个模板基因。1142243F和1160920F有189个模板基因未按原名匹配，使用15581个交集基因；没有猜测别名或补零。omitted_template_genes.csv保存明细。',
    'all_slices_T_atlas.png使用相同T权重色阶展示六张切片；每张独立的CAF_T_NK_expression图使用片内99%分位显示上限，CSV保留全部未截断数值。',
    '病理区域汇总使用作者标签，Missing和Uncertain不参与病理汇总，但仍可保留在表达与拟合中。没有进行跨切片配准，也没有将不同切片的spot当成独立患者。',
    '', '## 6. 核验记录','',
    f'本轮Python核验{validation["n_passed"]}/{validation["n_checks"]}通过，覆盖输入哈希、参考计数、每片矩阵与坐标、QC并集、模型排除记录、权重汇总、表达分数和相关计算。RDS、CSV与参数另外逐片核对。',
    f'此前{len(previous)}个计算核验脚本已重新运行，全部退出码为0。日志位于previous_result_review；计算通过不抹去此前已经发现的CD4/CD8混淆和低信息供者问题。',
    '原始数据和既有结果保持原路径，本轮输入、模型对象和导出结果分别位于data/processed/Wu2021_multislice_CAF_T_v1与results/Wu2021_multislice_CAF_T_v1。']
(OUT/'interpretation.md').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines))
