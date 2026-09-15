# CID44971 CAF/T空间分析

作者亚型：TNBC。质控保留1159个spot，RCTD完成1159个，模型额外排除0个。
参考为CID4471和CID4535，先按13类拟合，再将CD4、CD8、NKT、Cycling T相加为T cells；NK单列。
本片使用15770个公共基因，模板中未按原名匹配的基因有0个。

表达一致性：
- T cells与对应表达面板的Spearman相关为0.562。
- CAFs与对应表达面板的Spearman相关为0.892。
- NK cells与对应表达面板的Spearman相关为0.160。

T细胞最高10%的116个spot中，1个T表达分数不高于本片25%分位阈值，已在spot_review.csv标记；零值并列时低表达集合可超过25%。
表达面板与反卷积使用同一份RNA数据，这是表达一致性核对，不是独立细胞身份验证。NK相关面板也可由细胞毒性T细胞表达。
病理汇总排除Missing和Uncertain，但这些spot仍保留在表达与拟合中。空间图和逐spot表不用于CAF/T因果或空间排斥判断。
本片为TNBC，参考供者为ER+；本次按固定方案进行跨亚型应用检查，不作为同亚型独立验证。
上述分歧spot中，1个T面板UMI总数为零；平均预测T权重为30.24%。T子类平均贡献为：T cells CD4+ 0.01%，T cells CD8+ 2.13%，NKT cells 0.07%，Cycling T-cells 28.03%。
T_expression_discordant_spots.csv保留这些点的坐标、原病理标签、完整文库UMI和T子类分量，供进一步复核；未删除或重新标注。
新增敏感性诊断：仅汇总CD4、CD8、NKT而不纳入Cycling T时，与T面板的相关为0.662。这不是预先固定的主定义，不据此替换主结果，也不能把所有Cycling T判定为错误。
