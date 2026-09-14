# CID4535 参考重抽样：方法与运行记录

## 1. 本轮要检查什么

前面已经得到一次空间预测，这一步检查同一个供者的参考细胞重抽样后，CAF、CD8 的相对权重与高值区域是否保持一致。癌上皮和 PVL 作为辅助观察对象，所有13类都计算指标。

## 2. 怎么抽样

在作者参考标签内部有放回抽取原来的细胞数，共10次，种子1001—1010。每次预处理和拟合前设种子9；使用full模式、2核心、ref_UMI_min=100、ref_n_cells_min=15，其余设置沿用spacexr 1.2.0默认值。

重复抽到的细胞保存原条码，并使用新的唯一条码构建对象。每类抽取数量不变，但不同原始细胞的数量会减少，尤其小群体需要查看 unique_cells_by_type.csv，不能把重复细胞当成新增样本。

空间spot、参考标签和输入共同基因保持不变。每轮重新调用createRctd，因此模型内部筛选基因可能变化，保存基因列表与相对基线的Jaccard。这里评价整条参考构建到空间映射流程对重抽样的敏感性。

## 3. 文件和使用顺序

- `src/sc_sncell_github/rctd_reference_bootstrap.R`：运行10次重抽样和拟合，逐次保存抽样条码、参数、版本、原生RDS、权重、日志和运行时间。
- `src/sc_sncell_github/inspect_reference_bootstrap.ipynb`：使用Python读取完整结果，比较权重、排序、高值spot集合、病理区域分布，并生成空间稳定性图和数值摘要。
- 原始基线保持在 `data/processed/CID4535_rctd_v1`，新拟合输出位于 `data/processed/CID4535_reference_bootstrap_v1`，解读结果写入 `results/CID4535_reference_bootstrap_v1`。

从项目根目录运行：

```bash
Rscript src/sc_sncell_github/rctd_reference_bootstrap.R
```

拟合完成后，选择项目Python内核，从头运行解读notebook。完整运行有completed.rds标记；中断后重跑脚本会跳过完整运行，重新计算未完成运行。

## 4. 怎么判断

同时看Spearman排序、MAE百分点、偏移量、前5%/10%/20% spot集合的重合程度。逐spot保存重复均值、标准差、观察到的最小最大值和前10%出现频率。空间图中基线和重复均值共用色阶，标准差和高值出现频率单独显示。

10次重复是初步敏感性检查，最小最大值不是置信区间，高值出现频率不是细胞存在概率。结果稳定只能支持同一供者参考内部的稳定性，不能证明真实准确率。全切片CAF—CD8相关是描述性量，不能替代独立病理边界分析。

## 5. 本次实际运行状态

2026-09-15北京时间，10/10次拟合全部完成；每次均保留1102个spot和13种细胞类型。10次源文件MD5一致。每轮实际拟合基因数为1563—1885，基线为1354，基因集合Jaccard为0.566—0.702。

最初两次尝试在fitBulk读取Bioconductor缓存锁时被权限限制中断。用户调整执行权限并要求重试后，重新运行成功，没有将失败的尝试计入10次重复。

解读notebook已从头执行并保存输出；两张空间稳定性图和指标汇总图已逐张检查。生成的指标、区域表、逐spot波动表和PNG/PDF位于results/CID4535_reference_bootstrap_v1。

完整判读见 [CID4535_reference_bootstrap_interpretation.md](CID4535_reference_bootstrap_interpretation.md)。这些结果说明参考内部稳定性，不能替代跨供者准确性验证。
