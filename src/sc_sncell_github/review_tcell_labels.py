"""检查作者T细胞标签、表达信息、参考状态覆盖和已拟合特征；不改标签。"""
from pathlib import Path
import json
import hashlib
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/Wu2021_tcell_label_review_v1'
OUT.mkdir(parents=True,exist_ok=True)
DONORS=['CID4471','CID4535','CID4067','CID4463','CID3941','CID3948','CID4461']
TYPES=['T cells CD4+','T cells CD8+','NK cells','NKT cells','Cycling T-cells']
MARKERS=['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2','CD4','CD8A','CD8B','IL7R','CCR7','TCF7','FOXP3','IL2RA','CTLA4','NKG7','GNLY','PRF1','GZMB','GZMK','KLRD1','IFNG','ZFP36','IFIT1','ISG15','MKI67','TOP2A','EPCAM','KRT19','LYZ','MS4A1']
MARKERS += ['PTPRC','CD2','CD247','LCK','KRT8','KRT18','CD79A','LST1']
def source(donor):
    if donor in ['CID4471','CID4535','CID4067']:folder='Wu2021_ER_donor_qc_v1'
    elif donor=='CID4463':folder='Wu2021_independent_validation_v1'
    else:folder='Wu2021_tcell_validation_v1'
    return ROOT/f'data/processed/{folder}/{donor}_counts_before_qc.h5ad'

cells=[];summaries=[];means={};detection={};source_hashes={}
for donor in DONORS:
    path=source(donor)
    a=ad.read_h5ad(path)
    labels=a.obs.celltype_major.astype(str).where(~a.obs.celltype_major.eq('T-cells'),a.obs.celltype_minor.astype(str))
    keep=labels.isin(TYPES);a=a[keep].copy();labels=labels[keep]
    assert np.array_equal(np.asarray(a.X.sum(axis=1)).ravel(),a.obs.nCount_RNA)
    assert np.array_equal(np.asarray((a.X>0).sum(axis=1)).ravel(),a.obs.nFeature_RNA)
    genes=a.var_names
    raw=a[:,MARKERS].X.toarray()
    table=a.obs[['orig.ident','celltype_major','celltype_minor','celltype_subset','nCount_RNA','nFeature_RNA','percent.mito']].copy()
    table['reference_label']=labels;table['role']='reference' if donor in DONORS[:2] else 'previously_inspected_validation'
    for j,g in enumerate(MARKERS):table[g]=raw[:,j]
    table['T_marker_umi']=table[['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2']].sum(axis=1)
    table['CD8_marker_umi']=table[['CD8A','CD8B']].sum(axis=1)
    table['T_markers_detected']=table.T_marker_umi.gt(0)
    table['CD8_detected']=table.CD8_marker_umi.gt(0)
    table['CD4_detected']=table.CD4.gt(0)
    table['no_T_or_CD8_detected']=~table.T_markers_detected & ~table.CD8_detected
    table['epithelial_panel_umi']=table[['EPCAM','KRT8','KRT18','KRT19']].sum(axis=1)
    table['epithelial_panel_n_genes']=table[['EPCAM','KRT8','KRT18','KRT19']].gt(0).sum(axis=1)
    cells.append(table)
    normalized=a.X.astype(float).multiply((10000/a.obs.nCount_RNA.to_numpy())[:,None]).tocsr()
    for t in TYPES:
        mask=labels.eq(t).to_numpy()
        if not mask.any():continue
        group=table.loc[mask]
        means[(donor,t)]=np.asarray(normalized[mask].mean(axis=0)).ravel()
        detection[(donor,t)]=np.asarray((a.X[mask]>0).mean(axis=0)).ravel()
        summaries.append(dict(donor_id=donor,reference_label=t,n_cells=int(mask.sum()),umi_median=group.nCount_RNA.median(),genes_median=group.nFeature_RNA.median(),mito_median=group['percent.mito'].median(),T_marker_positive=int(group.T_markers_detected.sum()),CD8_positive=int(group.CD8_detected.sum()),CD4_positive=int(group.CD4_detected.sum()),no_T_or_CD8=int(group.no_T_or_CD8_detected.sum())))
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    source_hashes[str(path.relative_to(ROOT))]=h.hexdigest()
cells=pd.concat(cells);cells.index.name='barcode'
assert cells.index.is_unique
cells.to_csv(OUT/'tcell_label_expression_review.csv')
summary=pd.DataFrame(summaries);summary.to_csv(OUT/'label_marker_summary.csv',index=False)
pd.crosstab([cells['orig.ident'],cells.reference_label],cells.celltype_subset).to_csv(OUT/'author_state_counts.csv')
focus=cells[cells['orig.ident'].eq('CID4461')&cells.reference_label.eq('T cells CD8+')]
focus.to_csv(OUT/'CID4461_11_CD8_review.csv')
cells[cells['orig.ident'].eq('CID4535')&cells.reference_label.eq('T cells CD4+')].to_csv(OUT/'CID4535_CD4_review.csv')
marker_rows=[]
for (donor,t),profile in means.items():
    for g in MARKERS:
        i=genes.get_loc(g)
        marker_rows.append(dict(donor_id=donor,reference_label=t,gene=g,mean_cp10k=profile[i],detected_fraction=detection[(donor,t)][i]))
marker_summary=pd.DataFrame(marker_rows);marker_summary.to_csv(OUT/'marker_expression_summary.csv',index=False)

# 只用两位参考供者检查CD8相对CD4的表达方向，再查看其他供者。
contrasts=pd.DataFrame(index=genes)
for donor in DONORS:
    contrasts[donor]=np.log2((means[(donor,'T cells CD8+')]+.1)/(means[(donor,'T cells CD4+')]+.1))
weights=summary.set_index(['donor_id','reference_label']).n_cells
pooled={t:sum(means[(d,t)]*weights.loc[(d,t)] for d in DONORS[:2])/sum(weights.loc[(d,t)] for d in DONORS[:2]) for t in TYPES}
contrasts['pooled_reference_log2ratio']=np.log2((pooled['T cells CD8+']+.1)/(pooled['T cells CD4+']+.1))
contrasts['pooled_reference_mean_cp10k']=(pooled['T cells CD8+']+pooled['T cells CD4+'])/2
contrasts['reference_direction_agrees']=(contrasts.CID4471*contrasts.CID4535)>0
contrasts['both_reference_abs_log2ratio_ge_0_5']=(contrasts.CID4471.abs()>=.5)&(contrasts.CID4535.abs()>=.5)
contrasts.to_csv(OUT/'cd8_vs_cd4_expression_contrasts.csv',index_label='gene')
eligible=contrasts[contrasts.pooled_reference_mean_cp10k>=1].copy()
eligible['magnitude']=eligible.pooled_reference_log2ratio.abs()
top=eligible.sort_values(['magnitude'],ascending=False,kind='stable').head(100)
top.to_csv(OUT/'top100_reference_contrasts.csv',index_label='gene')

# 将零CD8配方的预测关联回来源细胞；这不是单细胞预测。
spot_records=[];exposures=[]
for donor in ['CID3941','CID3948','CID4461']:
    base=ROOT/'results/Wu2021_tcell_validation_v1'/donor
    design=pd.read_csv(base/'simulation_design.csv',index_col=0)
    members=pd.read_csv(base/'mixture_membership.csv')
    selected=design.index[design.scenario.eq('gradient')&design.n_cd8.eq(0)&design.n_caf.eq(0)&design.depth_fraction.eq(1)]
    for variant in ['baseline','donor_balanced','epithelial_split','balanced_split']:
        pred=pd.read_csv(base/variant/'rctd_relative_weights.csv',index_col=0)
        for spot in selected:
            group=members[members.spot_id.eq(spot)]
            barcodes=group.source_barcode.tolist();source_cells=cells.loc[barcodes]
            states=source_cells.celltype_subset.value_counts().to_dict()
            spot_records.append(dict(donor_id=donor,variant=variant,spot_id=spot,predicted_CD8_pct=pred.loc[spot,'T cells CD8+']*100,source_CD4_detected=int(source_cells.CD4_detected.sum()),source_CD8_detected=int(source_cells.CD8_detected.sum()),source_T_detected=int(source_cells.T_markers_detected.sum()),source_states=json.dumps(states,ensure_ascii=False)))
            if variant=='baseline':
                exposures.extend(dict(barcode=b,spot_id=spot,mixture_CD8_prediction_pct=pred.loc[spot,'T cells CD8+']*100) for b in barcodes)
pd.DataFrame(spot_records).to_csv(OUT/'zero_CD8_mixture_source_review.csv',index=False)
exposures=pd.DataFrame(exposures)
exposure_summary=exposures.groupby('barcode').agg(n_mixtures=('spot_id','nunique'),mean_mixture_CD8_prediction_pct=('mixture_CD8_prediction_pct','mean'))
cells.join(exposure_summary,how='inner').to_csv(OUT/'CD4_cells_in_zero_CD8_mixtures.csv')

# 检查实际输入的参考表达和基因筛选是否保留了所查看的标志。
feature_rows=[];cosines=[]
for variant in ['baseline','donor_balanced','epithelial_split','balanced_split']:
    profile=pd.read_csv(OUT/f'{variant}_rctd_profiles.csv',index_col=0)
    for donor in ['CID3941','CID3948','CID4461']:
        used=set((ROOT/'results/Wu2021_tcell_validation_v1'/donor/variant/'rctd_selected_genes.txt').read_text().splitlines())
        ref=ROOT/'data/processed/Wu2021_known_composition_v1' if variant=='baseline' else ROOT/'data/processed/Wu2021_reference_comparison_v1'/variant
        nnls_genes=set(pd.read_csv(ref/'nnls_reference_profiles.csv',index_col=0).index)
        for g in MARKERS:
            feature_rows.append(dict(variant=variant,validation_donor=donor,gene=g,in_RCTD= g in used,in_NNLS=g in nnls_genes,reference_CD4_profile=profile.loc[g,'T cells CD4+'],reference_CD8_profile=profile.loc[g,'T cells CD8+']))
        shared=[g for g in profile.index if g in used]
        for other in ['T cells CD4+','NK cells','NKT cells','Cycling T-cells']:
            x=profile.loc[shared,'T cells CD8+'].to_numpy();y=profile.loc[shared,other].to_numpy()
            cosines.append(dict(variant=variant,validation_donor=donor,comparison=other,n_genes=len(shared),cosine=float(x@y/(np.linalg.norm(x)*np.linalg.norm(y)))))
pd.DataFrame(feature_rows).to_csv(OUT/'marker_feature_coverage.csv',index=False)
pd.DataFrame(cosines).to_csv(OUT/'reference_profile_cosines.csv',index=False)
reference_qc=[]
for variant in ['baseline','donor_balanced','epithelial_split','balanced_split']:
    ref=ROOT/'data/processed/Wu2021_known_composition_v1' if variant=='baseline' else ROOT/'data/processed/Wu2021_reference_comparison_v1'/variant
    meta=pd.read_csv(ref/'reference_metadata.csv',index_col=0)
    selected=cells.loc[cells.index.intersection(meta.index)]
    for t,group in selected.groupby('reference_label'):
        reference_qc.append(dict(variant=variant,reference_label=t,n_cells=len(group),CID4535_cells=int(group['orig.ident'].eq('CID4535').sum()),no_T_or_CD8=int(group.no_T_or_CD8_detected.sum()),no_T_or_CD8_fraction=group.no_T_or_CD8_detected.mean()))
pd.DataFrame(reference_qc).to_csv(OUT/'reference_marker_coverage_by_variant.csv',index=False)
(OUT/'provenance.json').write_text(json.dumps(dict(source_hashes=source_hashes,markers=MARKERS,normalization='Counts per cell / full-library UMI * 10000; detection means raw count >=1.',contrast='Descriptive log2 ratio of mean CP10k with pseudocount 0.1; not a significance test.',top100='Reference-only mean CP10k >=1, ranked by absolute pooled CD8/CD4 log2 ratio.',labels_changed=False,cells_removed=False),ensure_ascii=False,indent=2)+'\n')

fig,ax=plt.subplots(figsize=(12,6),layout='constrained')
display_genes=MARKERS[:26]
im=ax.imshow(np.log1p(focus[display_genes].to_numpy()),aspect='auto',cmap='YlGnBu')
ax.set_xticks(range(len(display_genes)),display_genes,rotation=65,ha='right',fontsize=8)
ax.set_yticks(range(len(focus)),[b.split('_',1)[1] for b in focus.index],fontsize=8)
ax.set_title('CID4461 | 11 author-labeled CD8 cells');fig.colorbar(im,ax=ax,label='log(1 + raw UMI)',shrink=.8)
fig.savefig(OUT/'CID4461_marker_counts.png',dpi=200);fig.savefig(OUT/'CID4461_marker_counts.pdf');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(13,6),layout='constrained')
show=['CD3D','CD3E','TRAC','CD4','CD8A','CD8B','IL7R','FOXP3','NKG7','GNLY','GZMK','GZMB','IFIT1','MKI67']
for ax,t in zip(axes,['T cells CD4+','T cells CD8+']):
    table=marker_summary[marker_summary.reference_label.eq(t)].pivot(index='donor_id',columns='gene',values='detected_fraction').loc[DONORS,show]
    im=ax.imshow(table,vmin=0,vmax=1,cmap='YlGnBu',aspect='auto')
    ax.set_xticks(range(len(show)),show,rotation=65,ha='right',fontsize=8);ax.set_yticks(range(7),DONORS)
    ax.set_title(t)
fig.colorbar(im,ax=axes,label='Fraction of cells with raw UMI >= 1',shrink=.8)
fig.savefig(OUT/'marker_detection_by_donor.png',dpi=200);fig.savefig(OUT/'marker_detection_by_donor.pdf');plt.close(fig)
print(summary[summary.reference_label.isin(['T cells CD4+','T cells CD8+'])].to_string(index=False))
print('Top100 direction agreement:',int(top.reference_direction_agrees.sum()))
print(pd.crosstab([cells['orig.ident'],cells.reference_label],cells.celltype_subset).loc[(slice(None),['T cells CD4+','T cells CD8+']),:].to_string())

lines=['# T细胞标签与参考表达核对','',
    '## 1. 本次检查','',
    f'检查7位供者的{len(cells)}个作者标注CD4、CD8、NK、NKT及Cycling T细胞。读取原始UMI，逐细胞保存39个查看基因；没有删除细胞、修改标签或重新拟合。',
    'T标志检测指CD3D、CD3E、CD3G、TRAC、TRBC1、TRBC2中任一个原始计数至少1；CD8检测指CD8A或CD8B至少1。未检出只描述这份RNA数据，不等于证明细胞不属于该类型。',
    '', '## 2. CID4461的11个CD8','',
    'UMI中位数400，检测基因中位数287。9/11检出至少一个本次定义的T标志，只有2/11检出CD8A或CD8B，0/11检出CD4。',
    '作者状态为5个CD8+ ZFP36、5个CD8+ IFNG及1个IFIT1。参考侧也有这些状态，不能说它们在参考中完全缺失；但这11个细胞提供的CD8直接表达信息很少。',
    '其中部分细胞检出GNLY或KLRD1，不能仅凭这些共享表达把作者CD8标签改为NK或NKT。应保留原标签，同时记录为低信息验证组。',
    '', '## 3. CID4535的CD4参考需要重点复核','',
    'CID4471的412个CD4中，408个检出T标志；CID4535的239个CD4中只有105个。CID4535有132个同时未检出上述T标志及CD8标志，CID4471对应为4个。',
    '问题集中在CID4535的CD4+ IL7R状态：166个细胞，UMI中位数377，只有34个检出T标志。其余CD4状态UMI中位数为2594至3438。',
    'CID4535的CD4中，145个检出KRT19，63个在EPCAM、KRT8、KRT18、KRT19中检出至少两个；33个同时没有检出上述T/CD8标志。这需要区分环境RNA、低信息细胞和标签混杂，当前证据不足以直接改成上皮标签。',
    '供者均衡把CID4535在CD4参考中的占比从239/651（36.7%）提高到239/478（50.0%）；无T/CD8标志的细胞比例从136/651（20.9%）提高到135/478（28.2%）。均衡确实增加了这部分低信息参考的占比，是否导致预测恶化还需单独对照。',
    '', '## 4. 实际模型特征','',
    '本轮12次RCTD拟合均保留CD8A、CD8B，但没有保留CD4和FOXP3。原参考CD4基因在全部类型中的最高平均相对表达约0.000163，低于当前回归基因阈值0.0002；FOXP3最高约0.0000666，也低于阈值。安装版本getDeGenes要求类型表达超过此阈值，因此这两者未能进入回归特征。',
    '这符合当前筛选规则，不是代码报错。降低阈值或定向保留标志可以作为后续对照，但本次未修改参数，也没有证明只加入CD4就能解决混淆。NNLS保留了CD4，仍然存在误差。',
    f'按两位参考的合并CD8/CD4表达差异挑选前100个基因，其中{int(top.reference_direction_agrees.sum())}个在两位供者的方向一致。该统计为描述性比值，不是显著性检验。部分偏向CD4的差异主要由CID4535贡献，例如KRT19、AGR2、SLC39A6，不能直接当作稳定的CD4区分特征。',
    '实际筛选基因上的CD4/CD8参考余弦相似度较高，但相似度受共享表达影响，不能单独证明类型不可区分。完整数值保存在reference_profile_cosines.csv。',
    '', '## 5. 零CD8配方里的来源CD4','',
    'CID3941误分配最高的原参考配方预测CD8为31.07%，来源为8个作者标注CD4+ IL7R及2个CD4+ CCR7。CID3948最高配方为52.11%，来源为6个CD4+ IL7R、2个Treg和2个CD4+ CCR7。两组都没有混入作者标注CD8。',
    '这些配方分别只有1个来源细胞检出CD8A或CD8B。来源组成支持检查CD4状态与参考表达之间的差异，但spot级误差不能归责为某一个来源细胞预测错误。',
    'CD4_cells_in_zero_CD8_mixtures.csv列出全部来源细胞与其出现配方的平均预测；该列是关联配方的预测，不是单细胞预测。',
    '', '## 6. 当前修订方向','',
    '优先复核CID4535的CD4+ IL7R参考，设计保留全量与仅使用表达支持较充分细胞的对照，避免直接按预测好坏删细胞。CID4461的11个CD8继续保留，单列为低信息组。',
    '另做基因特征对照，检查CD4等标志保留前后的影响，并同时查看CD4→CD8误分配和真实CD8的检出。两项调整应分别进行，才能判断各自作用。',
    '原标签来源为本地Wu作者metadata；原研究：[Wu et al., 2021](https://pubmed.ncbi.nlm.nih.gov/34493872/)。本次基因查看和表达对照不能代替蛋白或TCR证据。']
(OUT/'interpretation.md').write_text('\n'.join(lines)+'\n')
