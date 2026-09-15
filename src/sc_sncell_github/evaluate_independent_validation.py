"""评估冻结参考在CID4463上的预测，保留全部方案与类型。"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import anndata as ad
from scipy.optimize import nnls
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data/processed/Wu2021_independent_validation_v1'
OUT=ROOT/'results/Wu2021_independent_validation_v1'
VARIANTS=['baseline','donor_balanced','epithelial_split','balanced_split']
NAMES=dict(zip(VARIANTS,['原参考','供者均衡','上皮拆分','均衡并拆分']))
ENGLISH=['Baseline','Donor balanced','Epithelial split','Balanced + split']

def main(nnls_only=False):
    design=pd.read_csv(OUT/'simulation_design.csv',index_col=0)
    truths={n:pd.read_csv(OUT/f'truth_{f}_fraction.csv',index_col=0) for n,f in [('rna','rna'),('cells','cell')]}
    types=truths['rna'].columns
    spots=ad.read_h5ad(DATA/'simulated_spots.h5ad')
    assert spots.obs_names.equals(design.index)
    members=pd.read_csv(OUT/'mixture_membership.csv')
    source_sets=members.groupby('spot_id').source_barcode.agg(lambda x:'|'.join(sorted(x)))
    diversity=design[['scenario','recipe','depth_fraction']].copy()
    diversity['source_cell_set']=source_sets
    diversity.groupby(['scenario','recipe','depth_fraction']).agg(n_spots=('source_cell_set','size'),distinct_source_sets=('source_cell_set','nunique')).to_csv(OUT/'source_mixture_diversity.csv')
    predictions={};rows=[];pure_rows=[];error_rows=[];titration_rows=[]
    for variant in VARIANTS:
        folder=OUT/variant;folder.mkdir(exist_ok=True)
        ref=ROOT/'data/processed/Wu2021_known_composition_v1' if variant=='baseline' else ROOT/'data/processed/Wu2021_reference_comparison_v1'/variant
        profiles=pd.read_csv(ref/'nnls_reference_profiles.csv',index_col=0)
        y=spots[:,profiles.index].X.toarray()/design.nUMI.to_numpy()[:,None]
        raw=np.vstack([nnls(profiles.to_numpy(),row)[0] for row in y])
        assert (raw.sum(axis=1)>0).all()
        frame=pd.DataFrame(raw,index=design.index,columns=profiles.columns)
        frame.to_csv(folder/'nnls_raw_coefficients.csv')
        frame.div(frame.sum(axis=1),axis=0).to_csv(folder/'nnls_relative_weights.csv')
        if nnls_only:continue
        for method in ['RCTD','NNLS']:
            pred=pd.read_csv(folder/f'{method.lower()}_relative_weights.csv',index_col=0)
            assert pred.index.is_unique and set(pred.index)==set(design.index)
            assert np.isfinite(pred.values).all() and (pred.values>=0).all() and np.allclose(pred.sum(axis=1),1)
            if variant!='baseline':
                mapping=pd.read_csv(ROOT/'results/Wu2021_reference_comparison_v1'/variant/'label_mapping.csv',index_col=0).parent_label
                assert set(mapping.index)==set(pred.columns)
                pred=pred.T.groupby(mapping).sum().T
            assert set(pred.columns)==set(types)
            pred=pred.loc[design.index,types]
            pred.to_csv(folder/f'{method.lower()}_parent_weights.csv')
            predictions[(variant,method)]=pred
            for truth_name,truth in truths.items():
                for (scenario,depth),group in design.groupby(['scenario','depth_fraction']):
                    ids=group.index
                    for t in list(types)+['Epithelial total']:
                        cols=['Cancer Epithelial','Normal Epithelial'] if t=='Epithelial total' else [t]
                        actual=truth.loc[ids,cols].sum(axis=1)
                        estimated=pred.loc[ids,cols].sum(axis=1)
                        error=(estimated-actual)*100
                        present=truths['cells'].loc[ids,cols].sum(axis=1)>0
                        rows.append(dict(donor_id='CID4463',variant=variant,method=method,truth=truth_name,scenario=scenario,depth_fraction=depth,cell_type=t,n_spots=len(ids),mae_pp=error.abs().mean(),rmse_pp=np.sqrt(np.square(error).mean()),bias_pp=error.mean(),present_mae_pp=error[present].abs().mean(),n_present=int(present.sum()),n_absent=int((~present).sum()),absent_mean_prediction_pct=estimated[~present].mean()*100,absent_ge1pct_fraction=(estimated[~present]>=.01).mean() if (~present).any() else np.nan))
            for t in types:
                errors=design[['donor_id','base_mixture','scenario','depth_fraction']].copy()
                errors['variant']=variant;errors['method']=method;errors['cell_type']=t
                errors['truth_rna']=truths['rna'][t];errors['prediction']=pred[t]
                errors['rna_error_pp']=(pred[t]-truths['rna'][t])*100
                error_rows.append(errors)
            for depth in [1.,.25]:
                ids=design.index[design.scenario.eq('pure')&design.depth_fraction.eq(depth)]
                actual=truths['cells'].loc[ids].idxmax(axis=1)
                for t in sorted(actual.unique()):
                    keep=actual.index[actual.eq(t)]
                    pure_rows.append(dict(donor_id='CID4463',variant=variant,method=method,depth_fraction=depth,cell_type=t,n_spots=len(keep),dominant_matches=int(pred.loc[keep].idxmax(axis=1).eq(t).sum()),mean_weight_on_true_type=pred.loc[keep,t].mean(),mean_normal_weight=pred.loc[keep,'Normal Epithelial'].mean(),mean_cancer_weight=pred.loc[keep,'Cancer Epithelial'].mean()))
            ids=design.index[design.scenario.eq('epithelial_titration')]
            table=design.loc[ids,['depth_fraction','recipe','repeat']].copy()
            table['variant']=variant;table['method']=method
            table['normal_cell_fraction']=truths['cells'].loc[ids,'Normal Epithelial']
            table['normal_rna_fraction']=truths['rna'].loc[ids,'Normal Epithelial']
            table['normal_prediction']=pred.loc[ids,'Normal Epithelial']
            table['cancer_prediction']=pred.loc[ids,'Cancer Epithelial']
            titration_rows.append(table)
    if nnls_only:return
    metrics=pd.DataFrame(rows);pure=pd.DataFrame(pure_rows)
    metrics.to_csv(OUT/'validation_metrics.csv',index=False)
    pure.to_csv(OUT/'pure_control_summary.csv',index=False)
    pd.concat(error_rows).to_csv(OUT/'spot_prediction_errors.csv',index_label='spot_id')
    titration=pd.concat(titration_rows);titration.to_csv(OUT/'epithelial_titration_predictions.csv',index_label='spot_id')
    key=metrics[metrics.truth.eq('rna')&metrics.scenario.eq('factorial')]
    colors=['#636363','#437b9e','#ca8a3d','#5a976c']
    fig,axes=plt.subplots(2,3,figsize=(13,7),layout='constrained')
    for i,method in enumerate(['RCTD','NNLS']):
        for j,t in enumerate(['CAFs','T cells CD8+','Cancer Epithelial']):
            ax=axes[i,j]
            for k,v in enumerate(VARIANTS):
                sub=key[key.method.eq(method)&key.variant.eq(v)&key.cell_type.eq(t)].set_index('depth_fraction')
                ax.bar(np.arange(2)+(k-1.5)*.19,sub.loc[[1.,.25],'mae_pp'],width=.18,color=colors[k],label=ENGLISH[k])
            ax.set_xticks([0,1],['Original','25% depth']);ax.set_title(f'{method} | {t}')
            ax.set_ylabel('RNA-fraction MAE (pp)');ax.set_ylim(0,key.loc[key.cell_type.eq(t),'mae_pp'].max()*1.1)
            ax.spines[['top','right']].set_visible(False)
    fig.legend(*axes[0,0].get_legend_handles_labels(),loc='outside lower center',ncol=4,frameon=False)
    fig.savefig(OUT/'independent_validation.png',dpi=200);fig.savefig(OUT/'independent_validation.pdf');plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(9,8),layout='constrained')
    for i,method in enumerate(['RCTD','NNLS']):
        for j,depth in enumerate([1.,.25]):
            ax=axes[i,j]
            for v,color,label in zip(VARIANTS,colors,ENGLISH):
                sub=titration[titration.method.eq(method)&titration.variant.eq(v)&titration.depth_fraction.eq(depth)]
                ax.scatter(sub.normal_rna_fraction*100,sub.normal_prediction*100,s=24,alpha=.65,color=color,label=label)
            ax.plot([0,100],[0,100],color='#444444',lw=1)
            ax.set(xlim=(-2,102),ylim=(-2,102),xlabel='Known normal epithelial RNA (%)',ylabel='Predicted normal epithelial weight (%)',title=f'{method} | depth {depth:g}')
    fig.legend(*axes[0,0].get_legend_handles_labels(),loc='outside lower center',ncol=2,frameon=False)
    fig.savefig(OUT/'epithelial_titration.png',dpi=200);fig.savefig(OUT/'epithelial_titration.pdf');plt.close(fig)
    lines=['# CID4463独立验证','',
        '## 1. 供者与固定设计','',
        'CID4463在临床补充表1第19行精确对应Case 4463，为未治疗ER+、IDC。作者病理描述包含小叶样生长区域，E-cadherin阳性；本轮沿用IDC记录。',
        '参考固定为CID4471、CID4535，CID4067已用于开发。本轮1138个验证细胞包括CAF 25个、CD8 75个、癌上皮659个和正常上皮26个。',
        '四套参考及NNLS特征均复用原文件，拟合前保存哈希和配方。237个基础混合生成原深度和25%深度两个版本，共474个模拟spot。',
        '每个深度包含162个CAF/CD8混合、50个纯群对照，以及25个正常/癌上皮比例梯度。B细胞缺失，未构造B背景；NK只有3个，保留1细胞NK背景，跳过5细胞NK纯群。',
        'Cycling T只有5个来源细胞，5次纯群使用同一组细胞。原深度计数相同，不能把5/5匹配理解为5组独立细胞的重复验证。source_mixture_diversity.csv记录各配方的不同来源细胞组合数。',
        '## 2. 混合配方RNA份额误差','',
        'MAE单位为百分点，两个深度分别列出，不把同一供者的模拟spot当作独立患者。','',
        '| 方法 | 参考 | 深度 | CAF | CD8 | 癌上皮 | 正常上皮 | PVL |',
        '|---|---|---|---:|---:|---:|---:|---:|']
    for method in ['RCTD','NNLS']:
        for v in VARIANTS:
            for depth in [1.,.25]:
                sub=key[key.method.eq(method)&key.variant.eq(v)&key.depth_fraction.eq(depth)].set_index('cell_type')
                values=' | '.join(f'{sub.loc[t,"mae_pp"]:.2f}' for t in ['CAFs','T cells CD8+','Cancer Epithelial','Normal Epithelial','PVL'])
                lines.append(f'| {method} | {NAMES[v]} | {depth:g} | {values} |')
    lines+=['','## 3. 正常上皮与癌上皮检查','',
        '比例梯度每个spot含10个上皮细胞，正常上皮个数为0、2、5、8、10，各重复5次。评价使用实际RNA份额；细胞数量比例单独保存。','',
        '| 方法 | 参考 | 深度 | 梯度正常上皮MAE | 纯正常主导匹配 | 纯正常平均正确权重 | 纯癌平均正常权重 |',
        '|---|---|---|---:|---:|---:|---:|']
    for method in ['RCTD','NNLS']:
        for v in VARIANTS:
            for depth in [1.,.25]:
                m=metrics[metrics.truth.eq('rna')&metrics.scenario.eq('epithelial_titration')&metrics.method.eq(method)&metrics.variant.eq(v)&metrics.depth_fraction.eq(depth)&metrics.cell_type.eq('Normal Epithelial')].iloc[0]
                p=pure[pure.method.eq(method)&pure.variant.eq(v)&pure.depth_fraction.eq(depth)].set_index('cell_type')
                normal=p.loc['Normal Epithelial'];cancer=p.loc['Cancer Epithelial']
                lines.append(f'| {method} | {NAMES[v]} | {depth:g} | {m.mae_pp:.2f} | {int(normal.dominant_matches)}/5 | {normal.mean_weight_on_true_type:.1%} | {cancer.mean_normal_weight:.1%} |')
    lines+=['','## 4. 本轮观察','',
        '原深度CAF/CD8混合中，RCTD上皮拆分的CAF、CD8和癌上皮MAE分别为1.19、1.08和2.63个百分点，原参考分别为1.28、1.13和10.26。25%深度下上皮拆分的对应误差为1.30、1.24和3.90。',
        '供者均衡的癌上皮MAE也降到2.95，但组合方案为4.37。开发供者CID4067上表现最好的方案，没有在CID4463的所有指标上继续排第一。',
        '上皮梯度中，RCTD单独拆分后的正常上皮MAE从16.39升到21.33个百分点；组合方案为14.53。纯正常上皮的平均正确权重分别为原参考66.3%、供者均衡68.6%、上皮拆分59.6%、组合71.6%。纯群主导类型虽然都是5/5匹配，仍有较多正常上皮RNA被分配给癌上皮。',
        '原深度纯CD4对照的主导类型匹配，原参考和上皮拆分均为3/5，供者均衡为1/5，组合方案为0/5。混合配方里的CD8误差较小，不能代替对CD4/CD8混淆的检查。',
        'Cycling T的主导类型为5/5匹配，但正确类型的平均权重只有31.5%至34.6%，且五次纯群使用相同来源细胞。这一结果不足以说明Cycling T定量已经准确。',
        '本轮支持保留上皮拆分作为CAF/CD8混合分析的候选，同时保留原参考作对照。正常/癌上皮梯度和CD4纯群显示不同方案存在取舍，目前不把任何一套替换为所有类型的最终参考。',
        '', '## 5. 供者使用与文件','',
        'CID3941、CID3948、CID4461保留为后续验证候选，本轮未提取其计数或查看预测。CID4290A和CID4530N的后缀对应关系尚未确认；CID4040缺少癌上皮；CID4398接受过治疗。',
        '本轮是一个新供者的独立验证，正常上皮来自26个作者标注细胞。四套参考的结果均报告，本轮没有按结果修改参考或参数。',
        'donor_candidates.csv保存全队列候选及原因；protocol.json和frozen_recipes.csv保存验证设计；validation_metrics.csv保存全部类型、配方和两种真值的指标。',
        '临床来源：[Wu等人的补充表](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-021-00911-1/MediaObjects/41588_2021_911_MOESM4_ESM.xlsx)，本地原件位于data/provenance/wu2021_supplementary_tables.xlsx。']
    (OUT/'interpretation.md').write_text('\n'.join(lines)+'\n')
    print(key[key.cell_type.isin(['CAFs','T cells CD8+','Cancer Epithelial','Normal Epithelial'])][['variant','method','depth_fraction','cell_type','mae_pp']].to_string(index=False))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--nnls-only',action='store_true')
    main(parser.parse_args().nnls_only)
