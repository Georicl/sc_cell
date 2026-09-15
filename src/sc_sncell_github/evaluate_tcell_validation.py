"""分别评价三位供者的T细胞梯度、零CD8配方及CAF配对变化。"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import anndata as ad
from scipy.optimize import nnls
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data/processed/Wu2021_tcell_validation_v1'
OUT=ROOT/'results/Wu2021_tcell_validation_v1'
DONORS=['CID3941','CID3948','CID4461']
VARIANTS=['baseline','donor_balanced','epithelial_split','balanced_split']
NAMES=dict(zip(VARIANTS,['原参考','供者均衡','上皮拆分','均衡并拆分']))
ENGLISH=['Baseline','Donor balanced','Epithelial split','Balanced + split']
COLORS=['#636363','#437b9e','#ca8a3d','#5a976c']

def main(nnls_only=False):
    metrics=[];pure_rows=[];spot_rows=[];pair_rows=[];diversity_rows=[]
    for donor in DONORS:
        out=OUT/donor
        design=pd.read_csv(out/'simulation_design.csv',index_col=0)
        truths={n:pd.read_csv(out/f'truth_{f}_fraction.csv',index_col=0) for n,f in [('rna','rna'),('cells','cell')]}
        types=truths['rna'].columns
        spots=ad.read_h5ad(DATA/donor/'simulated_spots.h5ad')
        assert spots.obs_names.equals(design.index)
        members=pd.read_csv(out/'mixture_membership.csv')
        source_sets=members.groupby('spot_id').source_barcode.agg(lambda s:'|'.join(sorted(s)))
        diversity=design[['scenario','n_cd8','n_caf','depth_fraction']].copy();diversity['source_set']=source_sets
        distinct=diversity.groupby(['scenario','n_cd8','n_caf','depth_fraction']).agg(n_spots=('source_set','size'),distinct_source_sets=('source_set','nunique')).reset_index();distinct['donor_id']=donor;diversity_rows.append(distinct)
        for variant in VARIANTS:
            folder=out/variant;folder.mkdir(exist_ok=True)
            ref=ROOT/'data/processed/Wu2021_known_composition_v1' if variant=='baseline' else ROOT/'data/processed/Wu2021_reference_comparison_v1'/variant
            profiles=pd.read_csv(ref/'nnls_reference_profiles.csv',index_col=0)
            y=spots[:,profiles.index].X.toarray()/design.nUMI.to_numpy()[:,None]
            raw=np.vstack([nnls(profiles.to_numpy(),row)[0] for row in y])
            assert (raw.sum(axis=1)>0).all()
            coef=pd.DataFrame(raw,index=design.index,columns=profiles.columns)
            coef.to_csv(folder/'nnls_raw_coefficients.csv');coef.div(coef.sum(axis=1),axis=0).to_csv(folder/'nnls_relative_weights.csv')
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
                pred=pred.loc[design.index,types];pred.to_csv(folder/f'{method.lower()}_parent_weights.csv')
                for truth_name,truth in truths.items():
                    for (depth,caf),group in design[design.scenario.eq('gradient')].groupby(['depth_fraction','n_caf']):
                        ids=group.index
                        for t in types:
                            actual=truth.loc[ids,t];estimated=pred.loc[ids,t];error=(estimated-actual)*100
                            present=truths['cells'].loc[ids,t]>0
                            metrics.append(dict(donor_id=donor,variant=variant,method=method,truth=truth_name,depth_fraction=depth,n_caf=caf,cell_type=t,n_spots=len(ids),mae_pp=error.abs().mean(),rmse_pp=np.sqrt(np.square(error).mean()),bias_pp=error.mean(),present_mae_pp=error[present].abs().mean(),n_absent=int((~present).sum()),absent_mean_prediction_pct=estimated[~present].mean()*100,absent_ge1pct_fraction=(estimated[~present]>=.01).mean() if (~present).any() else np.nan,absent_ge5pct_fraction=(estimated[~present]>=.05).mean() if (~present).any() else np.nan))
                for depth in [1.,.25]:
                    ids=design.index[design.scenario.eq('pure')&design.depth_fraction.eq(depth)]
                    actual=truths['cells'].loc[ids].idxmax(axis=1)
                    for t in sorted(actual.unique()):
                        keep=actual.index[actual.eq(t)]
                        pure_rows.append(dict(donor_id=donor,variant=variant,method=method,depth_fraction=depth,cell_type=t,n_spots=len(keep),dominant_matches=int(pred.loc[keep].idxmax(axis=1).eq(t).sum()),mean_weight_on_true_type=pred.loc[keep,t].mean(),mean_cd4_weight=pred.loc[keep,'T cells CD4+'].mean(),mean_cd8_weight=pred.loc[keep,'T cells CD8+'].mean()))
                for t in types:
                    frame=design.copy();frame['variant']=variant;frame['method']=method;frame['cell_type']=t
                    frame['truth_rna']=truths['rna'][t];frame['truth_cells']=truths['cells'][t];frame['prediction']=pred[t]
                    frame['rna_error_pp']=(pred[t]-truths['rna'][t])*100;spot_rows.append(frame)
                for (depth,pair_id),group in design[design.scenario.eq('gradient')].groupby(['depth_fraction','pair_id']):
                    ids=group.set_index('n_caf')
                    base=group.index[group.n_caf.eq(0)][0]
                    for caf in [1,3]:
                        target=group.index[group.n_caf.eq(caf)][0]
                        for t in ['CAFs','T cells CD4+','T cells CD8+']:
                            e0=(pred.loc[base,t]-truths['rna'].loc[base,t])*100
                            e1=(pred.loc[target,t]-truths['rna'].loc[target,t])*100
                            pair_rows.append(dict(donor_id=donor,variant=variant,method=method,depth_fraction=depth,pair_id=pair_id,n_caf=caf,cell_type=t,base_spot=base,caf_spot=target,error_change_pp=e1-e0,absolute_error_change_pp=abs(e1)-abs(e0)))
    if nnls_only:return
    metrics=pd.DataFrame(metrics);pure=pd.DataFrame(pure_rows);pairs=pd.DataFrame(pair_rows);spots=pd.concat(spot_rows)
    metrics.to_csv(OUT/'gradient_metrics.csv',index=False);pure.to_csv(OUT/'pure_control_summary.csv',index=False)
    spots.to_csv(OUT/'spot_predictions.csv',index_label='spot_id');pairs.to_csv(OUT/'paired_caf_effects.csv',index=False)
    pairs.groupby(['donor_id','variant','method','depth_fraction','n_caf','cell_type'])[['error_change_pp','absolute_error_change_pp']].mean().to_csv(OUT/'paired_caf_effect_summary.csv')
    pd.concat(diversity_rows).to_csv(OUT/'source_mixture_diversity.csv',index=False)
    key=metrics[metrics.truth.eq('rna')&metrics.cell_type.eq('T cells CD8+')]
    key.to_csv(OUT/'cd8_key_metrics.csv',index=False)
    for depth in [1.,.25]:
        fig,axes=plt.subplots(3,3,figsize=(13,11),layout='constrained')
        for i,donor in enumerate(DONORS):
            for j,caf in enumerate([0,1,3]):
                ax=axes[i,j]
                for variant,color,label in zip(VARIANTS,COLORS,ENGLISH):
                    sub=spots[spots.donor_id.eq(donor)&spots.variant.eq(variant)&spots.method.eq('RCTD')&spots.scenario.eq('gradient')&spots.depth_fraction.eq(depth)&spots.n_caf.eq(caf)&spots.cell_type.eq('T cells CD8+')]
                    means=sub.groupby('n_cd8')[['truth_rna','prediction']].mean()*100
                    ax.scatter(sub.truth_rna*100,sub.prediction*100,color=color,alpha=.25,s=12)
                    ax.plot(means.truth_rna,means.prediction,marker='o',markersize=4,color=color,label=label,lw=1.3)
                ax.plot([0,100],[0,100],color='#333333',ls='--',lw=.8)
                ax.set(xlim=(-2,102),ylim=(-2,102),xlabel='Known CD8 RNA (%)',ylabel='Predicted CD8 weight (%)',title=f'{donor} | CAF cells: {caf}')
        fig.legend(*axes[0,0].get_legend_handles_labels(),loc='outside lower center',ncol=4,frameon=False)
        fig.savefig(OUT/f'rctd_cd8_gradient_depth{int(depth*100)}.png',dpi=180)
        fig.savefig(OUT/f'rctd_cd8_gradient_depth{int(depth*100)}.pdf');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    for ax,method in zip(axes,['RCTD','NNLS']):
        table=key[key.method.eq(method)&key.depth_fraction.eq(1)].pivot(index=['donor_id','n_caf'],columns='variant',values='absent_mean_prediction_pct').reindex(columns=VARIANTS)
        im=ax.imshow(table,vmin=0,vmax=100,cmap='YlOrRd',aspect='auto')
        ax.set_xticks(range(4),ENGLISH,rotation=30,ha='right',fontsize=9)
        ax.set_yticks(range(len(table)),[f'{d} | CAF {c}' for d,c in table.index],fontsize=9)
        for i in range(len(table)):
            for j in range(4):
                value=table.iloc[i,j];ax.text(j,i,f'{value:.1f}',ha='center',va='center',color='white' if value>55 else 'black',fontsize=9)
        ax.set_title(f'{method} | CD8 absent, original depth')
    fig.colorbar(im,ax=axes,label='Mean predicted CD8 weight (%)',shrink=.8)
    fig.savefig(OUT/'cd8_zero_predictions.png',dpi=180);fig.savefig(OUT/'cd8_zero_predictions.pdf');plt.close(fig)
    lines=['# 三位供者的CD4/CD8与CAF背景验证','',
        '## 1. 固定设计','',
        'CID3941、CID3948、CID4461均为未治疗ER+、IDC，临床补充表1分别对应第8、10、18行。三位供者此前保留为验证候选，本轮未加入参考。',
        '每个梯度spot含10个T细胞，CD8个数为0、2、5、8、10，其余为CD4；分别加入0、1、3个CAF。每种配方重复5次，另为CAF、CD4和CD8各生成5个五细胞纯群。',
        '每位供者90个基础混合，在原深度和25%深度各生成一个版本，共180个spot；三位合计540个。相同T细胞及其抽稀计数跨CAF背景配对，CAF 1使用CAF 3中的第一个细胞。',
        'CID3941只有8个CAF，CID4461只有11个CD8，配方之间会复用来源细胞。不同来源细胞组合数保存在source_mixture_diversity.csv。',
        '四套参考、NNLS特征和RCTD参数均在预测前固定；RCTD按供者分别进行平台校正和拟合。本轮是T细胞与CAF组成的针对性模拟，不混入其他细胞类型。',
        '', '## 2. RCTD原深度CD8预测','',
        '每行MAE来自25个梯度spot；CD8实际缺失的配方为其中5个。MAE单位为百分点，缺失预测均值单位为百分比。','',
        '| 供者 | 参考 | CAF个数 | CD8 MAE | CD8缺失时均值 | 缺失时≥1% | 缺失时≥5% |',
        '|---|---|---:|---:|---:|---:|---:|']
    for donor in DONORS:
        for v in VARIANTS:
            for caf in [0,1,3]:
                r=key[key.donor_id.eq(donor)&key.variant.eq(v)&key.method.eq('RCTD')&key.depth_fraction.eq(1)&key.n_caf.eq(caf)].iloc[0]
                lines.append(f'| {donor} | {NAMES[v]} | {caf} | {r.mae_pp:.2f} | {r.absent_mean_prediction_pct:.2f} | {r.absent_ge1pct_fraction:.0%} | {r.absent_ge5pct_fraction:.0%} |')
    lines+=['','## 3. 本轮观察','',
        '原深度、没有CAF的梯度中，原参考RCTD的CD8 MAE在CID3941、CID3948和CID4461分别为20.02、24.29和34.48个百分点。单独上皮拆分分别为20.46、23.91和34.10，没有消除T细胞定量误差。',
        'CD8实际为零时，CID3941的原参考平均预测为7.40%，供者均衡为37.72%，组合方案为43.40%；CID3948分别为32.15%、41.85%和46.68%。这些数值来自10个CD4、没有CAF的配方，每组5个spot。',
        'CID4461的零CD8预测较低，但真实CD8存在时严重低估。原参考五细胞纯CD8对照中，正确类型平均权重为20.3%，主导类型匹配2/5；另外三套参考也只有2/5匹配。低假阳性不能替代对真实CD8的识别。',
        '已知组成沿用作者细胞标签。CID4461的11个CD8细胞UMI中位数为400，CID3941和CID3948分别为1641和1160；CID4461的低估可能同时受来源细胞表达信息较少和参考覆盖影响，本轮没有把两者分开。',
        '相同CD4细胞、CD8始终为零，CAF从0个增加到3个时，CID3948原参考的预测CD8从32.15%降到5.38%，CID3941从7.40%降到2.22%。模型在这些零CD8配方中也能产生CAF增加、预测CD8下降的表象。因此真实切片上的权重负相关本身不足以说明CAF与CD8发生空间排斥。',
        '25%深度下，没有CAF的原参考CD8 MAE分别为16.83、20.31和34.29个百分点，误差仍然存在。3个CAF背景下，各供者及方案的RCTD CAF MAE为6.47至16.89个百分点，CAF定量也随供者和配方变化。',
        'NNLS原参考在没有CAF的原深度梯度中，三位供者CD8 MAE分别为14.14、16.81和49.40个百分点；它在前两位低于RCTD，在CID4461更差，不能据此直接替换为统一方法。',
        '此前CID4463的癌上皮背景配方中，CD8真实RNA份额均值只有2.83%，中位数1.20%；本轮覆盖0%至100%的T细胞RNA组成。此前较小的绝对误差不能直接推广到完整比例范围，也不能推广到其他供者。',
        '本轮没有一套参考解决了三个供者的CD8识别。原参考和上皮拆分仍可作为后续诊断对照；供者均衡及组合方案不能仅凭之前的上皮改善就用于正式CAF/CD8关系推断。',
        '', '## 4. 纯CD4和纯CD8对照','',
        '下表为原深度五细胞纯群，每类5个。主导匹配只检查权重最高的标签，同时列出纯CD4中错误分配给CD8的平均权重。','',
        '| 供者 | 方法 | 参考 | 纯CD4主导匹配 | 纯CD4的CD8权重 | 纯CD8主导匹配 | 纯CD8正确权重 |',
        '|---|---|---|---:|---:|---:|---:|']
    for donor in DONORS:
        for method in ['RCTD','NNLS']:
            for v in VARIANTS:
                sub=pure[pure.donor_id.eq(donor)&pure.method.eq(method)&pure.variant.eq(v)&pure.depth_fraction.eq(1)].set_index('cell_type')
                cd4,cd8=sub.loc['T cells CD4+'],sub.loc['T cells CD8+']
                lines.append(f'| {donor} | {method} | {NAMES[v]} | {int(cd4.dominant_matches)}/5 | {cd4.mean_cd8_weight:.1%} | {int(cd8.dominant_matches)}/5 | {cd8.mean_weight_on_true_type:.1%} |')
    lines+=['','## 5. 配对变化与完整文件','',
        '加入CAF会改变真实CD8 RNA份额。配对评价使用“预测减真实”的误差变化，以及绝对误差变化，避免把正常的RNA稀释当作错误。正的绝对误差变化表示加入CAF后误差增大。',
        'gradient_metrics.csv包含全部13类、两种真值和两个深度；cd8_key_metrics.csv汇总CD8；paired_caf_effects.csv保存逐配对变化；pure_control_summary.csv保存纯群表现。',
        '本轮未根据验证结果修改参考或参数，也没有把540个模拟spot当作540位独立供者。']
    (OUT/'interpretation.md').write_text('\n'.join(lines)+'\n')
    print(key[key.method.eq('RCTD')&key.depth_fraction.eq(1)][['donor_id','variant','n_caf','mae_pp','absent_mean_prediction_pct']].to_string(index=False))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--nnls-only',action='store_true')
    main(parser.parse_args().nnls_only)
