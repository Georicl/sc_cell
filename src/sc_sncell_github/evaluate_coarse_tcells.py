"""比较细分拟合后汇总与粗分类重新拟合，按原配方分别评价。"""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy.optimize import nnls
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data/processed/Wu2021_coarse_T_v1'
OUT=ROOT/'results/Wu2021_coarse_T_v1'
APPROACHES=['fine_fit_then_sum','coarse_reference_refit']
NAMES={'fine_fit_then_sum':'细分拟合后汇总','coarse_reference_refit':'粗分类重新拟合'}
def main(nnls_only=False):
    protocol=json.loads((OUT/'protocol.json').read_text())
    mapping=pd.read_csv(OUT/'label_mapping.csv',index_col=0).reference_label
    profiles=pd.read_csv(DATA/'nnls_reference_profiles.csv',index_col=0)
    types=sorted(mapping.unique());rows=[];pure_rows=[];errors=[];contrast_rows=[]
    for entry in protocol['datasets']:
        donor=entry['donor_id'];inp=ROOT/entry['input_dir'];res=ROOT/entry['result_dir'];out=OUT/donor;out.mkdir(exist_ok=True)
        design=pd.read_csv(res/'simulation_design.csv',index_col=0)
        fine_truths={n:pd.read_csv(res/f'truth_{f}_fraction.csv',index_col=0) for n,f in [('rna','rna'),('cells','cell')]}
        truths={n:t.T.groupby(mapping).sum().T.reindex(columns=types) for n,t in fine_truths.items()}
        for n,t in truths.items():t.to_csv(out/f'truth_{n}_coarse.csv')
        spots=ad.read_h5ad(inp/'simulated_spots.h5ad')
        assert spots.obs_names.equals(design.index)
        x=spots[:,profiles.index].X.toarray()/design.nUMI.to_numpy()[:,None]
        coefficients=np.vstack([nnls(profiles.to_numpy(),v)[0] for v in x])
        assert (coefficients.sum(axis=1)>0).all()
        coef=pd.DataFrame(coefficients,index=design.index,columns=profiles.columns)
        coef.to_csv(out/'nnls_raw_coefficients.csv');coef.div(coef.sum(axis=1),axis=0).to_csv(out/'nnls_relative_weights.csv')
        if nnls_only:continue
        for method in ['RCTD','NNLS']:
            fine=pd.read_csv(ROOT/entry['fine_prediction_dir']/f'{method.lower()}_relative_weights.csv',index_col=0)
            aggregated=fine.T.groupby(mapping).sum().T
            coarse=pd.read_csv(out/f'{method.lower()}_relative_weights.csv',index_col=0)
            predictions=dict(zip(APPROACHES,[aggregated,coarse]))
            for approach,pred in predictions.items():
                assert set(pred.index)==set(design.index) and set(pred.columns)==set(types)
                pred=pred.loc[design.index,types]
                assert np.isfinite(pred.values).all() and (pred.values>=0).all() and np.allclose(pred.sum(axis=1),1)
                pred.to_csv(out/f'{method.lower()}_{approach}_weights.csv')
                for truth_name,truth in truths.items():
                    for (scenario,depth),group in design.groupby(['scenario','depth_fraction']):
                        ids=group.index
                        for t in types:
                            estimated=pred.loc[ids,t];actual=truth.loc[ids,t];error=(estimated-actual)*100
                            present=truths['cells'].loc[ids,t]>0
                            rows.append(dict(donor_id=donor,method=method,approach=approach,truth=truth_name,scenario=scenario,depth_fraction=depth,cell_type=t,n_spots=len(ids),mae_pp=error.abs().mean(),rmse_pp=np.sqrt((error**2).mean()),bias_pp=error.mean(),present_mae_pp=error[present].abs().mean(),n_absent=int((~present).sum()),absent_mean_prediction_pct=estimated[~present].mean()*100,absent_ge1pct_fraction=(estimated[~present]>=.01).mean() if (~present).any() else np.nan))
                for t in types:
                    e=design.copy();e['donor_id']=donor;e['method']=method;e['approach']=approach;e['cell_type']=t
                    e['truth_rna']=truths['rna'][t];e['truth_cells']=truths['cells'][t];e['prediction']=pred[t];e['rna_error_pp']=(pred[t]-truths['rna'][t])*100
                    errors.append(e)
                for depth in [1.,.25]:
                    ids=design.index[design.scenario.eq('pure')&design.depth_fraction.eq(depth)]
                    actual=fine_truths['cells'].loc[ids].idxmax(axis=1)
                    for original in sorted(actual.unique()):
                        keep=actual.index[actual.eq(original)];parent=mapping.loc[original]
                        pure_rows.append(dict(donor_id=donor,method=method,approach=approach,depth_fraction=depth,original_type=original,coarse_type=parent,n_spots=len(keep),dominant_matches=int(pred.loc[keep].idxmax(axis=1).eq(parent).sum()),mean_correct_weight=pred.loc[keep,parent].mean(),mean_T_weight=pred.loc[keep,'T cells'].mean(),mean_NK_weight=pred.loc[keep,'NK cells'].mean()))
                if 'n_caf' in design:
                    for (depth,caf),group in design[design.scenario.eq('gradient')].groupby(['depth_fraction','n_caf']):
                        for t in ['T cells','CAFs']:
                            error=(pred.loc[group.index,t]-truths['rna'].loc[group.index,t])*100
                            contrast_rows.append(dict(donor_id=donor,method=method,approach=approach,depth_fraction=depth,n_caf=caf,cell_type=t,mae_pp=error.abs().mean(),bias_pp=error.mean()))
    if nnls_only:return
    metrics=pd.DataFrame(rows);pure=pd.DataFrame(pure_rows);errors=pd.concat(errors)
    metrics.to_csv(OUT/'comparison_metrics.csv',index=False);pure.to_csv(OUT/'pure_control_summary.csv',index=False)
    errors.to_csv(OUT/'spot_errors.csv',index_label='spot_id');pd.DataFrame(contrast_rows).to_csv(OUT/'caf_background_metrics.csv',index=False)
    target=metrics[metrics.truth.eq('rna')&metrics.cell_type.eq('T cells')&metrics.scenario.isin(['factorial','gradient'])]
    target.to_csv(OUT/'T_cell_key_metrics.csv',index=False)
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    donors=[e['donor_id'] for e in protocol['datasets']]
    for ax,method in zip(axes,['RCTD','NNLS']):
        for offset,approach,color,label in [(-.18,APPROACHES[0],'#617d91','Fine fit, then sum'),(.18,APPROACHES[1],'#c08a45','Coarse refit')]:
            sub=target[target.method.eq(method)&target.approach.eq(approach)&target.depth_fraction.eq(1)].set_index('donor_id').loc[donors]
            ax.bar(np.arange(5)+offset,sub.mae_pp,width=.35,color=color,label=label)
        ax.set_xticks(range(5),donors,rotation=25);ax.set_ylabel('T-cell RNA-fraction MAE (pp)');ax.set_title(f'{method} | original depth')
        ax.set_ylim(0,target[target.depth_fraction.eq(1)].mae_pp.max()*1.15);ax.spines[['top','right']].set_visible(False)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='outside lower center',ncol=2,frameon=False)
    fig.savefig(OUT/'coarse_T_comparison.png',dpi=200);fig.savefig(OUT/'coarse_T_comparison.pdf');plt.close(fig)
    # 原始细分类纯群分别展示，保留NK负对照。
    display_pure=pure[pure.method.eq('RCTD')&pure.depth_fraction.eq(1)&pure.original_type.isin(protocol['merge_to_T_cells']+['NK cells'])]
    table=display_pure.pivot(index=['donor_id','original_type'],columns='approach',values='mean_T_weight').reindex(columns=APPROACHES)
    fig,ax=plt.subplots(figsize=(7,max(5,len(table)*.32)),layout='constrained')
    im=ax.imshow(table,vmin=0,vmax=1,cmap='YlGnBu',aspect='auto')
    ax.set_xticks([0,1],['Fine fit, then sum','Coarse refit']);ax.set_yticks(range(len(table)),[f'{d} | {t}' for d,t in table.index],fontsize=8)
    for i in range(len(table)):
        for j in range(2):ax.text(j,i,f'{table.iloc[i,j]:.2f}',ha='center',va='center',color='white' if table.iloc[i,j]>.6 else 'black',fontsize=8)
    fig.colorbar(im,ax=ax,label='Mean predicted T-cell weight',shrink=.6)
    fig.savefig(OUT/'pure_T_and_NK_controls.png',dpi=180);fig.savefig(OUT/'pure_T_and_NK_controls.pdf');plt.close(fig)
    lines=['# T细胞粗分类比较','',
        '## 1. 分类与比较方式','',
        '将作者CD4、CD8、NKT和Cycling T合并为T cells，NK单独保留，其余类型不变，共10类。这里的T cells是这四个作者标签的操作性汇总，NKT标签未进一步解释为经TCR确认的特定亚群。',
        '参考使用CID4471、CID4535的原12570个细胞，其中合并后的T细胞982个。本轮没有均衡供者，也没有删除低信息细胞。',
        '比较两种方式：先用原13类模型拟合再相加；直接使用10类参考重新拟合。两个方法都以相加后的已知RNA份额和细胞份额评价。',
        '复用此前五位供者的1502个模拟spot，未生成新的模拟组成。这些供者已经参与前期检查，本轮用于参考方案修订比较。',
        '', '## 2. T细胞RNA份额误差','',
        'CID4067和CID4463使用原CAF/CD8混合配方；另外三位使用T细胞梯度及CAF背景配方。逐供者列出，不跨不同配方合并一个总体误差。单位为百分点。','',
        '| 供者 | 方法 | 深度 | 细分后汇总MAE | 粗分类重拟合MAE |',
        '|---|---|---:|---:|---:|']
    for donor in donors:
        for method in ['RCTD','NNLS']:
            for depth in [1.,.25]:
                sub=target[target.donor_id.eq(donor)&target.method.eq(method)&target.depth_fraction.eq(depth)].set_index('approach')
                lines.append(f'| {donor} | {method} | {depth:g} | {sub.loc[APPROACHES[0],"mae_pp"]:.2f} | {sub.loc[APPROACHES[1],"mae_pp"]:.2f} |')
    lines+=['','## 3. 本轮采用的方式','',
        '当前采用细分拟合后汇总作为T细胞大类的输出方式。RCTD原深度下，五位供者中四位的T细胞MAE低于粗分类重拟合；CID4067是例外，重拟合为1.62个百分点，汇总为1.92。两个深度下这项排序一致。',
        'CID3941、CID3948、CID4461的细分后汇总MAE分别为4.68、7.18、15.10个百分点；粗分类重拟合分别为10.60、14.27、23.77。直接用一个平均T细胞参考替代不同状态，并没有带来更好的定量。',
        '纯CD8对照中，CID3948细分后汇总的平均T权重为87.8%，直接重拟合为61.6%；后者还将平均28.4%分配到NK。CID4461对应的T权重为68.2%和55.3%，低信息来源细胞的问题仍然存在。',
        'CID4067纯NK对照中，细分后汇总仍平均给T细胞16.6%的权重，直接重拟合为6.2%。汇总解决的是T细胞内部标签分配对总量的影响，不能自动解决T与NK之间的误分配。',
        'CAF指标在细分后汇总中保持原值；例如三位T细胞梯度供者的原深度CAF MAE分别为6.90、6.30和7.69个百分点。T细胞大类与原CD8是不同的评价对象，不能把两者MAE之差当成同一个目标准确性提高。',
        '', '## 4. 纯群及缺失类型','',
        'pure_control_summary.csv保留原始纯群标签，分别查看CD4、CD8、NKT、Cycling T合并后的表现，以及NK是否被错误分配到T细胞。不能只看原有CD4/CD8互相混淆是否消失。',
        'comparison_metrics.csv保存全部10类、各配方、两种真值和两个深度；caf_background_metrics.csv按CAF背景分别评价。',
        '粗分类改变了研究对象：现在评价T细胞整体RNA贡献，不再把这项结果解释为CD8数量或CD8特异空间分布。',
        '', '## 5. CID4535空间结果','',
        'spatial_CID4535目录保存原空间拟合的同定义汇总及CAF/T图。该图使用原CID4535单供者参考的1102个spot，与本轮双供者参考比较分开记录；没有据此宣称新参考的空间准确性已得到验证。',
        '同时汇总既有10次参考重抽样，用于查看T细胞大类的分布稳定性。没有开展CAF/T因果或空间排斥判断。']
    (OUT/'interpretation.md').write_text('\n'.join(lines)+'\n')
    print(target[['donor_id','method','approach','depth_fraction','mae_pp']].round(3).to_string(index=False))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--nnls-only',action='store_true')
    main(parser.parse_args().nnls_only)
