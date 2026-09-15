"""在相同模拟spot上比较参考方案，所有预测汇总回原13类。"""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy.optimize import nnls
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'results/Wu2021_known_composition_v1'
DATA = ROOT / 'data/processed/Wu2021_reference_comparison_v1'
OUT = ROOT / 'results/Wu2021_reference_comparison_v1'
VARIANTS = ['baseline', 'donor_balanced', 'epithelial_split', 'balanced_split']
LABELS = ['Baseline', 'Donor balanced', 'Epithelial split', 'Balanced + split']

def main(nnls_only=False):
    design = pd.read_csv(BASE/'simulation_design.csv', index_col=0)
    truths = {name: pd.read_csv(BASE/f'truth_{filename}_fraction.csv', index_col=0)
              for name, filename in [('rna', 'rna'), ('cells', 'cell')]}
    types = truths['rna'].columns
    spots = ad.read_h5ad(ROOT/'data/processed/Wu2021_known_composition_v1/simulated_spots.h5ad')
    assert spots.obs_names.equals(design.index)
    for variant in VARIANTS[1:]:
        profiles = pd.read_csv(DATA/variant/'nnls_reference_profiles.csv', index_col=0)
        y = spots[:, profiles.index].X.toarray()/design.nUMI.to_numpy()[:, None]
        coefficients = np.vstack([nnls(profiles.to_numpy(), row)[0] for row in y])
        assert (coefficients.sum(axis=1)>0).all()
        raw = pd.DataFrame(coefficients,index=design.index,columns=profiles.columns)
        raw.to_csv(OUT/variant/'nnls_raw_coefficients.csv')
        raw.div(raw.sum(axis=1),axis=0).to_csv(OUT/variant/'nnls_relative_weights.csv')
    if nnls_only:
        return
    predictions, rows, pure_rows = {}, [], []
    for variant in VARIANTS:
        folder = BASE if variant == 'baseline' else OUT/variant
        for method in ['RCTD', 'NNLS']:
            pred = pd.read_csv(folder/f'{method.lower()}_relative_weights.csv', index_col=0)
            assert pred.index.is_unique and set(pred.index)==set(design.index)
            assert np.isfinite(pred.values).all() and (pred.values>=0).all() and np.allclose(pred.sum(axis=1),1)
            if variant != 'baseline':
                mapping = pd.read_csv(folder/'label_mapping.csv', index_col=0).parent_label
                assert set(pred.columns)==set(mapping.index)
                pred = pred.T.groupby(mapping).sum().T
            assert set(pred.columns)==set(types)
            pred = pred.loc[design.index, types]
            predictions[(variant,method)] = pred
            if variant != 'baseline': pred.to_csv(folder/f'{method.lower()}_parent_weights.csv')
            for truth_name, truth in truths.items():
                for (scenario,depth), group in design.groupby(['scenario','depth_fraction']):
                    ids = group.index
                    for cell_type in list(types)+['Epithelial total']:
                        cols = ['Cancer Epithelial','Normal Epithelial'] if cell_type=='Epithelial total' else [cell_type]
                        actual = truth.loc[ids,cols].sum(axis=1)
                        estimated = pred.loc[ids,cols].sum(axis=1)
                        error = (estimated-actual)*100
                        present = truths['cells'].loc[ids,cols].sum(axis=1)>0
                        rows.append(dict(variant=variant,method=method,truth=truth_name,scenario=scenario,
                            depth_fraction=depth,cell_type=cell_type,n_spots=len(ids),
                            mae_pp=error.abs().mean(),rmse_pp=np.sqrt(np.square(error).mean()),bias_pp=error.mean(),
                            present_mae_pp=error[present].abs().mean(),
                            absent_mean_prediction_pct=estimated[~present].mean()*100,
                            absent_ge1pct_fraction=(estimated[~present]>=.01).mean() if (~present).any() else np.nan))
            for depth in [1.,.25]:
                ids = design.index[design.scenario.eq('pure') & design.depth_fraction.eq(depth)]
                actual = truths['cells'].loc[ids].idxmax(axis=1)
                for cell_type in sorted(actual.unique()):
                    keep = actual.index[actual.eq(cell_type)]
                    pure_rows.append(dict(variant=variant,method=method,depth_fraction=depth,cell_type=cell_type,
                        n_spots=len(keep),dominant_matches=pred.loc[keep].idxmax(axis=1).eq(cell_type).sum(),
                        mean_weight_on_true_type=pred.loc[keep,cell_type].mean(),
                        mean_normal_epithelial_weight=pred.loc[keep,'Normal Epithelial'].mean()))
    metrics = pd.DataFrame(rows)
    metrics.to_csv(OUT/'comparison_metrics.csv',index=False)
    pd.DataFrame(pure_rows).to_csv(OUT/'pure_control_summary.csv',index=False)
    main_metrics = metrics[metrics.truth.eq('rna') & metrics.scenario.eq('factorial')]
    selected = ['CAFs','T cells CD8+','Cancer Epithelial','Normal Epithelial','Epithelial total']
    main_metrics[main_metrics.cell_type.isin(selected)].to_csv(OUT/'key_metrics.csv',index=False)
    deltas=[]
    for (variant,method), pred in predictions.items():
        if variant=='baseline': continue
        baseline=predictions[('baseline',method)]
        for cell_type in types:
            delta=(pred[cell_type]-truths['rna'][cell_type]).abs()-(baseline[cell_type]-truths['rna'][cell_type]).abs()
            frame=design[['base_mixture','scenario','depth_fraction']].copy()
            frame['variant']=variant;frame['method']=method;frame['cell_type']=cell_type
            frame['absolute_error_change_pp']=delta*100
            deltas.append(frame)
    pd.concat(deltas).to_csv(OUT/'paired_error_changes.csv',index_label='spot_id')
    colors=['#636363','#437b9e','#ca8a3d','#5a976c']
    fig,axes=plt.subplots(2,3,figsize=(13,7),layout='constrained')
    for i,method in enumerate(['RCTD','NNLS']):
        for j,cell_type in enumerate(['CAFs','T cells CD8+','Cancer Epithelial']):
            ax=axes[i,j]
            for k,variant in enumerate(VARIANTS):
                sub=main_metrics[main_metrics.method.eq(method)&main_metrics.variant.eq(variant)&main_metrics.cell_type.eq(cell_type)].set_index('depth_fraction')
                ax.bar(np.arange(2)+(k-1.5)*.19,sub.loc[[1.,.25],'mae_pp'],width=.18,color=colors[k],label=LABELS[k])
            ax.set_xticks([0,1],['Original','25% depth'])
            ax.set_title(f'{method} | {cell_type}')
            ax.set_ylabel('RNA-fraction MAE (pp)')
            ax.spines[['top','right']].set_visible(False)
    for j, cell_type in enumerate(['CAFs','T cells CD8+','Cancer Epithelial']):
        upper=main_metrics.loc[main_metrics.cell_type.eq(cell_type),'mae_pp'].max()*1.1
        for ax in axes[:,j]: ax.set_ylim(0,upper)
    handles, labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=4,frameon=False)
    fig.savefig(OUT/'reference_comparison.png',dpi=200)
    fig.savefig(OUT/'reference_comparison.pdf')
    plt.close(fig)
    lines=['# 参考方案比较','',
        '参考为CID4471和CID4535，开发测试使用CID4067的同一批488个模拟spot。原深度和25%深度共享244个基础混合。', '',
        '## 1. 本轮参考设置','',
        '原参考含12570个细胞。供者均衡后保留3799个细胞，每个原类型按较少供者的数量抽样；Cycling T均衡后不足25个，因此保留原27个。抽样使用固定种子20260916。',
        '上皮拆分保留作者LumA和LumB标签，其余癌上皮归为Cancer Other。单独拆分使用全部参考细胞；组合方案使用与供者均衡方案完全相同的细胞。预测后将三类癌上皮相加，按原13类比较。',
        '正常上皮没有删除；它在开发测试中缺失，仍然用于检查误分配。', '',
        '## 2. RNA份额误差','',
        '指标为混合配方MAE，单位为百分点；原深度和25%深度各189个spot。', '',
        '| 方法 | 参考方案 | 深度 | CAF | CD8 | 癌上皮 | 正常上皮 | 上皮合计 |',
        '|---|---|---|---:|---:|---:|---:|---:|']
    names=dict(zip(VARIANTS,['原参考','供者均衡','上皮拆分','均衡并拆分']))
    for method in ['RCTD','NNLS']:
        for variant in VARIANTS:
            for depth in [1.,.25]:
                sub=main_metrics[main_metrics.method.eq(method)&main_metrics.variant.eq(variant)&main_metrics.depth_fraction.eq(depth)].set_index('cell_type')
                values=' | '.join(f'{sub.loc[t,"mae_pp"]:.2f}' for t in selected)
                lines.append(f'| {method} | {names[variant]} | {depth:g} | {values} |')
    lines += ['', '## 3. 本轮观察','',
        'RCTD原深度下，供者均衡将癌上皮MAE从31.39降至15.41个百分点；正常上皮平均错误权重从36.32%降至18.24%。单独上皮拆分的癌上皮MAE为16.72，组合方案为17.37，两项调整没有叠加改善癌上皮识别。',
        '组合方案的CAF MAE为2.30个百分点，原参考为3.60；CD8为1.74，原参考为1.73。单独上皮拆分的CD8 MAE为1.65。25%深度下，各方案的CAF、CD8及癌上皮MAE排序与原深度一致。',
        '供者均衡也使RCTD的PVL MAE从3.49升至3.83，CD4从2.61升至2.86个百分点。没有一套方案在所有类型上都更好。',
        '纯癌上皮对照中，供者均衡后的平均正常上皮权重仍为20.63%，原参考为40.45%。原深度下四套参考的癌上皮主导类型均为5/5匹配，单看主导类型会漏掉这部分误分配。Cycling T仍为0/5匹配。',
        'NNLS的上皮拆分将癌上皮MAE从42.86降至23.71；供者均衡则使CAF MAE从4.21升至6.51。参考调整在两种方法上的影响不同。',
        '当前可以保留供者均衡和上皮拆分作为后续比较方案，优先补充供者与亚型覆盖。正常上皮在本轮开发测试中全部缺失，均衡方案又将其参考从1988个减至44个；下一轮还需要包含正常上皮的独立供者，检查正常上皮识别是否保持。',
        '', '## 4. 如何阅读这次比较','',
        '正常上皮的真值全部为零，所以该列MAE也等于平均错误分配的权重百分比。上皮合计把癌上皮和正常上皮相加，用于查看总量是否接近已知组成。',
        '供者均衡同时减少了参考细胞数量，例如正常上皮从1988个降为44个，CAF从1394个降为204个；结果变化不能全部归因于供者权重。',
        'LumA参考全部来自CID4471，LumB参考全部来自CID4535。亚型与供者相互对应，本轮不能单独分离两者的影响。',
        '本轮复用已经查看过误差的CID4067，只用于比较参考方案；尚未执行新供者的独立检验。', '',
        '## 5. 参数和文件','',
        'RCTD使用full模式、2个核心、种子9，参考UMI最低100，每个模型类型至少25个细胞。NNLS继续按参考侧类型均值选择2000个基因。拆分类别后，各方法会按新的参考重新选择基因。',
        '参数用法核对：[Bioconductor spacexr文档](https://www.bioconductor.org/packages/release/bioc/manuals/spacexr/man/spacexr.pdf)。本机版本及实际参数另存于各方案的sessionInfo和RDS文件。',
        'comparison_metrics.csv包含全部类型、两种真值和两类配方；pure_control_summary.csv保存纯群识别；paired_error_changes.csv保存相对于原参考的逐spot误差变化。']
    (OUT/'interpretation.md').write_text('\n'.join(lines)+'\n')
    print(main_metrics[main_metrics.cell_type.isin(selected)][['variant','method','depth_fraction','cell_type','mae_pp']].to_string(index=False))

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--nnls-only',action='store_true')
    main(parser.parse_args().nnls_only)
