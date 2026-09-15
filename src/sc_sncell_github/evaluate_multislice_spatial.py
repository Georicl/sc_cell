"""逐切片核对CAF/T空间表达，比较CID4535新旧参考并汇总扩展结果。"""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data/processed/Wu2021_multislice_CAF_T_v1';OUT=ROOT/'results/Wu2021_multislice_CAF_T_v1'
RAW=ROOT/'data/raw/wu2021_visium/spatial'
MAPPING=pd.read_csv(ROOT/'results/Wu2021_coarse_T_v1/label_mapping.csv',index_col=0).reference_label
PROTOCOL=json.loads((OUT/'protocol.json').read_text())
PANELS=PROTOCOL['marker_checks']

def high(series):
    return series.sort_values(ascending=False,kind='stable').head(int(np.ceil(len(series)*.1))).index

def context(sample,obs):
    folder=RAW/f'{sample}_spatial';scale=json.loads((folder/'scalefactors_json.json').read_text())
    image=plt.imread(folder/'tissue_hires_image.png')
    xy=obs[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy()*scale['tissue_hires_scalef']
    radius=scale['spot_diameter_fullres']*scale['tissue_hires_scalef']/2
    return image,xy,radius

def draw(fig,ax,ctx,values,title,label,vmin=0,vmax=None,cmap='viridis'):
    image,xy,radius=ctx;ax.imshow(image,alpha=.32)
    p=PatchCollection([Circle((x,y),radius) for x,y in xy],linewidth=0,cmap=cmap,rasterized=True)
    p.set_array(np.asarray(values));p.set_clim(vmin,vmax);ax.add_collection(p)
    ax.set_xlim(xy[:,0].min()-50,xy[:,0].max()+50);ax.set_ylim(xy[:,1].max()+50,xy[:,1].min()-50)
    ax.set_aspect('equal');ax.axis('off');ax.set_title(title,fontsize=10)
    fig.colorbar(p,ax=ax,orientation='horizontal',label=label,shrink=.8)

def main(samples):
    inventory=pd.read_csv(OUT/'slice_inventory.csv').set_index('sample_id')
    all_summary=[];all_marker=[];comparison=[]
    for sample in samples:
        out=OUT/sample;a=ad.read_h5ad(DATA/sample/'spatial_counts_full.h5ad')
        raw=pd.read_csv(out/'rctd_raw_weights.csv',index_col=0)
        assert raw.index.is_unique and set(raw.index).issubset(a.obs_names)
        assert np.isfinite(raw.values).all() and (raw.values>=0).all() and (raw.sum(axis=1)>0).all()
        fine=raw.div(raw.sum(axis=1),axis=0)
        assert set(fine.columns)==set(MAPPING.index)
        weights=fine.T.groupby(MAPPING).sum().T
        weights.to_csv(out/'coarse_relative_weights.csv')
        obs=a.obs.loc[weights.index].copy()
        library=obs.total_counts.to_numpy()
        scores=pd.DataFrame(index=weights.index);marker_frames=[]
        for name,genes in PANELS.items():
            assert set(genes).issubset(a.var_names)
            counts=a[weights.index,genes].X.toarray()
            norm=counts/library[:,None]*10000
            scores[name+'_score']=np.log1p(norm).mean(axis=1)
            scores[name+'_genes_detected']=(counts>0).sum(axis=1)
            scores[name+'_panel_umi']=counts.sum(axis=1)
            for i,gene in enumerate(genes):
                marker_frames.append(pd.DataFrame(dict(gene=gene,raw_umi=counts[:,i],cp10k=norm[:,i]),index=weights.index))
        scores.to_csv(out/'marker_scores.csv');pd.concat(marker_frames).to_csv(out/'marker_gene_values.csv',index_label='spot_id')
        # 扩展后新增的诊断：不改主结果，只分开查看Cycling T分量。
        noncycling=fine[['T cells CD4+','T cells CD8+','NKT cells']].sum(axis=1)
        sensitivity=pd.DataFrame({'T_all':weights['T cells'],'T_without_cycling':noncycling,'Cycling_T':fine['Cycling T-cells'],'T_score':scores.T_score},index=weights.index)
        sensitivity.to_csv(out/'T_cycling_sensitivity_by_spot.csv')
        pd.DataFrame([dict(sample_id=sample,definition=definition,spearman=float(spearmanr(values,scores.T_score).statistic)) for definition,values in [('main_T_all',weights['T cells']),('diagnostic_without_cycling',noncycling),('Cycling_T_component',fine['Cycling T-cells'])]]).to_csv(out/'T_cycling_sensitivity_summary.csv',index=False)
        for t,group in [('T cells','T'),('CAFs','CAF'),('NK cells','NK')]:
            for score in ['T','CAF','NK']:
                value=float(spearmanr(weights[t],scores[score+'_score']).statistic)
                all_marker.append(dict(sample_id=sample,subtype=inventory.loc[sample,'subtype'],cell_type=t,marker_panel=score,spearman=value,n_spots=len(weights)))
        classified=obs.pathology_eligible.astype(bool)
        region=weights.loc[classified].assign(pathology=obs.loc[classified,'Classification']).groupby('pathology').median()
        region.insert(0,'n_spots',obs.loc[classified].groupby('Classification').size())
        region.to_csv(out/'pathology_region_medians.csv')
        review=obs.join(weights[['CAFs','T cells','NK cells']]).join(scores)
        for label in ['T cells CD4+','T cells CD8+','NKT cells','Cycling T-cells']:
            review[label+'_component_pct']=fine[label]*100
        review['high_T']=review.index.isin(high(weights['T cells']))
        review['high_CAF']=review.index.isin(high(weights.CAFs))
        review['high_T_low_T_expression']=review.high_T & review.T_score.le(review.T_score.quantile(.25))
        review['high_CAF_low_CAF_expression']=review.high_CAF & review.CAF_score.le(review.CAF_score.quantile(.25))
        review.to_csv(out/'spot_review.csv')
        review.loc[review.high_T_low_T_expression].to_csv(out/'T_expression_discordant_spots.csv')
        ctx=context(sample,obs)
        fig,axes=plt.subplots(2,3,figsize=(13,9),layout='constrained')
        caps={}
        for j,(t,panel) in enumerate([('CAFs','CAF'),('T cells','T'),('NK cells','NK')]):
            values=weights[t]*100;cap=float(values.quantile(.99));caps[t]=cap
            draw(fig,axes[0,j],ctx,values,f'{sample} | {t}','Relative weight (%) | P99 cap',vmax=cap)
            values=scores[panel+'_score'];cap=float(values.quantile(.99))
            draw(fig,axes[1,j],ctx,values,f'{panel}-associated expression','Mean log1p(CP10k) | P99 cap',vmax=cap)
        fig.savefig(out/'CAF_T_NK_expression.png',dpi=170);fig.savefig(out/'CAF_T_NK_expression.pdf');plt.close(fig)
        (out/'display_limits.json').write_text(json.dumps(caps,indent=2)+'\n')
        all_summary.append(dict(sample_id=sample,subtype=inventory.loc[sample,'subtype'],n_qc_retained=len(a),n_fitted=len(weights),n_model_removed=len(a)-len(weights),T_median_pct=weights['T cells'].median()*100,CAF_median_pct=weights.CAFs.median()*100,NK_median_pct=weights['NK cells'].median()*100,T_high_count=int(review.high_T.sum()),T_high_low_marker=int(review.high_T_low_T_expression.sum()),CAF_high_low_marker=int(review.high_CAF_low_CAF_expression.sum()),gene_count=int(inventory.loc[sample,'model_genes'])))
        if sample=='CID4535':
            old_raw=pd.read_csv(ROOT/'data/processed/CID4535_rctd_v1/rctd_weights_spot_by_celltype.csv',index_col=0).loc[weights.index]
            old=old_raw.div(old_raw.sum(axis=1),axis=0).T.groupby(MAPPING).sum().T
            status=pd.DataFrame(index=weights.index)
            legacy_bootstrap=pd.read_csv(ROOT/'results/Wu2021_coarse_T_v1/spatial_CID4535/spot_stability.csv',index_col=0)
            for t in weights.columns:
                e=(weights[t]-old[t])*100;b=set(high(old[t]));n=set(high(weights[t]))
                comparison.append(dict(cell_type=t,spearman=float(spearmanr(old[t],weights[t]).statistic),mean_abs_change_pp=e.abs().mean(),rmse_change_pp=np.sqrt((e**2).mean()),mean_change_pp=e.mean(),old_mean_pct=old[t].mean()*100,new_mean_pct=weights[t].mean()*100,top_k=len(b),top_overlap=len(b&n),top_retention=len(b&n)/len(b)))
                status[t+'_old_pct']=old[t]*100;status[t+'_new_pct']=weights[t]*100;status[t+'_change_pp']=e
                status[t+'_high_status']=['both_high' if s in b&n else 'old_only' if s in b else 'new_only' if s in n else 'other' for s in status.index]
                if t in ['CAFs','T cells']:
                    legacy=legacy_bootstrap[legacy_bootstrap.cell_type.eq(t)].loc[status.index]
                    status[t+'_old_reference_bootstrap_high_frequency']=legacy.top10_frequency
                    status[t+'_old_reference_bootstrap_sd_pp']=legacy.bootstrap_sd*100
            status.to_csv(out/'reference_comparison_by_spot.csv')
            pd.DataFrame(comparison).to_csv(out/'reference_comparison_metrics.csv',index=False)
            expression_comparison=[]
            for t,panel in [('CAFs','CAF'),('T cells','T'),('NK cells','NK')]:
                expression_comparison.append(dict(cell_type=t,marker_panel=panel,old_spearman=float(spearmanr(old[t],scores[panel+'_score']).statistic),new_spearman=float(spearmanr(weights[t],scores[panel+'_score']).statistic)))
            pd.DataFrame(expression_comparison).to_csv(out/'reference_expression_concordance.csv',index=False)
            fig,axes=plt.subplots(2,3,figsize=(13,9),layout='constrained')
            for i,t in enumerate(['CAFs','T cells']):
                cap=float(pd.concat([old[t],weights[t]]).quantile(.99)*100)
                draw(fig,axes[i,0],ctx,old[t]*100,f'{t} | old single donor','Relative weight (%)',vmax=cap)
                draw(fig,axes[i,1],ctx,weights[t]*100,f'{t} | new two donors','Relative weight (%)',vmax=cap)
                delta=(weights[t]-old[t])*100;limit=float(delta.abs().quantile(.99))
                draw(fig,axes[i,2],ctx,delta,'New minus old','Difference (pp) | P99 cap',vmin=-limit,vmax=limit,cmap='RdBu_r')
            fig.savefig(out/'old_new_CAF_T.png',dpi=180);fig.savefig(out/'old_new_CAF_T.pdf');plt.close(fig)
            # 空间高权重区域的参考一致性；颜色是集合成员关系，不是准确性标签。
            fig,axes=plt.subplots(1,2,figsize=(10,6),layout='constrained')
            palette={'other':'#d0d0d0','both_high':'#31688e','old_only':'#d98e3a','new_only':'#6dba73'}
            from matplotlib.lines import Line2D
            for ax,t in zip(axes,['CAFs','T cells']):
                image,xy,radius=ctx;ax.imshow(image,alpha=.32)
                p=PatchCollection([Circle((x,y),radius) for x,y in xy],linewidth=0,rasterized=True)
                p.set_facecolor(status[t+'_high_status'].map(palette));ax.add_collection(p)
                ax.set_xlim(xy[:,0].min()-50,xy[:,0].max()+50);ax.set_ylim(xy[:,1].max()+50,xy[:,1].min()-50);ax.set_aspect('equal');ax.axis('off');ax.set_title(t)
            fig.legend(handles=[Line2D([0],[0],marker='o',ls='',color=c,label=k) for k,c in palette.items()],loc='outside lower center',ncol=4,frameon=False)
            fig.savefig(out/'reference_high_spot_agreement.png',dpi=180);fig.savefig(out/'reference_high_spot_agreement.pdf');plt.close(fig)
        # 保存每张切片的简要记录。
        local_marker=pd.DataFrame(all_marker);local_marker=local_marker[local_marker.sample_id.eq(sample)]
        local_marker.to_csv(out/'marker_concordance.csv',index=False)
        lines=[f'# {sample} CAF/T空间分析','',f'作者亚型：{inventory.loc[sample,"subtype"]}。质控保留{len(a)}个spot，RCTD完成{len(weights)}个，模型额外排除{len(a)-len(weights)}个。',
            '参考为CID4471和CID4535，先按13类拟合，再将CD4、CD8、NKT、Cycling T相加为T cells；NK单列。',
            f'本片使用{inventory.loc[sample,"model_genes"]}个公共基因，模板中未按原名匹配的基因有{inventory.loc[sample,"omitted_genes"]}个。', '', '表达一致性：']
        for t,panel in [('T cells','T'),('CAFs','CAF'),('NK cells','NK')]:
            r=local_marker[local_marker.cell_type.eq(t)&local_marker.marker_panel.eq(panel)].iloc[0]
            lines.append(f'- {t}与对应表达面板的Spearman相关为{r.spearman:.3f}。')
        lines += ['',f'T细胞最高10%的{int(review.high_T.sum())}个spot中，{int(review.high_T_low_T_expression.sum())}个T表达分数不高于本片25%分位阈值，已在spot_review.csv标记；零值并列时低表达集合可超过25%。',
            '表达面板与反卷积使用同一份RNA数据，这是表达一致性核对，不是独立细胞身份验证。NK相关面板也可由细胞毒性T细胞表达。',
            '病理汇总排除Missing和Uncertain，但这些spot仍保留在表达与拟合中。空间图和逐spot表不用于CAF/T因果或空间排斥判断。']
        if inventory.loc[sample,'subtype']=='TNBC':lines+=['本片为TNBC，参考供者为ER+；本次按固定方案进行跨亚型应用检查，不作为同亚型独立验证。']
        if review.high_T_low_T_expression.any():
            flagged=review.loc[review.high_T_low_T_expression]
            zero=int(flagged.T_panel_umi.eq(0).sum())
            means=fine.loc[flagged.index,['T cells CD4+','T cells CD8+','NKT cells','Cycling T-cells']].mean()*100
            lines += [f'上述分歧spot中，{zero}个T面板UMI总数为零；平均预测T权重为{flagged["T cells"].mean()*100:.2f}%。T子类平均贡献为：'+ '，'.join(f'{t} {value:.2f}%' for t,value in means.items())+'。',
                'T_expression_discordant_spots.csv保留这些点的坐标、原病理标签、完整文库UMI和T子类分量，供进一步复核；未删除或重新标注。']
        lines += [f'新增敏感性诊断：仅汇总CD4、CD8、NKT而不纳入Cycling T时，与T面板的相关为{spearmanr(noncycling,scores.T_score).statistic:.3f}。这不是预先固定的主定义，不据此替换主结果，也不能把所有Cycling T判定为错误。']
        (out/'interpretation.md').write_text('\n'.join(lines)+'\n')
    pd.DataFrame(all_summary).to_csv(OUT/'completed_slice_summary.csv',index=False)
    pd.DataFrame(all_marker).to_csv(OUT/'all_marker_concordance.csv',index=False)
    if len(samples)>1:
        fig,axes=plt.subplots(2,3,figsize=(13,9),layout='constrained')
        caps=[]
        for sample in samples:
            w=pd.read_csv(OUT/sample/'coarse_relative_weights.csv',index_col=0);caps.append(w['T cells'].quantile(.99)*100)
        vmax=max(caps)
        for ax,sample in zip(axes.ravel(),samples):
            a=ad.read_h5ad(DATA/sample/'spatial_counts_full.h5ad');w=pd.read_csv(OUT/sample/'coarse_relative_weights.csv',index_col=0)
            draw(fig,ax,context(sample,a.obs.loc[w.index]),w['T cells']*100,f'{sample} | {inventory.loc[sample,"subtype"]}','T relative weight (%) | shared scale',vmax=vmax)
        for ax in axes.ravel()[len(samples):]:ax.axis('off')
        fig.savefig(OUT/'all_slices_T_atlas.png',dpi=180);fig.savefig(OUT/'all_slices_T_atlas.pdf');plt.close(fig)
    print(pd.DataFrame(all_summary).to_string(index=False))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('samples',nargs='*',default=PROTOCOL['samples']);args=parser.parse_args()
    main(args.samples)
