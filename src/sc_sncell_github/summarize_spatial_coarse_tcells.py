"""将CID4535既有细分拟合汇总为T细胞大类，重绘分布与参考重抽样稳定性。"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/Wu2021_coarse_T_v1/spatial_CID4535'
OUT.mkdir(parents=True,exist_ok=True)
mapping=pd.read_csv(ROOT/'results/Wu2021_coarse_T_v1/label_mapping.csv',index_col=0).reference_label
source=ROOT/'data/processed/CID4535_rctd_v1/rctd_weights_spot_by_celltype.csv'
spatial_path=ROOT/'data/processed/CID4535_spatial/spatial_qc_filtering.h5ad'
a=ad.read_h5ad(spatial_path)
def read_coarse(path):
    raw=pd.read_csv(path,index_col=0).loc[a.obs_names]
    assert set(raw.columns)==set(mapping.index) and np.isfinite(raw.values).all() and (raw.values>=0).all() and (raw.sum(axis=1)>0).all()
    fine=raw.div(raw.sum(axis=1),axis=0)
    coarse=fine.T.groupby(mapping).sum().T
    assert np.allclose(coarse.sum(axis=1),1)
    return coarse
weights=read_coarse(source);weights.to_csv(OUT/'coarse_relative_weights.csv')
topk=int(np.ceil(len(weights)*.1));targets=['CAFs','T cells']
runs=[];arrays=[];paths=sorted((ROOT/'data/processed/CID4535_reference_bootstrap_v1').glob('seed_*/weights.csv'))
assert len(paths)==10
for path in paths:
    w=read_coarse(path);arrays.append(w[targets].to_numpy())
    for t in targets:
        b=set(weights[t].sort_values(ascending=False,kind='stable').head(topk).index)
        h=set(w[t].sort_values(ascending=False,kind='stable').head(topk).index)
        runs.append(dict(run=path.parent.name,cell_type=t,spearman=float(spearmanr(weights[t],w[t]).statistic),top_k=topk,top_retention=len(b&h)/topk))
array=np.stack(arrays);metrics=pd.DataFrame(runs);metrics.to_csv(OUT/'bootstrap_metrics.csv',index=False)
spot_rows=[]
for j,t in enumerate(targets):
    frequency=np.zeros(len(weights))
    for replicate in array[:,:,j]:frequency[np.argsort(-replicate,kind='stable')[:topk]]+=1/len(paths)
    frame=pd.DataFrame(dict(cell_type=t,weight=weights[t],bootstrap_mean=array[:,:,j].mean(axis=0),bootstrap_sd=array[:,:,j].std(axis=0,ddof=1),top10_frequency=frequency),index=weights.index)
    spot_rows.append(frame)
stats=pd.concat(spot_rows);stats.to_csv(OUT/'spot_stability.csv',index_label='spot_id')
image=a.uns['spatial']['CID4535']['images']['hires'];scale=a.uns['spatial']['CID4535']['scalefactors']['tissue_hires_scalef']
xy=a.obsm['spatial']*scale;radius=a.uns['spatial']['CID4535']['scalefactors']['spot_diameter_fullres']*scale/2
fig,axes=plt.subplots(2,2,figsize=(10,9),layout='constrained')
caps={}
for j,t in enumerate(targets):
    for i in [0,1]:
        ax=axes[i,j];ax.imshow(image,alpha=.35)
        collection=PatchCollection([Circle((x,y),radius) for x,y in xy],linewidth=0,cmap='viridis')
        if i==0:
            vals=weights[t]*100;cap=float(vals.quantile(.99));caps[t]=cap
            collection.set_array(vals.to_numpy());collection.set_clim(0,cap)
            ax.set_title('CAF' if t=='CAFs' else 'T cells (CD4/CD8/NKT/cycling)')
            label='Relative weight (%) | P99 display cap'
        else:
            vals=stats[stats.cell_type.eq(t)].loc[weights.index,'top10_frequency']
            collection.set_array(vals.to_numpy());collection.set_clim(0,1)
            ax.set_title('Reference-bootstrap high-spot frequency');label='Fraction of 10 bootstrap runs'
        ax.add_collection(collection);ax.set_xlim(xy[:,0].min()-70,xy[:,0].max()+70);ax.set_ylim(xy[:,1].max()+70,xy[:,1].min()-70)
        ax.set_aspect('equal');ax.axis('off');fig.colorbar(collection,ax=ax,orientation='horizontal',label=label,shrink=.8)
fig.savefig(OUT/'CAF_T_spatial.png',dpi=220);fig.savefig(OUT/'CAF_T_spatial.pdf');plt.close(fig)
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()
provenance=dict(source_fit_reference='Original CID4535-only reference; existing 13-type fit, not the two-donor reference used in the coarse benchmark.',operation='Normalize existing raw weights by spot, then sum CD4/CD8/NKT/cycling. No new spatial fitting.',n_spots=len(weights),n_bootstrap=len(paths),top_k=topk,display_caps_percent=caps,source_hashes={str(p.relative_to(ROOT)):sha(p) for p in [source,spatial_path]+paths})
(OUT/'provenance.json').write_text(json.dumps(provenance,ensure_ascii=False,indent=2)+'\n')
summary=metrics.groupby('cell_type')[['spearman','top_retention']].mean()
lines=['# CID4535的CAF与T细胞大类分布','',
    '空间图将原CID4535单供者参考的13类拟合结果相加，沿用1102个保留spot；没有重新拟合。它用于与原图逐点对照，和本轮两供者参考验证不是同一套拟合。',
    'T cells由CD4、CD8、NKT和Cycling T组成，NK单独保留。上排为相对RNA权重，颜色上限使用各类型99%分位数；CSV保存未截断数值。下排为既有10次参考重抽样中进入该类型最高10% spot的频率。','']
for t,row in summary.iterrows():lines.append(f'- {t}：平均空间排序相关{row.spearman:.3f}，前{topk}个高权重spot平均保留{row.top_retention:.1%}。')
lines+=['','重抽样描述参考扰动下的稳定性，不是空间真值准确性。本图不用于CD8特异分布或CAF与CD8空间排斥的判断。']
(OUT/'interpretation.md').write_text('\n'.join(lines)+'\n')
print(summary)
