"""Redraw the stage-report figures from saved results; no model is rerun."""
from pathlib import Path
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/Wu2021_stage_paper_v1'
OUT.mkdir(parents=True, exist_ok=True)
SOURCES = {}
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'savefig.dpi': 300})

def read(path):
    p = ROOT / 'results' / path
    SOURCES[str(p.relative_to(ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
    return pd.read_csv(p)

def save(fig, name):
    fig.savefig(OUT / f'{name}.png', bbox_inches='tight', facecolor='white')
    fig.savefig(OUT / f'{name}.pdf', bbox_inches='tight', facecolor='white')
    plt.close(fig)

def title(ax, text):
    ax.set_title(text, loc='left', weight='bold', fontsize=9, pad=10)

# 1. Reference composition retains the author parent labels.
inv = read('Wu2021_reference_comparison_v1/reference_inventory.csv')
ref = inv.query("variant == 'epithelial_split'").groupby(
    ['parent_label', 'donor_id']).n_cells.sum().unstack(fill_value=0)
assert ref.values.sum() == 12570
ref.to_csv(OUT / 'figure1_reference_counts.csv')
boot = read('CID4535_reference_bootstrap_v1/stability_summary.csv')
fig = plt.figure(figsize=(9, 6.5), layout='constrained')
gs = fig.add_gridspec(2, 2, width_ratios=[1.45, 1], height_ratios=[.55, 1])
ax = fig.add_subplot(gs[0, :]); ax.axis('off')
title(ax, 'A  Study sequence and donor roles')
steps = [('Baseline', 'CID4535\n3,961 cells; 1,102 spots\n10 reference resamples'),
         ('Reference tests', 'CID4471 + CID4535\n12,570 reference cells\nCID4067: 488 mixtures'),
         ('Donor / label tests', 'CID4463: 474 mixtures\nCID3941 / 3948 / 4461\n540 T-gradient mixtures'),
         ('Spatial application', '2 ER + 4 TNBC slices\n15,396 fitted spots\nCAF and aggregate T')]
for i, (heading, body) in enumerate(steps):
    x = .125 + i*.25
    ax.text(x, .52, heading+'\n\n'+body, ha='center', va='center', fontsize=8,
            bbox=dict(boxstyle='round,pad=.65', fc='#edf3f7', ec='#c6d6e1'))
    if i < 3: ax.annotate('', xy=(x+.145,.52), xytext=(x+.105,.52), arrowprops=dict(arrowstyle='->', color='#61717d'))
ax = fig.add_subplot(gs[1,0]); ref.plot.barh(stacked=True, ax=ax, color=['#4682a9','#e5a15c'], width=.75)
title(ax, 'B  Two-donor reference'); ax.set_xlabel('Cells'); ax.set_ylabel(''); ax.legend(frameon=False, fontsize=7)
ax = fig.add_subplot(gs[1,1]); selected = boot.set_index('cell_type').loc[['CAFs','T cells CD8+','Cancer Epithelial','PVL']]
ax.barh(['CAF','CD8 T','Cancer epithelial','PVL'], selected.spearman_median, color='#4682a9')
for i,v in enumerate(selected.spearman_median): ax.text(v-.02,i,f'{v:.3f}',va='center',ha='right',color='white')
ax.set_xlim(0,1.04); ax.set_xlabel('Median Spearman correlation'); title(ax,'C  CID4535 single-donor resampling')
save(fig,'figure1_study_reference')

# 2. Keep donor and depth separate, use one common error scale.
dev = read('Wu2021_reference_comparison_v1/key_metrics.csv'); val = read('Wu2021_independent_validation_v1/validation_metrics.csv')
variants = ['baseline','donor_balanced','epithelial_split','balanced_split']
types = ['CAFs','T cells CD8+','Cancer Epithelial','Normal Epithelial']
fig, axes = plt.subplots(2,2,figsize=(9,5.7),layout='constrained')
exports=[]
for j,(donor,data) in enumerate([('CID4067',dev),('CID4463',val)]):
    for i,depth in enumerate([1.,.25]):
        d=data.query("method == 'RCTD' and truth == 'rna' and scenario == 'factorial' and depth_fraction == @depth")
        matrix=d.pivot(index='variant',columns='cell_type',values='mae_pp').loc[variants,types]
        ax=axes[i,j]; im=ax.imshow(matrix, cmap='YlOrRd',vmin=0,vmax=37,aspect='auto')
        ax.set_xticks(range(4),['CAF','CD8 T','Cancer epi.','Normal epi.']); ax.set_yticks(range(4),['Baseline','Balanced','Epi. split','Both'])
        for r in range(4):
            for c in range(4): ax.text(c,r,f'{matrix.iloc[r,c]:.2f}',ha='center',va='center',color='white' if matrix.iloc[r,c]>22 else '#222222')
        title(ax,f'{"ABCD"[i*2+j]}  {donor} | depth {depth:g}')
        exports.append(d.assign(report_donor=donor))
fig.colorbar(im,ax=axes,shrink=.7,label='RNA-fraction MAE (percentage points)')
pd.concat(exports).to_csv(OUT/'figure2_reference_metrics.csv',index=False)
save(fig,'figure2_reference_tests')

# 3. T-specific tests and the choice of post-fit aggregation.
g=read('Wu2021_tcell_validation_v1/gradient_metrics.csv').query("variant == 'baseline' and method == 'RCTD' and truth == 'rna' and depth_fraction == 1 and cell_type == 'T cells CD8+'").copy()
t=read('Wu2021_coarse_T_v1/T_cell_key_metrics.csv').query("method == 'RCTD' and truth == 'rna' and depth_fraction == 1").copy()
g.to_csv(OUT/'figure3_cd8_gradient_metrics.csv',index=False);t.to_csv(OUT/'figure3_aggregate_T_metrics.csv',index=False)
fig, axes=plt.subplots(1,3,figsize=(10,3.6),layout='constrained')
for donor,color in zip(['CID3941','CID3948','CID4461'],['#4682a9','#d77941','#65874e']):
    d=g.query('donor_id == @donor').sort_values('n_caf')
    axes[0].plot(d.n_caf,d.mae_pp,'o-',label=donor,color=color)
    axes[1].plot(d.n_caf,d.absent_mean_prediction_pct,'o-',color=color)
for ax in axes[:2]: ax.set_xticks([0,1,3]);ax.set_xlabel('CAF cells added to 10 T cells');ax.set_ylim(bottom=0)
title(axes[0],'A  CD8 gradient error');axes[0].set_ylabel('RNA-fraction MAE (pp)');axes[0].legend(frameon=False,fontsize=7)
title(axes[1],'B  Predicted CD8 when absent');axes[1].set_ylabel('Mean predicted CD8 weight (%)')
donors=['CID4067','CID4463','CID3941','CID3948','CID4461']; ax=axes[2]
for approach,offset,color,label in [('fine_fit_then_sum',-.18,'#4682a9','Fine fit, then sum'),('coarse_reference_refit',.18,'#e5a15c','Coarse reference refit')]:
    d=t.query('approach == @approach').set_index('donor_id').loc[donors]
    ax.barh(np.arange(5)+offset,d.mae_pp,height=.34,color=color,label=label)
ax.set_yticks(range(5),donors);ax.invert_yaxis();ax.axhline(1.5,color='#aaaaaa',ls=':',lw=1)
ax.set_xlabel('Aggregate T MAE (pp)');ax.legend(frameon=False,fontsize=7,loc='upper right');title(ax,'C  T-label aggregation')
save(fig,'figure3_T_label_tests')

# 4. A common scale per cell type permits between-slice visual comparison.
summary=read('Wu2021_multislice_CAF_T_v1/completed_slice_summary.csv')
corr=read('Wu2021_multislice_CAF_T_v1/all_marker_concordance.csv')
samples=summary.sample_id.tolist(); spots={s:read(f'Wu2021_multislice_CAF_T_v1/{s}/spot_review.csv') for s in samples}
fig=plt.figure(figsize=(12,7.8),layout='constrained');gs=fig.add_gridspec(3,6,height_ratios=[1,1,.9])
limits={}
for row,(ct,cmap,label) in enumerate([('CAFs','YlOrBr','CAF'),('T cells','viridis','T')]):
    vmax=float(np.quantile(np.concatenate([d[ct].to_numpy()*100 for d in spots.values()]),.99));limits[ct]=vmax
    row_axes=[]
    for col,s in enumerate(samples):
        ax=fig.add_subplot(gs[row,col]);row_axes.append(ax);d=spots[s]
        im=ax.scatter(d.pxl_col_in_fullres,d.pxl_row_in_fullres,c=d[ct]*100,s=2.4,cmap=cmap,vmin=0,vmax=vmax,rasterized=True,linewidths=0)
        ax.set_aspect('equal');ax.invert_yaxis();ax.axis('off')
        title(ax,f'{chr(65+row*6+col)}  {s}\n{summary.set_index("sample_id").loc[s,"subtype"]} | {label}')
    fig.colorbar(im,ax=row_axes,shrink=.65,label=f'{label} relative weight (%)')
ax=fig.add_subplot(gs[2,:3]);x=np.arange(6)
for ct,panel,offset,color in [('CAFs','CAF',-.16,'#b57b35'),('T cells','T',.16,'#397e86')]:
    d=corr.query('cell_type == @ct and marker_panel == @panel').set_index('sample_id').loc[samples]
    ax.bar(x+offset,d.spearman,width=.3,label=panel,color=color)
ax.set_xticks(x,samples,rotation=25,ha='right');ax.set_ylim(0,1);ax.set_ylabel('Spearman correlation');ax.legend(frameon=False);title(ax,'M  Concordance with marker expression')
ax=fig.add_subplot(gs[2,3:]);v=summary.T_high_low_marker/summary.T_high_count*100
ax.bar(x,v,color='#758694');ax.set_xticks(x,samples,rotation=25,ha='right');ax.set_ylim(0,85);ax.set_ylabel('High-T spots with low T score (%)')
for i,(a,b) in enumerate(zip(summary.T_high_low_marker,summary.T_high_count)):ax.text(i,v.iloc[i]+2,f'{a}/{b}',ha='center',fontsize=7)
title(ax,'N  Expression support in high-T spots')
summary.to_csv(OUT/'figure4_slice_summary.csv',index=False);corr.to_csv(OUT/'figure4_marker_correlations.csv',index=False)
save(fig,'figure4_spatial_atlas')
(OUT/'manifest.json').write_text(json.dumps({'script':str(Path(__file__).relative_to(ROOT)),'source_sha256':SOURCES,'spatial_shared_p99_percent':limits,'figure2_3_truth':'captured RNA fraction; errors in percentage points','figure3_depth_fraction':1,'models_rerun':False},indent=2)+'\n')
print(OUT)
