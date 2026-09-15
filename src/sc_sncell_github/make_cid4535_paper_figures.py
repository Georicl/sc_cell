"""从已有分析对象重新绘制README论文图"""
# 1. 读取已有结果，统一图形样式和来源记录
from pathlib import Path
import hashlib
import json

import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/CID4535_paper_v1"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7,
    "axes.titlesize": 8, "axes.labelsize": 7, "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": .6, "pdf.fonttype": 42, "svg.fonttype": "none"})
SC_PATH = ROOT / "data/processed/CID4535_baseline_v1/single_cell_analysis.h5ad"
SP_PATH = ROOT / "data/processed/CID4535_spatial/spatial_qc_filtering.h5ad"
W_PATH = ROOT / "data/processed/CID4535_rctd_v1/rctd_weights_spot_by_celltype.csv"
BOOT = ROOT / "results/CID4535_reference_bootstrap_v1"
sc = ad.read_h5ad(SC_PATH)
sp = ad.read_h5ad(SP_PATH)
raw = pd.read_csv(W_PATH, index_col=0).loc[sp.obs_names]
assert raw.index.is_unique and np.isfinite(raw.to_numpy()).all()
assert (raw.to_numpy() >= 0).all() and (raw.sum(axis=1) > 0).all()
weights = raw.div(raw.sum(axis=1), axis=0)
targets = ["CAFs", "T cells CD8+", "Cancer Epithelial", "PVL"]
short = {"CAFs": "CAF", "T cells CD8+": "CD8 T", "Cancer Epithelial": "Cancer", "PVL": "PVL"}
metrics = pd.read_csv(BOOT / "metrics_by_run.csv")
overlap = pd.read_csv(BOOT / "high_spot_overlap.csv")
spot_stats = pd.read_csv(BOOT / "spot_stability.csv")
spatial = sp.uns["spatial"]["CID4535"]
image = spatial["images"]["hires"]
scale = spatial["scalefactors"]["tissue_hires_scalef"]
xy = sp.obsm["spatial"] * scale
radius = spatial["scalefactors"]["spot_diameter_fullres"] * scale / 2
extent = (xy[:, 0].min()-70, xy[:, 0].max()+70, xy[:, 1].max()+70, xy[:, 1].min()-70)
major = sc.obs["celltype_major"].cat.categories.tolist()
palette = dict(zip(major, sc.uns["celltype_major_colors"]))

def panel(ax, letter, title):
    ax.set_title(title, loc="left", pad=7)
    ax.text(-.07, 1.055, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")

def tissue(ax, alpha=.27):
    ax.imshow(image, alpha=alpha)
    ax.set_xlim(extent[:2]); ax.set_ylim(extent[2:])
    ax.set_aspect("equal"); ax.set_axis_off()

def spot_collection(ax, positions, values=None, colors=None, cmap="magma", vmax=None):
    patches = [Circle((x, y), radius) for x, y in positions]
    coll = PatchCollection(patches, linewidths=0, rasterized=True)
    if colors is not None:
        coll.set_facecolor(colors)
    else:
        coll.set_array(np.asarray(values)); coll.set_cmap(cmap); coll.set_clim(0, vmax)
    ax.add_collection(coll)
    return coll

def spatial_panel(fig, rect, values, letter, title, vmax=None, unit="Relative weight", cmap="magma"):
    ax = fig.add_axes(rect); tissue(ax)
    coll = spot_collection(ax, xy, values=values, vmax=vmax, cmap=cmap)
    panel(ax, letter, title)
    box = ax.get_position()
    cax = fig.add_axes([box.x0, box.y0-.034, box.width, .011])
    cb = fig.colorbar(coll, cax=cax, orientation="horizontal")
    cb.ax.tick_params(labelsize=6, pad=1, length=2); cb.set_label(unit, fontsize=6, labelpad=1)
    return ax

def save(fig, name):
    # 使用各图设定的固定画布，不通过tight裁剪改变实际纸面尺寸
    fig.savefig(OUT / f"{name}.png", dpi=300, facecolor="white")
    fig.savefig(OUT / f"{name}.pdf", facecolor="white")
    plt.close(fig)

def canvas(title):
    fig = plt.figure(figsize=(210/25.4, 148.5/25.4))
    fig.text(.045, .975, title, va="top", fontsize=10, fontweight="bold")
    return fig

# 2. Figure 1：单细胞图谱、细胞组成及作者大类marker
fig = canvas("CID4535 | Single-cell reference")
ax = fig.add_axes([.07,.53,.34,.37])
umap = sc.obsm["X_umap"]
for cell_type in major:
    mask = sc.obs["celltype_major"].eq(cell_type).to_numpy()
    ax.scatter(umap[mask,0],umap[mask,1],s=1.8,c=palette[cell_type],linewidths=0,rasterized=True)
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2"); ax.set_xticks([]); ax.set_yticks([])
panel(ax,"A","Author annotations | 3,961 cells")
ax = fig.add_axes([.64,.54,.30,.35])
composition = sc.obs["celltype_major"].value_counts().reindex(major)
ax.barh(range(len(major)), composition, color=[palette[t] for t in major],height=.68)
ax.set_yticks(range(len(major)),major,fontsize=6.5); ax.invert_yaxis()
ax.set_xlim(0,2600); ax.set_xlabel("Number of cells")
for y,v in enumerate(composition): ax.text(v+25,y,f"{v:,}",va="center",fontsize=6.5)
panel(ax,"B","Captured cell composition")
composition.rename("n_cells").to_csv(OUT / "figure1_cell_counts.csv")

genes = ["EPCAM","KRT19","COL1A1","DCN","RGS5","PDGFRB","PECAM1","VWF","CD3D","CD8A","LYZ","CD79A","JCHAIN"]
x = sc[:,genes].X.toarray()
means = pd.DataFrame(x,index=sc.obs_names,columns=genes).groupby(sc.obs["celltype_major"],observed=True).mean().reindex(major)
fracs = pd.DataFrame(x>0,index=sc.obs_names,columns=genes).groupby(sc.obs["celltype_major"],observed=True).mean().reindex(major)
scaled = (means-means.min())/(means.max()-means.min())
means.to_csv(OUT / "figure1_marker_log_means.csv"); fracs.to_csv(OUT / "figure1_marker_detection_fraction.csv")
ax = fig.add_axes([.22,.16,.62,.245])
gx,gy=np.meshgrid(np.arange(len(genes)),np.arange(len(major)))
dots=ax.scatter(gx.ravel(),gy.ravel(),s=fracs.to_numpy().ravel()*30,c=scaled.to_numpy().ravel(),cmap="Blues",vmin=0,vmax=1,linewidths=.15,edgecolors="#526575")
ax.set_xticks(range(len(genes)),genes,rotation=50,ha="right",fontsize=6.5)
ax.set_yticks(range(len(major)),major,fontsize=6.5); ax.invert_yaxis()
ax.set_xlim(-.6,len(genes)-.4); ax.set_ylim(len(major)-.5,-.5)
ax.spines[["left","bottom"]].set_visible(False);ax.tick_params(length=0)
panel(ax,"C","Marker expression by author major type")
cax=fig.add_axes([.89,.23,.012,.15]);cb=fig.colorbar(dots,cax=cax);cb.set_label("Scaled mean",fontsize=6)
handles=[plt.scatter([],[],s=v*30,c="#8395a7") for v in [.25,.5,1]]
fig.legend(handles,["25%","50%","100%"],title="Cells expressing",loc="lower right",bbox_to_anchor=(.99,.055),frameon=False,fontsize=6,title_fontsize=6)
save(fig,"figure1_single_cell")

# 3. Figure 2：空间对齐、过滤记录与QC分布
fig=canvas("CID4535 | Spatial registration and quality control")
fig.set_size_inches(210/25.4, 110/25.4)
rects=[[.04,.29,.205,.60],[.28,.29,.205,.60],[.52,.29,.205,.60],[.76,.29,.205,.60]]
ax=fig.add_axes(rects[0]);tissue(ax,alpha=.65)
spot_collection(ax,xy,colors="#98a4b0")
excluded=pd.read_csv(ROOT / "data/processed/CID4535_spatial/excluded_spots_qc.csv",index_col=0)
pos=pd.read_csv(ROOT / "data/raw/wu2021_visium/spatial/CID4535_spatial/tissue_positions_list.csv",header=None,index_col=0)
exxy=pos.loc[excluded.index,[5,4]].to_numpy()*scale
ex=ax.scatter(exxy[:,0],exxy[:,1],s=6,facecolors="none",edgecolors="#b93637",linewidths=.5)
panel(ax,"A","Spot review")
fig.text(.043,.22,"1,127 input spots\n25 excluded; 1,102 retained",fontsize=7,va="top")
pathology=sp.obs["Classification"].astype("string").fillna("Missing")
path_names=["Invasive cancer","Invasive cancer + lymphocytes","Stroma","Lymphocytes","Uncertain","Adipose tissue","Invasive cancer + adipose tissue + lymphocytes","Missing"]
path_colors=["#c76e55","#dca89a","#588f83","#517cb1","#c6b86a","#ac94b3","#7c5165","#c9cdd1"]
pathmap=dict(zip(path_names,path_colors))
ax=fig.add_axes(rects[1]);tissue(ax);spot_collection(ax,xy,colors=pathology.map(pathmap).tolist());panel(ax,"B","Author pathology")
spatial_panel(fig,rects[2],sp.obs.total_counts/1000,"C","UMI counts",vmax=sp.obs.total_counts.max()/1000,unit="UMIs (×1,000)",cmap="viridis")
spatial_panel(fig,rects[3],sp.obs.pct_counts_mt,"D","Mitochondrial RNA",vmax=sp.obs.pct_counts_mt.max(),unit="Percent of counts",cmap="viridis")
handles=[Line2D([0],[0],marker="o",color="none",markerfacecolor=pathmap[t],markeredgecolor="none",markersize=4) for t in path_names]
leg_labels=["Invasive cancer","Cancer + lymphocytes","Stroma","Lymphocytes","Uncertain","Adipose tissue","Cancer + adipose + lymphocytes","Missing"]
fig.legend(handles,leg_labels,loc="lower center",bbox_to_anchor=(.5,.075),ncol=4,frameon=False,fontsize=6.3,columnspacing=1.6)
fig.text(.045,.04,"Grey spots: retained   |   Red rings: excluded   |   Coordinates and H&E: Wu et al., 2021",fontsize=6.5,color="#505963")
pathology.value_counts().rename("n_spots").to_csv(OUT / "figure2_pathology_counts.csv")
sp.obs[["total_counts","n_genes_by_counts","pct_counts_mt"]].describe().to_csv(OUT / "figure2_qc_summary.csv")
save(fig,"figure2_spatial_qc")

# 4. Figure 3：反卷积空间分布、区域中位数与marker表达对应
fig=canvas("CID4535 | Spatial deconvolution")
limits={}
for j,t in enumerate(targets):
    upper=float(weights[t].quantile(.99));limits[t]=upper
    spatial_panel(fig,[.04+j*.24,.47,.205,.43],weights[t],"ABCD"[j],t,vmax=upper,unit="Relative weight (P99 cap)")
regions=["Invasive cancer","Stroma","Lymphocytes"]
region=weights.assign(pathology=pathology).groupby("pathology")[targets].median()
region.to_csv(OUT / "figure3_pathology_medians.csv")
ax=fig.add_axes([.15,.13,.35,.19]);h=region.loc[regions].to_numpy()*100
im=ax.imshow(h,cmap="Blues",vmin=0,vmax=100,aspect="auto")
ax.set_yticks(range(3),["Invasive cancer (418)","Stroma (169)","Lymphocytes (69)"],fontsize=6.5)
ax.set_xticks(range(4),[short[t] for t in targets]);ax.tick_params(length=0)
for i in range(3):
    for j in range(4):ax.text(j,i,f"{h[i,j]:.2f}",ha="center",va="center",color="white" if h[i,j]>65 else "#243846",fontsize=7)
panel(ax,"E","Regional median relative weight (%)")
cax=fig.add_axes([.515,.13,.012,.19]);fig.colorbar(im,cax=cax).ax.tick_params(labelsize=6)
marker_sets={"CAFs":["COL1A1","COL1A2","DCN","LUM"],"T cells CD8+":["CD3D","CD3E","TRAC","CD8A","CD8B"],"Cancer Epithelial":["EPCAM","KRT8","KRT18","KRT19"],"PVL":["RGS5","PDGFRB","MCAM","CSPG4"]}
total=np.asarray(sp.X.sum(axis=1)).ravel()
marker_rows=[]
for t,gs in marker_sets.items():
    e=np.log1p(sp[:,gs].X.toarray()/total[:,None]*1e4).mean(axis=1)
    marker_rows.append({"cell_type":t,"markers":";".join(gs),"spearman":spearmanr(weights[t],e).statistic})
marker_table=pd.DataFrame(marker_rows).set_index("cell_type");marker_table.to_csv(OUT / "figure3_marker_correlations.csv")
ax=fig.add_axes([.67,.13,.27,.19]);v=marker_table.loc[targets,"spearman"]
ax.hlines(range(4),0,v,color="#b6c6d0",lw=2);ax.scatter(v,range(4),color="#356886",s=20)
ax.set_yticks(range(4),[short[t] for t in targets]);ax.invert_yaxis();ax.set_xlim(0,1.05);ax.set_xlabel("Spearman correlation")
for y,n in enumerate(v):ax.text(n+.025,y,f"{n:.3f}",va="center",fontsize=6.5)
panel(ax,"F","Marker-expression concordance")
save(fig,"figure3_deconvolution")

# 5. Figure 4：重抽样指标和CAF/CD8的高值位置稳定性
fig=canvas("CID4535 | Reference bootstrap stability (10 runs)")
top10=overlap[overlap.top_fraction.eq(.1)]
metric_panels=[(metrics,"spearman","Rank concordance","Spearman vs baseline"),(metrics,"mae_pp","Weight difference","MAE (percentage points)"),(top10,"retention","High-spot retention","Fraction retained")]
order=["CAFs","T cells CD8+","PVL","Cancer Epithelial"]
for j,(table,key,title,ylabel) in enumerate(metric_panels):
    ax=fig.add_axes([.09+j*.305,.68,.235,.21])
    for i,t in enumerate(order):
        val=table.loc[table.cell_type.eq(t),key].to_numpy()
        ax.scatter(i+np.linspace(-.11,.11,len(val)),val,s=9,c="#356886",alpha=.8)
        ax.plot([i-.18,i+.18],[np.median(val)]*2,color="#b76e35",lw=1.4)
    ax.set_xticks(range(4),[short[t] for t in order]);ax.set_ylabel(ylabel,fontsize=6.5)
    ax.set_ylim(0,1.04 if key!="mae_pp" else 2.9);panel(ax,"ABC"[j],title)
for j,t in enumerate(["CAFs","T cells CD8+"]):
    values=spot_stats[spot_stats.cell_type.eq(t)].set_index("barcode").loc[sp.obs_names]
    spatial_panel(fig,[.07+j*.315,.15,.245,.40],values.top10_frequency,"DE"[j],f"{short[t]} high-spot frequency",vmax=1,unit="Fraction of 10 runs",cmap="viridis")
v=spot_stats[spot_stats.cell_type.eq("T cells CD8+")].set_index("barcode").loc[sp.obs_names]
spatial_panel(fig,[.70,.15,.245,.40],v.bootstrap_sd*100,"F","CD8 weight variability",vmax=v.bootstrap_sd.max()*100,unit="SD (percentage points)",cmap="viridis")
fig.text(.065,.025,"High spots: top 10% per run (111 spots). Orange lines in A–C: median across 10 runs.",fontsize=6.5)
save(fig,"figure4_bootstrap")

# 6. 保存正文数值及成图来源，供README逐项核对
summary={"n_sc":sc.n_obs,"n_sc_genes":sc.n_vars,"n_spots":sp.n_obs,"n_spatial_genes":sp.n_vars,
    "n_leiden":sc.obs.leiden.nunique(),"relative_weight_mean":weights.mean().to_dict(),
    "relative_weight_p99":limits,"marker_correlations":marker_table.spearman.to_dict(),
    "weight_sum_range":[float(raw.sum(axis=1).min()),float(raw.sum(axis=1).max())],
    "caf_cd8_spearman":float(spearmanr(weights.CAFs,weights["T cells CD8+"]).statistic),
    "sources":[{"path":str(p.relative_to(ROOT)),"sha256":hashlib.sha256(p.read_bytes()).hexdigest()} for p in [SC_PATH,SP_PATH,W_PATH,BOOT/"metrics_by_run.csv",BOOT/"spot_stability.csv"]],
    "figure_dimensions_mm":{"figure1":[210,148.5],"figure2":[210,110],"figure3":[210,148.5],"figure4":[210,148.5]},"png_dpi":300,"expression_panel_normalization":"log1p(counts / total_counts * 10000), mean across listed genes"}
(OUT / "figure_data_manifest.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
print(json.dumps({k:v for k,v in summary.items() if k!="sources"},ensure_ascii=False,indent=2))
