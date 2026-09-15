"""固定三位保留供者的CD4/CD8梯度和CAF背景配方。"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.io import mmwrite
import extract_er_qc_donors as extraction

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data/processed/Wu2021_tcell_validation_v1'
OUT=ROOT/'results/Wu2021_tcell_validation_v1'
DONORS=['CID3941','CID3948','CID4461']
VARIANTS=['baseline','donor_balanced','epithelial_split','balanced_split']
DATA.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

old=ROOT/'results/Wu2021_independent_validation_v1'
previous=json.loads((old/'frozen_reference_hashes.json').read_text())
for name in ['data/provenance/wu2021_supplementary_tables.xlsx','data/raw/wu2021_scrna/Wu_etal_2021_BRCA_scRNASeq/metadata.csv']:
    assert sha(ROOT/name)==previous[name]
candidates=pd.read_csv(old/'donor_candidates.csv').set_index('donor_id').loc[DONORS]
assert candidates.role.eq('reserve_validation').all()
assert candidates.clinical_subtype.eq('ER+').all() and candidates.treatment_status.eq('Naïve').all()
clinical=pd.read_csv(old/'clinical_source_records.csv').set_index('case_id')
for donor in DONORS:
    row=clinical.loc[int(donor[3:])]
    assert row.clinical_subtype=='ER+' and row.histology=='IDC' and row.treatment_status=='Naïve'
current=candidates.copy()
current['previous_role']=current.role
current['role']='tcell_validation'
current.to_csv(OUT/'validation_donors.csv')
protocol=dict(frozen_at_utc=datetime.now(timezone.utc).isoformat(),validation_donors=DONORS,
    reference_donors=['CID4471','CID4535'],previously_inspected=['CID4067','CID4463'],
    variants=VARIANTS,seeds=dict(zip(DONORS,[20260918,20260919,20260920])),rctd_seed=9,
    gradient=dict(t_cells=10,cd8_cells=[0,2,5,8,10],caf_cells=[0,1,3],repeats=5),
    pure=dict(types=['CAFs','T cells CD4+','T cells CD8+'],cells=5,repeats=5),
    depths=[1.,.25],n_base_mixtures_per_donor=90,n_spots_per_donor=180,
    pairing='Each CD4/CD8 draw is shared across CAF backgrounds; CAF 1 uses first cell of the same CAF 3 draw. Thinning is per component, paired across CAF backgrounds and depths.',
    sampling='No repeated cells within spot; reuse across spots; report distinct unordered source-cell sets.',
    cell_policy='Retain author cells; verify raw counts and detected genes; no new QC exclusion.',
    reference_policy='Freeze all reference counts, labels and NNLS profiles; no fitting-result-based tuning.',
    rctd=dict(mode='full',max_cores=2,ref_UMI_min=100,ref_n_cells_min=25,platform_calibration='Separate run per validation donor using its unlabeled synthetic counts.'),
    metrics=['RNA-fraction MAE','RMSE','bias','cell-fraction MAE','CD8-zero mean prediction','CD8-zero weight >=1%','CD8-zero weight >=5%','pure dominant matches','paired CAF effect'],
    thresholds=dict(primary_positive_weight=.01,secondary_positive_weight=.05),
    primary_truth='captured RNA after thinning',secondary_truth='cell counts',
    pooling='Report each donor separately; do not treat simulated spots as independent donors.')
(OUT/'protocol.json').write_text(json.dumps(protocol,ensure_ascii=False,indent=2)+'\n')
frozen={}
for path in previous:
    frozen[path]=sha(ROOT/path)
    assert frozen[path]==previous[path]
frozen[str((OUT/'protocol.json').relative_to(ROOT))]=sha(OUT/'protocol.json')
(OUT/'frozen_inputs.json').write_text(json.dumps(frozen,indent=2)+'\n')
extraction.DONORS=DONORS;extraction.OUTPUT=DATA
if not all((DATA/f'{d}_counts_before_qc.h5ad').exists() for d in DONORS):extraction.extract_donors()
types=sorted(pd.read_csv(ROOT/'data/processed/Wu2021_known_composition_v1/reference_metadata.csv').reference_label.unique())
all_qc=[];source_hashes={}
for donor in DONORS:
    inp=DATA/donor;out=OUT/donor;inp.mkdir(exist_ok=True);out.mkdir(exist_ok=True)
    a=ad.read_h5ad(DATA/f'{donor}_counts_before_qc.h5ad')
    a.obs['reference_label']=a.obs.celltype_major.astype(str).where(~a.obs.celltype_major.eq('T-cells'),a.obs.celltype_minor.astype(str))
    assert set(a.obs['orig.ident'])=={donor}
    assert np.array_equal(np.asarray(a.X.sum(axis=1)).ravel(),a.obs.nCount_RNA)
    assert np.array_equal(np.asarray((a.X>0).sum(axis=1)).ravel(),a.obs.nFeature_RNA)
    for t in types:assert a.obs.reference_label.eq(t).sum()==candidates.loc[donor,t]
    a.obs.to_csv(out/'cell_inventory.csv')
    qc=a.obs.groupby('reference_label')[['nCount_RNA','nFeature_RNA','percent.mito']].agg(['count','min','median','max'])
    qc.to_csv(out/'qc_summary.csv')
    rng=np.random.default_rng(protocol['seeds'][donor])
    pools={t:np.flatnonzero(a.obs.reference_label.eq(t)) for t in types}
    matrices=[];design=[];members=[];rna=[];cells=[];recipes=[]
    def save_recipe(selected,thinned,scenario,n_cd8,n_caf,repeat,pair_id):
        number=len(recipes);labels=a.obs.reference_label.iloc[selected].to_numpy()
        assert len(selected)==len(set(selected))
        composition={t:int((labels==t).sum()) for t in types}
        recipes.append(dict(base_mixture=number,scenario=scenario,pair_id=pair_id,n_cd8=n_cd8,n_caf=n_caf,repeat=repeat,**composition))
        for depth,component in [(1.,a.X[selected].tocsr()),(.25,thinned)]:
            spot_id=f'{donor}_mix{number:03d}_depth{int(depth*100)}'
            umis=np.asarray(component.sum(axis=1)).ravel();assert umis.sum()>0
            matrices.append(sparse.csr_matrix(component.sum(axis=0)))
            design.append(dict(spot_id=spot_id,donor_id=donor,base_mixture=number,pair_id=pair_id,scenario=scenario,n_cd8=n_cd8,n_caf=n_caf,repeat=repeat,depth_fraction=depth,n_cells=len(selected),nUMI=int(umis.sum()),x=len(design),y=0))
            rna.append([umis[labels==t].sum()/umis.sum() for t in types]);cells.append([(labels==t).sum()/len(labels) for t in types])
            members.extend(dict(spot_id=spot_id,base_mixture=number,source_barcode=a.obs_names[i],reference_label=t,component_umi=int(u)) for i,t,u in zip(selected,labels,umis))
    for n_cd8 in [0,2,5,8,10]:
        for repeat in range(5):
            selected=np.concatenate([rng.choice(pools['T cells CD4+'],10-n_cd8,replace=False),rng.choice(pools['T cells CD8+'],n_cd8,replace=False),rng.choice(pools['CAFs'],3,replace=False)])
            low=a.X[selected].tocsr().copy();low.data=rng.binomial(low.data,.25);low.eliminate_zeros()
            for n_caf in [0,1,3]:save_recipe(selected[:10+n_caf],low[:10+n_caf],'gradient',n_cd8,n_caf,repeat,f'cd8{n_cd8}_r{repeat}')
    for t in ['CAFs','T cells CD4+','T cells CD8+']:
        for repeat in range(5):
            selected=rng.choice(pools[t],5,replace=False)
            low=a.X[selected].tocsr().copy();low.data=rng.binomial(low.data,.25);low.eliminate_zeros()
            save_recipe(selected,low,'pure',5 if t=='T cells CD8+' else 0,5 if t=='CAFs' else 0,repeat,f'pure_{t}_r{repeat}')
    d=pd.DataFrame(design).set_index('spot_id');x=sparse.vstack(matrices,format='csr')
    assert len(d)==180 and len(recipes)==90
    d.to_csv(out/'simulation_design.csv');pd.DataFrame(recipes).to_csv(out/'frozen_recipes.csv',index=False)
    pd.DataFrame(members).to_csv(out/'mixture_membership.csv',index=False)
    pd.DataFrame(rna,index=d.index,columns=types).to_csv(out/'truth_rna_fraction.csv')
    pd.DataFrame(cells,index=d.index,columns=types).to_csv(out/'truth_cell_fraction.csv')
    mmwrite(inp/'spatial_counts.mtx',x.T)
    d[['x','y','nUMI']].to_csv(inp/'spatial_metadata.csv')
    pd.Series(a.var_names).to_csv(inp/'genes.tsv',header=False,index=False)
    ad.AnnData(X=x,obs=d,var=a.var.copy()).write_h5ad(inp/'simulated_spots.h5ad',compression='gzip')
    for p in [DATA/f'{donor}_counts_before_qc.h5ad',inp/'spatial_counts.mtx',inp/'spatial_metadata.csv',inp/'genes.tsv',out/'frozen_recipes.csv']:
        source_hashes[str(p.relative_to(ROOT))]=sha(p)
    print(donor,a.shape,d.groupby(['scenario','depth_fraction']).nUMI.agg(['size','min','median','max']).to_string(),flush=True)
(OUT/'simulation_hashes.json').write_text(json.dumps(source_hashes,indent=2)+'\n')
