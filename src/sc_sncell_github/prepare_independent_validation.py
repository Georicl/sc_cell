"""核对候选供者，冻结CID4463验证设计并从作者矩阵提取计数。"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import xml.etree.ElementTree as ET
import zipfile
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.io import mmwrite

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'data/raw/wu2021_scrna/Wu_etal_2021_BRCA_scRNASeq'
OUT=ROOT/'results/Wu2021_independent_validation_v1'
DATA=ROOT/'data/processed/Wu2021_independent_validation_v1'
OUT.mkdir(parents=True,exist_ok=True);DATA.mkdir(parents=True,exist_ok=True)
DONOR='CID4463'
BASE=ROOT/'data/processed/Wu2021_known_composition_v1'
COMPARE=ROOT/'data/processed/Wu2021_reference_comparison_v1'
VARIANTS=['baseline','donor_balanced','epithelial_split','balanced_split']

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def clinical_rows():
    # 补充表1为工作簿第一张表；直接读取单元格值，保留Excel行号。
    ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(ROOT/'data/provenance/wu2021_supplementary_tables.xlsx') as z:
        workbook=ET.fromstring(z.read('xl/workbook.xml'))
        assert workbook.find('s:sheets',ns)[0].get('name')=='Supplementary Table 1'
        strings=[''.join(x.itertext()) for x in ET.fromstring(z.read('xl/sharedStrings.xml'))]
        rows=[]
        for row in ET.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('.//s:row',ns):
            values={}
            for cell in row.findall('s:c',ns):
                value=cell.find('s:v',ns)
                if value is None:continue
                value=value.text
                if cell.get('t')=='s':value=strings[int(value)]
                col=''.join(c for c in cell.get('r') if c.isalpha())
                values[col]=value
            if values.get('A','').isdigit():
                rows.append(dict(case_id=values['A'],clinical_excel_row=int(row.get('r')),
                    clinical_subtype=values.get('K'),histology=values.get('E'),
                    treatment_status=values.get('L'),treatment_details=values.get('M'),
                    pathological_features=values.get('N')))
    return pd.DataFrame(rows)

metadata=pd.read_csv(SOURCE/'metadata.csv',index_col=0)
assert metadata.index.is_unique
metadata['reference_label']=metadata.celltype_major.where(~metadata.celltype_major.eq('T-cells'),metadata.celltype_minor)
clinical=clinical_rows()
clinical.to_csv(OUT/'clinical_source_records.csv',index=False)
counts=pd.crosstab(metadata['orig.ident'],metadata.reference_label)
subtype=metadata.groupby('orig.ident').subtype.agg(lambda x: '|'.join(sorted(x.unique())))
candidates=counts.copy()
candidates.insert(0,'subtype',subtype)
candidates.insert(1,'total_cells',counts.sum(axis=1))
candidates.index.name='donor_id'
candidates=candidates.reset_index()
candidates['case_id']=candidates.donor_id.str.extract(r'^CID(\d+)$',expand=False)
candidates=candidates.merge(clinical,on='case_id',how='left',validate='many_to_one')
roles=[]
for r in candidates.itertuples(index=False):
    donor=r.donor_id
    if donor in ['CID4471','CID4535']:role,reason='reference','已用于参考构建'
    elif donor=='CID4067':role,reason='development','已用于误差检查和方案开发'
    elif pd.isna(r.case_id):role,reason='mapping_unresolved','带后缀样本未确认与临床Case ID的对应'
    elif r.subtype!=r.clinical_subtype:role,reason='subtype_conflict','作者metadata与临床亚型字段不一致'
    elif r.subtype!='ER+':role,reason='outside_ER_scope','不属于本轮ER+验证范围'
    elif r.treatment_status!='Naïve':role,reason='treated_or_unknown','不是已确认未治疗供者'
    elif donor==DONOR:role,reason='independent_validation','临床精确对应，未治疗ER+，CAF/CD8/癌上皮/正常上皮均有覆盖'
    elif counts.loc[donor,'Cancer Epithelial']==0:role,reason='no_cancer_control','缺少作者标注癌上皮，不能完成本轮上皮验证'
    else:role,reason='reserve_validation','未治疗ER+，保留为后续验证；本轮未使用'
    roles.append((role,reason))
candidates[['role','reason']]=roles
candidates.to_csv(OUT/'donor_candidates.csv',index=False)
row=candidates.set_index('donor_id').loc[DONOR]
assert row.clinical_subtype=='ER+' and row.treatment_status=='Naïve' and row.histology=='IDC'
assert (counts.loc[DONOR,['CAFs','T cells CD8+','Cancer Epithelial','Normal Epithelial']]>=10).all()
old=pd.read_csv(ROOT/'results/Wu2021_donor_reference_inventory/donor_by_reference_celltype.csv',index_col=0)
for t in counts.columns:assert np.array_equal(counts.loc[old.index,t],old[t])
backgrounds=['T cells CD4+','Myeloid','Endothelial','B-cells','NK cells','NKT cells','Cycling T-cells']
available=[t for t in backgrounds if counts.loc[DONOR,t]>=1]
pure=[t for t in sorted(counts.columns) if counts.loc[DONOR,t]>=5]
protocol=dict(frozen_at_utc=datetime.now(timezone.utc).isoformat(),validation_donors=[DONOR],
    reference_donors=['CID4471','CID4535'],development_donors=['CID4067'],variants=VARIANTS,
    seed=20260917,rctd_seed=9,reference_adjustments='none; reuse frozen reference matrices, labels and NNLS profiles',
    donor_selection='Exact clinical mapping, untreated ER+, contains CAF/CD8/cancer/normal epithelium; no validation predictions inspected before selection.',
    cell_policy='Retain all author cells as in previous benchmark; verify raw UMI and detected genes against metadata; QC summaries are descriptive.',
    factorial=dict(caf=[0,1,3],cd8=[0,1,3],fixed_pvl=1,backgrounds=available,background_cells=1,total_cells=10,repeats=3),
    pure=dict(types=pure,cells=5,repeats=5,skip_below_cells=5),
    epithelial_titration=dict(normal_cells=[0,2,5,8,10],total_cells=10,repeats=5),
    depths=[1.,.25],sampling='Without replacement within spot; reuse across spots; paired depths share source cells; binomial thinning per component cell.',
    primary_truth='post-thinning captured RNA UMI fraction',secondary_truth='cell-number fraction',
    metrics=['MAE','RMSE','bias','absent mean prediction','absent >=1% frequency','pure dominant type match'],
    rctd=dict(mode='full',max_cores=2,ref_UMI_min=100,ref_n_cells_min=25),
    nnls=dict(features='reuse each frozen reference profile with its original 2000 reference-selected genes',normalization='full-library UMI; coefficients normalized per spot'),
    interpretation='Report all four references and both methods; no tuning on this validation donor.')
(OUT/'protocol.json').write_text(json.dumps(protocol,ensure_ascii=False,indent=2)+'\n')
paths=[SOURCE/'metadata.csv',ROOT/'data/provenance/wu2021_supplementary_tables.xlsx']
for v in VARIANTS:
    folder=BASE if v=='baseline' else COMPARE/v
    paths.extend(folder/x for x in ['reference_counts.mtx','reference_metadata.csv','genes.tsv','nnls_reference_profiles.csv'])
frozen={str(p.relative_to(ROOT)):digest(p) for p in paths}
(OUT/'frozen_reference_hashes.json').write_text(json.dumps(frozen,indent=2)+'\n')
print('Frozen donor:',DONOR,'backgrounds:',available,'pure types:',pure,flush=True)

# 提取新供者计数，不加载其他供者的表达矩阵。
import extract_er_qc_donors as extraction
extraction.DONORS=[DONOR]
extraction.OUTPUT=DATA
target=DATA/f'{DONOR}_counts_before_qc.h5ad'
if not target.exists():extraction.extract_donors()
test=ad.read_h5ad(target)
test.obs['reference_label']=metadata.loc[test.obs_names,'reference_label'].astype(str)
assert set(test.obs['orig.ident'])=={DONOR}
assert np.array_equal(np.asarray(test.X.sum(axis=1)).ravel(),test.obs.nCount_RNA)
assert np.array_equal(np.asarray((test.X>0).sum(axis=1)).ravel(),test.obs.nFeature_RNA)
pd.crosstab(test.obs.reference_label,test.obs.celltype_minor).to_csv(OUT/'validation_subtype_counts.csv')
test.obs.groupby('reference_label')[['nCount_RNA','nFeature_RNA','percent.mito']].agg(['count','min','median','max']).to_csv(OUT/'validation_qc_summary.csv')
test.obs[['orig.ident','reference_label','celltype_minor','nCount_RNA','nFeature_RNA','percent.mito']].to_csv(OUT/'validation_cell_inventory.csv')

types=sorted(counts.columns)
recipes=[]
for caf in [0,1,3]:
    for cd8 in [0,1,3]:
        for bg in available:
            for repeat in range(3):recipes.append(('factorial',f'caf{caf}_cd8{cd8}_{bg}',repeat,{'CAFs':caf,'T cells CD8+':cd8,'PVL':1,bg:1,'Cancer Epithelial':8-caf-cd8}))
for t in pure:
    for repeat in range(5):recipes.append(('pure',t,repeat,{t:5}))
for normal in [0,2,5,8,10]:
    for repeat in range(5):recipes.append(('epithelial_titration',f'normal{normal}_cancer{10-normal}',repeat,{'Normal Epithelial':normal,'Cancer Epithelial':10-normal}))
rng=np.random.default_rng(protocol['seed'])
pools={t:np.flatnonzero(test.obs.reference_label.eq(t)) for t in types}
matrices=[];design=[];members=[];rna=[];cells=[];recipe_rows=[]
for number,(scenario,recipe,repeat,composition) in enumerate(recipes):
    selected=np.concatenate([rng.choice(pools[t],n,replace=False) for t,n in composition.items() if n])
    labels=test.obs.reference_label.iloc[selected].to_numpy()
    source=test.X[selected].tocsr()
    recipe_rows.append(dict(base_mixture=number,scenario=scenario,recipe=recipe,repeat=repeat,**{t:composition.get(t,0) for t in types}))
    for depth in protocol['depths']:
        spot_id=f'{DONOR}_mix{number:04d}_depth{int(depth*100)}'
        component=source.copy()
        if depth<1:
            component.data=rng.binomial(component.data,depth);component.eliminate_zeros()
        umis=np.asarray(component.sum(axis=1)).ravel()
        assert umis.sum()>0
        rna.append([umis[labels==t].sum()/umis.sum() for t in types])
        cells.append([(labels==t).sum()/len(labels) for t in types])
        matrices.append(sparse.csr_matrix(component.sum(axis=0)))
        design.append(dict(spot_id=spot_id,donor_id=DONOR,base_mixture=number,scenario=scenario,recipe=recipe,repeat=repeat,depth_fraction=depth,n_cells=len(selected),nUMI=int(umis.sum()),x=len(design),y=0))
        members.extend(dict(spot_id=spot_id,base_mixture=number,source_barcode=test.obs_names[i],reference_label=t,component_umi=int(u)) for i,t,u in zip(selected,labels,umis))
design=pd.DataFrame(design).set_index('spot_id')
matrix=sparse.vstack(matrices,format='csr')
assert len(design)==474 and len(recipes)==237
pd.DataFrame(recipe_rows).to_csv(OUT/'frozen_recipes.csv',index=False)
design.to_csv(OUT/'simulation_design.csv')
pd.DataFrame(members).to_csv(OUT/'mixture_membership.csv',index=False)
pd.DataFrame(rna,index=design.index,columns=types).to_csv(OUT/'truth_rna_fraction.csv')
pd.DataFrame(cells,index=design.index,columns=types).to_csv(OUT/'truth_cell_fraction.csv')
mmwrite(DATA/'spatial_counts.mtx',matrix.T)
pd.Series(test.var_names).to_csv(DATA/'genes.tsv',header=False,index=False)
design[['x','y','nUMI']].to_csv(DATA/'spatial_metadata.csv')
ad.AnnData(X=matrix,obs=design,var=test.var.copy()).write_h5ad(DATA/'simulated_spots.h5ad',compression='gzip')
(OUT/'simulation_hashes.json').write_text(json.dumps({str(p.relative_to(ROOT)):digest(p) for p in [target,DATA/'spatial_counts.mtx',DATA/'spatial_metadata.csv',DATA/'genes.tsv',OUT/'frozen_recipes.csv',OUT/'protocol.json']},indent=2)+'\n')
print(design.groupby(['scenario','depth_fraction']).nUMI.agg(['size','min','median','max']).to_string(),flush=True)
