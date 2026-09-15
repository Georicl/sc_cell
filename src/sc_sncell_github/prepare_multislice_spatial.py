"""固定双供者参考和旧CID4535基因范围，逐切片读取作者计数与质控。"""
from pathlib import Path
from datetime import datetime, timezone
import gzip
import hashlib
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy.io import mmread, mmwrite

ROOT=Path(__file__).resolve().parents[2]
RAW=ROOT/'data/raw/wu2021_visium'
DATA=ROOT/'data/processed/Wu2021_multislice_CAF_T_v1'
OUT=ROOT/'results/Wu2021_multislice_CAF_T_v1'
DATA.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
SAMPLES=['CID4535','CID4290','CID4465','CID44971','1142243F','1160920F']
protocol=dict(frozen_at_utc=datetime.now(timezone.utc).isoformat(),samples=SAMPLES,reference_donors=['CID4471','CID4535'],
    labels='Original 13 reference labels; sum CD4, CD8, NKT and Cycling T after fitting. NK separate.',
    gene_scope='Use original CID4535 15770-gene input list as template. Intersect with each slice feature names, without alias guessing or zero padding; record omitted genes. CID4535 keeps its exact original gene input.',
    qc='Exclude author Artefact OR in_tissue==0; missing pathology and Uncertain retained in expression analysis, excluded from pathology summaries. No new UMI/mitochondrial filter.',
    rctd=dict(seed=9,mode='full',max_cores=2,ref_UMI_min=100,ref_n_cells_min=25),
    comparison='CID4535 old single-donor fit vs new two-donor fit; both use identical spatial counts and gene list.',
    expansion='Other local Wu slices; ER and TNBC displayed separately; ER reference applied to TNBC as transfer check.',
    marker_checks=dict(T=['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2'],CAF=['COL1A1','COL1A2','DCN','LUM'],NK=['NKG7','GNLY','KLRD1']),
    marker_score='Mean log1p(CP10k) using full measured spatial-library UMI; expression concordance, not independent cell identity validation.',
    high_spots='Top ceil(10% * n_spots), stable barcode order breaks ties',
    inference='Describe CAF and aggregate T distributions and reference sensitivity; no causal/exclusion inference from relative-weight correlations.')
(OUT/'protocol.json').write_text(json.dumps(protocol,ensure_ascii=False,indent=2)+'\n')
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()
def stream(path):
    with path.open('rb') as f:magic=f.read(2)
    return gzip.open(path,'rb') if magic==b'\x1f\x8b' else path.open('rb')
old_input=ROOT/'data/processed/CID4535_combine_v1/rctd_input'
genes=pd.read_csv(old_input/'genes.tsv',header=None)[0].astype(str)
assert len(genes)==15770 and genes.is_unique
reference=ad.read_h5ad(ROOT/'data/processed/Wu2021_known_composition_v1/reference_counts.h5ad')
assert set(genes).issubset(reference.var_names)
ref=reference[:,genes].copy()
meta=ref.obs[['donor_id','reference_label']].copy();meta['nUMI']=np.asarray(ref.X.sum(axis=1)).ravel()
assert meta.loc[meta.nUMI>=100,'reference_label'].value_counts().min()>=25
meta.to_csv(DATA/'reference_metadata.csv',index_label='barcode')
mmwrite(DATA/'reference_counts.mtx',ref.X.T)
genes.to_csv(DATA/'genes.tsv',header=False,index=False)
meta.groupby(['donor_id','reference_label']).size().rename('n_cells').to_csv(OUT/'reference_inventory.csv')
hashes={str(p.relative_to(ROOT)):sha(p) for p in [DATA/'reference_counts.mtx',DATA/'reference_metadata.csv',DATA/'genes.tsv',OUT/'protocol.json',old_input/'spatial_counts.mtx']}
inventory=[];checks={}
for sample in SAMPLES:
    inp=DATA/sample;out=OUT/sample;inp.mkdir(exist_ok=True);out.mkdir(exist_ok=True)
    counts_dir=RAW/'filtered_count_matrices'/f'{sample}_filtered_count_matrix'
    spatial_dir=RAW/'spatial'/f'{sample}_spatial'
    paths={n:counts_dir/n for n in ['matrix.mtx.gz','barcodes.tsv.gz','features.tsv.gz']}
    with stream(paths['matrix.mtx.gz']) as f:x=mmread(f).T.tocsr()
    with stream(paths['barcodes.tsv.gz']) as f:barcodes=pd.read_csv(f,sep='\t',header=None)[0].astype(str)
    with stream(paths['features.tsv.gz']) as f:features=pd.read_csv(f,sep='\t',header=None)[0].astype(str)
    assert barcodes.is_unique and features.is_unique and x.shape==(len(barcodes),len(features))
    metadata_path=RAW/'metadata'/f'{sample}_metadata.csv'
    metadata=pd.read_csv(metadata_path,index_col=0)
    position_path=spatial_dir/'tissue_positions_list.csv'
    position=pd.read_csv(position_path,header=None,index_col=0,names=['barcode','in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres'])
    assert metadata.index.is_unique and position.index.is_unique
    assert set(barcodes)==set(metadata.index) and set(barcodes).issubset(position.index)
    obs=metadata.loc[barcodes].join(position.loc[barcodes]);obs.index.name='barcode'
    total=np.asarray(x.sum(axis=1)).ravel();detected=np.asarray((x>0).sum(axis=1)).ravel()
    checks[sample+'_author_umi_exact']=bool(np.array_equal(total,obs.nCount_RNA))
    checks[sample+'_author_genes_exact']=bool(np.array_equal(detected,obs.nFeature_RNA))
    assert checks[sample+'_author_umi_exact'] and checks[sample+'_author_genes_exact']
    obs['total_counts']=total;obs['detected_genes']=detected
    mt=features.str.startswith('MT-').to_numpy()
    obs['mito_pct']=np.asarray(x[:,mt].sum(axis=1)).ravel()/total*100
    obs['outside_tissue']=obs.in_tissue.eq(0)
    obs['author_artefact']=obs.Classification.eq('Artefact').fillna(False)
    obs['excluded']=obs.outside_tissue|obs.author_artefact
    obs['pathology_eligible']=obs.Classification.notna()&~obs.Classification.isin(['Artefact','Uncertain'])&~obs.excluded
    obs.to_csv(out/'all_spot_qc.csv')
    keep=~obs.excluded.to_numpy();retained=obs.loc[keep].copy();full=x[keep]
    retained['pathology_display']=retained.Classification.fillna('Missing')
    available=genes.isin(features)
    sample_genes=genes[available]
    pd.DataFrame({'gene':genes[~available],'reason':'not present under exact name in slice features'}).to_csv(out/'omitted_template_genes.csv',index=False)
    sample_genes.to_csv(inp/'genes.tsv',header=False,index=False)
    idx=pd.Index(features).get_indexer(sample_genes)
    assert (idx>=0).all()
    shared=full[:,idx].tocsr()
    umi=np.asarray(shared.sum(axis=1)).ravel()
    spatial_meta=pd.DataFrame(dict(x=retained.pxl_col_in_fullres,y=retained.pxl_row_in_fullres,nUMI=umi),index=retained.index)
    spatial_meta.to_csv(inp/'spatial_metadata.csv');mmwrite(inp/'spatial_counts.mtx',shared.T)
    ad.AnnData(X=full,obs=retained,var=pd.DataFrame(index=pd.Index(features,name='gene'))).write_h5ad(inp/'spatial_counts_full.h5ad',compression='gzip')
    if sample=='CID4535':
        old=mmread(old_input/'spatial_counts.mtx').T.tocsr()
        old_meta=pd.read_csv(old_input/'spatial_metadata.csv',index_col=0)
        checks['CID4535_same_spot_order']=retained.index.equals(old_meta.index)
        checks['CID4535_same_spatial_matrix']=(shared!=old).nnz==0
        assert checks['CID4535_same_spot_order'] and checks['CID4535_same_spatial_matrix']
    inventory.append(dict(sample_id=sample,subtype='|'.join(sorted(retained.subtype.unique())),n_before=len(obs),n_outside=int(obs.outside_tissue.sum()),n_artefact=int(obs.author_artefact.sum()),n_excluded=int(obs.excluded.sum()),n_retained=len(retained),missing_pathology=int(retained.Classification.isna().sum()),uncertain_pathology=int(retained.Classification.eq('Uncertain').sum()),model_genes=len(sample_genes),omitted_genes=int((~available).sum()),model_umi_min=int(umi.min()),model_umi_median=float(np.median(umi)),full_umi_median=float(np.median(retained.total_counts)),reference_overlap=sample=='CID4535'))
    for p in list(paths.values())+[metadata_path,position_path,inp/'spatial_counts.mtx',inp/'spatial_metadata.csv',inp/'genes.tsv']:
        hashes[str(p.relative_to(ROOT))]=sha(p)
    print(sample,inventory[-1],flush=True)
pd.DataFrame(inventory).to_csv(OUT/'slice_inventory.csv',index=False)
(OUT/'input_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
(OUT/'input_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
