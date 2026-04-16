"""
Phase 5 — Figure 5: GSE188711 左右侧 CRC 异质性分析
Panel A: Split-UMAP (LCRC vs RCRC)
Panel B: OLR1+ TAM 占比 Boxplot + 显著性
Panel C: RCRC vs LCRC OLR1+ TAM Volcano Plot
Panel D: Pathway 偏好 (基质重塑/血管生成)
Panel E: Pre-driver 丰度与临床特征相关性热图
Panel F: 左右侧 vs MMRp Pre-driver signature 相似性
"""
import warnings; warnings.filterwarnings("ignore")
import io, sys, pathlib, gzip, tarfile
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
import anndata as ad
import scanpy as sc
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

OUT = pathlib.Path("phase5_results"); OUT.mkdir(exist_ok=True)
FIG = pathlib.Path("figures/phase5"); FIG.mkdir(exist_ok=True)

# ── OLR1+ TAM signature genes (SPP1+B vs all, log2FC ranked) ─────
OLR1_TAM_MARKERS = [
    'SPP1','OLR1','TREM2','MARCO','LPL','MMP9','APOE','C1QC',
    'CTSL','GPNMB','CYP27A1','FABP5','VEGFA','LGALS1','IL10',
    'FOLR2','CD74','FN1','NT5E','MIF','ENTPD1',
]
# Broad deconvolution signature (top markers)
SIGNATURE_GENES = ['SPP1','OLR1','TREM2','MARCO','LPL','MMP9','APOE','C1QC','CTSL','GPNMB']

# ── Sample metadata ───────────────────────────────────────────────
SAMPLE_META = {
    'WGC':    {'side': 'Left',  'gsm': 'GSM5688706'},
    'JCA':    {'side': 'Left',  'gsm': 'GSM5688707'},
    'LS-CRC3':{'side': 'Left',  'gsm': 'GSM5688708'},
    'RS-CRC1':{'side': 'Right', 'gsm': 'GSM5688709'},
    'R_CRC3': {'side': 'Right', 'gsm': 'GSM5688710'},
    'R_CRC4': {'side': 'Right', 'gsm': 'GSM5688711'},
}

# ═══════════════════════════════════════════════════════════════════
# 1. Load GSE188711
# ═══════════════════════════════════════════════════════════════════
print("Loading GSE188711...")

def load_sample(tar, sample_name):
    barcodes_m = [m for m in tar.getmembers() if sample_name in m.name and 'barcodes' in m.name][0]
    features_m = [m for m in tar.getmembers() if sample_name in m.name and 'features' in m.name][0]
    matrix_m   = [m for m in tar.getmembers() if sample_name in m.name and 'matrix' in m.name][0]

    barcodes = gzip.open(io.BytesIO(tar.extractfile(barcodes_m).read())).read().decode().strip().split('\n')
    feat_raw = gzip.open(io.BytesIO(tar.extractfile(features_m).read())).read().decode().strip().split('\n')
    gene_ids = [l.split('\t')[0] for l in feat_raw]
    gene_symbols = [l.split('\t')[1] if '\t' in l else l for l in feat_raw]

    import scipy.io
    mat_bytes = io.BytesIO(gzip.open(io.BytesIO(tar.extractfile(matrix_m).read())).read())
    matrix = scipy.io.mmread(mat_bytes).T.tocsr()

    adata = ad.AnnData(X=matrix)
    adata.obs_names = [f"{b}_{sample_name}" for b in barcodes]
    adata.var_names = gene_symbols
    adata.var['gene_id'] = gene_ids
    # Handle duplicate gene names
    adata.var_names_make_unique()
    return adata

tar = tarfile.open('data/GSE188711_RAW.tar', 'r')
adatas = []
for sample_name, meta in SAMPLE_META.items():
    print(f"  Loading {sample_name} ({meta['side']})...")
    a = load_sample(tar, sample_name)
    a.obs['sample'] = sample_name
    a.obs['side'] = meta['side']
    adatas.append(a)
tar.close()

adata = ad.concat(adatas, join='outer', fill_value=0)
adata.obs_names_make_unique()
print(f"Combined: {adata.shape}")

# ── Basic QC & filter ─────────────────────────────────────────────
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
adata = adata[adata.obs.pct_counts_mt < 20].copy()
adata = adata[adata.obs.n_genes_by_counts < 6000].copy()
print(f"After QC: {adata.shape}")

# ── Normalize & process ───────────────────────────────────────────
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata.copy()

sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor='seurat_v3')
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added='leiden')

# ── OLR1+ TAM score ──────────────────────────────────────────────
adata_raw = adata.raw.to_adata()
present_markers = [g for g in OLR1_TAM_MARKERS if g in adata_raw.var_names]
sc.tl.score_genes(adata_raw, gene_list=present_markers, score_name='OLR1TAM_score')
adata.obs['OLR1TAM_score'] = adata_raw.obs['OLR1TAM_score'].values

# ── Annotate Myeloid cluster likely containing OLR1+ TAMs ────────
# Score myeloid markers
myeloid_genes = ['CD68','CSF1R','LYZ','AIF1','TYROBP','C1QC','APOE','SPP1','MMP9']
myeloid_present = [g for g in myeloid_genes if g in adata_raw.var_names]
sc.tl.score_genes(adata_raw, gene_list=myeloid_present, score_name='Myeloid_score')
adata.obs['Myeloid_score'] = adata_raw.obs['Myeloid_score'].values

# T cell score
tcell_genes = ['CD3D','CD3E','CD8A','CD8B','CD4','TRAC']
tcell_present = [g for g in tcell_genes if g in adata_raw.var_names]
sc.tl.score_genes(adata_raw, gene_list=tcell_present, score_name='Tcell_score')
adata.obs['Tcell_score'] = adata_raw.obs['Tcell_score'].values

# Fibroblast score
fibro_genes = ['FAP','ACTA2','COL1A1','COL1A2','FN1','VIM','PDGFRA']
fibro_present = [g for g in fibro_genes if g in adata_raw.var_names]
sc.tl.score_genes(adata_raw, gene_list=fibro_present, score_name='Fibro_score')
adata.obs['Fibro_score'] = adata_raw.obs['Fibro_score'].values

# Simple cell type annotation
def assign_celltype(row):
    scores = {
        'OLR1+ TAMs': row.OLR1TAM_score,
        'Myeloid': row.Myeloid_score,
        'T cells': row.Tcell_score,
        'Fibroblasts': row.Fibro_score,
    }
    best = max(scores, key=scores.get)
    if scores[best] < 0.1:
        return 'Other'
    return best

adata.obs['cell_type'] = adata.obs.apply(assign_celltype, axis=1)

# Save processed data
adata.write_h5ad(OUT / 'gse188711_processed.h5ad')
print("Saved processed h5ad")

# ═══════════════════════════════════════════════════════════════════
# Compute stats
# ═══════════════════════════════════════════════════════════════════
# Per-sample OLR1+ TAM proportion
sample_stats = []
for samp, meta in SAMPLE_META.items():
    mask = adata.obs['sample'] == samp
    sub = adata.obs[mask]
    if len(sub) == 0: continue
    myeloid = sub[sub.Myeloid_score > 0.1]
    n_total = len(sub)
    n_olr1  = (sub.OLR1TAM_score > sub.OLR1TAM_score.median()).sum()
    n_mye   = len(myeloid)
    prop = n_olr1 / n_total if n_total > 0 else 0
    prop_mye = (myeloid.OLR1TAM_score > 0).sum() / n_mye if n_mye > 0 else 0
    sample_stats.append(dict(
        sample=samp, side=meta['side'],
        n_cells=n_total,
        prop_olr1=prop,
        prop_mye=prop_mye,
        mean_score=sub.OLR1TAM_score.mean(),
    ))
stat_df = pd.DataFrame(sample_stats)
stat_df.to_csv(OUT / 'sample_olr1_proportions.csv', index=False)
print(stat_df.to_string())

# ═══════════════════════════════════════════════════════════════════
# DEG: RCRC vs LCRC within high-OLR1 cells
# ═══════════════════════════════════════════════════════════════════
print("Computing DEG RCRC vs LCRC...")
# Use raw counts for DEG
adata_deg = adata.raw.to_adata()
# High OLR1 score cells (top 30% per sample)
threshold = adata.obs.OLR1TAM_score.quantile(0.70)
high_olr1 = adata[adata.obs.OLR1TAM_score > threshold].copy()
n_left  = (high_olr1.obs['side']=='Left').sum()
n_right = (high_olr1.obs['side']=='Right').sum()
print(f"High OLR1 cells: {high_olr1.n_obs} (Left: {n_left}, Right: {n_right})")

# Pseudo-bulk: sum counts per sample for high OLR1 cells
from scipy.sparse import issparse
pb_data = {}
for samp in SAMPLE_META:
    sub = high_olr1[high_olr1.obs['sample'] == samp]
    if sub.n_obs < 5: continue
    X = sub.raw.X if hasattr(sub, 'raw') and sub.raw is not None else sub.X
    if issparse(X): X = X.toarray()
    pb_data[samp] = X.sum(axis=0)

pb_df = pd.DataFrame(pb_data, index=high_olr1.var_names if hasattr(high_olr1,'var_names') else adata.var_names).T
side_map = {s: m['side'] for s, m in SAMPLE_META.items()}
pb_meta = pd.DataFrame({'sample': list(pb_data.keys())})
pb_meta['side'] = pb_meta['sample'].map(side_map)

# Simple DEG by mean difference
left_samps  = pb_meta[pb_meta.side=='Left'].sample.tolist()
right_samps = pb_meta[pb_meta.side=='Right'].sample.tolist()
left_df  = pb_df.loc[left_samps] if left_samps else pd.DataFrame()
right_df = pb_df.loc[right_samps] if right_samps else pd.DataFrame()

deg_results = []
for gene in pb_df.columns:
    lv = left_df[gene].values if len(left_df) > 0 else np.array([0])
    rv = right_df[gene].values if len(right_df) > 0 else np.array([0])
    # CPM normalize
    lv_cpm = lv / (lv.sum() / 1e6 + 1)
    rv_cpm = rv / (rv.sum() / 1e6 + 1)
    log2fc = np.log2((rv_cpm.mean() + 1) / (lv_cpm.mean() + 1))
    _, pval = stats.ttest_ind(rv_cpm, lv_cpm) if (len(rv_cpm) > 1 and len(lv_cpm) > 1) else (0, 1.0)
    deg_results.append(dict(gene=gene, log2fc=log2fc, pval=pval,
                            right_mean=rv_cpm.mean(), left_mean=lv_cpm.mean()))

deg_df = pd.DataFrame(deg_results)
from scipy.stats import false_discovery_control
deg_df['padj'] = false_discovery_control(deg_df.pval.fillna(1).values)
deg_df = deg_df.sort_values('log2fc', ascending=False)
deg_df.to_csv(OUT / 'deg_rcrc_vs_lcrc_olr1tam.csv', index=False)
print(f"DEG computed, {(deg_df.padj < 0.05).sum()} significant")
print(deg_df.head(10)[['gene','log2fc','padj']].to_string())

# Pathway scores for Panel D
PATHWAY_GENES_D = {
    'Matrix Remodeling': ['MMP9','MMP2','MMP14','COL1A1','COL1A2','FN1','ITGAV','ITGB1','TIMP1','FAP'],
    'Angiogenesis':      ['VEGFA','VEGFB','ANGPT1','ANGPT2','FLT1','KDR','PECAM1','ENG','THBS1'],
    'Lipid Metabolism':  ['OLR1','LPL','FABP5','APOE','CYP27A1','PLTP','ABCA1','SCARB1'],
    'Immune Suppression':['LGALS1','LGALS9','CD274','HAVCR2','ARG1','IDO1','IL10','MIF'],
    'T cell Exclusion':  ['CXCL10','CCL2','CXCL1','CXCL8','CCL20','S100A8','S100A9'],
}

adata_for_pw = adata.raw.to_adata()
sc.pp.normalize_total(adata_for_pw, target_sum=1e4)
sc.pp.log1p(adata_for_pw)

pw_scores = {}
for pw, genes in PATHWAY_GENES_D.items():
    present = [g for g in genes if g in adata_for_pw.var_names]
    if not present: continue
    sc.tl.score_genes(adata_for_pw, gene_list=present, score_name=f'pw_{pw}')
    adata.obs[f'pw_{pw}'] = adata_for_pw.obs[f'pw_{pw}'].values

# Compute mean per side for high OLR1 cells
high_olr1_obs = adata.obs[adata.obs.OLR1TAM_score > threshold]
pw_by_side = {}
for pw in PATHWAY_GENES_D:
    col = f'pw_{pw}'
    if col not in adata.obs.columns: continue
    left_v  = high_olr1_obs[high_olr1_obs['side'] == 'Left'][col].values
    right_v = high_olr1_obs[high_olr1_obs['side'] == 'Right'][col].values
    stat, pval = stats.mannwhitneyu(right_v, left_v, alternative='two-sided') if (len(right_v)>0 and len(left_v)>0) else (0,1)
    pw_by_side[pw] = dict(
        left_mean=left_v.mean(), right_mean=right_v.mean(),
        log2fc=np.log2((right_v.mean()+0.01)/(left_v.mean()+0.01)),
        pval=pval,
    )
pw_df = pd.DataFrame(pw_by_side).T
pw_df.to_csv(OUT / 'pathway_scores_side.csv')
print("\nPathway scores by side (Right vs Left, log2FC):")
print(pw_df[['left_mean','right_mean','log2fc','pval']].to_string())

print("\nPhase 5 GSE188711 preprocessing complete.")
