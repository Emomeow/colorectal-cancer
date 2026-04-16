"""
Phase 5 Panel C — Sensitivity Analysis: OLR1+ TAM score threshold comparison
Tests 4 thresholds (top 10% / 20% / 30% / 50%) for RCRC vs LCRC volcano plots
Output: phase5_results/panelC_sensitivity_4panel.{png,pdf}
"""
import warnings; warnings.filterwarnings("ignore")
import io, sys, pathlib
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
import pandas as pd
import scipy.stats as stats
import anndata as ad
import scanpy as sc
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

OUT = pathlib.Path("phase5_results"); OUT.mkdir(exist_ok=True)
DEL = pathlib.Path("Project_Delivery/results/Figure5_Anatomical_Heterogeneity/Sensitivity_Analysis/figures")
DEL.mkdir(parents=True, exist_ok=True)

C_LEFT  = "#3498DB"
C_RIGHT = "#E74C3C"

# ── Load processed data ───────────────────────────────────────────
print("Loading processed h5ad...")
adata = ad.read_h5ad(OUT / "gse188711_processed.h5ad")

# Re-score (same as phase5_figure5.py)
adata_rn = adata.raw.to_adata()
sc.pp.normalize_total(adata_rn, target_sum=1e4)
sc.pp.log1p(adata_rn)

OLR1_MARKERS = ['SPP1','OLR1','TREM2','MARCO','LPL','MMP9','APOE','C1QC',
                'CTSL','GPNMB','CYP27A1','FABP5','VEGFA','LGALS1','FOLR2']
present = [g for g in OLR1_MARKERS if g in adata_rn.var_names]
sc.tl.score_genes(adata_rn, gene_list=present, score_name='OLR1TAM_score')
adata.obs['OLR1TAM_score'] = adata_rn.obs['OLR1TAM_score'].values
print(f"Scored {adata.n_obs} cells")

# Sample metadata (derived from gse188711_processed.h5ad obs)
SAMPLE_META = {
    'RS-CRC1': {'side': 'Right'},
    'R_CRC3':  {'side': 'Right'},
    'R_CRC4':  {'side': 'Right'},
    'JCA':     {'side': 'Left'},
    'LS-CRC3': {'side': 'Left'},
    'WGC':     {'side': 'Left'},
}

def compute_deg_at_threshold(adata, quantile_cutoff):
    """Compute RCRC vs LCRC DEG for cells above the quantile_cutoff OLR1 score."""
    threshold = adata.obs.OLR1TAM_score.quantile(quantile_cutoff)
    high_olr1 = adata[adata.obs.OLR1TAM_score > threshold].copy()
    n_left  = (high_olr1.obs['side'] == 'Left').sum()
    n_right = (high_olr1.obs['side'] == 'Right').sum()
    print(f"  quantile={quantile_cutoff:.0%}: {high_olr1.n_obs} cells (L={n_left}, R={n_right})")

    # Pseudo-bulk per sample
    from scipy.sparse import issparse
    pb_data = {}
    for samp in SAMPLE_META:
        sub = high_olr1[high_olr1.obs['sample'] == samp]
        if sub.n_obs < 3:
            continue
        # Use raw counts if available
        if sub.raw is not None:
            X = sub.raw.X
        else:
            X = sub.X
        if issparse(X): X = X.toarray()
        pb_data[samp] = X.sum(axis=0)

    if len(pb_data) < 2:
        return pd.DataFrame(), n_left, n_right

    var_names = high_olr1.raw.var_names if (high_olr1.raw is not None) else high_olr1.var_names
    pb_df = pd.DataFrame(pb_data, index=var_names).T
    side_map = {s: m['side'] for s, m in SAMPLE_META.items()}
    pb_meta = pd.DataFrame({'sample': list(pb_data.keys())})
    pb_meta['side'] = pb_meta['sample'].map(side_map)

    left_samps  = pb_meta[pb_meta.side == 'Left']['sample'].tolist()
    right_samps = pb_meta[pb_meta.side == 'Right']['sample'].tolist()
    left_pb  = pb_df.loc[[s for s in left_samps  if s in pb_df.index]]
    right_pb = pb_df.loc[[s for s in right_samps if s in pb_df.index]]

    deg_results = []
    for gene in pb_df.columns:
        lv = left_pb[gene].values  if len(left_pb)  > 0 else np.array([0])
        rv = right_pb[gene].values if len(right_pb) > 0 else np.array([0])
        lv_cpm = lv / (lv.sum() / 1e6 + 1)
        rv_cpm = rv / (rv.sum() / 1e6 + 1)
        log2fc = np.log2((rv_cpm.mean() + 1) / (lv_cpm.mean() + 1))
        _, pval = (stats.ttest_ind(rv_cpm, lv_cpm)
                   if (len(rv_cpm) > 1 and len(lv_cpm) > 1)
                   else (0, 1.0))
        deg_results.append(dict(gene=gene, log2fc=log2fc, pval=pval))

    deg_df = pd.DataFrame(deg_results)
    deg_df['pval'] = deg_df['pval'].fillna(1.0)
    from scipy.stats import false_discovery_control
    deg_df['padj'] = false_discovery_control(deg_df['pval'].values)
    deg_df['-log10p'] = -np.log10(deg_df['pval'].clip(1e-300) + 1e-300)
    # Filter uninformative genes
    deg_df = deg_df[deg_df.gene.apply(
        lambda g: not g.startswith('MT-') and not g.startswith('RP') and
                  'LINC' not in g and not g[:3].count('-') > 0)]
    return deg_df.reset_index(drop=True), n_left, n_right


# ── Run 4 thresholds ─────────────────────────────────────────────
THRESHOLDS = [
    (0.90, 'Top 10%'),
    (0.80, 'Top 20%'),
    (0.70, 'Top 30%'),   # current default
    (0.50, 'Top 50%\n(median split)'),
]

print("Computing DEG for 4 thresholds...")
results = []
for q, label in THRESHOLDS:
    deg, n_l, n_r = compute_deg_at_threshold(adata, q)
    results.append((q, label, deg, n_l, n_r))

# ── 4-panel Volcano figure ────────────────────────────────────────
fig, axes = plt.subplots(1, 4, figsize=(18, 5.5))
fig.subplots_adjust(wspace=0.35, left=0.05, right=0.97, top=0.85, bottom=0.15)

def volcano_color(row):
    if abs(row.log2fc) > 1.5 and row.pval < 0.05:
        return C_RIGHT if row.log2fc > 0 else C_LEFT
    return '#CCCCCC'

for ax, (q, label, deg, n_l, n_r) in zip(axes, results):
    if len(deg) == 0:
        ax.text(0.5, 0.5, 'Insufficient cells', ha='center', va='center',
                transform=ax.transAxes, fontsize=10)
        ax.set_title(f'{label}', fontsize=10, fontweight='bold')
        continue

    deg['vcolor'] = deg.apply(volcano_color, axis=1)
    n_sig_r = ((deg.log2fc > 1.5) & (deg.pval < 0.05)).sum()
    n_sig_l = ((deg.log2fc < -1.5) & (deg.pval < 0.05)).sum()

    ax.scatter(deg.log2fc, deg['-log10p'], c=deg.vcolor,
               s=4, alpha=0.6, rasterized=True)
    ax.axvline(1.5,  ls='--', color='#AAAAAA', lw=0.8)
    ax.axvline(-1.5, ls='--', color='#AAAAAA', lw=0.8)
    ax.axhline(-np.log10(0.05), ls='--', color='#AAAAAA', lw=0.8)

    # Label top genes
    top_r = deg[(deg.log2fc > 1.5) & (deg.pval < 0.1)].nlargest(6, 'log2fc')
    top_l = deg[(deg.log2fc < -1.5) & (deg.pval < 0.1)].nsmallest(6, 'log2fc')
    for _, row in pd.concat([top_r, top_l]).iterrows():
        ax.text(row.log2fc, row['-log10p'], row.gene,
                fontsize=6, ha='center', va='bottom',
                color=C_RIGHT if row.log2fc > 0 else C_LEFT)

    # Sig count annotation
    ax.text(0.97, 0.97, f'RCRC: {n_sig_r} DEGs\nLCRC: {n_sig_l} DEGs',
            transform=ax.transAxes, ha='right', va='top', fontsize=7.5,
            bbox=dict(boxstyle='round', fc='white', ec='#CCCCCC', lw=0.8))

    pct = int((1 - q) * 100)
    n_total = n_l + n_r
    ax.set_title(f'{label}\n(L={n_l}, R={n_r}, total={n_total})',
                 fontsize=9.5, fontweight='bold')
    # Mark current default
    if q == 0.70:
        ax.set_title(f'{label} [default]\n(L={n_l}, R={n_r}, total={n_total})',
                     fontsize=9.5, fontweight='bold', color='#1A1A2E')
        for spine in ax.spines.values():
            spine.set_edgecolor('#E74C3C')
            spine.set_linewidth(1.8)

    ax.set_xlabel('log$_2$FC (RCRC / LCRC)', fontsize=9)
    if ax is axes[0]:
        ax.set_ylabel('-log$_{10}$(p-value)', fontsize=9)
    ax.spines[['top','right']].set_visible(False)

patches = [mpatches.Patch(color=C_RIGHT, label='RCRC enriched (|FC|>1.5, p<0.05)'),
           mpatches.Patch(color=C_LEFT,  label='LCRC enriched (|FC|>1.5, p<0.05)'),
           mpatches.Patch(color='#CCCCCC', label='Not significant')]
fig.legend(handles=patches, fontsize=8, frameon=False,
           loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.01))

fig.suptitle(
    'Panel C Sensitivity Analysis — OLR1+ TAM Score Threshold Comparison\n'
    'RCRC vs LCRC DEG (GSE188711, pseudo-bulk t-test, |log2FC|>1.5, p<0.05)',
    fontsize=11, fontweight='bold', y=0.97)

for dest in [OUT, DEL]:
    dest.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest / 'panelC_sensitivity_4panel.png', dpi=200, bbox_inches='tight')
    fig.savefig(dest / 'panelC_sensitivity_4panel.pdf', bbox_inches='tight')
plt.close()
print("Saved panelC_sensitivity_4panel.{png,pdf}")
