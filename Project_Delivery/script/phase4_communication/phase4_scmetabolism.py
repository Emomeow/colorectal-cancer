"""
Phase 4 Panel D/E — AUCell-based metabolic pathway scoring
Equivalent to scMetabolism (AUCell mode) using decoupler-py + KEGG 2021
Replaces ssGSEA-like mean expression scoring.
"""
import sys, io, warnings, pathlib
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import gseapy
import decoupler as dc
import scipy.sparse as sp

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

OUT = pathlib.Path("phase4_results")
OUT.mkdir(exist_ok=True)

# ── Load data ──────────────────────────────────────────────────────────
print("Loading GSE144735...")
adata = ad.read_h5ad("processed_data/GSE144735.h5ad")
print(f"  {adata.n_obs} cells, {adata.n_vars} genes")
print(f"  Cell subtypes: {adata.obs.Cell_subtype.value_counts().head(8).to_dict()}")

# ── Build KEGG metabolic gene sets ────────────────────────────────────
print("Fetching KEGG 2021 Human metabolic gene sets...")
kegg_all = gseapy.get_library('KEGG_2021_Human')

# Select metabolic pathways (same logic as scMetabolism)
METAB_KEYWORDS = [
    'Metabolism', 'Biosynthesis', 'Degradation', 'Fatty acid',
    'Oxidative phosphorylation', 'Glycolysis', 'Steroid', 'Sphingolipid',
    'Purine', 'Pyrimidine', 'Arachidonic', 'Linoleic', 'Butanoate',
    'Pyruvate', 'Citrate', 'Carbon', 'Nitrogen', 'Glyoxylate',
    'Ascorbate', 'Folate', 'Thiamine', 'Porphyrin', 'Terpenoid',
    'Propanoate', 'Inositol', 'Amino',
]
metab_sets = {k: v for k, v in kegg_all.items()
              if any(kw in k for kw in METAB_KEYWORDS)}

# Focus on the most relevant pathways for TAM biology (keep manageable)
FOCAL_PATHWAYS = {
    'Lipid Uptake\n& FAO':             ['OLR1','CD36','FABP4','FABP5','LPL','APOE',
                                         'CPT1A','ACSL4','HADHA','ACADM','ACADVL'],
    'Cholesterol\nBiosynthesis':        metab_sets.get('Steroid biosynthesis', [])[:20],
    'Glycolysis':                       metab_sets.get('Glycolysis / Gluconeogenesis', [])[:20],
    'Oxidative\nPhosphorylation':       metab_sets.get('Oxidative phosphorylation', [])[:20],
    'Fatty Acid\nBiosynthesis':         metab_sets.get('Fatty acid biosynthesis', [])[:18],
    'Fatty Acid\nDegradation':          metab_sets.get('Fatty acid degradation', [])[:20],
    'Purine\nMetabolism':               metab_sets.get('Purine metabolism', [])[:20],
    'Amino Sugar\nMetabolism':          metab_sets.get('Amino sugar and nucleotide sugar metabolism', [])[:20],
    'TCA Cycle':                        metab_sets.get('Citrate cycle (TCA cycle)', [])[:20],
    'Arachidonic\nAcid Metabolism':     metab_sets.get('Arachidonic acid metabolism', [])[:20],
}

# Filter to genes present in adata
focal_clean = {}
for pw, genes in FOCAL_PATHWAYS.items():
    present = [g for g in genes if g in adata.var_names]
    if len(present) >= 3:
        focal_clean[pw] = present
        print(f"  {pw.replace(chr(10),' ')}: {len(present)}/{len(genes)} genes present")

# ── Convert to decoupler long-format network ─────────────────────────
net_rows = []
for pw, genes in focal_clean.items():
    for g in genes:
        net_rows.append({'source': pw, 'target': g, 'weight': 1.0})
net_df = pd.DataFrame(net_rows)
print(f"\nNetwork: {len(net_df)} gene-pathway links, {net_df.source.nunique()} pathways")

# ── Subset to raw counts for AUCell ──────────────────────────────────
# Use raw layer if available
if adata.raw is not None:
    adata_raw = adata.raw.to_adata()
else:
    adata_raw = adata.copy()

# Ensure integer counts
if sp.issparse(adata_raw.X):
    adata_raw.X = adata_raw.X.toarray()
adata_raw.X = np.round(adata_raw.X).astype(np.float32)

# ── Run AUCell (decoupler) ────────────────────────────────────────────
print("\nRunning AUCell scoring (scMetabolism-equivalent)...")
dc.run_aucell(
    adata_raw,
    net=net_df,
    source='source',
    target='target',
    use_raw=False,
    verbose=True,
)
# AUCell scores are stored in adata_raw.obsm['aucell_estimate']
aucell_df = adata_raw.obsm['aucell_estimate'].copy()
aucell_df.index = adata_raw.obs_names
print(f"AUCell done: {aucell_df.shape} (cells x pathways)")

# Copy scores back to main adata
for col in aucell_df.columns:
    adata.obs[f'aucell_{col}'] = aucell_df[col].values

# ── Cell type groups for Panel D/E ───────────────────────────────────
CELL_TYPES = {
    'OLR1+ TAMs\n(SPP1+B)':     'SPP1+B',
    'Anti-inflam\nTAMs':         'Anti-inflammatory',
    'CD8+ T cells':              'CD8+ T cells',
    'Myofibroblasts':            'Myofibroblasts',
    'Reg T cells':               'Regulatory T cells',
}
CT_COLORS = {
    'OLR1+ TAMs\n(SPP1+B)':  '#E74C3C',
    'Anti-inflam\nTAMs':      '#F1948A',
    'CD8+ T cells':           '#2ECC71',
    'Myofibroblasts':         '#3498DB',
    'Reg T cells':            '#9B59B6',
}

pw_names = list(focal_clean.keys())

# Compute mean AUCell score per cell type per pathway
score_data = {}
for ct_label, subtype in CELL_TYPES.items():
    mask = adata.obs['Cell_subtype'] == subtype
    sub_scores = aucell_df[mask] if mask.any() else pd.DataFrame()
    row = {}
    for pw in pw_names:
        if pw in sub_scores.columns and len(sub_scores) > 0:
            row[pw] = sub_scores[pw].mean()
        else:
            row[pw] = 0.0
    score_data[ct_label] = row

score_df = pd.DataFrame(score_data).T   # cell types x pathways

# Z-score across cell types per pathway
score_z = score_df.copy()
for col in score_z.columns:
    m, s = score_z[col].mean(), score_z[col].std()
    score_z[col] = (score_z[col] - m) / (s + 1e-9)

print("\nAUCell mean scores (raw):")
print(score_df.round(4).to_string())
print("\nZ-scores:")
print(score_z.round(3).to_string())

# ═══════════════════════════════════════════════════════════════════════
# PANEL D: Grouped bar chart — AUCell z-scores by cell type
# ═══════════════════════════════════════════════════════════════════════
print("\nDrawing Panel D...")

n_pw = len(pw_names)
n_ct = len(CELL_TYPES)
bar_w = 0.14
x = np.arange(n_pw)

fig_d, ax_d = plt.subplots(figsize=(13, 5))

ct_labels = list(CELL_TYPES.keys())
ct_colors = [CT_COLORS[ct] for ct in ct_labels]

for i, (ct, color) in enumerate(zip(ct_labels, ct_colors)):
    vals = [score_z.loc[ct, pw] for pw in pw_names]
    offset = (i - n_ct / 2 + 0.5) * bar_w
    ax_d.bar(x + offset, vals, bar_w, color=color, alpha=0.85,
             label=ct.replace('\n', ' '), edgecolor='white', lw=0.5)

# Highlight OLR1+ TAM bars
for i_pw in range(n_pw):
    val = score_z.loc[ct_labels[0], pw_names[i_pw]]
    if val > 0.5:
        ax_d.axvspan(i_pw - 0.45, i_pw + 0.45, alpha=0.05, color='#E74C3C')

ax_d.set_xticks(x)
ax_d.set_xticklabels([p.replace('\n', ' ') for p in pw_names],
                      fontsize=8.5, ha='right', rotation=35)
ax_d.set_ylabel('AUCell score (z-score across cell types)', fontsize=9)
ax_d.set_title('D  Metabolic Pathway Activity by Cell Type (GSE144735)\n'
               'AUCell scoring on KEGG 2021 gene sets (scMetabolism-equivalent)',
               fontsize=11, fontweight='bold', loc='left')
ax_d.axhline(0, color='black', lw=0.8, ls='--')
ax_d.spines[['top', 'right']].set_visible(False)
ax_d.legend(fontsize=7.5, frameon=False, ncol=3, loc='upper right')

fig_d.tight_layout()
fig_d.savefig(OUT / 'panelD_metabolism.png', dpi=200, bbox_inches='tight')
fig_d.savefig(OUT / 'panelD_metabolism.pdf', bbox_inches='tight')
plt.close()
print("Panel D saved.")

# ═══════════════════════════════════════════════════════════════════════
# PANEL E: Heatmap — genes x cell types (transposed, with annotation bar)
# ═══════════════════════════════════════════════════════════════════════
print("Drawing Panel E (heatmap)...")

# For Panel E: use raw AUCell pathway scores (not z-scored) across all pathways
# Rows = pathways, Cols = cell types  (cleaner than gene-level for AUCell)
# Z-score across cell types per pathway for colour scale
hm_df = score_z.T   # pathways x cell types

grp_colors_e = {
    'Lipid Uptake\n& FAO':         '#E74C3C',
    'Fatty Acid\nDegradation':     '#E67E22',
    'Fatty Acid\nBiosynthesis':    '#F39C12',
    'Cholesterol\nBiosynthesis':   '#8E44AD',
    'TCA Cycle':                   '#2980B9',
    'Oxidative\nPhosphorylation':  '#1ABC9C',
    'Glycolysis':                  '#27AE60',
    'Purine\nMetabolism':          '#16A085',
    'Arachidonic\nAcid Metabolism':'#7F8C8D',
    'Amino Sugar\nMetabolism':     '#95A5A6',
}

ct_order  = list(CELL_TYPES.keys())
pw_order  = list(hm_df.index)
hm_vals   = hm_df[ct_order].values   # pathways x cell_types

n_pw_e  = len(pw_order)
n_ct_e  = len(ct_order)

fig_e = plt.figure(figsize=(7, 8))
gs_e  = fig_e.add_gridspec(1, 3, width_ratios=[0.04, 1, 0.04],
                             wspace=0.02, left=0.38, right=0.97,
                             top=0.92, bottom=0.12)
ax_grp = fig_e.add_subplot(gs_e[0])
ax_hm  = fig_e.add_subplot(gs_e[1])
ax_cb  = fig_e.add_subplot(gs_e[2])

im = ax_hm.imshow(hm_vals, cmap='RdBu_r', aspect='auto', vmin=-2, vmax=2)

ax_hm.set_xticks(range(n_ct_e))
ax_hm.set_xticklabels([c.replace('\n', ' ') for c in ct_order],
                       fontsize=8.5, rotation=30, ha='right')
ax_hm.set_yticks([])   # pathway names in ax_grp

ax_hm.set_xticks(np.arange(-0.5, n_ct_e, 1), minor=True)
ax_hm.set_yticks(np.arange(-0.5, n_pw_e, 1), minor=True)
ax_hm.grid(which='minor', color='white', linewidth=0.6)
ax_hm.tick_params(which='minor', length=0)
ax_hm.set_ylim(-0.5, n_pw_e - 0.5)

# Highlight OLR1+ TAMs column
olr1_col = ct_order.index('OLR1+ TAMs\n(SPP1+B)')
ax_hm.add_patch(mpatches.Rectangle(
    (olr1_col - 0.5, -0.5), 1, n_pw_e,
    fill=False, ec='#E74C3C', lw=2.5, zorder=5))

# ── Left annotation bar ───────────────────────────────────────────────
ax_grp.set_xlim(0, 1)
ax_grp.set_ylim(-0.5, n_pw_e - 0.5)
ax_grp.axis('off')

for i, pw in enumerate(pw_order):
    color = grp_colors_e.get(pw, '#888888')
    # Colored strip
    ax_grp.add_patch(mpatches.Rectangle(
        (0.55, i - 0.45), 0.4, 0.9,
        fc=color, ec='white', lw=0.5, zorder=2))
    # Pathway name — right of strip, flush against heatmap
    ax_grp.text(0.50, i, pw.replace('\n', ' '), ha='right', va='center',
                fontsize=8, color='#222222')

fig_e.colorbar(im, cax=ax_cb)
ax_cb.set_ylabel('AUCell z-score', fontsize=8)
ax_cb.tick_params(labelsize=7)

fig_e.suptitle(
    'E  Metabolic Pathway Activity Heatmap: OLR1+ TAMs vs Other Cell Types\n'
    '(GSE144735  |  AUCell on KEGG 2021  |  z-score across cell types)',
    fontsize=10, fontweight='bold', x=0.02, ha='left')

fig_e.savefig(OUT / 'panelE_metabolic_heatmap.png', dpi=200, bbox_inches='tight')
fig_e.savefig(OUT / 'panelE_metabolic_heatmap.pdf', bbox_inches='tight')
plt.close()
print("Panel E saved.")

# Save scores
score_df.to_csv(OUT / 'aucell_scores_raw.csv')
score_z.to_csv(OUT / 'aucell_scores_zscore.csv')
print("\nDone. AUCell scores saved to phase4_results/aucell_scores_*.csv")
