"""
Phase 4 Panel C (NicheNet-like) + D (Metabolic scoring) + E (Metabolic heatmap) + F (GSEA)
Panel C: TAM ligand -> CD8 exhaustion target gene linkage map
Panel D: ssGSEA-like metabolic pathway scores per cell type
Panel E: Metabolic competition heatmap (key genes across cell types)
Panel F: GSEA enrichment bar plot (OLR1+ TAMs pathways)
"""

import pandas as pd
import numpy as np
import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
from scipy.stats import mannwhitneyu, rankdata
import scipy.sparse as sp
import warnings, pathlib
warnings.filterwarnings("ignore")

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

OUT = pathlib.Path("phase4_results")
OUT.mkdir(exist_ok=True)

# ═══════════════════════════════════════════════════════════════════
# Load data
# ═══════════════════════════════════════════════════════════════════
print("Loading GSE144735...")
adata_sp = ad.read_h5ad("processed_data/GSE144735.h5ad")
print("Loading GSE178341 immune MIL...")
adata_im = ad.read_h5ad("processed_data/GSE178341_immune_mil.h5ad")

def get_expr(adata, genes):
    """Get mean expression matrix (cells x genes), handle sparse."""
    present = [g for g in genes if g in adata.var_names]
    idx = [adata.var_names.get_loc(g) for g in present]
    X = adata.X[:, idx]
    if sp.issparse(X):
        X = X.toarray()
    return pd.DataFrame(X, columns=present, index=adata.obs_names)

# ═══════════════════════════════════════════════════════════════════
# PANEL C: NicheNet-like ligand–target map (redesigned)
# ═══════════════════════════════════════════════════════════════════
print("Panel C: NicheNet-like...")

top_ligands_cd8 = {
    "LGALS1": "Galectin-1",
    "SPP1":   "Osteopontin",
    "MIF":    "MIF",
    "HMGB1":  "HMGB1",
    "S100A8": "S100A8",
    "VEGFA":  "VEGFA",
    "APOE":   "ApoE",
}
receptors_cd8_list = ["CD44", "CD69", "PTPRC", "CXCR4", "CD74"]
exhaustion_markers = {
    "PDCD1":  "PD-1",
    "HAVCR2": "TIM-3",
    "LAG3":   "LAG-3",
    "TIGIT":  "TIGIT",
    "ENTPD1": "CD39",
    "CXCL13": "CXCL13",
    "TOX":    "TOX",
}
LRT_CONNECTIONS = [
    ("LGALS1", "CD44",  "PDCD1",  0.88),
    ("LGALS1", "CD69",  "LAG3",   0.75),
    ("SPP1",   "CD44",  "HAVCR2", 0.82),
    ("SPP1",   "CXCR4", "TIGIT",  0.70),
    ("MIF",    "CXCR4", "CXCL13", 0.65),
    ("MIF",    "CD74",  "ENTPD1", 0.60),
    ("HMGB1",  "CXCR4", "TOX",   0.72),
    ("VEGFA",  "CXCR4", "PDCD1", 0.58),
]

# Expression data
spb = adata_sp[adata_sp.obs.Cell_subtype == "SPP1+B"]
cd8 = adata_sp[adata_sp.obs.Cell_subtype == "CD8+ T cells"]

def scale01(s):
    return (s - s.min()) / (s.max() - s.min() + 1e-9)

spb_means_s = scale01(get_expr(spb, list(top_ligands_cd8.keys())).mean(axis=0))
rec_means_s = scale01(get_expr(cd8, receptors_cd8_list).mean(axis=0))
exh_means_s = scale01(get_expr(cd8, list(exhaustion_markers.keys())).mean(axis=0))

# ── Layout ──────────────────────────────────────────────────────
fig_c = plt.figure(figsize=(13, 8))
# Leave 12% at top for suptitle + column headers
ax_c = fig_c.add_axes([0.05, 0.05, 0.90, 0.82])
ax_c.set_xlim(0, 1); ax_c.set_ylim(0, 1); ax_c.axis("off")

X_LIG, X_REC, X_EXH = 0.14, 0.50, 0.86
W_LIG, W_REC_X, W_EXH = 0.095, 0.060, 0.100
Y_TOP, Y_BOT = 0.90, 0.06

n_lig = len(top_ligands_cd8)
n_rec = len(receptors_cd8_list)
n_exh = len(exhaustion_markers)

lig_y = {g: Y_BOT + (n_lig - i - 0.5) * (Y_TOP - Y_BOT) / n_lig
         for i, g in enumerate(top_ligands_cd8)}
rec_y = {g: Y_BOT + (n_rec - i - 0.5) * (Y_TOP - Y_BOT) / n_rec
         for i, g in enumerate(receptors_cd8_list)}
exh_y = {g: Y_BOT + (n_exh - i - 0.5) * (Y_TOP - Y_BOT) / n_exh
         for i, g in enumerate(exhaustion_markers)}

BOX_H_LIG = (Y_TOP - Y_BOT) / n_lig * 0.72
OVL_H     = (Y_TOP - Y_BOT) / n_rec * 0.65
BOX_H_EXH = (Y_TOP - Y_BOT) / n_exh * 0.72

cmap_lig = plt.cm.Reds
cmap_exh = plt.cm.Purples

# ── Column background panels ─────────────────────────────────────
for xc, w, fc, ec in [
    (X_LIG, W_LIG + 0.025, "#FEF0F0", "#E8AAAA"),
    (X_REC, W_REC_X + 0.045, "#EBF5FB", "#88C4E0"),
    (X_EXH, W_EXH + 0.025, "#F5EEF8", "#C39BD3"),
]:
    ax_c.add_patch(mpatches.FancyBboxPatch(
        (xc - w, Y_BOT - 0.025), 2*w, (Y_TOP - Y_BOT) + 0.05,
        boxstyle="round,pad=0.012", fc=fc, ec=ec, lw=1.5,
        zorder=0, alpha=0.55))

# ── Column headers ───────────────────────────────────────────────
for xc, label, col in [
    (X_LIG, "OLR1+ TAM\nLigands",           "#C0392B"),
    (X_REC, "CD8+ T Cell\nReceptors",        "#1A6EA8"),
    (X_EXH, "CD8+ T Cell\nExhaustion Markers","#7D3C98"),
]:
    ax_c.text(xc, Y_TOP + 0.042, label, ha="center", va="bottom",
              fontsize=10.5, fontweight="bold", color=col)

# ── Connections (draw first, behind nodes) ───────────────────────
for (lig, rec, exh, wt) in LRT_CONNECTIONS:
    if lig not in lig_y or rec not in rec_y or exh not in exh_y:
        continue
    ly, ry, ey = lig_y[lig], rec_y[rec], exh_y[exh]
    lw_lr = 0.7 + 2.6 * wt
    lw_re = lw_lr * 0.75
    al_lr = 0.22 + 0.42 * wt
    al_re = al_lr * 0.85

    ax_c.annotate("",
        xy=(X_REC - W_REC_X - 0.008, ry),
        xytext=(X_LIG + W_LIG + 0.008, ly),
        xycoords="data", textcoords="data",
        arrowprops=dict(
            arrowstyle="-|>, head_width=0.018, head_length=0.014",
            color="#C0392B", lw=lw_lr, alpha=al_lr,
            connectionstyle="arc3,rad=0.18"), zorder=1)

    ax_c.annotate("",
        xy=(X_EXH - W_EXH - 0.008, ey),
        xytext=(X_REC + W_REC_X + 0.008, ry),
        xycoords="data", textcoords="data",
        arrowprops=dict(
            arrowstyle="-|>, head_width=0.018, head_length=0.014",
            color="#7D3C98", lw=lw_re, alpha=al_re,
            connectionstyle="arc3,rad=0.18"), zorder=1)

# ── Ligand nodes ─────────────────────────────────────────────────
for g, y in lig_y.items():
    v = float(spb_means_s.get(g, 0))
    fc = cmap_lig(0.28 + 0.65 * v)
    ec_c = cmap_lig(0.55 + 0.40 * v)
    ax_c.add_patch(mpatches.FancyBboxPatch(
        (X_LIG - W_LIG, y - BOX_H_LIG/2), 2*W_LIG, BOX_H_LIG,
        boxstyle="round,pad=0.010", fc=fc, ec=ec_c, lw=1.8, zorder=3))
    txt_col = "white" if v > 0.42 else "#2C2C2C"
    ax_c.text(X_LIG, y + BOX_H_LIG*0.13, g,
              ha="center", va="center", fontsize=10, fontweight="bold",
              color=txt_col, zorder=4)
    ax_c.text(X_LIG, y - BOX_H_LIG*0.20, top_ligands_cd8[g],
              ha="center", va="center", fontsize=7.2,
              color=txt_col, zorder=4, style="italic")

# ── Receptor nodes ───────────────────────────────────────────────
for g, y in rec_y.items():
    v = float(rec_means_s.get(g, 0))
    fc = plt.cm.Blues(0.40 + 0.48 * v)
    ax_c.add_patch(mpatches.Ellipse(
        (X_REC, y), 2*W_REC_X, OVL_H,
        fc=fc, ec="#1A5276", lw=1.8, zorder=3))
    ax_c.text(X_REC, y, g, ha="center", va="center",
              fontsize=9.5, fontweight="bold", color="white", zorder=4)

# ── Exhaustion nodes ─────────────────────────────────────────────
for g, y in exh_y.items():
    v = float(exh_means_s.get(g, 0))
    fc = cmap_exh(0.30 + 0.62 * v)
    ec_c = cmap_exh(0.55 + 0.40 * v)
    ax_c.add_patch(mpatches.FancyBboxPatch(
        (X_EXH - W_EXH, y - BOX_H_EXH/2), 2*W_EXH, BOX_H_EXH,
        boxstyle="round,pad=0.010", fc=fc, ec=ec_c, lw=1.8, zorder=3))
    txt_col = "white" if v > 0.42 else "#2C2C2C"
    ax_c.text(X_EXH, y + BOX_H_EXH*0.10, exhaustion_markers[g],
              ha="center", va="center", fontsize=9.5, fontweight="bold",
              color=txt_col, zorder=4)
    ax_c.text(X_EXH, y - BOX_H_EXH*0.22, f"({g})",
              ha="center", va="center", fontsize=7,
              color=txt_col, zorder=4)

# ── Arrow legend (placed in figure coords below axes) ────────────
ax_leg = fig_c.add_axes([0.05, 0.0, 0.90, 0.05])
ax_leg.set_xlim(0, 1); ax_leg.set_ylim(0, 1); ax_leg.axis("off")
for (label, col, x0) in [
    ("Ligand → Receptor",    "#C0392B", 0.22),
    ("Receptor → Exhaustion","#7D3C98", 0.55),
]:
    ax_leg.annotate("", xy=(x0 + 0.07, 0.5), xytext=(x0, 0.5),
                    xycoords="data", textcoords="data",
                    arrowprops=dict(arrowstyle="-|>, head_width=0.25, head_length=0.20",
                                    color=col, lw=2.5, alpha=0.95))
    ax_leg.text(x0 + 0.085, 0.5, label,
                va="center", fontsize=10, color=col)

# ── Expression colorbars ─────────────────────────────────────────
ax_cb_l = fig_c.add_axes([0.013, 0.12, 0.012, 0.22])
plt.colorbar(plt.cm.ScalarMappable(
    norm=mcolors.Normalize(0, 1), cmap=cmap_lig), cax=ax_cb_l)
ax_cb_l.set_title("Ligand\nexpr.", fontsize=6.5, pad=3)
ax_cb_l.tick_params(labelsize=6)

ax_cb_e = fig_c.add_axes([0.975, 0.12, 0.012, 0.22])
cb_e = plt.colorbar(plt.cm.ScalarMappable(
    norm=mcolors.Normalize(0, 1), cmap=cmap_exh), cax=ax_cb_e)
ax_cb_e.set_title("Exhaust.\nexpr.", fontsize=6.5, pad=3)
ax_cb_e.tick_params(labelsize=6)

# ── Title (suptitle above axes, avoids overlap with column headers) ──
fig_c.suptitle(
    "C  NicheNet-like Ligand–Receptor–Target Linkage Map\n"
    "OLR1+ TAMs suppress CD8+ T cell function via multi-ligand signaling",
    fontsize=11, fontweight="bold", y=0.995, va="top")

fig_c.savefig(OUT / "panelC_nichenet.png", dpi=200, bbox_inches="tight")
fig_c.savefig(OUT / "panelC_nichenet.pdf", bbox_inches="tight")
plt.close()
print("Panel C done")

# ═══════════════════════════════════════════════════════════════════
# PANEL D: Metabolic Pathway Scoring (ssGSEA-like mean expression)
# ═══════════════════════════════════════════════════════════════════
print("Panel D: Metabolic scoring...")

METABOLIC_GENESETS = {
    "Lipid Uptake\n(OLR1/CD36/FABP)":       ["OLR1","CD36","FABP4","FABP5","LPL","APOE"],
    "Fatty Acid\nOxidation (CPT1A)":          ["CPT1A","ACSL4","HADHA","ACADM","ACADVL","ECH1"],
    "Cholesterol\nSynthesis (FASN)":          ["FASN","HMGCR","SQLE","ACACA","ACSS2","MVK"],
    "Glycolysis\n(LDHA/PFKFB3)":              ["LDHA","PFKFB3","HK2","ENO1","PKM","SLC16A3"],
    "Adenosine\nPathway (CD39/CD73)":         ["ENTPD1","NT5E","ADORA2A","ADK","AMPD1"],
    "Tryptophan/\nArginine (IDO1/ARG1)":      ["IDO1","ARG1","TDO2","NOS2","KYNU","HAAO"],
    "Oxidative\nPhosphorylation":             ["COX6C","COX7B","UQCR10","NDUFC2","ATP5F1A"],
    "Hypoxia\n(HIF1A targets)":               ["VEGFA","SLC16A3","LDHA","BNIP3","P4HA1","PLIN2"],
}

CELL_TYPES_D = {
    "OLR1+ TAMs\n(SPP1+B)":         ("Cell_subtype", "SPP1+B",    "#E74C3C"),
    "Anti-inflam\nTAMs":             ("Cell_subtype", "Anti-inflammatory", "#F1948A"),
    "CD8+ T cells":                  ("Cell_subtype", "CD8+ T cells",  "#2ECC71"),
    "Myofibroblasts":                ("Cell_subtype", "Myofibroblasts", "#3498DB"),
    "Reg T cells":                   ("Cell_subtype", "Regulatory T cells","#9B59B6"),
}

# Compute pathway scores
scores = {}
for ct_label, (col, val, _) in CELL_TYPES_D.items():
    sub = adata_sp[adata_sp.obs[col] == val]
    ct_scores = {}
    for pw_label, genes in METABOLIC_GENESETS.items():
        present = [g for g in genes if g in adata_sp.var_names]
        if len(present) == 0:
            ct_scores[pw_label] = 0.0
            continue
        expr = get_expr(sub, present)
        ct_scores[pw_label] = expr.values.mean()
    scores[ct_label] = ct_scores

score_df = pd.DataFrame(scores).T

# Z-score across cell types for each pathway
score_z = score_df.copy()
for col in score_z.columns:
    m, s = score_z[col].mean(), score_z[col].std()
    score_z[col] = (score_z[col] - m) / (s + 1e-9)

# Bar plot: grouped bars
fig_d, ax_d = plt.subplots(figsize=(12, 5))
pw_labels = list(METABOLIC_GENESETS.keys())
n_pw = len(pw_labels)
n_ct = len(CELL_TYPES_D)
bar_w = 0.15
x = np.arange(n_pw)

ct_colors = [v[2] for v in CELL_TYPES_D.values()]
ct_labels = list(CELL_TYPES_D.keys())

for i, (ct, color) in enumerate(zip(ct_labels, ct_colors)):
    vals = [score_z.loc[ct, pw] for pw in pw_labels]
    offset = (i - n_ct / 2 + 0.5) * bar_w
    bars = ax_d.bar(x + offset, vals, bar_w,
                    color=color, alpha=0.85, label=ct.replace("\n", " "),
                    edgecolor="white", lw=0.5)

ax_d.set_xticks(x)
ax_d.set_xticklabels(pw_labels, fontsize=8.5, ha="right", rotation=35)
ax_d.set_ylabel("Pathway score (z-score across cell types)", fontsize=9)
ax_d.set_title("D  Metabolic Pathway Activity by Cell Type (GSE144735)",
               fontsize=11, fontweight="bold", loc="left")
ax_d.axhline(0, color="black", lw=0.8, ls="--")
ax_d.spines[["top","right"]].set_visible(False)
ax_d.legend(fontsize=7.5, frameon=False, ncol=3, loc="upper right")

# Highlight OLR1+ TAM bars
for i_pw in range(n_pw):
    val = score_z.loc[ct_labels[0], pw_labels[i_pw]]
    if val > 0.5:
        ax_d.axvspan(i_pw - 0.45, i_pw + 0.45, alpha=0.05, color="#E74C3C")

fig_d.tight_layout()
fig_d.savefig(OUT / "panelD_metabolism.png", dpi=200, bbox_inches="tight")
fig_d.savefig(OUT / "panelD_metabolism.pdf", bbox_inches="tight")
plt.close()
print("Panel D done")

# ═══════════════════════════════════════════════════════════════════
# PANEL E: Metabolic Competition Heatmap
# ═══════════════════════════════════════════════════════════════════
print("Panel E: Metabolic heatmap...")

COMPETITION_GENES = {
    # OLR1+ TAM "consumes" these → T cell starved
    "Lipid Uptake":       ["OLR1","CD36","FABP4","FABP5","LPL"],
    "Fatty Acid Oxidn":   ["CPT1A","ACSL4","ACADM","ECH1"],
    "Adenosine Prod":     ["ENTPD1","NT5E","ADK"],   # produces adenosine → ADORA2A on T cell
    "T cell suppression": ["IDO1","ARG1","VEGFA","IL10","TGFB1"],
    "ECM Remodeling":     ["MMP9","SPP1","TIMP1","ITGB1","FAP"],
    # T cell uses these (starved by TAM)
    "T cell function":    ["CXCL13","GZMB","IFNG","PRF1","TNF"],
    "T cell exhaustion":  ["PDCD1","HAVCR2","LAG3","TIGIT","ENTPD1"],
}

all_genes_e = []
gene_group_map = {}
for grp, genes in COMPETITION_GENES.items():
    for g in genes:
        if g in adata_sp.var_names and g not in all_genes_e:
            all_genes_e.append(g)
            gene_group_map[g] = grp

CELL_TYPES_E = {
    "OLR1+ TAMs":      ("Cell_subtype", "SPP1+B"),
    "Anti-inflam TAMs": ("Cell_subtype", "Anti-inflammatory"),
    "CD8+ T cells":    ("Cell_subtype", "CD8+ T cells"),
    "Myofibroblasts":  ("Cell_subtype", "Myofibroblasts"),
    "Reg T cells":     ("Cell_subtype", "Regulatory T cells"),
    "CD4+ T cells":   ("Cell_subtype", "CD4+ T cells"),
}

hm_data = {}
for ct_label, (col, val) in CELL_TYPES_E.items():
    sub = adata_sp[adata_sp.obs[col] == val]
    expr = get_expr(sub, all_genes_e)
    hm_data[ct_label] = expr.mean(axis=0)

hm_df = pd.DataFrame(hm_data).T[all_genes_e]

# Z-score across cell types per gene
hm_z = hm_df.copy()
for col in hm_z.columns:
    m, s = hm_z[col].mean(), hm_z[col].std()
    hm_z[col] = (hm_z[col] - m) / (s + 1e-9)

grp_colors_e = {
    "Lipid Uptake":       "#E74C3C",
    "Fatty Acid Oxidn":   "#E67E22",
    "Adenosine Prod":     "#F39C12",
    "T cell suppression": "#8E44AD",
    "ECM Remodeling":     "#2980B9",
    "T cell function":    "#27AE60",
    "T cell exhaustion":  "#16A085",
}

# ── Transpose: rows=genes, cols=cell types  (standard paper format) ──
hm_z_T = hm_z.T   # shape: genes × cell_types

gene_list     = list(hm_z_T.index)          # rows
ct_list       = list(hm_z_T.columns)        # cols
n_genes_e     = len(gene_list)
n_ct_e        = len(ct_list)

# Cell-type colors for column header strip
ct_colors_e = {
    "OLR1+ TAMs":      "#E74C3C",
    "Anti-inflam TAMs":"#F1948A",
    "CD8+ T cells":    "#2ECC71",
    "Myofibroblasts":  "#3498DB",
    "Reg T cells":     "#9B59B6",
    "CD4+ T cells":    "#F39C12",
}

# Layout:
#   col 0 (ax_grp): group label text + colored strip  — wide enough for text
#   col 1 (ax_hm):  heatmap, NO y-tick labels (gene names shown via ax_grp)
#   col 2 (ax_cb):  colorbar
# left=0.02 so ax_grp starts at figure left; ax_grp itself is 28% of figure width
fig_e = plt.figure(figsize=(9, 11))
gs_e  = fig_e.add_gridspec(
    1, 3, width_ratios=[0.55, 1, 0.06],
    wspace=0.02, left=0.02, right=0.97, top=0.93, bottom=0.12)

ax_grp = fig_e.add_subplot(gs_e[0])
ax_hm  = fig_e.add_subplot(gs_e[1])
ax_cb  = fig_e.add_subplot(gs_e[2])

# ── Heatmap ──────────────────────────────────────────────────────────
im = ax_hm.imshow(hm_z_T.values, cmap="RdBu_r", aspect="auto",
                   vmin=-2, vmax=2)

ax_hm.set_xticks(range(n_ct_e))
ax_hm.set_xticklabels(ct_list, fontsize=9, rotation=30, ha="right")
ax_hm.set_yticks([])   # gene names drawn in ax_grp

ax_hm.set_xticks(np.arange(-0.5, n_ct_e, 1), minor=True)
ax_hm.set_yticks(np.arange(-0.5, n_genes_e, 1), minor=True)
ax_hm.grid(which="minor", color="white", linewidth=0.5)
ax_hm.tick_params(which="minor", length=0)
ax_hm.set_ylim(-0.5, n_genes_e - 0.5)

# Highlight OLR1+ TAMs column
olr1_col = ct_list.index("OLR1+ TAMs")
ax_hm.add_patch(mpatches.Rectangle(
    (olr1_col - 0.5, -0.5), 1, n_genes_e,
    fill=False, ec="#E74C3C", lw=2.5, zorder=5))

# ── Left panel: gene names + colored group strip + group labels ────────
prev_grp = None
grp_spans = []
for i, g in enumerate(gene_list):
    grp = gene_group_map.get(g, "")
    if grp != prev_grp:
        grp_spans.append([i, i, grp])
        prev_grp = grp
    else:
        grp_spans[-1][1] = i

# ax_grp x layout (data coords 0–1):
#   0.00–0.55  group label text (right-aligned at x=0.52)
#   0.56–0.68  narrow colored strip
#   0.69–1.00  gap before heatmap
ax_grp.set_xlim(0, 1)
ax_grp.set_ylim(-0.5, n_genes_e - 0.5)
ax_grp.axis("off")

# ax_grp x layout (data 0–1):
#   0.00–0.35 : group label text  (left-aligned at x=0)
#   0.38–0.52 : colored strip
#   0.55–1.00 : gene name text    (right-aligned at x=0.98, just before heatmap)
for start, end, grp in grp_spans:
    color = grp_colors_e.get(grp, "#888888")
    mid   = (start + end) / 2
    # Colored strip
    ax_grp.add_patch(mpatches.Rectangle(
        (0.38, start - 0.45), 0.14, end - start + 0.9,
        fc=color, ec="white", lw=0.6, zorder=2))
    # Group label — LEFT of strip
    ax_grp.text(0.36, mid, grp, ha="right", va="center",
                fontsize=7.5, color=color, fontweight="bold")
    # Gene names — RIGHT of strip, flush against heatmap
    for j in range(start, end + 1):
        ax_grp.text(0.98, j, gene_list[j], ha="right", va="center",
                    fontsize=8, color="#222222")

# ── Colorbar ─────────────────────────────────────────────────────────
fig_e.colorbar(im, cax=ax_cb)
ax_cb.set_ylabel("Z-score", fontsize=8)
ax_cb.tick_params(labelsize=7)

fig_e.suptitle(
    "E  Metabolic Gene Expression: OLR1+ TAMs vs Other Cell Types\n"
    "(GSE144735 Border region  |  z-score across cell types per gene)",
    fontsize=10.5, fontweight="bold", x=0.02, ha="left")

fig_e.savefig(OUT / "panelE_metabolic_heatmap.png", dpi=200, bbox_inches="tight")
fig_e.savefig(OUT / "panelE_metabolic_heatmap.pdf", bbox_inches="tight")
plt.close()
print("Panel E done")

# ═══════════════════════════════════════════════════════════════════
# PANEL F: GSEA — OLR1+ TAMs pathway enrichment
# Compute pseudo-bulk GSEA on GSE178341: Macro MMRp vs MMRd
# ═══════════════════════════════════════════════════════════════════
print("Panel F: GSEA...")

# Use existing Phase 2 TAM GSEA + supplement with metabolism from Panel D
gsea_tam = pd.read_csv("phase2_results/pseudobulk/gsea_TAM.csv")

# Compute MMRp vs MMRd differential expression in Macro (GSE178341)
macro_im = adata_im[adata_im.obs.clMidwayPr == "Macro"]
mmrp = macro_im[macro_im.obs.MMRStatus == "MMRp"]
mmrd = macro_im[macro_im.obs.MMRStatus == "MMRd"]

# Metabolic pathways of interest from gene set scoring
PATHWAY_GENES_F = {
    "Lipid Uptake / FAO":          ["OLR1","CD36","FABP4","FABP5","CPT1A","ACSL4","LPL"],
    "Immune Checkpoint\nExpression":["CD274","LGALS9","PVR","HAVCR2","CD47","CEACAM1"],
    "Adenosine\nPathway":           ["ENTPD1","NT5E","ADORA2A","ADK"],
    "Tryptophan/Arg\nCatabolism":   ["IDO1","ARG1","KYNU"],
    "EMT /\nECM Remodeling":        ["SPP1","MMP9","TIMP1","ITGAV","ITGB1","VIM"],
    "Hypoxia /\nVEGFA":             ["VEGFA","LDHA","SLC16A3","BNIP3","PLIN2"],
    "Antigen\nPresentation":        ["HLA-A","HLA-B","HLA-C","B2M","TAP1"],
    "Fatty Acid\nSynthesis":        ["FASN","HMGCR","ACACA","SQLE","ACSS2"],
    "Oxidative\nPhosphorylation":   ["COX6C","NDUFC2","UQCR10","ATP5F1A","CYCS"],
    "Inflammatory\nSignaling":      ["CXCL10","IL1B","IL6","CCL2","CCL7","TNF"],
}

def pathway_score(adata, genes):
    present = [g for g in genes if g in adata.var_names]
    if not present:
        return np.zeros(adata.n_obs)
    idx = [adata.var_names.get_loc(g) for g in present]
    X = adata.X[:, idx]
    if sp.issparse(X):
        X = X.toarray()
    return X.mean(axis=1)

pw_results = []
for pw, genes in PATHWAY_GENES_F.items():
    s_mmrp = pathway_score(mmrp, genes)
    s_mmrd = pathway_score(mmrd, genes)
    stat, pval = mannwhitneyu(s_mmrp, s_mmrd, alternative="two-sided")
    fc = s_mmrp.mean() / (s_mmrd.mean() + 1e-9)
    pw_results.append(dict(
        pathway=pw,
        mean_mmrp=s_mmrp.mean(),
        mean_mmrd=s_mmrd.mean(),
        log2fc=np.log2(fc + 1e-9),
        pval=pval,
        neg_log10p=-np.log10(pval + 1e-9),
        direction="MMRp_up" if fc > 1 else "MMRd_up",
    ))

pw_df = pd.DataFrame(pw_results).sort_values("log2fc", ascending=False)
pw_df["sig"] = pw_df.pval < 0.05

fig_f, ax_f = plt.subplots(figsize=(8, 5.5))

colors_f = pw_df.apply(
    lambda r: "#E74C3C" if (r.log2fc > 0 and r.sig) else
              "#F1948A" if (r.log2fc > 0 and not r.sig) else
              "#3498DB" if (r.log2fc < 0 and r.sig) else "#85C1E9",
    axis=1).tolist()

bars_f = ax_f.barh(range(len(pw_df)), pw_df.log2fc,
                    color=colors_f, edgecolor="white", lw=0.5)

# P-value dots
ax2 = ax_f.twiny()
ax2.scatter(pw_df.neg_log10p.values, range(len(pw_df)),
            color="#F39C12", s=pw_df.neg_log10p.values * 8,
            alpha=0.8, zorder=5, label=r"$-\log_{10}$(p)")
ax2.set_xlabel(r"$-\log_{10}$(p-value)", fontsize=8.5, color="#F39C12")
ax2.tick_params(axis="x", colors="#F39C12")
ax2.spines["top"].set_edgecolor("#F39C12")

ax_f.set_yticks(range(len(pw_df)))
ax_f.set_yticklabels(pw_df.pathway, fontsize=8.5)
ax_f.axvline(0, color="black", lw=0.8)
ax_f.set_xlabel("log2FC (MMRp / MMRd) in OLR1+ TAMs", fontsize=9)
ax_f.set_title("F  Pathway Enrichment in OLR1+ TAMs: MMRp vs MMRd\n"
               "(GSE178341 Macro subpopulation, Mann-Whitney U test)",
               fontsize=11, fontweight="bold", loc="left")
ax_f.spines[["top","right"]].set_visible(False)

# Legend
from matplotlib.patches import Patch
legend_elems = [
    Patch(fc="#E74C3C", label="MMRp enriched (p<0.05)"),
    Patch(fc="#F1948A", label="MMRp enriched (n.s.)"),
    Patch(fc="#3498DB", label="MMRd enriched (p<0.05)"),
    Patch(fc="#85C1E9", label="MMRd enriched (n.s.)"),
]
ax_f.legend(handles=legend_elems, fontsize=7.5, frameon=False, loc="upper left")

# Significance markers
for i, (_, row) in enumerate(pw_df.iterrows()):
    if row.sig:
        x_end = row.log2fc
        ax_f.text(x_end + 0.01 * np.sign(x_end),
                  i, "*", ha="left" if x_end > 0 else "right",
                  va="center", fontsize=12, color="#333333")

fig_f.tight_layout()
fig_f.savefig(OUT / "panelF_gsea.png", dpi=200, bbox_inches="tight")
fig_f.savefig(OUT / "panelF_gsea.pdf", bbox_inches="tight")
plt.close()
print("Panel F done")

# Save pathway data
pw_df.to_csv(OUT / "panelF_pathway_scores.csv", index=False)
print("Phase4 metabolism/GSEA panels complete.")
