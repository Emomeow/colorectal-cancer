"""
Phase 4 Panel A + B: CellChat / LIANA 通讯网络
Panel A: OLR1+ TAM 与关键亚群的通讯强度 Circle Plot
Panel B: Inhibitory + SPP1→Integrin LR pairs Bubble Plot
数据源: phase3_results/cellchat/liana_border_all.csv (GSE144735 Border cells)
"""

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import matplotlib.patheffects as pe
from scipy.stats import mannwhitneyu
import warnings, pathlib
warnings.filterwarnings("ignore")

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

OUT = pathlib.Path("phase4_results")
OUT.mkdir(exist_ok=True)

# ── 读 LIANA ─────────────────────────────────────────────────────
print("Loading LIANA...")
df = pd.read_csv("phase3_results/cellchat/liana_border_all.csv")
# magnitude_rank: 0=strongest, 1=weakest; cellphone_pvals < 0.05 = significant
sig = df[df.cellphone_pvals < 0.05].copy()
# communication strength = -log10(magnitude_rank + 1e-6)
sig["strength"] = -np.log10(sig["magnitude_rank"] + 1e-6)

# ── Cell type mapping ─────────────────────────────────────────────
KEY_TYPES = {
    "OLR1+ TAMs": ["SPP1+B"],
    "CD8+ T cells": ["CD8+ T cells"],
    "Myofibroblasts": ["Myofibroblasts"],
    "Reg T cells": ["Regulatory T cells"],
    "CD4+ T cells": ["CD4+ T cells"],
    "Anti-inflamMAC": ["Anti-inflammatory", "Pro-inflammatory"],
    "Stromal/Other": ["Stromal 1", "Stromal 2", "Stromal 3",
                       "Pericytes", "Tip-like ECs", "Stalk-like ECs",
                       "Smooth muscle cells", "Lymphatic ECs"],
}

def map_ct(ct):
    for k, vs in KEY_TYPES.items():
        if ct in vs:
            return k
    return None

sig["src_group"] = sig["source"].apply(map_ct)
sig["tgt_group"] = sig["target"].apply(map_ct)
sig_filtered = sig.dropna(subset=["src_group","tgt_group"]).copy()
sig_filtered = sig_filtered[sig_filtered.src_group != sig_filtered.tgt_group]

# Aggregate: mean strength per src-tgt pair
comm = (sig_filtered.groupby(["src_group","tgt_group"])
        .agg(mean_strength=("strength","mean"),
             n_pairs=("strength","count"))
        .reset_index())

print("Communication matrix computed, n_rows:", len(comm))

# ═══════════════════════════════════════════════════════════════════
# PANEL A: Circle Plot
# ═══════════════════════════════════════════════════════════════════

GROUPS = list(KEY_TYPES.keys())
N = len(GROUPS)
COLORS = {
    "OLR1+ TAMs":     "#E74C3C",
    "CD8+ T cells":   "#2ECC71",
    "Myofibroblasts": "#3498DB",
    "Reg T cells":     "#9B59B6",
    "CD4+ T cells":   "#F39C12",
    "Anti-inflamMAC": "#1ABC9C",
    "Stromal/Other":  "#95A5A6",
}

def polar_pos(i, n, r=1.0):
    angle = 2 * np.pi * i / n - np.pi / 2
    return r * np.cos(angle), r * np.sin(angle)

fig_a, ax_a = plt.subplots(figsize=(7, 7))
ax_a.set_aspect("equal")
ax_a.axis("off")

# Draw nodes
node_r = 0.12
for i, g in enumerate(GROUPS):
    x, y = polar_pos(i, N, r=1.0)
    circ = plt.Circle((x, y), node_r, color=COLORS[g], zorder=4, ec="white", lw=1.5)
    ax_a.add_patch(circ)
    # Label
    angle_deg = np.degrees(2 * np.pi * i / N - np.pi / 2)
    ha = "left" if -90 < angle_deg < 90 else "right"
    lx = x * 1.22
    ly = y * 1.22
    ax_a.text(lx, ly, g, ha=ha, va="center", fontsize=8.5, fontweight="bold",
              color=COLORS[g])

# Node size by total outgoing strength
out_strength = {}
for g in GROUPS:
    sub = comm[comm.src_group == g]
    out_strength[g] = sub.mean_strength.sum() if len(sub) > 0 else 0.1

# Draw edges
max_s = comm.mean_strength.max()
for _, row in comm.iterrows():
    if row.mean_strength < 0.30:
        continue
    src = row.src_group
    tgt = row.tgt_group
    si = GROUPS.index(src)
    ti = GROUPS.index(tgt)
    x0, y0 = polar_pos(si, N)
    x1, y1 = polar_pos(ti, N)
    lw = 0.5 + 3.5 * row.mean_strength / max_s
    alpha = 0.25 + 0.55 * row.mean_strength / max_s
    # Draw arc (bezier via ConnectionPatch)
    ax_a.annotate("", xy=(x1, y1), xytext=(x0, y0),
                  arrowprops=dict(
                      arrowstyle="-|>, head_width=0.25, head_length=0.15",
                      color=COLORS.get(src, "#888888"),
                      lw=lw, alpha=max(alpha, 0.7),
                      connectionstyle="arc3,rad=0.25",
                      shrinkA=22, shrinkB=22,
                  ),
                  zorder=3)

# OLR1+ TAM highlight ring
olr1_i = GROUPS.index("OLR1+ TAMs")
x0, y0 = polar_pos(olr1_i, N)
ring = plt.Circle((x0, y0), node_r + 0.04, fill=False, ec="#E74C3C",
                   lw=2.5, linestyle="--", zorder=3, alpha=0.7)
ax_a.add_patch(ring)

ax_a.set_xlim(-1.55, 1.55)
ax_a.set_ylim(-1.55, 1.55)
ax_a.set_title("A  Cell-Cell Communication Network (GSE144735 Border)",
               fontsize=11, fontweight="bold", pad=8)

# Legend: edge width
for lw_val, label in [(1.0, "Weak"), (2.5, "Medium"), (4.0, "Strong")]:
    ax_a.plot([], [], lw=lw_val, color="#888888", alpha=0.6, label=label)
ax_a.legend(title="Signal strength", fontsize=7.5, title_fontsize=8,
            loc="lower right", frameon=False)

fig_a.tight_layout()
fig_a.savefig(OUT / "panelA_circle_plot.png", dpi=200, bbox_inches="tight")
fig_a.savefig(OUT / "panelA_circle_plot.pdf", bbox_inches="tight")
plt.close()
print("Panel A done")

# ═══════════════════════════════════════════════════════════════════
# PANEL B: Bubble Plot - Key LR pairs from OLR1+ TAMs
# ═══════════════════════════════════════════════════════════════════

# Define key LR pairs to highlight
KEY_PAIRS = [
    # Inhibitory: TAM → CD8 (immune suppression)
    ("SPP1+B", "CD8+ T cells",  "LGALS1",  "CD69",       "Immune suppression"),
    ("SPP1+B", "CD8+ T cells",  "LGALS1",  "PTPRC",      "Immune suppression"),
    ("SPP1+B", "CD8+ T cells",  "SPP1",    "CD44",       "Immune suppression"),
    ("SPP1+B", "CD8+ T cells",  "MIF",     "CD74_CXCR4", "Immune suppression"),
    ("SPP1+B", "CD8+ T cells",  "HMGB1",   "CXCR4",      "Immune suppression"),
    ("SPP1+B", "CD8+ T cells",  "S100A8",  "CD69",       "Immune suppression"),
    # TAM → Treg (recruitment/activation)
    ("SPP1+B", "Regulatory T cells", "LGALS1", "CD69",   "Treg crosstalk"),
    ("SPP1+B", "Regulatory T cells", "LGALS1", "PTPRC",  "Treg crosstalk"),
    ("SPP1+B", "Regulatory T cells", "SPP1",   "CD44",   "Treg crosstalk"),
    ("SPP1+B", "Regulatory T cells", "MIF",    "CD74_CXCR4", "Treg crosstalk"),
    # TAM → Myofibroblast (structural barrier)
    ("SPP1+B", "Myofibroblasts", "SPP1",   "ITGAV_ITGB1","Matrix remodeling"),
    ("SPP1+B", "Myofibroblasts", "SPP1",   "CD44",       "Matrix remodeling"),
    ("SPP1+B", "Myofibroblasts", "SPP1",   "ITGAV_ITGB5","Matrix remodeling"),
    ("SPP1+B", "Myofibroblasts", "APOE",   "LRP1",       "Matrix remodeling"),
    ("SPP1+B", "Myofibroblasts", "LGALS1", "ITGB1",      "Matrix remodeling"),
    ("SPP1+B", "Myofibroblasts", "TIMP1",  "CD63",       "Matrix remodeling"),
]

cat_colors = {
    "Immune suppression": "#E74C3C",
    "Treg crosstalk":     "#9B59B6",
    "Matrix remodeling":  "#3498DB",
}

rows_b = []
for (src, tgt, lig, rec, cat) in KEY_PAIRS:
    sub = df[(df.source == src) & (df.target == tgt) &
             (df.ligand_complex == lig) & (df.receptor_complex == rec)]
    if len(sub) == 0:
        continue
    r = sub.iloc[0]
    rows_b.append(dict(
        label=f"{lig} -> {rec}",
        target_group=tgt.replace("Regulatory T cells","Reg T cells"),
        category=cat,
        lr_means=r.lr_means,
        pval=r.cellphone_pvals,
        magnitude_rank=r.magnitude_rank,
        specificity_rank=r.specificity_rank,
        strength=-np.log10(r.magnitude_rank + 1e-6),
    ))

df_b = pd.DataFrame(rows_b)

# Add SPP1+A as comparison for SPP1-related pairs (specificity check)
# Filter to remove very low-strength
df_b = df_b[df_b.lr_means > 1.3].copy()
df_b = df_b.sort_values(["category","strength"], ascending=[True, False])

fig_b, ax_b = plt.subplots(figsize=(9, 6))

cat_order = ["Immune suppression", "Treg crosstalk", "Matrix remodeling"]
y_labels = []
y_pos = []
cur_y = 0
cat_boundaries = []

for cat in cat_order:
    sub = df_b[df_b.category == cat].reset_index(drop=True)
    if len(sub) == 0:
        continue
    cat_start = cur_y
    for idx, row in sub.iterrows():
        label = f"{row.label}  [{row.target_group}]"
        y_labels.append(label)
        y_pos.append(cur_y)
        # Bubble: x = lr_means, size = -log10(magnitude_rank)
        size = max(30, min(300, row.strength * 40))
        ax_b.scatter(row.lr_means, cur_y,
                     s=size, c=cat_colors[cat], alpha=0.85,
                     edgecolors="white", lw=0.8, zorder=3)
        cur_y += 1
    cat_boundaries.append((cat, cat_start, cur_y - 1))
    cur_y += 0.6  # gap between categories

# Category background bands
for cat, y0, y1 in cat_boundaries:
    ax_b.axhspan(y0 - 0.4, y1 + 0.4, alpha=0.06,
                 color=cat_colors[cat], zorder=0)
    ax_b.text(-0.18, (y0 + y1) / 2, cat,
              ha="right", va="center", fontsize=8.5,
              color=cat_colors[cat], fontweight="bold",
              transform=ax_b.get_yaxis_transform())

ax_b.set_yticks(y_pos)
ax_b.set_yticklabels(y_labels, fontsize=8)
ax_b.set_xlabel("Mean LR expression (LIANA lr_means)", fontsize=9)
ax_b.set_title("B  Key Ligand–Receptor Pairs from OLR1+ TAMs\n"
               "(GSE144735 Border region, SPP1+B as sender)",
               fontsize=11, fontweight="bold", loc="left")
ax_b.spines[["top","right","left"]].set_visible(False)
ax_b.axvline(1.5, ls="--", color="#AAAAAA", lw=0.8, label="LR mean = 1.5")

# Bubble size legend
for s_val, label in [(2, "Low"), (5, "Medium"), (8, "High")]:
    ax_b.scatter([], [], s=s_val*40, c="#999999", alpha=0.7,
                 label=f"Strength: {label} (-log10 rank)")
ax_b.legend(fontsize=7.5, frameon=False, loc="upper right")

ax_b.set_xlim(ax_b.get_xlim()[0] - 0.3, ax_b.get_xlim()[1] + 0.3)
fig_b.tight_layout()
fig_b.savefig(OUT / "panelB_bubble_plot.png", dpi=200, bbox_inches="tight")
fig_b.savefig(OUT / "panelB_bubble_plot.pdf", bbox_inches="tight")
plt.close()
print("Panel B done")
print("Phase4 CellChat panels complete.")
