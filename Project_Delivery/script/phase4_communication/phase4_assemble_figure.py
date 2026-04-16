"""
Phase 4: Assemble Figure 4 (Panels A-F) and Supp Fig 5
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pathlib, numpy as np
import matplotlib.patches as mpatches

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

SRC = pathlib.Path("phase4_results")
OUT_DIRS = [
    pathlib.Path("figures/phase4"),
    pathlib.Path("Project_Delivery/results/Figure4_ImmunoSuppression_Mechanism/figures"),
    pathlib.Path("Project_Delivery/results/Supplementary/Supp5_Metabolism_Supp/figures"),
]
for d in OUT_DIRS:
    d.mkdir(parents=True, exist_ok=True)

def load(name):
    p = SRC / name
    if not p.exists():
        print(f"[WARN] missing: {p}")
        return None
    return mpimg.imread(str(p))

# ── Figure 4 (2x3 grid) ──────────────────────────────────────────
fig4, axes = plt.subplots(2, 3, figsize=(18, 12))
fig4.subplots_adjust(wspace=0.06, hspace=0.08,
                     left=0.02, right=0.98, top=0.94, bottom=0.02)

panels = [
    ("panelA_circle_plot.png",     "A"),
    ("panelB_bubble_plot.png",     "B"),
    ("panelC_nichenet.png",        "C"),
    ("panelD_metabolism.png",      "D"),
    ("panelE_metabolic_heatmap.png","E"),
    ("panelF_gsea.png",            "F"),
]

for ax, (fname, label) in zip(axes.flat, panels):
    img = load(fname)
    ax.axis("off")
    if img is not None:
        ax.imshow(img)
    ax.text(0.01, 0.99, label, transform=ax.transAxes,
            fontsize=16, fontweight="bold", va="top", ha="left",
            color="black")

fig4.suptitle(
    "Figure 4  —  OLR1⁺ TAMs Induce Immune Evasion: Communication & Metabolic Mechanisms",
    fontsize=14, fontweight="bold", y=0.98)

for out in [OUT_DIRS[0], OUT_DIRS[1]]:
    fig4.savefig(out / "figure4_immune_evasion.png", dpi=200, bbox_inches="tight")
    fig4.savefig(out / "figure4_immune_evasion.pdf", bbox_inches="tight")
    print(f"Saved figure4 -> {out}")
plt.close()

# ── Supp Fig 5: Supplementary panels ─────────────────────────────
# Panel A: Stromal/Endothelial communication (from LIANA, non-TAM senders)
# Panel B: Extended metabolic gene list table (saved as figure)

import pandas as pd
import matplotlib.gridspec as gridspec

fig5, axes5 = plt.subplots(1, 2, figsize=(16, 7))
fig5.subplots_adjust(wspace=0.1, left=0.02, right=0.98, top=0.92, bottom=0.02)

# Panel A: Myofibroblast -> SPP1+B communication (complement reverse direction)
print("Supp Fig 5A: Myofibroblast -> OLR1+ TAM reverse comm...")
df_liana = pd.read_csv("phase3_results/cellchat/liana_border_all.csv")
sig5 = df_liana[df_liana.cellphone_pvals < 0.05].copy()
sig5["strength"] = -np.log10(sig5["magnitude_rank"] + 1e-6)

# Extended senders -> SPP1+B
extended_senders = ["Myofibroblasts", "Stromal 1", "Stromal 2",
                    "Tip-like ECs", "Stalk-like ECs", "Pericytes"]
rev = sig5[(sig5.source.isin(extended_senders)) & (sig5.target == "SPP1+B")]
top_rev = rev.nsmallest(25, "magnitude_rank")[
    ["source","ligand_complex","receptor_complex","lr_means","magnitude_rank","strength"]]
top_rev = top_rev.sort_values(["source","strength"], ascending=[True,False])

ax5a = axes5[0]
ax5a.axis("off")

# Color by source
src_colors5 = {
    "Myofibroblasts": "#3498DB",
    "Stromal 1":      "#5DADE2",
    "Stromal 2":      "#85C1E9",
    "Tip-like ECs":   "#1A5276",
    "Stalk-like ECs": "#2471A3",
    "Pericytes":      "#2E86C1",
}

rows_plot = []
for _, r in top_rev.iterrows():
    rows_plot.append((
        r.source,
        f"{r.ligand_complex} → {r.receptor_complex}",
        r.lr_means,
        r.strength,
        src_colors5.get(r.source, "#999999"),
    ))

y_pos = np.arange(len(rows_plot))
for i, (src, lr_label, lr_mean, strength, color) in enumerate(rows_plot):
    ax5a.barh(i, lr_mean, color=color, alpha=0.8, height=0.75)
    ax5a.text(-0.05, i, f"{src}:\n{lr_label}",
              ha="right", va="center", fontsize=6.5, color="#333333")
    ax5a.text(lr_mean + 0.01, i, f"s={strength:.1f}",
              ha="left", va="center", fontsize=6, color="#777777")

ax5a.set_yticks([])
ax5a.set_xlabel("Mean LR expression", fontsize=9)
ax5a.set_title("A  Reverse Communication:\nStromal/Myofibroblast → OLR1⁺ TAMs",
               fontsize=10, fontweight="bold", loc="left")
ax5a.spines[["top","right","left"]].set_visible(False)

# Legend
legend_patches = [mpatches.Patch(color=c, label=s)
                  for s, c in src_colors5.items()]
ax5a.legend(handles=legend_patches, fontsize=7, frameon=False,
            loc="lower right")

# Panel B: Metabolic pathway summary table
ax5b = axes5[1]
ax5b.axis("off")

pw_df = pd.read_csv("phase4_results/panelF_pathway_scores.csv")
pw_df["log2fc"] = pw_df["log2fc"].round(3)
pw_df["pval"]   = pw_df["pval"].apply(lambda p: f"{p:.3e}")
pw_df["dir"]    = pw_df.apply(
    lambda r: "MMRp↑" if r["direction"] == "MMRp↑" else "MMRd↑", axis=1)
pw_df["sig"]    = pw_df["sig"].map({True:"*", False:"n.s."})

cols = ["pathway","log2fc","pval","sig","dir"]
col_labels = ["Pathway","log2FC","p-value","Sig","Direction"]
tbl_data = pw_df[cols].values.tolist()

table = ax5b.table(
    cellText=tbl_data,
    colLabels=col_labels,
    cellLoc="center",
    loc="center",
    bbox=[0, 0, 1, 1],
)
table.auto_set_font_size(False)
table.set_fontsize(8)

# Color header
for j in range(len(col_labels)):
    table[(0, j)].set_facecolor("#2C3E50")
    table[(0, j)].set_text_props(color="white", fontweight="bold")

# Color rows by direction
for i, row in enumerate(pw_df.itertuples()):
    fc_color = "#FDECEA" if row.direction == "MMRp↑" else "#EBF5FB"
    for j in range(len(cols)):
        table[(i+1, j)].set_facecolor(fc_color)

ax5b.set_title("B  Metabolic Pathway Scores: OLR1⁺ TAMs\nMMRp vs MMRd (GSE178341)",
               fontsize=10, fontweight="bold", loc="left")

fig5.suptitle("Supp Fig 5  —  Supplementary: Reverse Communication & Metabolic Summary",
              fontsize=12, fontweight="bold", y=0.97)

for out in [OUT_DIRS[0], OUT_DIRS[2]]:
    fig5.savefig(out / "suppFig5_metabolism_supp.png", dpi=200, bbox_inches="tight")
    fig5.savefig(out / "suppFig5_metabolism_supp.pdf", bbox_inches="tight")
    print(f"Saved suppFig5 -> {out}")
plt.close()

print("Assembly complete.")
