import sys, io, numpy as np, pathlib
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

OUT = pathlib.Path("phase4_results")

CONNECTIONS = [
    dict(src="OLR1+ TAMs", tgt="CD8+ T cells",
         n=168, ligands="LGALS1, SPP1, HMGB1",
         color="#C0392B", lw=4.5, rad=0.28),
    dict(src="OLR1+ TAMs", tgt="CD4+ T cells",
         n=167, ligands="LGALS1, SPP1",
         color="#E67E22", lw=2.8, rad=0.22),
    dict(src="OLR1+ TAMs", tgt="Reg T cells",
         n=192, ligands="SPP1, MIF",
         color="#7D3C98", lw=2.8, rad=0.22),
    dict(src="OLR1+ TAMs", tgt="Myofibroblasts",
         n=379, ligands="SPP1->ITGAV:ITGB1\nLGALS1->ITGB1",
         color="#2980B9", lw=3.5, rad=-0.32),
    dict(src="Myofibroblasts", tgt="OLR1+ TAMs",
         n=471, ligands="TIMP1->CD63\nCOL1A1->CD44",
         color="#1A5276", lw=3.0, rad=-0.32),
]

NODE_POS = {
    "OLR1+ TAMs":    (0.50, 0.44),
    "CD8+ T cells":  (0.85, 0.72),
    "CD4+ T cells":  (0.85, 0.20),
    "Reg T cells":   (0.50, 0.90),
    "Myofibroblasts":(0.12, 0.44),
}
NODE_COLORS = {
    "OLR1+ TAMs":    "#C0392B",
    "CD8+ T cells":  "#1E8449",
    "CD4+ T cells":  "#D68910",
    "Reg T cells":   "#7D3C98",
    "Myofibroblasts":"#1A5276",
}
NODE_R = {
    "OLR1+ TAMs": 0.088,
    "CD8+ T cells": 0.060,
    "CD4+ T cells": 0.060,
    "Reg T cells":  0.060,
    "Myofibroblasts": 0.072,
}

fig = plt.figure(figsize=(11, 7))
ax  = fig.add_axes([0.00, 0.0, 0.62, 1.0])
ax2 = fig.add_axes([0.63, 0.05, 0.36, 0.90])

ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_aspect('equal'); ax.axis('off')
ax2.set_xlim(0, 1); ax2.set_ylim(0, 1)
ax2.axis('off')

# ── Arrows ───────────────────────────────────────────────────────
for conn in CONNECTIONS:
    x0, y0 = NODE_POS[conn['src']]
    x1, y1 = NODE_POS[conn['tgt']]
    shrink_s = NODE_R[conn['src']] * 530
    shrink_t = NODE_R[conn['tgt']] * 530

    ax.annotate("",
        xy=(x1, y1), xytext=(x0, y0),
        xycoords='data', textcoords='data',
        arrowprops=dict(
            arrowstyle="-|>, head_width=0.30, head_length=0.18",
            color=conn['color'], lw=conn['lw'], alpha=0.88,
            connectionstyle="arc3,rad=%.2f" % conn['rad'],
            shrinkA=shrink_s, shrinkB=shrink_t,
        ), zorder=2)

    # Ligand label: proper quadratic Bezier arc3 midpoint
    # matplotlib arc3 control point: ctrl = mid + rad*(dy, -dx)
    dx = x1-x0; dy = y1-y0
    ctrl_x = (x0+x1)/2 + conn['rad']*dy
    ctrl_y = (y0+y1)/2 - conn['rad']*dx
    arc_mx = 0.25*x0 + 0.5*ctrl_x + 0.25*x1
    arc_my = 0.25*y0 + 0.5*ctrl_y + 0.25*y1
    # offset in bulge direction (away from chord, toward arc exterior)
    bulge_x = conn['rad']*dy; bulge_y = -conn['rad']*dx
    bulge_norm = np.sqrt(bulge_x**2+bulge_y**2)+1e-9
    off_x = bulge_x/bulge_norm * 0.055
    off_y = bulge_y/bulge_norm * 0.055
    ax.text(arc_mx+off_x, arc_my+off_y, conn['ligands'],
            ha='center', va='center', fontsize=7.8, color=conn['color'],
            bbox=dict(boxstyle='round,pad=0.25', fc='white',
                      ec=conn['color'], alpha=0.92, lw=0.8), zorder=5)

    # n= label
    frac = 0.20
    tx = x1+(x0-x1)*frac; ty = y1+(y0-y1)*frac
    ax.text(tx, ty, "n=%d" % conn['n'],
            fontsize=6.5, color=conn['color'],
            ha='center', va='center', alpha=0.85, zorder=5)

# ── Nodes ────────────────────────────────────────────────────────
for name, (x, y) in NODE_POS.items():
    r = NODE_R[name]
    c = NODE_COLORS[name]
    ax.add_patch(plt.Circle((x, y), r, color=c, zorder=4, ec='white', lw=2.2))
    if name == "OLR1+ TAMs":
        ax.add_patch(plt.Circle((x, y), r+0.022, fill=False,
                                ec=c, lw=2.0, ls='--', zorder=3, alpha=0.55))
        ax.text(x, y, "Border", ha='center', va='center',
                fontsize=7.5, color='white', fontweight='bold', zorder=6)
    elif name == "Myofibroblasts":
        ax.text(x, y+0.010, "FAP+", ha='center', va='center',
                fontsize=7.5, color='white', fontweight='bold', zorder=6)
        ax.text(x, y-0.016, "CAF", ha='center', va='center',
                fontsize=6.5, color='white', zorder=6)

    if name == "Reg T cells":
        ax.text(x, y+r+0.040, name, ha='center', va='bottom',
                fontsize=9.5, fontweight='bold', color=c, zorder=6)
    elif name == "Myofibroblasts":
        ax.text(x-r-0.015, y, name, ha='right', va='center',
                fontsize=9.5, fontweight='bold', color=c, zorder=6)
    elif name == "OLR1+ TAMs":
        ax.text(x, y-r-0.042, name, ha='center', va='top',
                fontsize=10, fontweight='bold', color=c, zorder=6)
    else:
        ax.text(x+r+0.015, y, name, ha='left', va='center',
                fontsize=9.5, fontweight='bold', color=c, zorder=6)

ax.text(0.02, 0.98, "A   OLR1+ TAM Communication Network",
        fontsize=12, fontweight='bold', va='top', transform=ax.transAxes)
ax.text(0.02, 0.93, "GSE144735 Border | LIANA rank_aggregate | p < 0.05",
        fontsize=8.5, color='#555555', va='top', transform=ax.transAxes)
ax.text(0.02, 0.03, "Arrow width = communication strength  |  n = LR pairs",
        fontsize=7.5, color='#888888', va='bottom', transform=ax.transAxes)

# ── Legend ────────────────────────────────────────────────────────
ax2.text(0.0, 0.97, "Signal type", fontsize=10.5,
         fontweight='bold', va='top')

legend_items = [
    ("#C0392B", "Immune suppression",
     "OLR1+ TAMs -> CD8+ T cells", "n=168 | LGALS1, SPP1, HMGB1"),
    ("#E67E22", "T helper modulation",
     "OLR1+ TAMs -> CD4+ T cells", "n=167 | LGALS1, SPP1"),
    ("#7D3C98", "Treg recruitment",
     "OLR1+ TAMs -> Reg T cells",  "n=192 | SPP1, MIF"),
    ("#2980B9", "ECM remodeling (out)",
     "OLR1+ TAMs -> Myofibroblasts","n=379 | SPP1->ITGAV:ITGB1"),
    ("#1A5276", "CAF feedback (in)",
     "Myofibroblasts -> OLR1+ TAMs","n=471 | TIMP1->CD63, COL1A1->CD44"),
]

y0 = 0.88
for col, title, line1, line2 in legend_items:
    ax2.annotate("", xy=(0.10, y0), xytext=(0.0, y0),
                 arrowprops=dict(arrowstyle="-|>", color=col, lw=2.5))
    ax2.text(0.14, y0+0.016, title, fontsize=9, fontweight='bold', color=col, va='bottom')
    ax2.text(0.14, y0-0.008, line1, fontsize=7.8, color='#333333', va='top')
    ax2.text(0.14, y0-0.030, line2, fontsize=7.5, color='#666666', va='top')
    y0 -= 0.165

# Key finding box
ax2.add_patch(plt.Rectangle((0.0, 0.02), 0.98, 0.12,
              fc='#FEF9E7', ec='#E6AC00', lw=1.2, zorder=1,
              transform=ax2.transData))
ax2.text(0.49, 0.135, "Key finding", ha='center', fontsize=9,
         fontweight='bold', color='#9A6700')
ax2.text(0.49, 0.090,
         "OLR1+ TAMs act as the dominant immunosuppressive\nhub, targeting CD8+ T cells via LGALS1/SPP1/HMGB1\nwhile forming a positive feedback loop with FAP+ CAFs",
         ha='center', va='center', fontsize=7.8, color='#444444', linespacing=1.5)

for ext in ['png', 'pdf']:
    fig.savefig(OUT / ("panelA_circle_plot.%s" % ext), dpi=200, bbox_inches='tight')
plt.close()
print("Panel A redesign done.")
