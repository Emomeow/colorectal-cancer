import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pathlib

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

OUT = pathlib.Path("phase4_results")
DEL = pathlib.Path("Project_Delivery/results/Figure4_Communication_Metabolism"
                   "/Cell_Communication/figures")
DEL.mkdir(parents=True, exist_ok=True)

W, H = 14, 15
fig, ax = plt.subplots(figsize=(W, H))
ax.set_xlim(0, W); ax.set_ylim(0, H)
ax.axis("off")
fig.patch.set_facecolor("white")

# ── Palette ───────────────────────────────────────────────────────
DARK  = "#1C2833"
RED   = "#C0392B"; RED_D  = "#922B21"; RED_BG  = "#FEF5F5"
BLUE  = "#1F618D"; BLUE_D = "#154360"; BLUE_BG = "#F0F7FD"
GRN   = "#1E8449"; GRN_D  = "#155235"; GRN_BG  = "#F0FAF4"
GOLD  = "#9A7D0A"; GOLD_D = "#7D6608"; GOLD_BG = "#FDFBEE"
GRAY  = "#555555"

LX, MX, RX = 2.6, 6.0, 9.4   # column centres
AW = 3.1                        # column width

# ── Helpers ───────────────────────────────────────────────────────
def clean_card(x, y, w, h, hdr_col, hdr_txt, body_fc, body_lines):
    """
    Clean card: colored top header strip + light body.
    No accent bar — minimal design.
    """
    HDR_H = 0.48

    # Outer border rectangle (whole card)
    ax.add_patch(mpatches.Rectangle(
        (x - w/2, y - h/2), w, h,
        fc=body_fc, ec=hdr_col, lw=1.8, zorder=3))

    # Top header strip
    ax.add_patch(mpatches.Rectangle(
        (x - w/2, y + h/2 - HDR_H), w, HDR_H,
        fc=hdr_col, ec="none", zorder=4))

    # Header label — centered
    ax.text(x, y + h/2 - HDR_H/2, hdr_txt,
            ha="center", va="center", fontsize=9.5, fontweight="bold",
            color="white", zorder=5)

    # Body lines
    TOP = y + h/2 - HDR_H - 0.10
    row_h = (h - HDR_H - 0.15) / max(len(body_lines), 1)
    for i, (line, bold, col, fs, ital) in enumerate(body_lines):
        ax.text(x - w/2 + 0.18,
                TOP - (i + 0.5) * row_h,
                line, ha="left", va="center",
                fontsize=fs, fontweight="bold" if bold else "normal",
                color=col, fontstyle="italic" if ital else "normal", zorder=5)

def rect(x, y, w, h, fc, ec="#DDDDDD", lw=1.2, zorder=3):
    ax.add_patch(mpatches.Rectangle(
        (x - w/2, y - h/2), w, h,
        fc=fc, ec=ec, lw=lw, zorder=zorder))

def t(x, y, s, fs=9.5, fw="normal", col="black",
      va="center", ha="center", style="normal"):
    ax.text(x, y, s, fontsize=fs, fontweight=fw, color=col,
            va=va, ha=ha, fontstyle=style, zorder=4,
            multialignment="center")

def arrow_v(x0, y0, x1, y1, col="#888888", lw=1.8, head=0.22):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
        arrowprops=dict(
            arrowstyle=f"-|>, head_width={head}, head_length={head*0.85}",
            color=col, lw=lw, zorder=2))

def hline(y, x0, x1, col="#888888", lw=1.6):
    ax.plot([x0, x1], [y, y], color=col, lw=lw, zorder=2)

def vline(x, y0, y1, col="#888888", lw=1.6):
    ax.plot([x, x], [y0, y1], color=col, lw=lw, zorder=2)

# ═══════════════════════════════════════════════════════════════════
# TITLE
t(6.0, 14.60, "Integrated Mechanism Model", fs=16, fw="bold", col=DARK)
t(6.0, 14.10,
  "OLR1+ lipid-associated TAMs drive immune exclusion in MMRp CRC",
  fs=10.5, col=GRAY, style="italic")

# LIPID ENVIRONMENT
t(6.0, 13.50, "Lipid microenvironment  (oxLDL enrichment)", fs=9.5, col=GOLD_D)
arrow_v(6.0, 13.22, 6.0, 12.55, col=GOLD_D, lw=2.0, head=0.24)
# label beside the arrow, not floating between sections
t(6.22, 12.88, "OLR1-mediated lipid uptake", fs=9, col=GOLD_D, style="italic", ha="left")

# ── OLR1+ TAM HUB ────────────────────────────────────────────────
HY = 11.85
rect(6.0, HY, 7.4, 1.10, fc=RED_BG, ec=RED, lw=2.5)
# Colored left stripe
ax.add_patch(mpatches.Rectangle((6.0 - 7.4/2, HY - 0.55), 0.20, 1.10,
    fc=RED, ec="none", zorder=4))
t(6.0 + 0.10, HY + 0.26, "OLR1+ lipid-associated TAMs", fs=13, fw="bold", col=RED_D)
t(6.0 + 0.10, HY - 0.25,
  "SPP1+B   |   Border zone   |   35.5% of Border cells",
  fs=9, col=RED_D, style="italic")

# fork
FORK_Y = 11.01
vline(6.0, HY - 0.58, FORK_Y, col=DARK, lw=2.0)
hline(FORK_Y, LX, RX, col=DARK, lw=2.0)
for xc in [LX, MX, RX]:
    arrow_v(xc, FORK_Y, xc, 10.47, col=DARK, lw=2.0, head=0.20)

# ── ARM BANNERS ───────────────────────────────────────────────────
AHY, AHH = 10.16, 0.48
for xc, col, label in [
    (LX, BLUE, "[1]  Ligand-Receptor Axis"),
    (MX, GRN,  "[2]  Metabolic Competition"),
    (RX, GOLD, "[3]  ECM Remodeling Axis"),
]:
    rect(xc, AHY, AW, AHH, fc=col, ec=col, lw=0, zorder=3)
    t(xc, AHY, label, fs=10.5, fw="bold", col="white")

# ── CONTENT CARDS ─────────────────────────────────────────────────
CY, CH = 8.20, 3.10

for xc in [LX, MX, RX]:
    arrow_v(xc, AHY - AHH/2, xc, CY + CH/2, col="#CCCCCC", lw=1.2, head=0.14)

clean_card(LX, CY, AW, CH, BLUE, "Key Ligand-Receptor Pairs", BLUE_BG, [
    ("LGALS1  →  CD69 / PTPRC", True,  BLUE_D, 9.0, False),
    ("SPP1      →  CD44",        True,  BLUE_D, 9.0, False),
    ("HMGB1  →  CXCR4",         False, BLUE_D, 9.0, False),
    ("MIF  →  CD74 / CXCR4",    False, BLUE_D, 9.0, False),
    ("S100A8  →  CD69",          False, BLUE_D, 9.0, False),
    ("VEGFA  →  CXCR4",          False, BLUE_D, 9.0, False),
])

clean_card(MX, CY, AW, CH, GRN, "Metabolic Axes", GRN_BG, [
    ("Adenosine: ENTPD1 / NT5E",        True,  GRN_D,   9.0, False),
    ("MMRp-specific  p=1.9e-10",         False, GRAY,    8.0, True),
    ("IDO1/ARG1: Trp & Arg depletion",  True,  GRN_D,   9.0, False),
    ("T cell amino acid starvation",    False, GRAY,    8.5, True),
    ("Lipid/FAO    z-score = +1.70",    False, GRN,     9.0, False),
    ("OXPHOS       z-score = +1.70",    False, GRN,     9.0, False),
])

clean_card(RX, CY, AW, CH, GOLD, "ECM Remodeling", GOLD_BG, [
    ("SPP1 → ITGAV:ITGB1  (TAM→CAF)", True,  GOLD_D, 9.0, False),
    ("APOE → LRP1",                    False, GOLD_D, 9.0, False),
    ("COL1A1 → CD44  (feedback)",      False, GOLD,   9.0, False),
    ("TIMP1   → CD63   (feedback)",     False, GOLD,   9.0, False),
    ("→  Collagen deposition",          False, GOLD_D, 9.0, True),
    ("→  Physical ECM barrier",         False, GOLD_D, 9.0, True),
])

# ── OUTCOME CARDS ─────────────────────────────────────────────────
OY, OH = 5.76, 1.60

for xc in [LX, MX, RX]:
    arrow_v(xc, CY - CH/2, xc, OY + OH/2, col="#AAAAAA", lw=1.6, head=0.20)

clean_card(LX, OY, AW, OH, RED, "CD8+ T Cell Exhaustion", RED_BG, [
    ("PD-1 / TIM-3 / LAG3  up",      True,  RED_D, 9.0, False),
    ("TOX up  (terminal exhaustion)", False, RED,   9.0, False),
])
clean_card(MX, OY, AW, OH, GRN, "T Cell Metabolic Deprivation", GRN_BG, [
    ("Adenosine → ADORA2A → cAMP up",   True,  "#1A5276", 9.0, False),
    ("Amino acid depletion → mTOR down", False, GRN,       9.0, False),
])
clean_card(RX, OY, AW, OH, GOLD, "FAP+ CAF Activation", GOLD_BG, [
    ("Collagen secretion → Dense ECM",  True,  GOLD_D, 9.0, False),
    ("TAM-CAF positive feedback loop",  False, GOLD,   9.0, False),
])

# ── POSITIVE FEEDBACK — clean L-shaped path through right margin ──
# Route: FAP+ CAF right edge → right  →  up  → TAM hub right edge
MARGIN_X = 12.10           # stays in right margin, clear of all content
FB_SX    = RX + AW/2       # ~10.95  right edge of FAP+ CAF box
FB_SY    = OY              # ~5.76
FB_EX    = 6.0 + 7.4/2    # ~9.70   right edge of TAM hub
FB_EY    = HY              # ~11.85

# Segment 1: horizontal right from FAP+ CAF box
ax.plot([FB_SX, MARGIN_X], [FB_SY, FB_SY], color=GOLD_D, lw=2.0, zorder=2)
# Segment 2: vertical up through margin
ax.plot([MARGIN_X, MARGIN_X], [FB_SY, FB_EY], color=GOLD_D, lw=2.0, zorder=2)
# Segment 3: horizontal arrow into TAM hub right edge
ax.annotate("",
    xy=(FB_EX, FB_EY), xytext=(MARGIN_X, FB_EY),
    arrowprops=dict(
        arrowstyle="-|>, head_width=0.24, head_length=0.20",
        color=GOLD_D, lw=2.0, zorder=2))

# Label to the right of the vertical segment
t(13.10, (FB_SY + FB_EY) / 2,
  "Positive\nfeedback\n(TAM-CAF\nloop)",
  fs=9, col=GOLD_D, style="italic", ha="center")

# ── CONVERGENCE → FINAL ───────────────────────────────────────────
CONV_Y = 4.94
hline(CONV_Y, LX, RX, col=DARK, lw=2.0)
for xc in [LX, MX, RX]:
    vline(xc, OY - OH/2, CONV_Y, col=DARK, lw=1.8)
arrow_v(6.0, CONV_Y, 6.0, 4.32, col=DARK, lw=2.4, head=0.28)

rect(6.0, 3.86, 10.0, 1.40, fc=DARK, ec=DARK, lw=0, zorder=3)
t(6.0, 4.16, "MMRp CRC Immune Exclusion", fs=15, fw="bold", col="white")
t(6.0, 3.66,
  '"Cold tumor" phenotype   |   T cell-excluded microenvironment',
  fs=10, col="#AAAAAA", style="italic")

# ── CAPTION ───────────────────────────────────────────────────────
ax.add_patch(mpatches.Rectangle(
    (0.5, 0.10), 13.0, 1.00,
    fc="#F8F9FA", ec="#DDDDDD", lw=1.0, zorder=1))
t(7.0, 0.62,
  "[1] Direct ligand-receptor inhibition   "
  "[2] Metabolic competition & deprivation   "
  "[3] Physical ECM barrier via TAM-CAF crosstalk",
  fs=9, col="#444444", style="italic")

fig.tight_layout(pad=0.2)
for dest in [OUT, DEL]:
    fig.savefig(dest / "mechanism_model_EN.png", dpi=200, bbox_inches="tight")
    fig.savefig(dest / "mechanism_model_EN.pdf", bbox_inches="tight")
plt.close()
print("Done.")
