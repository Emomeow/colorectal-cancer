#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 3 Panel A — ABMIL 模型架构示意图（中文版）
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch
import warnings, io, sys, os
warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

plt.rcParams['font.family']  = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42

FIG_DIR = 'figures/phase3'
os.makedirs(FIG_DIR, exist_ok=True)

# ── 颜色 ─────────────────────────────────────────────────────────────
C = {
    'macro':  '#E64B35',
    'tcd8':   '#3C5488',
    'tcd4':   '#7EB4D4',
    'mono':   '#F4A460',
    'nk':     '#00A087',
    'other':  '#CCCCCC',
    'enc':    '#2C3E6A',
    'att_v':  '#C0392B',
    'att_u':  '#E67E22',
    'agg':    '#1A5C3A',
    'cls':    '#006080',
    'mmrp':   '#E64B35',
    'mmrd':   '#3A8FC0',
    'zone':   '#F4F6FA',
    'arrow':  '#444444',
}

fig = plt.figure(figsize=(20, 8.5))
fig.patch.set_facecolor('white')
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 20)
ax.set_ylim(0, 8.5)
ax.axis('off')

# ── 工具函数 ─────────────────────────────────────────────────────────
def rbox(x, y, w, h, fc='#E8F0FE', ec='#335599', lw=1.2,
         radius=0.15, alpha=1.0, zorder=2):
    p = FancyBboxPatch((x, y), w, h,
                        boxstyle=f'round,pad=0,rounding_size={radius}',
                        fc=fc, ec=ec, lw=lw, alpha=alpha, zorder=zorder)
    ax.add_patch(p)

def dot(x, y, r, fc, ec='white', lw=0.7, alpha=0.92, zorder=4):
    ax.add_patch(Circle((x, y), r, fc=fc, ec=ec, lw=lw, alpha=alpha, zorder=zorder))

def txt(x, y, s, fs=8.5, ha='center', va='center', c='#222222',
        fw='normal', zorder=6, style='normal'):
    ax.text(x, y, s, fontsize=fs, ha=ha, va=va, color=c,
            fontweight=fw, zorder=zorder, fontstyle=style)

def arr(x1, y1, x2, y2, c=C['arrow'], lw=1.6, hs=0.35, hw=0.18, zorder=5):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                 arrowprops=dict(
                     arrowstyle=f'->,head_length={hs},head_width={hw}',
                     color=c, lw=lw),
                 zorder=zorder)

def arr_curve(x1, y1, x2, y2, c=C['arrow'], lw=1.4, rad=0.3, zorder=5):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                 arrowprops=dict(
                     arrowstyle='->,head_length=0.3,head_width=0.16',
                     color=c, lw=lw,
                     connectionstyle=f'arc3,rad={rad}'),
                 zorder=zorder)

# ══════════════════════════════════════════════════════════════════════
#  Zone backgrounds (5 zones)
# ══════════════════════════════════════════════════════════════════════
ZONES = [
    (0.15, 0.4, 3.55, 7.7),   # A: 患者样本袋
    (3.9,  0.4, 3.6,  7.7),   # B: 细胞编码器
    (7.65, 0.4, 3.0,  7.7),   # C: 门控注意力
    (10.8, 0.4, 3.7,  7.7),   # D: 聚合与分类
    (14.65,0.4, 5.2,  7.7),   # E: 注意力结果
]
ZONE_COLORS = ['#FFF3F1','#F0F3FF','#F0FFF5','#FFFFF0','#F8F0FF']
ZONE_EDGES  = ['#E8B8B0','#B8C0E8','#A8D8B8','#D8D8A0','#C8B0E0']
for (x,y,w,h), fc, ec in zip(ZONES, ZONE_COLORS, ZONE_EDGES):
    rbox(x, y, w, h, fc=fc, ec=ec, lw=0.8, radius=0.3, zorder=1)

# Zone step labels
step_labels = ['(1) Patient Bag', '(2) Cell Encoder',
               '(3) Gated Attention', '(4) Aggregation & Classifier',
               '(5) Attention Result']
step_xs = [1.925, 5.7, 9.15, 12.65, 17.25]
for lbl, sx in zip(step_labels, step_xs):
    txt(sx, 7.85, lbl, fs=8.5, fw='bold', c='#333355')

# ══════════════════════════════════════════════════════════════════════
#  A: 患者样本袋
# ══════════════════════════════════════════════════════════════════════
txt(1.925, 7.45, 'MMRp Patient', fs=9.5, fw='bold', c=C['macro'])
txt(1.925, 7.08, '1 bag = 1 patient', fs=7.5, c='#666666', style='italic')

# Bag dashed border
rbox(0.35, 1.5, 3.1, 5.35, fc='white', ec='#CC3300',
     lw=1.8, radius=0.25, zorder=2)
ax.add_patch(FancyBboxPatch((0.35, 1.5), 3.1, 5.35,
     boxstyle='round,pad=0,rounding_size=0.25',
     fc='none', ec='#CC3300', lw=1.8, linestyle='--', zorder=3))

txt(0.80, 6.7, 'Bag', fs=8, c='#CC3300', fw='bold')

# Cell type legend (2 columns x 3 rows, below bag)
leg = [('TAM/Macro*', C['macro']), ('CD8_Tex/TCD8*', C['tcd8']),
       ('TCD4', C['tcd4']), ('Mono', C['mono']),
       ('NK', C['nk']),    ('Other', C['other'])]
for i, (lbl, col) in enumerate(leg):
    lx = 0.6 + (i % 2) * 1.55
    ly = 1.42 - (i // 2) * 0.25
    dot(lx, ly, 0.07, fc=col, zorder=5)
    fw = 'bold' if '*' in lbl else 'normal'
    col_t = col if '*' in lbl else '#444444'
    txt(lx + 0.12, ly, lbl.replace('*',''), fs=6.5, ha='left',
        c=col_t, fw=fw)

# Cells scattered
np.random.seed(7)
positions = [
    (0.75,6.2),(1.2,6.4),(1.7,6.1),(2.2,6.3),(2.7,6.15),(3.1,6.35),
    (0.65,5.5),(1.1,5.65),(1.6,5.4),(2.1,5.55),(2.6,5.45),(3.15,5.6),
    (0.75,4.75),(1.25,4.9),(1.75,4.65),(2.25,4.8),(2.75,4.7),(3.1,4.85),
    (0.65,4.0),(1.15,4.15),(1.65,3.95),(2.15,4.1),(2.65,4.0),(3.1,4.1),
    (0.75,3.25),(1.25,3.4),(1.75,3.2),(2.25,3.35),(2.75,3.25),(3.15,3.35),
    (0.85,2.55),(1.4,2.65),(1.9,2.5),(2.4,2.6),(2.9,2.55),
    (0.7,1.85),(1.3,1.95),(1.9,1.8),(2.5,1.9),(3.0,1.82),
]
cell_colors = (
    [C['macro']]*6 + [C['tcd8']]*6 + [C['tcd4']]*6 + [C['mono']]*5 +
    [C['nk']]*5 + [C['other']]*7 + [C['other']]*3
)[:len(positions)]

for (cx, cy), col in zip(positions, cell_colors):
    dot(cx, cy, 0.115, fc=col, alpha=0.88)

# Highlight one TAM cell
hx, hy = 1.2, 6.4
dot(hx, hy, 0.14, fc=C['macro'], ec='#880000', lw=2.0, zorder=7)
arr(hx + 0.14, hy, 3.6, hy, c='#880000', lw=1.2, hs=0.28, hw=0.14)

# Feature vector mini bars
bx0 = 3.62
for bi in range(6):
    bh = [0.28,0.18,0.32,0.12,0.22,0.08][bi]
    bc = C['macro'] if bi < 2 else '#BBBBBB'
    ax.bar(bx0 + bi*0.065, bh, bottom=hy - 0.16,
           width=0.05, color=bc, alpha=0.85, zorder=6)
txt(3.72, hy - 0.38, 'PCA-50\nvector', fs=6.5, c='#555555')
txt(1.925, 0.65, 'N cells / patient\n(median ~1,300)', fs=7, c='#777777')

# ══════════════════════════════════════════════════════════════════════
#  B: 细胞编码器
# ══════════════════════════════════════════════════════════════════════
txt(5.7, 7.45, 'Cell-level features', fs=8.5, fw='bold', c='#334488')

# Input dots column
in_ys = np.linspace(6.0, 3.0, 8)
for iy, col in zip(in_ys, [C['macro'],C['macro'],C['tcd8'],C['tcd8'],
                             C['tcd4'],C['mono'],C['nk'],C['other']]):
    dot(4.2, iy, 0.10, fc=col, alpha=0.85)

txt(4.2, 2.65, 'Input\nN cells\nPCA-50', fs=7, c='#555555')

# Arrows: input → encoder
enc_x, enc_y, enc_w, enc_h = 4.6, 3.0, 1.7, 3.2
mid_enc = enc_x + enc_w/2
for iy in in_ys:
    target_y = max(enc_y + 0.1, min(enc_y + enc_h - 0.1, iy))
    ax.annotate('', xy=(enc_x, target_y),
                 xytext=(4.32, iy),
                 arrowprops=dict(arrowstyle='->', color='#AAAACC',
                                 lw=0.7), zorder=4)

# Encoder MLP box
rbox(enc_x, enc_y, enc_w, enc_h, fc=C['enc'], ec='#1A2860',
     lw=1.3, radius=0.12, zorder=4)
txt(mid_enc, enc_y+enc_h-0.42, 'Encoder', fs=9, fw='bold', c='white')
txt(mid_enc, enc_y+enc_h-0.78, '(2-layer MLP)', fs=7.5, c='#BBCCFF')
txt(mid_enc, enc_y+enc_h-1.2, '50 → 256', fs=8, c='#DDDDFF', fw='bold')
txt(mid_enc, enc_y+enc_h-1.55, 'ReLU + Dropout', fs=7, c='#AABBDD')
txt(mid_enc, enc_y+enc_h-1.9, '256 → 256', fs=8, c='#DDDDFF', fw='bold')
txt(mid_enc, enc_y+enc_h-2.25, 'ReLU', fs=7, c='#AABBDD')
txt(mid_enc, enc_y+0.3, r'$H \in \mathbb{R}^{N \times 256}$', fs=8, c='#CCFFCC', fw='bold')

# Output dots column
out_ys = np.linspace(5.7, 3.3, 8)
for iy in out_ys:
    dot(6.65, iy, 0.10, fc='#778899', ec='white', lw=0.5)
txt(6.65, 2.85, 'H', fs=9, c='#334466', fw='bold')
txt(6.65, 2.55, 'N × 256', fs=7, c='#556677')

# Arrows: encoder → output
for iy in out_ys:
    ax.annotate('', xy=(6.55, iy),
                 xytext=(enc_x + enc_w, max(enc_y+0.1, min(enc_y+enc_h-0.1, iy))),
                 arrowprops=dict(arrowstyle='->', color='#AAAACC', lw=0.7), zorder=4)

# Zone A → Zone B main arrow
arr(3.88, 4.6, 4.12, 4.6, c='#884422', lw=2.2, hs=0.38, hw=0.20)
txt(4.0, 4.3, 'cells', fs=7, c='#884422')

# Zone B → Zone C main arrow
arr(6.78, 4.5, 7.6, 4.5, c='#334488', lw=2.2, hs=0.38, hw=0.20)
txt(7.18, 4.2, 'H', fs=9, c='#334488', fw='bold')

# ══════════════════════════════════════════════════════════════════════
#  C: 门控注意力
# ══════════════════════════════════════════════════════════════════════
txt(9.15, 7.45, 'Attention weights per cell', fs=8.5, fw='bold', c='#993322')

# V-path box (tanh)
rbox(7.8, 5.4, 1.7, 1.35, fc=C['att_v'], ec='#8B1A10', lw=1.2,
     radius=0.12, alpha=0.88, zorder=4)
txt(8.65, 6.35, 'V-path', fs=9, fw='bold', c='white')
txt(8.65, 6.0, '256 → 128', fs=8, c='#FFCCCC')
txt(8.65, 5.68, 'tanh  activation', fs=7.5, c='#FFE0DD')

# U-path box (sigmoid)
rbox(7.8, 3.8, 1.7, 1.35, fc=C['att_u'], ec='#8B4000', lw=1.2,
     radius=0.12, alpha=0.88, zorder=4)
txt(8.65, 4.78, 'U-path', fs=9, fw='bold', c='white')
txt(8.65, 4.43, '256 → 128', fs=8, c='#FFF0E0')
txt(8.65, 4.1, 'sigmoid  activation', fs=7.5, c='#FFE8CC')

# × multiply symbol
rbox(9.7, 4.55, 0.56, 0.56, fc='white', ec='#888888', lw=1.2,
     radius=0.12, zorder=5)
txt(9.98, 4.83, '×', fs=16, fw='bold', c='#555555')

# Linear(1) + Softmax
rbox(9.7, 3.2, 0.56, 1.0, fc='#6655AA', ec='#443388', lw=1.2,
     radius=0.10, alpha=0.9, zorder=4)
txt(9.98, 3.87, 'Linear', fs=7.5, fw='bold', c='white')
txt(9.98, 3.6, '(1)', fs=7.5, c='#DDCCFF')
txt(9.98, 3.33, 'Softmax', fs=7, c='#CCBBEE')

# Annotation
txt(9.98, 2.82, r'$A \in \mathbb{R}^{N}$', fs=8, fw='bold', c='#6655AA')
txt(9.98, 2.52, 'attention weights', fs=7, c='#888888', style='italic')

# Arrows within C
arr(7.78, 5.7,  7.6,  5.7,  c='#888888', lw=1.0, hs=0.22, hw=0.12)
arr(7.78, 4.47, 7.6,  4.47, c='#888888', lw=1.0, hs=0.22, hw=0.12)
ax.plot([7.6, 7.6], [4.47, 5.7], color='#AAAAAA', lw=1.0, ls='--', zorder=3)
arr(7.6, 5.1, 7.82, 6.07, c='#CCAAAA', lw=1.0, hs=0.22, hw=0.12)
arr(7.6, 5.1, 7.82, 4.47, c='#DDBB88', lw=1.0, hs=0.22, hw=0.12)
arr(9.52, 6.07, 9.72, 5.11, c=C['att_v'], lw=1.2, hs=0.25, hw=0.14)
arr(9.52, 4.47, 9.72, 4.83, c=C['att_u'], lw=1.2, hs=0.25, hw=0.14)
arr(9.98, 4.55, 9.98, 4.22, c='#6655AA', lw=1.3, hs=0.25, hw=0.14)

# Zone C → Zone D
arr(10.32, 3.0, 10.75, 3.0, c='#6655AA', lw=2.2, hs=0.38, hw=0.20)
txt(10.53, 2.7, 'A', fs=9, c='#6655AA', fw='bold')

# ══════════════════════════════════════════════════════════════════════
#  D: 聚合与分类
# ══════════════════════════════════════════════════════════════════════
txt(12.65, 7.45, 'Bag-level prediction', fs=8.5, fw='bold', c='#445500')

# Aggregation box
rbox(10.9, 5.4, 1.65, 1.7, fc='#1A5C3A', ec='#0D3A22', lw=1.3,
     radius=0.12, alpha=0.9, zorder=4)
txt(11.725, 6.65, 'Aggregation', fs=9, fw='bold', c='white')
txt(11.725, 6.3, r'$M = \sum_i a_i \cdot h_i$', fs=9.5, c='#AAFFCC', fw='bold')
txt(11.725, 5.95, 'weighted sum', fs=7.5, c='#88DDAA', style='italic')
txt(11.725, 5.6, r'$M \in \mathbb{R}^{256}$', fs=8, c='#CCFFEE', fw='bold')

# Classifier box
rbox(10.9, 3.5, 1.65, 1.55, fc='#005F7A', ec='#003C50', lw=1.3,
     radius=0.12, alpha=0.9, zorder=4)
txt(11.725, 4.7, 'Classifier', fs=9, fw='bold', c='white')
txt(11.725, 4.35, '256 → 1', fs=8, c='#AADDFF')
txt(11.725, 4.05, 'Dropout(0.25)', fs=7.5, c='#88BBDD')
txt(11.725, 3.72, 'σ(·) → P(MMRd)', fs=8, c='#CCEEFF', fw='bold')

arr(11.725, 5.4, 11.725, 5.07, c=C['agg'], lw=1.5, hs=0.28, hw=0.16)

# H feeds into aggregation
arr(7.6, 4.47, 10.87, 5.75, c='#334488', lw=1.2, hs=0.28, hw=0.14)
txt(9.3, 5.4, 'H (all cells)', fs=7, c='#334488', style='italic')

# Output boxes
rbox(12.75, 5.8, 2.1, 0.9, fc=C['mmrp'], ec='#990000', lw=1.4,
     radius=0.12, alpha=0.9, zorder=4)
txt(13.8, 6.42, 'P(MMRd) < 0.5', fs=8.5, fw='bold', c='white')
txt(13.8, 6.1, '→  Predicted: MMRp', fs=8, c='#FFCCCC')

rbox(12.75, 4.5, 2.1, 0.9, fc=C['mmrd'], ec='#004488', lw=1.4,
     radius=0.12, alpha=0.9, zorder=4)
txt(13.8, 5.12, 'P(MMRd) ≥ 0.5', fs=8.5, fw='bold', c='white')
txt(13.8, 4.8, '→  Predicted: MMRd', fs=8, c='#CCEEFF')

# Arrows → output
ax.plot([12.57, 12.75], [6.25, 6.25], color='#888888', lw=1.2, zorder=5)
ax.plot([12.57, 12.75], [4.95, 4.95], color='#888888', lw=1.2, zorder=5)
ax.plot([12.57, 12.57], [4.95, 6.25], color='#888888', lw=1.2, ls='--', zorder=4)
arr(12.57, 5.0, 12.57, 4.95, c='#888888', lw=1.0, hs=0.1, hw=0.08)
arr(11.725, 3.5, 11.725, 3.22, c=C['cls'], lw=1.5, hs=0.28, hw=0.16)

# BCE Loss box
rbox(10.9, 1.6, 1.65, 1.3, fc='#F4F4F4', ec='#AAAAAA', lw=1.0,
     radius=0.10, zorder=4)
txt(11.725, 2.58, 'BCE Loss', fs=8.5, fw='bold', c='#333333')
txt(11.725, 2.25, 'Adam + Weight Decay', fs=7, c='#555555')
txt(11.725, 1.95, 'Early Stopping (Pat.=20)', fs=6.5, c='#777777')
arr(11.725, 3.5, 11.725, 2.92, c='#AAAAAA', lw=1.0, hs=0.22, hw=0.12)

# CV result annotation
rbox(12.75, 2.8, 2.1, 1.6, fc='#FFFDE7', ec='#BBAA00', lw=1.3,
     radius=0.12, zorder=4)
txt(13.8, 4.1, '5-fold CV', fs=8.5, fw='bold', c='#554400')
txt(13.8, 3.75, 'AUC = 0.951', fs=10, fw='bold', c='#226600')
txt(13.8, 3.42, '± 0.067', fs=8.5, c='#444400')
txt(13.8, 3.08, 'Immune-only bags', fs=7, c='#888800', style='italic')

# Zone D → Zone E
arr(14.87, 3.5, 14.6, 3.5, c='#444444', lw=2.2, hs=0.38, hw=0.20)

# ══════════════════════════════════════════════════════════════════════
#  E: 注意力结果
# ══════════════════════════════════════════════════════════════════════
txt(17.25, 7.45, 'Attention → Discriminative Subtype', fs=8.5, fw='bold', c='#5533AA')

# Horizontal bar chart: mean attention by subtype (MMRp)
subtypes = ['Macro', 'Tgd', 'TCD8', 'NK', 'DC', 'Mono', 'TCD4']
att_raw  = [0.00230, 0.00111, 0.000275, 0.000156, 0.000141, 0.000053, 0.000004]
att_norm = np.array(att_raw) / max(att_raw)
bar_cols = [C['macro'], '#8888AA', C['tcd8'], C['nk'],
            '#CC6644', C['mono'], C['tcd4']]
pre_driver = {'Macro', 'TCD8'}

bar_h    = 0.52
bar_xoff = 15.3
bar_maxw = 2.8
y_top    = 6.85

for i, (st, val, col) in enumerate(zip(subtypes, att_norm, bar_cols)):
    by = y_top - i * (bar_h + 0.1)
    bw = val * bar_maxw
    # background track
    ax.barh(by, bar_maxw, height=bar_h, left=bar_xoff,
            color='#EEEEEE', alpha=0.6, zorder=3)
    # bar
    ax.barh(by, bw, height=bar_h, left=bar_xoff,
            color=col, alpha=0.88, zorder=4,
            edgecolor='white', linewidth=0.5)
    # subtype label (left)
    fw = 'bold' if st in pre_driver else 'normal'
    fc = col if st in pre_driver else '#333333'
    txt(bar_xoff - 0.12, by, st, fs=9, ha='right', va='center',
        c=fc, fw=fw)
    # value label (right of bar)
    txt(bar_xoff + bw + 0.1, by, f'{att_raw[i]:.5f}',
        fs=7, ha='left', va='center', c='#444444')
    # Pre-driver marker — right of the value label for each pre-driver row
    if st in pre_driver:
        txt(bar_xoff + bar_maxw + 0.72, by,
            '* Pre-driver', fs=7.5, ha='left', va='center',
            c='#AA0000', fw='bold')
    # significance — Macro: white text inside bar; TCD8: above its bar
    if st == 'Macro':
        txt(bar_xoff + bw - 0.12, by,
            'p≈0  FDR≈0', fs=6.5, c='white', fw='bold',
            ha='right', va='center')
    if st == 'TCD8':
        txt(bar_xoff + bw + 0.12, by + bar_h * 0.62,
            r'$p=1.7\times10^{-202}$', fs=6.5, c='#3344CC',
            fw='bold', ha='left')

# x-axis for bars
ax_y = y_top - len(subtypes) * (bar_h + 0.1) + bar_h/2 - 0.1
ax.plot([bar_xoff, bar_xoff + bar_maxw], [ax_y, ax_y],
        color='#AAAAAA', lw=0.8, zorder=3)
txt(bar_xoff + bar_maxw/2, ax_y - 0.28,
    'Mean Attention Score (MMRp patients)', fs=7.5, c='#555555')

# Curved arrow from attention result back to bag
arr_curve(15.0, 5.5, 1.9, 6.3, c='#9966CC', lw=1.4, rad=-0.45)
txt(8.4, 7.1, 'Attention back-maps to cell identity',
    fs=7.5, c='#7744AA', style='italic')

# ── Panel label & figure title ────────────────────────────────────────
txt(0.12, 8.4, 'A', fs=16, fw='bold', c='#111111', ha='left')
txt(10.0, 8.35,
    'Figure 3A — Attention-Based Multiple Instance Learning (ABMIL): '
    'MMRp vs MMRd Classification from Immune Microenvironment',
    fs=11, fw='bold', c='#1A1A2E', ha='center')
txt(10.0, 8.05,
    'Each patient = 1 bag  |  Each immune cell = 1 instance  |  '
    'Gated Attention identifies the most discriminative cell subtypes',
    fs=8, c='#555555', ha='center', style='italic')

# ── Save ──────────────────────────────────────────────────────────────
out = 'figure3A_schema'
fig.savefig(f'{FIG_DIR}/{out}.png', dpi=200, bbox_inches='tight',
            facecolor='white')
fig.savefig(f'{FIG_DIR}/{out}.pdf',           bbox_inches='tight',
            facecolor='white')
plt.close(fig)
print(f'Saved: {FIG_DIR}/{out}.png/.pdf')
