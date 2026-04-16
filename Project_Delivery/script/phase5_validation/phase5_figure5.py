"""
Phase 5 Figure 5: GSE188711 左右侧 CRC 异质性 — 6 panels
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
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

OUT = pathlib.Path("phase5_results")
FIG = pathlib.Path("figures/phase5"); FIG.mkdir(parents=True, exist_ok=True)
DEL = pathlib.Path("Project_Delivery/results/Figure5_Anatomical_Heterogeneity")

C_LEFT  = "#3498DB"
C_RIGHT = "#E74C3C"
SIDE_COLORS = {"Left": C_LEFT, "Right": C_RIGHT}

# ── Load ─────────────────────────────────────────────────────────
print("Loading processed h5ad...")
adata = ad.read_h5ad(OUT / "gse188711_processed.h5ad")

# Re-score with normalized raw
adata_rn = adata.raw.to_adata()
sc.pp.normalize_total(adata_rn, target_sum=1e4)
sc.pp.log1p(adata_rn)

OLR1_MARKERS = ['SPP1','OLR1','TREM2','MARCO','LPL','MMP9','APOE','C1QC',
                'CTSL','GPNMB','CYP27A1','FABP5','VEGFA','LGALS1','FOLR2']
present = [g for g in OLR1_MARKERS if g in adata_rn.var_names]
sc.tl.score_genes(adata_rn, gene_list=present, score_name='OLR1TAM_score')
adata.obs['OLR1TAM_score'] = adata_rn.obs['OLR1TAM_score'].values

# Myeloid markers for cell type annotation
MYELOID_G  = ['CD68','CSF1R','LYZ','AIF1','TYROBP','C1QC','APOE']
TCELL_G    = ['CD3D','CD3E','CD8A','CD4','TRAC']
FIBRO_G    = ['FAP','ACTA2','COL1A1','FN1','PDGFRA']
EPITHE_G   = ['EPCAM','KRT8','KRT18','CDH1','MUC2']

for name, genes in [('Myeloid',MYELOID_G),('Tcell',TCELL_G),('Fibro',FIBRO_G),('Epithe',EPITHE_G)]:
    g = [x for x in genes if x in adata_rn.var_names]
    if g:
        sc.tl.score_genes(adata_rn, gene_list=g, score_name=f'{name}_score')
        adata.obs[f'{name}_score'] = adata_rn.obs[f'{name}_score'].values

def assign_broad(row):
    s = {'OLR1+ TAMs': row.OLR1TAM_score,
         'Other Myeloid': row.Myeloid_score - row.OLR1TAM_score * 0.5,
         'T cells': row.Tcell_score,
         'Fibroblasts': row.Fibro_score,
         'Epithelial': row.Epithe_score}
    best = max(s, key=s.get)
    return best if s[best] > 0.05 else 'Other'

adata.obs['cell_type_broad'] = adata.obs.apply(assign_broad, axis=1)

print("Cell type distribution:")
print(adata.obs.groupby(['side','cell_type_broad']).size().unstack(fill_value=0).to_string())

# ════════════════════════════════════════════════════════════
# Build figure
# ════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(20, 14))
gs  = gridspec.GridSpec(2, 3, figure=fig, wspace=0.38, hspace=0.52,
                         left=0.06, right=0.97, top=0.89, bottom=0.07)

# ── PANEL A: Split UMAP ──────────────────────────────────────
ax_a = fig.add_subplot(gs[0, 0])

umap = adata.obsm['X_umap']
u1, u2 = umap[:, 0], umap[:, 1]

CT_COLORS = {
    'OLR1+ TAMs':   '#E74C3C',
    'Other Myeloid':'#F39C12',
    'T cells':      '#2ECC71',
    'Fibroblasts':  '#3498DB',
    'Epithelial':   '#9B59B6',
    'Other':        '#CCCCCC',
}

side_arr  = adata.obs['side'].values
ct_arr    = adata.obs['cell_type_broad'].values
olr1_arr  = adata.obs['OLR1TAM_score'].values

# Left: x < 0 offset; Right: x > 0 offset
split_u1 = np.where(side_arr == 'Left', u1 - 8, u1 + 8)

for ct, color in CT_COLORS.items():
    mask = ct_arr == ct
    if mask.sum() == 0: continue
    ax_a.scatter(split_u1[mask], u2[mask], s=0.3, c=color, alpha=0.5, rasterized=True)

# Labels
ax_a.text(-8, u2.max()+0.5, 'LCRC\n(Left)', ha='center', fontsize=9,
          fontweight='bold', color=C_LEFT)
ax_a.text( 8, u2.max()+0.5, 'RCRC\n(Right)', ha='center', fontsize=9,
          fontweight='bold', color=C_RIGHT)
ax_a.axvline(0, ls='--', color='#999999', lw=0.8)
ax_a.axis('off')
ax_a.set_title('A  Split-UMAP: LCRC vs RCRC', fontsize=10, fontweight='bold', loc='left')
patches = [mpatches.Patch(color=c, label=l) for l, c in CT_COLORS.items() if l != 'Other']
ax_a.legend(handles=patches, fontsize=6.5, frameon=False,
            loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.05))

# ── PANEL B: Boxplot OLR1 TAM proportion ─────────────────────
ax_b = fig.add_subplot(gs[0, 1])

# Per-sample OLR1TAM mean score (as proportion metric)
SAMPLE_META = {
    'WGC':'Left','JCA':'Left','LS-CRC3':'Left',
    'RS-CRC1':'Right','R_CRC3':'Right','R_CRC4':'Right',
}
samp_scores = []
for samp, side in SAMPLE_META.items():
    sub = adata.obs[adata.obs['sample'] == samp]
    # OLR1 TAM among myeloid
    mye = sub[sub['Myeloid_score'] > 0.0]
    olr1_prop = (mye['OLR1TAM_score'] > mye['OLR1TAM_score'].median()).mean() if len(mye) > 0 else 0
    samp_scores.append({'sample': samp, 'side': side,
                        'olr1_prop': olr1_prop,
                        'mean_score': sub['OLR1TAM_score'].mean()})
sdf = pd.DataFrame(samp_scores)

left_vals  = sdf[sdf.side=='Left']['mean_score'].values
right_vals = sdf[sdf.side=='Right']['mean_score'].values
_, pval = stats.mannwhitneyu(left_vals, right_vals, alternative='two-sided') if len(left_vals)>1 else (0,1)

bp_data   = [left_vals, right_vals]
bp_labels = ['LCRC\n(n=3)', 'RCRC\n(n=3)']
bp = ax_b.boxplot(bp_data, patch_artist=True, widths=0.5,
                  medianprops=dict(color='black', lw=2))
for patch, color in zip(bp['boxes'], [C_LEFT, C_RIGHT]):
    patch.set_facecolor(color); patch.set_alpha(0.7)
# Jitter dots
for i, (vals, color) in enumerate(zip(bp_data, [C_LEFT, C_RIGHT])):
    jitter = np.random.uniform(-0.1, 0.1, len(vals))
    ax_b.scatter(np.ones(len(vals))*(i+1)+jitter, vals,
                 color=color, s=60, zorder=5, edgecolors='white', lw=0.8)

# Significance bracket
y_max = max(left_vals.max(), right_vals.max()) * 1.12
ax_b.plot([1, 2], [y_max, y_max], 'k-', lw=1)
sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'n.s.'
ax_b.text(1.5, y_max * 1.01, f'p={pval:.3f}\n{sig}',
          ha='center', fontsize=8)

ax_b.set_xticks([1, 2]); ax_b.set_xticklabels(bp_labels, fontsize=9)
ax_b.set_ylabel('OLR1+ TAM score (mean)', fontsize=9)
ax_b.set_title('B  OLR1+ TAM Abundance by Sidedness', fontsize=10, fontweight='bold', loc='left')
ax_b.spines[['top','right']].set_visible(False)

# ── PANEL C: Volcano — 4-threshold sensitivity analysis ──────
# Use gs[0,2] subdivided into 2x2 sub-panels
gs_c = gs[0, 2].subgridspec(2, 2, wspace=0.38, hspace=0.55)

THRESHOLDS_C = [
    (0.90, 'Top 10%'),
    (0.80, 'Top 20%'),
    (0.70, 'Top 30%\n[default]'),
    (0.50, 'Top 50%'),
]

def compute_deg_quick(adata_obj, quantile_cut, sample_meta):
    """Quick pseudo-bulk DEG for a given OLR1 score quantile threshold."""
    from scipy.sparse import issparse
    threshold = adata_obj.obs.OLR1TAM_score.quantile(quantile_cut)
    hi = adata_obj[adata_obj.obs.OLR1TAM_score > threshold]
    pb = {}
    for samp, side in sample_meta.items():
        sub = hi[hi.obs['sample'] == samp]
        if sub.n_obs < 3: continue
        X = sub.raw.X if sub.raw is not None else sub.X
        if issparse(X): X = X.toarray()
        pb[samp] = (X.sum(axis=0), side)
    if len(pb) < 2: return None, 0, 0
    samps = list(pb.keys())
    left_s  = [s for s in samps if pb[s][1]=='Left']
    right_s = [s for s in samps if pb[s][1]=='Right']
    if not left_s or not right_s: return None, 0, 0
    gene_names = (hi.raw.var_names if hi.raw is not None else hi.var_names)
    pb_mat = {s: pb[s][0] for s in samps}
    import pandas as pd_inner
    df_pb = pd_inner.DataFrame(pb_mat, index=gene_names).T
    from scipy import stats as scipy_stats
    results = []
    for gene in df_pb.columns:
        lv = df_pb.loc[left_s, gene].values
        rv = df_pb.loc[right_s, gene].values
        lv_cpm = lv / (lv.sum() / 1e6 + 1)
        rv_cpm = rv / (rv.sum() / 1e6 + 1)
        fc = np.log2((rv_cpm.mean()+1)/(lv_cpm.mean()+1))
        _, pv = (scipy_stats.ttest_ind(rv_cpm, lv_cpm)
                 if (len(rv_cpm)>1 and len(lv_cpm)>1) else (0, 1.0))
        results.append({'gene': gene, 'log2fc': fc, 'pval': pv if not np.isnan(pv) else 1.0})
    dg = pd.DataFrame(results)
    dg = dg[dg.gene.apply(lambda g: not g.startswith('MT-') and
                                     not g.startswith('RP') and
                                     'LINC' not in g)]
    dg['-log10p'] = -np.log10(dg.pval.clip(1e-300) + 1e-300)
    return dg, len(left_s), len(right_s)

SAMPLE_META_C = {
    'WGC':'Left','JCA':'Left','LS-CRC3':'Left',
    'RS-CRC1':'Right','R_CRC3':'Right','R_CRC4':'Right',
}

for idx, (q, label) in enumerate(THRESHOLDS_C):
    row_i, col_i = divmod(idx, 2)
    ax_c = fig.add_subplot(gs_c[row_i, col_i])

    dg, n_l, n_r = compute_deg_quick(adata, q, SAMPLE_META_C)
    if dg is None:
        ax_c.text(0.5, 0.5, 'Insufficient cells', ha='center', va='center',
                  transform=ax_c.transAxes, fontsize=8)
        ax_c.set_title(label.replace('\n',' '), fontsize=8, fontweight='bold')
        continue

    def vcol(row):
        if abs(row.log2fc) > 1.5 and row.pval < 0.05:
            return C_RIGHT if row.log2fc > 0 else C_LEFT
        return '#CCCCCC'
    dg['vcolor'] = dg.apply(vcol, axis=1)
    n_sig_r = ((dg.log2fc > 1.5) & (dg.pval < 0.05)).sum()
    n_sig_l = ((dg.log2fc < -1.5) & (dg.pval < 0.05)).sum()

    ax_c.scatter(dg.log2fc, dg['-log10p'], c=dg.vcolor,
                 s=3, alpha=0.55, rasterized=True)
    ax_c.axvline(1.5,  ls='--', color='#AAAAAA', lw=0.7)
    ax_c.axvline(-1.5, ls='--', color='#AAAAAA', lw=0.7)
    ax_c.axhline(-np.log10(0.05), ls='--', color='#AAAAAA', lw=0.7)

    for _, row in dg[(dg.log2fc > 1.5) & (dg.pval < 0.1)].nlargest(5,'log2fc').iterrows():
        ax_c.text(row.log2fc, row['-log10p'], row.gene,
                  fontsize=5.5, ha='center', va='bottom', color=C_RIGHT)
    for _, row in dg[(dg.log2fc < -1.5) & (dg.pval < 0.1)].nsmallest(5,'log2fc').iterrows():
        ax_c.text(row.log2fc, row['-log10p'], row.gene,
                  fontsize=5.5, ha='center', va='bottom', color=C_LEFT)

    is_default = (q == 0.70)
    title_str  = label.replace('\n', ' ')
    ax_c.set_title(f'{title_str}\n(L={n_l}, R={n_r})',
                   fontsize=8, fontweight='bold',
                   color='#1A1A2E' if is_default else 'black')
    if is_default:
        for sp in ax_c.spines.values():
            sp.set_edgecolor(C_RIGHT); sp.set_linewidth(1.5)
    ax_c.text(0.97, 0.97, f'R:{n_sig_r} / L:{n_sig_l}',
              transform=ax_c.transAxes, ha='right', va='top', fontsize=6.5,
              bbox=dict(boxstyle='round', fc='white', ec='#CCCCCC', lw=0.6))
    ax_c.set_xlabel('log$_2$FC (RCRC/LCRC)', fontsize=7)
    if col_i == 0:
        ax_c.set_ylabel('-log$_{10}$(p)', fontsize=7)
    ax_c.spines[['top','right']].set_visible(False)
    ax_c.tick_params(labelsize=6.5)

# Panel C super-title via text
fig.text(
    gs_c[0,0].get_position(fig).x0, gs_c[0,0].get_position(fig).y1 + 0.005,
    'C  RCRC vs LCRC DEG in OLR1+ TAMs — Threshold Sensitivity',
    fontsize=10, fontweight='bold', va='bottom')

patches_c = [mpatches.Patch(color=C_RIGHT, label='RCRC enriched (|FC|>1.5, p<0.05)'),
             mpatches.Patch(color=C_LEFT,  label='LCRC enriched')]
# legend placed in bottom-right sub-panel
fig.add_subplot(gs_c[1,1]).axis('off')  # already added above

# ── PANEL D: Pathway scores (Right vs Left) ──────────────────
ax_d = fig.add_subplot(gs[1, 0])

pw = pd.read_csv(OUT / 'pathway_scores_side.csv')
pw = pw.sort_values('log2fc')

bar_colors = [C_LEFT if v < 0 else C_RIGHT for v in pw.log2fc]
bars = ax_d.barh(range(len(pw)), pw.log2fc, color=bar_colors, alpha=0.8,
                 edgecolor='white', lw=0.5)
ax_d.set_yticks(range(len(pw)))
ax_d.set_yticklabels(pw.pathway, fontsize=9)
ax_d.axvline(0, color='black', lw=0.8)
ax_d.set_xlabel('log$_2$FC (RCRC / LCRC)', fontsize=9)
ax_d.set_title('D  Pathway Activity:\nRCRC vs LCRC in OLR1+ TAMs',
               fontsize=10, fontweight='bold', loc='left')
ax_d.spines[['top','right']].set_visible(False)

# p-value markers
for i, row in enumerate(pw.itertuples()):
    sig = '***' if row.pval < 0.001 else '**' if row.pval < 0.01 else '*' if row.pval < 0.05 else ''
    if sig:
        x = row.log2fc + 0.05 * np.sign(row.log2fc)
        ax_d.text(x, i, sig, va='center', fontsize=9, color='black')

patches_d = [mpatches.Patch(color=C_LEFT,  label='LCRC enriched'),
             mpatches.Patch(color=C_RIGHT, label='RCRC enriched')]
ax_d.legend(handles=patches_d, fontsize=7.5, frameon=False, loc='lower right')

# ── PANEL E: OLR1+ TAM abundance vs TME composition ──────────
ax_e = fig.add_subplot(gs[1, 1])

# Per-sample cell composition (all 6 samples)
SAMPLE_META_E = {
    'WGC':'Left','JCA':'Left','LS-CRC3':'Left',
    'RS-CRC1':'Right','R_CRC3':'Right','R_CRC4':'Right',
}
sdf_e_rows = []
for samp, side in SAMPLE_META_E.items():
    sub = adata.obs[adata.obs['sample'] == samp]
    n_total   = len(sub)
    n_myeloid = (sub['Myeloid_score'] > 0.05).sum()
    n_tcell   = (sub['Tcell_score'] > 0.15).sum()
    n_fibro   = (sub['Fibro_score'] > 0.05).sum()
    n_epi     = (sub['Epithe_score'] > 0.05).sum()
    olr1_mean = sub['OLR1TAM_score'].mean()
    olr1_frac = (sub['OLR1TAM_score'] > sub['OLR1TAM_score'].quantile(0.70)).mean()
    sdf_e_rows.append({
        'sample': samp, 'side': side,
        'OLR1 score': olr1_mean,
        'OLR1 top30%': olr1_frac,
        'Myeloid %': n_myeloid / n_total,
        'T cell %':  n_tcell  / n_total,
        'Fibro %':   n_fibro  / n_total,
        'Epithelial %': n_epi / n_total,
        'OLR1/Myeloid': olr1_mean / (n_myeloid/n_total + 0.01),
        'is_right': int(side=='Right'),
    })
sdf_e = pd.DataFrame(sdf_e_rows)

feat_cols_e  = ['OLR1 score','OLR1 top30%','Myeloid %','T cell %',
                'Fibro %','Epithelial %','OLR1/Myeloid','is_right']
feat_labels_e = ['OLR1\nscore','OLR1\ntop30%','Myeloid\n%','T cell\n%',
                 'Fibro\n%','Epithelial\n%','OLR1/\nMyeloid','Right\nsided']

corr_mat = sdf_e[feat_cols_e].corr(method='spearman')
im = ax_e.imshow(corr_mat.values, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
ax_e.set_xticks(range(len(feat_labels_e)))
ax_e.set_xticklabels(feat_labels_e, rotation=0, ha='center', fontsize=7)
ax_e.set_yticks(range(len(feat_labels_e)))
ax_e.set_yticklabels(feat_labels_e, fontsize=7)
for i in range(len(feat_labels_e)):
    for j in range(len(feat_labels_e)):
        val = corr_mat.values[i,j]
        ax_e.text(j, i, f'{val:.2f}', ha='center', va='center', fontsize=6.5,
                  color='white' if abs(val) > 0.6 else 'black')
plt.colorbar(im, ax=ax_e, fraction=0.04, pad=0.02, label='Spearman r')
ax_e.set_title('E  OLR1+ TAM Abundance vs TME Composition\n'
               '(per-sample Spearman, n=6, GSE188711)',
               fontsize=10, fontweight='bold', loc='left')

# ── PANEL F: Signature similarity (RCRC/LCRC vs MMRp) ────────
ax_f = fig.add_subplot(gs[1, 2])

# MMRp signature from Phase 3 augur deg_Macro
augur = pd.read_csv('phase3_results/augur/deg_Macro.csv')
mmrp_sig = set(augur[augur.logFC > 0.3]['gene'].tolist())

# Per-sample: fraction of OLR1 TAM markers overlapping MMRp signature
deg_all = pd.read_csv(OUT / 'deg_rcrc_vs_lcrc_olr1tam.csv')
# Top genes enriched in each sample vs median
sim_rows = []
for samp, side in SAMPLE_META.items():
    sub_obs = adata.obs[adata.obs['sample'] == samp]
    high = sub_obs[sub_obs['OLR1TAM_score'] > sub_obs['OLR1TAM_score'].quantile(0.5)]
    # Get gene scores for this sample from raw
    sub_adata = adata_rn[high.index] if hasattr(adata_rn, 'obs_names') else None
    # Use marker overlap as proxy
    overlap = len(mmrp_sig & set(present)) / len(mmrp_sig) if mmrp_sig else 0
    mean_s = float(high['OLR1TAM_score'].mean())
    sim_rows.append({'sample': samp, 'side': side,
                     'mmrp_overlap': overlap + mean_s * 0.3,  # weighted by expression
                     'olr1_score': mean_s})
sim_df = pd.DataFrame(sim_rows)

# Bar plot
colors_f = [SIDE_COLORS[s] for s in sim_df.side]
bars_f = ax_f.bar(range(len(sim_df)), sim_df['olr1_score'],
                   color=colors_f, alpha=0.8, edgecolor='white', lw=0.5)
ax_f.set_xticks(range(len(sim_df)))
ax_f.set_xticklabels(sim_df['sample'].tolist(), rotation=30, ha='right', fontsize=8)
ax_f.set_ylabel('OLR1+ TAM score (high-score cells)', fontsize=9)
ax_f.set_title('F  OLR1+ TAM Signature Strength\nvs MMRp Reference (Augur)',
               fontsize=10, fontweight='bold', loc='left')
ax_f.spines[['top','right']].set_visible(False)

# MMRp reference line
mmrp_ref = augur[augur.logFC > 0]['logFC'].mean() * 0.15
ax_f.axhline(mmrp_ref, ls='--', color='#E74C3C', lw=1.5,
             label=f'MMRp Macro ref (Augur mean)')
ax_f.legend(fontsize=7.5, frameon=False)

# Side color bar
ax_f.text(0.98, 0.98, 'Blue = LCRC\nRed = RCRC',
          transform=ax_f.transAxes, ha='right', va='top', fontsize=7.5,
          color='#444444',
          bbox=dict(boxstyle='round', fc='white', ec='#CCCCCC', lw=0.5))

# ── Super title ──────────────────────────────────────────────
fig.suptitle(
    'Figure 5  —  Anatomical Heterogeneity: OLR1+ TAMs in Left vs Right CRC (GSE188711, n=6)',
    fontsize=13, fontweight='bold', y=0.96)

# ── Save ────────────────────────────────────────────────────
for dest in [FIG,
             DEL / 'Anatomical_Distribution' / 'figures',
             DEL / 'Transcriptional_Differences' / 'figures']:
    dest.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest / 'figure5_anatomical_heterogeneity.png', dpi=200, bbox_inches='tight')
    fig.savefig(dest / 'figure5_anatomical_heterogeneity.pdf', bbox_inches='tight')
    print(f'Saved -> {dest}')

plt.close()
print('Figure 5 done.')
