"""
Phase 5 — Figure 6: TCGA-COAD Bulk Deconvolution + Clinical Validation
Panel A: OLR1+ TAM signature score barcode heatmap (TCGA tumor samples)
Panel B: Kaplan-Meier OS curve — Right-sided vs Left-sided CRC
Panel C: Multivariable Cox regression Forest Plot
"""
import warnings; warnings.filterwarnings("ignore")
import io, sys, pathlib, gzip
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

OUT = pathlib.Path("phase5_results"); OUT.mkdir(exist_ok=True)
FIG = pathlib.Path("figures/phase5"); FIG.mkdir(parents=True, exist_ok=True)
DEL = pathlib.Path("Project_Delivery/results/Figure6_Bulk_Clinical_Validation")

# ── OLR1+ TAM signature genes ────────────────────────────────
OLR1_SIG = ['SPP1','OLR1','TREM2','MARCO','LPL','MMP9','APOE','C1QC',
             'CTSL','GPNMB','CYP27A1','FABP5','LGALS1','FOLR2','NT5E',
             'VEGFA','MIF','ENTPD1','PLTP','CD74']

# ── Load TCGA expression matrix ──────────────────────────────
print("Loading TCGA-COAD expression matrix...")
with gzip.open('data/TCGA-COAD.star_counts.tsv.gz', 'rt') as f:
    header = f.readline()
    samples = header.strip().split('\t')[1:]
    rows = []
    gene_ids = []
    for line in f:
        parts = line.strip().split('\t')
        gene_ids.append(parts[0])
        rows.append([float(x) for x in parts[1:]])

expr = pd.DataFrame(rows, index=gene_ids, columns=samples)
print(f"Expression: {expr.shape}")

# Strip Ensembl version and convert to symbol via mapping
# The TCGA matrix uses log2(count+1)-like values; check
print(f"Expression value range: {expr.values.min():.2f} - {expr.values.max():.2f}")

# Map Ensembl → gene symbol using a minimal mapping built from the features file
# We'll use the GSE188711 features (same gene universe) as the mapping table
import tarfile
tar = tarfile.open('data/GSE188711_RAW.tar', 'r')
members = tar.getmembers()
feat_m = [m for m in members if 'WGC' in m.name and 'features' in m.name][0]
feat_content = gzip.open(io.BytesIO(tar.extractfile(feat_m).read())).read().decode()
tar.close()

ensg_to_sym = {}
for line in feat_content.strip().split('\n'):
    parts = line.split('\t')
    if len(parts) >= 2:
        ensg_to_sym[parts[0]] = parts[1]

# Strip version from TCGA gene IDs
expr.index = [g.split('.')[0] for g in expr.index]
expr.index = [ensg_to_sym.get(g, g) for g in expr.index]
# Remove duplicate gene symbols (keep first)
expr = expr[~expr.index.duplicated(keep='first')]
print(f"After symbol mapping: {expr.shape}")

# Filter to tumor samples only (01A/01B)
tumor_cols = [c for c in samples if '-01' in c and c.endswith(('A','B'))]
expr_tumor = expr[tumor_cols].copy()
print(f"Tumor samples: {len(tumor_cols)}")

# ── Deconvolution: ssGSEA-like OLR1+ TAM score ───────────────
print("Computing OLR1+ TAM signature score...")
present_sig = [g for g in OLR1_SIG if g in expr_tumor.index]
print(f"Signature genes present: {len(present_sig)}/{len(OLR1_SIG)}: {present_sig}")

# ssGSEA-like: mean of signature gene expression per sample
olr1_score = expr_tumor.loc[present_sig].mean(axis=0)
olr1_score.name = 'OLR1TAM_score'
print(f"Score range: {olr1_score.min():.3f} - {olr1_score.max():.3f}")

# ── Load clinical data ────────────────────────────────────────
print("Loading clinical data...")
clin = pd.read_csv('data/TCGA-COAD_clinicalMatrix.tsv', sep='\t', low_memory=False)
surv_all = pd.read_csv('data/TCGA-COAD_survival_pancan.tsv', sep='\t')
coad_surv = surv_all[surv_all['cancer type abbreviation']=='COAD'].copy()

# Merge: sample IDs in TCGA expression = TCGA-XX-XXXX-01A
# Clinical sampleID = TCGA-XX-XXXX-01
clin['sample_expr'] = clin['sampleID'].apply(
    lambda x: x + '-01A' if not x.endswith('A') and not x.endswith('B') else x)

surv_score = pd.DataFrame({'sample': olr1_score.index, 'OLR1TAM_score': olr1_score.values})
# Extract TCGA patient barcode (first 12 chars)
surv_score['patient'] = surv_score['sample'].apply(lambda x: '-'.join(x.split('-')[:3]))

coad_surv['patient'] = coad_surv['sample'].apply(lambda x: '-'.join(x.split('-')[:3]))
merged = surv_score.merge(coad_surv[['patient','OS','OS.time','PFI','PFI.time',
                                     'age_at_initial_pathologic_diagnosis',
                                     'ajcc_pathologic_tumor_stage','gender']],
                          on='patient', how='inner')
# Also merge anatomic subdivision
clin['patient'] = clin['sampleID'].apply(lambda x: '-'.join(x.split('-')[:3]))
merged = merged.merge(clin[['patient','anatomic_neoplasm_subdivision','MSI_updated_Oct62011']],
                      on='patient', how='left')
print(f"Merged samples: {len(merged)}")
print(f"OS non-null: {merged['OS'].notna().sum()}")

# Clean OS
merged = merged.dropna(subset=['OS','OS.time']).copy()
merged['OS']      = pd.to_numeric(merged['OS'], errors='coerce')
merged['OS.time'] = pd.to_numeric(merged['OS.time'], errors='coerce')
merged = merged.dropna(subset=['OS','OS.time']).copy()
print(f"Final for KM: {len(merged)}")

# Stage encoding
def encode_stage(s):
    if pd.isna(s): return np.nan
    s = str(s).upper()
    if 'I' not in s and 'II' not in s and 'III' not in s and 'IV' not in s: return np.nan
    if 'IV' in s: return 4
    if 'III' in s: return 3
    if 'II' in s: return 2
    return 1
merged['stage_num'] = merged['ajcc_pathologic_tumor_stage'].apply(encode_stage)
merged['age'] = pd.to_numeric(merged['age_at_initial_pathologic_diagnosis'], errors='coerce')
merged['is_female'] = (merged['gender'] == 'FEMALE').astype(float)

# Sidedness
def assign_side(loc):
    if pd.isna(loc): return np.nan
    loc = str(loc)
    right_locs = ['Cecum','Ascending Colon','Transverse Colon','Hepatic Flexure']
    left_locs  = ['Sigmoid Colon','Descending Colon','Splenic Flexure','Rectosigmoid Junction']
    if any(r in loc for r in right_locs): return 'Right'
    if any(l in loc for l in left_locs):  return 'Left'
    return np.nan
merged['side'] = merged['anatomic_neoplasm_subdivision'].apply(assign_side)
merged['is_right'] = (merged['side'] == 'Right').astype(float)

# MSI
merged['is_msi_h'] = (merged['MSI_updated_Oct62011'] == 'MSI-H').astype(float)

# Save scores
merged.to_csv(OUT / 'tcga_olr1tam_scores.csv', index=False)

# ── High vs Low stratification (for Panel A) ─────────────────
median_score = merged['OLR1TAM_score'].median()
merged['olr1_group'] = np.where(merged['OLR1TAM_score'] >= median_score, 'High', 'Low')
print(f"High: {(merged.olr1_group=='High').sum()}, Low: {(merged.olr1_group=='Low').sum()}")

# ── Sidedness KM (Panel B) ────────────────────────────────────
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

m_side = merged.dropna(subset=['side']).copy()
right_df = m_side[m_side['side'] == 'Right']
left_df  = m_side[m_side['side'] == 'Left']
print(f"Right-sided: {len(right_df)}, Left-sided: {len(left_df)}")
lr = logrank_test(right_df['OS.time'], left_df['OS.time'],
                  event_observed_A=right_df['OS'], event_observed_B=left_df['OS'])
print(f"Log-rank p (Right vs Left) = {lr.p_value:.4f}")

# Cox regression — merge dtangle proportion if available
from lifelines import CoxPHFitter
dtangle_file = pathlib.Path('phase5_results/dtangle_with_clinical.csv')
if dtangle_file.exists():
    dt_all = pd.read_csv(dtangle_file)
    dt_all['patient'] = dt_all['patient'].astype(str)
    merged['patient_str'] = merged['patient'].astype(str)
    merged = merged.merge(dt_all[['patient','SPP1+B']].rename(
        columns={'patient':'patient_str','SPP1+B':'OLR1_dtangle'}),
        on='patient_str', how='left')
    # Scale to percent (0-100) to avoid extreme HR coefficients
    merged['OLR1_dtangle'] = merged['OLR1_dtangle'] * 100
    cox_cols = ['OS','OS.time','OLR1_dtangle','age','stage_num','is_female','is_right','is_msi_h']
    print("Using dtangle OLR1+ TAM proportion in Cox model")
else:
    cox_cols = ['OS','OS.time','OLR1TAM_score','age','stage_num','is_female','is_right','is_msi_h']
    print("Using ssGSEA OLR1+ TAM score in Cox model")

cox_df = merged[cox_cols].copy()
cox_df = cox_df.dropna()
print(f"Cox sample size: {len(cox_df)}")
cph = CoxPHFitter()
cph.fit(cox_df, duration_col='OS.time', event_col='OS')
print(cph.summary[['coef','exp(coef)','exp(coef) lower 95%','exp(coef) upper 95%','p']].to_string())
cph.summary.to_csv(OUT / 'cox_results.csv')

# ════════════════════════════════════════════════════════════════
# Figure 6
# ════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(18, 7.5))
gs  = gridspec.GridSpec(1, 3, figure=fig, wspace=0.45,
                         left=0.07, right=0.97, top=0.87, bottom=0.14)

C_HIGH = '#E74C3C'
C_LOW  = '#3498DB'

# ── Panel A: OLR1+ TAM proportion by sidedness (dtangle deconvolution) ──
# Load dtangle proportions
dtangle_file = pathlib.Path('phase5_results/dtangle_with_clinical.csv')
if dtangle_file.exists():
    dt = pd.read_csv(dtangle_file)
    dt['patient'] = dt['patient'].astype(str)
    # Use dtangle SPP1+B proportion
    dt_side = dt.dropna(subset=['side','SPP1+B']).copy()
    dt_right = dt_side[dt_side['side']=='Right']['SPP1+B'].values
    dt_left  = dt_side[dt_side['side']=='Left']['SPP1+B'].values
    from scipy.stats import mannwhitneyu
    _, p_side = mannwhitneyu(dt_right, dt_left, alternative='two-sided')
    use_dtangle = True
    print(f"dtangle: Right n={len(dt_right)}, Left n={len(dt_left)}, p={p_side:.4f}")
else:
    use_dtangle = False
    print("dtangle results not found, falling back to ssGSEA score barcode")

gs_a = gs[0, 0].subgridspec(2, 1, height_ratios=[3, 2], hspace=0.35)
ax_a_top = fig.add_subplot(gs_a[0])
ax_a_bot = fig.add_subplot(gs_a[1])

if use_dtangle:
    # Top: violin/box of OLR1+ TAM proportion by sidedness
    import matplotlib.patches as mpatches
    for i, (vals, color, label) in enumerate(zip(
            [dt_left, dt_right], [C_LOW, C_HIGH],
            ['Left-sided\nCRC', 'Right-sided\nCRC'])):
        parts = ax_a_top.violinplot(vals, positions=[i], widths=0.6,
                                    showmedians=True, showextrema=True)
        for pc in parts['bodies']:
            pc.set_facecolor(color); pc.set_alpha(0.6)
        parts['cmedians'].set_color('black'); parts['cmedians'].set_linewidth(2)
        jitter = np.random.uniform(-0.12, 0.12, len(vals))
        ax_a_top.scatter(np.full(len(vals), i) + jitter, vals,
                         s=8, color=color, alpha=0.4, zorder=3)

    sig_str = ('***' if p_side < 0.001 else '**' if p_side < 0.01
               else '*' if p_side < 0.05 else 'n.s.')
    y_max = max(dt_right.max(), dt_left.max())
    ax_a_top.plot([0, 1], [y_max*1.05, y_max*1.05], 'k-', lw=1)
    ax_a_top.text(0.5, y_max*1.07,
                  f'p={p_side:.3f}  {sig_str}', ha='center', fontsize=8)
    ax_a_top.set_xticks([0, 1])
    ax_a_top.set_xticklabels(['Left-sided\nCRC', 'Right-sided\nCRC'], fontsize=8.5)
    ax_a_top.set_ylabel('OLR1+ TAM proportion\n(dtangle)', fontsize=8)
    ax_a_top.set_ylim(0, y_max * 1.18)
    ax_a_top.spines[['top','right']].set_visible(False)
    ax_a_top.set_title('A  OLR1+ TAM Deconvolution (TCGA-COAD)\n'
                       'dtangle — scRNA-seq reference (GSE144735, 40 cell types)',
                       fontsize=10, fontweight='bold', loc='left')

    # Bottom: sorted barcode of OLR1+ TAM proportion with sidedness
    dt_sorted = dt_side.sort_values('SPP1+B').reset_index(drop=True)
    bar_colors = [C_HIGH if s=='Right' else C_LOW if s=='Left' else '#AAAAAA'
                  for s in dt_sorted['side']]
    ax_a_bot.bar(range(len(dt_sorted)), dt_sorted['SPP1+B'],
                 color=bar_colors, width=1.0, linewidth=0)
    ax_a_bot.set_xticks([]); ax_a_bot.set_yticks([0, 0.01, 0.02])
    ax_a_bot.set_ylabel('Proportion', fontsize=7)
    ax_a_bot.set_xlabel(f'TCGA-COAD samples (n={len(dt_sorted)}) sorted by OLR1+ TAM proportion',
                        fontsize=7.5)
    ax_a_bot.spines[['top','right']].set_visible(False)
    patch_r = mpatches.Patch(color=C_HIGH, label='Right-sided')
    patch_l = mpatches.Patch(color=C_LOW,  label='Left-sided')
    ax_a_bot.legend(handles=[patch_r, patch_l], fontsize=7, frameon=False, loc='upper left')
else:
    # Fallback: ssGSEA barcode
    sorted_merged = merged.sort_values('OLR1TAM_score').reset_index(drop=True)
    n = len(sorted_merged)
    im = ax_a_top.imshow(sorted_merged['OLR1TAM_score'].values.reshape(1,-1),
                         aspect='auto', cmap='RdYlBu_r', extent=[0,n,0,1])
    ax_a_top.set_yticks([]); ax_a_top.set_xticks([])
    fig.colorbar(im, ax=ax_a_top, fraction=0.04, pad=0.02).set_label('OLR1+ TAM score', fontsize=8)
    ax_a_top.set_title('A  OLR1+ TAM Signature Score (TCGA-COAD)',
                       fontsize=10, fontweight='bold', loc='left')
    ax_a_bot.axis('off')

# ── Panel B: Kaplan-Meier — Right vs Left sidedness ──────────
ax_b = fig.add_subplot(gs[0, 1])
n_right = len(right_df); n_left = len(left_df)

kmf_r = KaplanMeierFitter(label=f'Right-sided CRC (n={n_right})')
kmf_l_km = KaplanMeierFitter(label=f'Left-sided CRC (n={n_left})')

kmf_r.fit(right_df['OS.time'], event_observed=right_df['OS'])
kmf_l_km.fit(left_df['OS.time'],  event_observed=left_df['OS'])

kmf_r.plot_survival_function(ax=ax_b, color=C_HIGH, ci_show=True, ci_alpha=0.15, linewidth=2.0)
kmf_l_km.plot_survival_function(ax=ax_b, color=C_LOW,  ci_show=True, ci_alpha=0.15, linewidth=2.0)

# Median survival lines
med_r = kmf_r.median_survival_time_
med_l_val = kmf_l_km.median_survival_time_
if not np.isinf(med_r):
    ax_b.axvline(med_r, ls=':', color=C_HIGH, lw=1.2, alpha=0.7)
if not np.isinf(med_l_val):
    ax_b.axvline(med_l_val, ls=':', color=C_LOW, lw=1.2, alpha=0.7)

# Stats box
sig = '***' if lr.p_value < 0.001 else '**' if lr.p_value < 0.01 else '*' if lr.p_value < 0.05 else 'n.s.'
ax_b.text(0.97, 0.97, f'Log-rank p = {lr.p_value:.4f}\n{sig}',
          transform=ax_b.transAxes, ha='right', va='top', fontsize=8.5,
          bbox=dict(boxstyle='round', fc='white', ec='#CCCCCC', lw=0.8))

med_r_str = f'{med_r:.0f} d' if not np.isinf(med_r) else 'not reached'
med_l_str = f'{med_l_val:.0f} d' if not np.isinf(med_l_val) else 'not reached'
ax_b.text(0.03, 0.42,
          f'Right median OS: {med_r_str}\nLeft median OS:  {med_l_str}',
          transform=ax_b.transAxes, ha='left', va='top', fontsize=7.5,
          color='#444444',
          bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#DDDDDD',
                    alpha=0.9, lw=0.6))

ax_b.set_xlabel('Time (days)', fontsize=9)
ax_b.set_ylabel('Survival probability', fontsize=9)
ax_b.set_title('B  Overall Survival: Right vs Left\nAnatomic Sidedness (TCGA-COAD)',
               fontsize=10, fontweight='bold', loc='left')
ax_b.spines[['top','right']].set_visible(False)
ax_b.legend(fontsize=8, frameon=True, framealpha=0.9,
            loc='upper right', bbox_to_anchor=(1.0, 0.92))
ax_b.set_ylim(-0.05, 1.12)

# ── Panel C: Forest Plot (Cox) ────────────────────────────────
ax_c = fig.add_subplot(gs[0, 2])

cox_sum = cph.summary.copy()
var_labels = {
    'OLR1TAM_score': 'OLR1+ TAM score (ssGSEA)',
    'OLR1_dtangle':  'OLR1+ TAM % (dtangle)',
    'age':           'Age',
    'stage_num':     'Tumor stage',
    'is_female':     'Female sex',
    'is_right':      'Right-sided CRC',
    'is_msi_h':      'MSI-H',
}
cox_sum['label'] = cox_sum.index.map(lambda x: var_labels.get(x, x))
cox_sum = cox_sum.sort_values('exp(coef)', ascending=True)

y_pos = range(len(cox_sum))
colors_fp = [C_HIGH if row['p'] < 0.05 else '#AAAAAA'
             for _, row in cox_sum.iterrows()]

ax_c.errorbar(
    x=cox_sum['exp(coef)'],
    y=list(y_pos),
    xerr=[cox_sum['exp(coef)'] - cox_sum['exp(coef) lower 95%'],
          cox_sum['exp(coef) upper 95%'] - cox_sum['exp(coef)']],
    fmt='o', color='black', ecolor='#888888', elinewidth=1.5,
    markersize=7, capsize=3, zorder=3)

# Color dots
for i, (idx, row) in enumerate(cox_sum.iterrows()):
    c = C_HIGH if row['p'] < 0.05 else '#AAAAAA'
    ax_c.scatter(row['exp(coef)'], i, color=c, s=60, zorder=4)

ax_c.axvline(1.0, color='black', ls='--', lw=1)
ax_c.set_yticks(list(y_pos))
ax_c.set_yticklabels(cox_sum['label'], fontsize=9)
ax_c.set_xlabel('Hazard Ratio (95% CI)', fontsize=9)
ax_c.set_title('C  Multivariable Cox Regression\nIndependent Prognostic Value',
               fontsize=10, fontweight='bold', loc='left')
ax_c.spines[['top','right']].set_visible(False)

# Set xlim wide enough for CI bars + annotation text column
max_upper = cox_sum['exp(coef) upper 95%'].max()
ax_c.set_xlim(0.2, max_upper * 2.8)
text_x = max_upper * 1.15

# HR + p annotation in a clear text column beyond the CI bars
for i, (idx, row) in enumerate(cox_sum.iterrows()):
    sig_str = ' *' if row['p'] < 0.05 else ''
    p_str = 'p<0.001' if row['p'] < 0.001 else f'p={row["p"]:.3f}'
    ax_c.text(text_x, i,
              f"HR={row['exp(coef)']:.2f}, {p_str}{sig_str}",
              ha='left', va='center', fontsize=6.5)

patch_sig = mpatches.Patch(color=C_HIGH, label='p < 0.05')
patch_ns  = mpatches.Patch(color='#AAAAAA', label='n.s.')
ax_c.legend(handles=[patch_sig, patch_ns], fontsize=7.5, frameon=True,
            framealpha=0.9, loc='lower right')

# ── Super title ───────────────────────────────────────────────
fig.suptitle(
    'Figure 6  —  Bulk Clinical Validation: Right-sided CRC Has Worse Prognosis (TCGA-COAD)',
    fontsize=12, fontweight='bold', y=0.96)

# ── Save ─────────────────────────────────────────────────────
for dest in [FIG,
             DEL / 'Deconvolution' / 'figures',
             DEL / 'Survival_Analysis' / 'figures']:
    dest.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest / 'figure6_clinical_validation.png', dpi=200, bbox_inches='tight')
    fig.savefig(dest / 'figure6_clinical_validation.pdf', bbox_inches='tight')
    print(f'Saved -> {dest}')

plt.close()

# Save tables
cox_sum[['label','exp(coef)','exp(coef) lower 95%','exp(coef) upper 95%','p']].to_csv(
    OUT / 'cox_forest_summary.csv', index=False)
merged[['patient','OLR1TAM_score','olr1_group','OS','OS.time','stage_num','side']].to_csv(
    DEL / 'Survival_Analysis' / 'tables' / 'tcga_km_data.csv', index=False)
(DEL / 'Deconvolution' / 'tables').mkdir(parents=True, exist_ok=True)
olr1_score.sort_values(ascending=False).to_csv(
    DEL / 'Deconvolution' / 'tables' / 'tcga_olr1tam_scores.csv', header=True)

print('Figure 6 done.')
