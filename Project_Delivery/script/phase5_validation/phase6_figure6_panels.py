"""
Figure 6 — individual panel output (PNG + Illustrator-editable PDF).
Panel A -> Deconvolution/figures/
Panels B, C, D -> Survival_Analysis/figures/
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

OUT  = pathlib.Path("phase5_results"); OUT.mkdir(exist_ok=True)
DEL  = pathlib.Path("Project_Delivery/results/Figure6_Bulk_Clinical_Validation")
DIR_A = DEL / 'Deconvolution'   / 'figures'; DIR_A.mkdir(parents=True, exist_ok=True)
DIR_B = DEL / 'Survival_Analysis' / 'figures'; DIR_B.mkdir(parents=True, exist_ok=True)

C_HIGH = '#E74C3C'
C_LOW  = '#3498DB'
OLR1_SIG = ['SPP1','OLR1','TREM2','MARCO','LPL','MMP9','APOE','C1QC',
             'CTSL','GPNMB','CYP27A1','FABP5','LGALS1','FOLR2','NT5E',
             'VEGFA','MIF','ENTPD1','PLTP','CD74']

def save_fig(fig, directory, stem):
    fig.savefig(directory / f"{stem}.png", dpi=200, bbox_inches='tight')
    fig.savefig(directory / f"{stem}.pdf", bbox_inches='tight')
    print(f"Saved {stem}.png/pdf -> {directory}")
    plt.close(fig)

# ── Load TCGA ────────────────────────────────────────────────────────────────
print("Loading TCGA-COAD expression matrix...")
import tarfile
with gzip.open('data/TCGA-COAD.star_counts.tsv.gz', 'rt') as f:
    header  = f.readline()
    samples = header.strip().split('\t')[1:]
    rows = []; gene_ids = []
    for line in f:
        parts = line.strip().split('\t')
        gene_ids.append(parts[0])
        rows.append([float(x) for x in parts[1:]])
expr = pd.DataFrame(rows, index=gene_ids, columns=samples)
print(f"Expression: {expr.shape}")

tar = tarfile.open('data/GSE188711_RAW.tar', 'r')
members = tar.getmembers()
feat_m   = [m for m in members if 'WGC' in m.name and 'features' in m.name][0]
feat_content = gzip.open(io.BytesIO(tar.extractfile(feat_m).read())).read().decode()
tar.close()
ensg_to_sym = {}
for line in feat_content.strip().split('\n'):
    parts = line.split('\t')
    if len(parts) >= 2:
        ensg_to_sym[parts[0]] = parts[1]

expr.index = [g.split('.')[0] for g in expr.index]
expr.index = [ensg_to_sym.get(g, g) for g in expr.index]
expr = expr[~expr.index.duplicated(keep='first')]

tumor_cols  = [c for c in samples if '-01' in c and c.endswith(('A','B'))]
expr_tumor  = expr[tumor_cols].copy()
print(f"Tumor samples: {len(tumor_cols)}")

present_sig = [g for g in OLR1_SIG if g in expr_tumor.index]
olr1_score  = expr_tumor.loc[present_sig].mean(axis=0)
olr1_score.name = 'OLR1TAM_score'

# ── Clinical merge ────────────────────────────────────────────────────────────
print("Loading clinical data...")
clin     = pd.read_csv('data/TCGA-COAD_clinicalMatrix.tsv', sep='\t', low_memory=False)
surv_all = pd.read_csv('data/TCGA-COAD_survival_pancan.tsv', sep='\t')
coad_surv = surv_all[surv_all['cancer type abbreviation']=='COAD'].copy()

surv_score = pd.DataFrame({'sample': olr1_score.index, 'OLR1TAM_score': olr1_score.values})
surv_score['patient'] = surv_score['sample'].apply(lambda x: '-'.join(x.split('-')[:3]))
coad_surv['patient']  = coad_surv['sample'].apply(lambda x: '-'.join(x.split('-')[:3]))
merged = surv_score.merge(coad_surv[['patient','OS','OS.time',
                                     'age_at_initial_pathologic_diagnosis',
                                     'ajcc_pathologic_tumor_stage','gender']],
                          on='patient', how='inner')
clin['patient'] = clin['sampleID'].apply(lambda x: '-'.join(x.split('-')[:3]))
merged = merged.merge(clin[['patient','anatomic_neoplasm_subdivision','MSI_updated_Oct62011']],
                      on='patient', how='left')
merged = merged.dropna(subset=['OS','OS.time']).copy()
merged['OS']      = pd.to_numeric(merged['OS'], errors='coerce')
merged['OS.time'] = pd.to_numeric(merged['OS.time'], errors='coerce')
merged = merged.dropna(subset=['OS','OS.time']).copy()
print(f"Final for KM: {len(merged)}")

def encode_stage(s):
    if pd.isna(s): return np.nan
    s = str(s).upper()
    if 'IV'  in s: return 4
    if 'III' in s: return 3
    if 'II'  in s: return 2
    if 'I'   in s: return 1
    return np.nan

def assign_side(loc):
    if pd.isna(loc): return np.nan
    right_locs = ['Cecum','Ascending Colon','Transverse Colon','Hepatic Flexure']
    left_locs  = ['Sigmoid Colon','Descending Colon','Splenic Flexure','Rectosigmoid Junction']
    if any(r in str(loc) for r in right_locs): return 'Right'
    if any(l in str(loc) for l in left_locs):  return 'Left'
    return np.nan

merged['stage_num'] = merged['ajcc_pathologic_tumor_stage'].apply(encode_stage)
merged['age']       = pd.to_numeric(merged['age_at_initial_pathologic_diagnosis'], errors='coerce')
merged['is_female'] = (merged['gender'] == 'FEMALE').astype(float)
merged['side']      = merged['anatomic_neoplasm_subdivision'].apply(assign_side)
merged['is_right']  = (merged['side'] == 'Right').astype(float)
merged['is_msi_h']  = (merged['MSI_updated_Oct62011'] == 'MSI-H').astype(float)

# Sidedness KM data
m_side   = merged.dropna(subset=['side']).copy()
right_df = m_side[m_side['side'] == 'Right']
left_df  = m_side[m_side['side'] == 'Left']

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter

lr = logrank_test(right_df['OS.time'], left_df['OS.time'],
                  event_observed_A=right_df['OS'], event_observed_B=left_df['OS'])
print(f"Log-rank p (Right vs Left) = {lr.p_value:.4f}")

# Cox model
dtangle_file = pathlib.Path('phase5_results/dtangle_with_clinical.csv')
if dtangle_file.exists():
    dt_all = pd.read_csv(dtangle_file)
    dt_all['patient'] = dt_all['patient'].astype(str)
    merged['patient_str'] = merged['patient'].astype(str)
    merged = merged.merge(dt_all[['patient','SPP1+B']].rename(
        columns={'patient':'patient_str','SPP1+B':'OLR1_dtangle'}),
        on='patient_str', how='left')
    merged['OLR1_dtangle'] = merged['OLR1_dtangle'] * 100
    cox_cols = ['OS','OS.time','OLR1_dtangle','age','stage_num','is_female','is_right','is_msi_h']
else:
    cox_cols = ['OS','OS.time','OLR1TAM_score','age','stage_num','is_female','is_right','is_msi_h']

cox_df = merged[cox_cols].dropna()
cph = CoxPHFitter()
cph.fit(cox_df, duration_col='OS.time', event_col='OS')
print(cph.summary[['coef','exp(coef)','exp(coef) lower 95%','exp(coef) upper 95%','p']].to_string())

# dtangle violin data (Panel A + C)
use_dtangle = False
if dtangle_file.exists():
    dt = pd.read_csv(dtangle_file)
    dt['patient'] = dt['patient'].astype(str)
    dt_side = dt.dropna(subset=['side','SPP1+B']).copy()
    dt_right = dt_side[dt_side['side']=='Right']['SPP1+B'].values
    dt_left2 = dt_side[dt_side['side']=='Left']['SPP1+B'].values
    from scipy.stats import mannwhitneyu
    _, p_side = mannwhitneyu(dt_right, dt_left2, alternative='two-sided')
    use_dtangle = True
    print(f"dtangle: Right n={len(dt_right)}, Left n={len(dt_left2)}, p={p_side:.4f}")

# ══════════════════════════════════════════════════════════════════════════════
# PANEL A: deconvolution violin + barcode
# ══════════════════════════════════════════════════════════════════════════════
fig_a = plt.figure(figsize=(6.5, 7))
gs_a  = gridspec.GridSpec(2, 1, figure=fig_a, height_ratios=[3,2],
                           hspace=0.42, left=0.14, right=0.96, top=0.89, bottom=0.10)
ax_a_top = fig_a.add_subplot(gs_a[0])
ax_a_bot = fig_a.add_subplot(gs_a[1])

if use_dtangle:
    np.random.seed(0)
    for i, (vals, color) in enumerate(zip([dt_left2, dt_right], [C_LOW, C_HIGH])):
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
    y_max = max(dt_right.max(), dt_left2.max())
    ax_a_top.plot([0, 1], [y_max*1.05]*2, 'k-', lw=1)
    ax_a_top.text(0.5, y_max*1.07, f'p={p_side:.3f}  {sig_str}', ha='center', fontsize=9)
    ax_a_top.set_xticks([0, 1])
    ax_a_top.set_xticklabels(['Left-sided\nCRC', 'Right-sided\nCRC'], fontsize=9)
    ax_a_top.set_ylabel('OLR1+ TAM proportion\n(dtangle)', fontsize=9)
    ax_a_top.set_ylim(0, y_max * 1.18)
    ax_a_top.spines[['top','right']].set_visible(False)
    ax_a_top.set_title('A  OLR1+ TAM Deconvolution (TCGA-COAD)\n'
                       'dtangle — scRNA-seq reference (GSE144735, 40 cell types)',
                       fontsize=11, fontweight='bold', loc='left', pad=5)

    dt_sorted   = dt_side.sort_values('SPP1+B').reset_index(drop=True)
    bar_colors  = [C_HIGH if s=='Right' else C_LOW if s=='Left' else '#AAAAAA'
                   for s in dt_sorted['side']]
    ax_a_bot.bar(range(len(dt_sorted)), dt_sorted['SPP1+B'],
                 color=bar_colors, width=1.0, linewidth=0)
    ax_a_bot.set_xticks([])
    ax_a_bot.set_ylabel('Proportion', fontsize=8)
    ax_a_bot.set_xlabel(f'TCGA-COAD samples (n={len(dt_sorted)}) sorted by OLR1+ TAM proportion',
                        fontsize=8)
    ax_a_bot.spines[['top','right']].set_visible(False)
    patch_r = mpatches.Patch(color=C_HIGH, label='Right-sided')
    patch_l = mpatches.Patch(color=C_LOW,  label='Left-sided')
    ax_a_bot.legend(handles=[patch_r, patch_l], fontsize=8, frameon=False, loc='upper left')
else:
    sorted_merged = merged.sort_values('OLR1TAM_score').reset_index(drop=True)
    n = len(sorted_merged)
    im = ax_a_top.imshow(sorted_merged['OLR1TAM_score'].values.reshape(1,-1),
                         aspect='auto', cmap='RdYlBu_r', extent=[0,n,0,1])
    ax_a_top.set_yticks([]); ax_a_top.set_xticks([])
    fig_a.colorbar(im, ax=ax_a_top, fraction=0.04, pad=0.02).set_label('OLR1+ TAM score', fontsize=8)
    ax_a_top.set_title('A  OLR1+ TAM Signature Score (TCGA-COAD)',
                       fontsize=11, fontweight='bold', loc='left', pad=5)
    ax_a_bot.axis('off')

save_fig(fig_a, DIR_A, 'figure6A_deconvolution')

# ══════════════════════════════════════════════════════════════════════════════
# PANEL B: KM — Right vs Left sidedness
# ══════════════════════════════════════════════════════════════════════════════
fig_b, ax_b = plt.subplots(figsize=(6, 5.5))
fig_b.subplots_adjust(left=0.14, right=0.96, top=0.88, bottom=0.12)

n_right = len(right_df); n_left = len(left_df)
kmf_r  = KaplanMeierFitter(label=f'Right-sided CRC (n={n_right})')
kmf_l2 = KaplanMeierFitter(label=f'Left-sided CRC (n={n_left})')
kmf_r.fit(right_df['OS.time'],  event_observed=right_df['OS'])
kmf_l2.fit(left_df['OS.time'],  event_observed=left_df['OS'])

kmf_r.plot_survival_function(ax=ax_b, color=C_HIGH, ci_show=True, ci_alpha=0.15, linewidth=2.0)
kmf_l2.plot_survival_function(ax=ax_b, color=C_LOW,  ci_show=True, ci_alpha=0.15, linewidth=2.0)

med_r   = kmf_r.median_survival_time_
med_l_v = kmf_l2.median_survival_time_
if not np.isinf(med_r):   ax_b.axvline(med_r,   ls=':', color=C_HIGH, lw=1.2, alpha=0.7)
if not np.isinf(med_l_v): ax_b.axvline(med_l_v, ls=':', color=C_LOW,  lw=1.2, alpha=0.7)

sig_b = ('***' if lr.p_value<0.001 else '**' if lr.p_value<0.01
         else '*' if lr.p_value<0.05 else 'n.s.')
ax_b.text(0.97, 0.97, f'Log-rank p = {lr.p_value:.4f}\n{sig_b}',
          transform=ax_b.transAxes, ha='right', va='top', fontsize=9,
          bbox=dict(boxstyle='round', fc='white', ec='#CCCCCC', lw=0.8))

med_r_str = f'{med_r:.0f} d'   if not np.isinf(med_r)   else 'not reached'
med_l_str = f'{med_l_v:.0f} d' if not np.isinf(med_l_v) else 'not reached'
ax_b.text(0.03, 0.42,
          f'Right median OS: {med_r_str}\nLeft  median OS: {med_l_str}',
          transform=ax_b.transAxes, ha='left', va='top', fontsize=8,
          color='#444444',
          bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#DDDDDD', alpha=0.9, lw=0.6))

ax_b.set_xlabel('Time (days)', fontsize=10)
ax_b.set_ylabel('Survival probability', fontsize=10)
ax_b.set_title('B  Overall Survival: Right vs Left\nAnatomic Sidedness (TCGA-COAD)',
               fontsize=11, fontweight='bold', loc='left', pad=5)
ax_b.spines[['top','right']].set_visible(False)
ax_b.legend(fontsize=8.5, frameon=True, framealpha=0.9,
            loc='upper right', bbox_to_anchor=(1.0, 0.92))
ax_b.set_ylim(-0.05, 1.12)

save_fig(fig_b, DIR_B, 'figure6B_km_sidedness')

# ══════════════════════════════════════════════════════════════════════════════
# PANEL C: KM — dtangle OLR1+ TAM High vs Low
# ══════════════════════════════════════════════════════════════════════════════
fig_c, ax_c = plt.subplots(figsize=(6, 5.5))
fig_c.subplots_adjust(left=0.14, right=0.96, top=0.88, bottom=0.12)

if dtangle_file.exists():
    dt_km = pd.read_csv(dtangle_file)
    dt_km = dt_km.dropna(subset=['SPP1+B','OS','OS.time']).copy()
    dt_km['OS']      = pd.to_numeric(dt_km['OS'],      errors='coerce')
    dt_km['OS.time'] = pd.to_numeric(dt_km['OS.time'], errors='coerce')
    dt_km = dt_km.dropna(subset=['OS','OS.time']).copy()

    median_prop = dt_km['SPP1+B'].median()
    dt_km['dtangle_group'] = np.where(dt_km['SPP1+B'] >= median_prop, 'High', 'Low')
    dt_high = dt_km[dt_km['dtangle_group'] == 'High']
    dt_low  = dt_km[dt_km['dtangle_group'] == 'Low']

    lr_dt = logrank_test(dt_high['OS.time'], dt_low['OS.time'],
                         event_observed_A=dt_high['OS'], event_observed_B=dt_low['OS'])
    print(f"dtangle KM log-rank p = {lr_dt.p_value:.4f}")

    kmf_dh = KaplanMeierFitter(label=f'High OLR1+ TAM % (n={len(dt_high)})')
    kmf_dl = KaplanMeierFitter(label=f'Low OLR1+ TAM % (n={len(dt_low)})')
    kmf_dh.fit(dt_high['OS.time'], event_observed=dt_high['OS'])
    kmf_dl.fit(dt_low['OS.time'],  event_observed=dt_low['OS'])

    kmf_dh.plot_survival_function(ax=ax_c, color=C_HIGH, ci_show=True, ci_alpha=0.15, linewidth=2.0)
    kmf_dl.plot_survival_function(ax=ax_c, color=C_LOW,  ci_show=True, ci_alpha=0.15, linewidth=2.0)

    for kmf_obj, col in [(kmf_dh, C_HIGH), (kmf_dl, C_LOW)]:
        med = kmf_obj.median_survival_time_
        if not np.isinf(med):
            ax_c.axvline(med, ls=':', color=col, lw=1.2, alpha=0.7)

    sig_c = ('***' if lr_dt.p_value<0.001 else '**' if lr_dt.p_value<0.01
             else '*' if lr_dt.p_value<0.05 else 'n.s.')
    ax_c.text(0.97, 0.97, f'Log-rank p = {lr_dt.p_value:.4f}\n{sig_c}',
              transform=ax_c.transAxes, ha='right', va='top', fontsize=9,
              bbox=dict(boxstyle='round', fc='white', ec='#CCCCCC', lw=0.8))

    med_h_str  = f'{kmf_dh.median_survival_time_:.0f} d' if not np.isinf(kmf_dh.median_survival_time_) else 'not reached'
    med_l2_str = f'{kmf_dl.median_survival_time_:.0f} d' if not np.isinf(kmf_dl.median_survival_time_) else 'not reached'
    ax_c.text(0.03, 0.42,
              f'High median OS: {med_h_str}\nLow  median OS: {med_l2_str}',
              transform=ax_c.transAxes, ha='left', va='top', fontsize=8,
              color='#444444',
              bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#DDDDDD', alpha=0.9, lw=0.6))

    ax_c.set_xlabel('Time (days)', fontsize=10)
    ax_c.set_ylabel('Survival probability', fontsize=10)
    ax_c.set_title('C  Overall Survival: High vs Low OLR1+ TAM\ndtangle Deconvolution Median Split (TCGA-COAD)',
                   fontsize=11, fontweight='bold', loc='left', pad=5)
    ax_c.spines[['top','right']].set_visible(False)
    ax_c.legend(fontsize=8.5, frameon=True, framealpha=0.9,
                loc='upper right', bbox_to_anchor=(1.0, 0.92))
    ax_c.set_ylim(-0.05, 1.12)
else:
    ax_c.text(0.5, 0.5, 'dtangle_with_clinical.csv not found',
              ha='center', va='center', transform=ax_c.transAxes, fontsize=10)
    ax_c.set_title('C  dtangle KM (data missing)', fontsize=11, fontweight='bold', loc='left')

save_fig(fig_c, DIR_B, 'figure6C_km_olr1tam')

# ══════════════════════════════════════════════════════════════════════════════
# PANEL D: Cox forest plot
# ══════════════════════════════════════════════════════════════════════════════
var_labels = {
    'OLR1TAM_score': 'OLR1+ TAM score (ssGSEA)',
    'OLR1_dtangle':  'OLR1+ TAM % (dtangle)',
    'age':           'Age',
    'stage_num':     'Tumor stage',
    'is_female':     'Female sex',
    'is_right':      'Right-sided CRC',
    'is_msi_h':      'MSI-H',
}
cox_sum = cph.summary.copy()
cox_sum['label'] = cox_sum.index.map(lambda x: var_labels.get(x, x))
cox_sum = cox_sum.sort_values('exp(coef)', ascending=True)

n_vars = len(cox_sum)
fig_d, ax_d = plt.subplots(figsize=(7.5, max(4, n_vars * 0.65 + 1.5)))
fig_d.subplots_adjust(left=0.26, right=0.72, top=0.91, bottom=0.10)

y_pos = range(n_vars)
ax_d.errorbar(
    x=cox_sum['exp(coef)'],
    y=list(y_pos),
    xerr=[cox_sum['exp(coef)'] - cox_sum['exp(coef) lower 95%'],
          cox_sum['exp(coef) upper 95%'] - cox_sum['exp(coef)']],
    fmt='o', color='black', ecolor='#888888', elinewidth=1.5,
    markersize=7, capsize=3, zorder=3)

for i, (_, row) in enumerate(cox_sum.iterrows()):
    c = C_HIGH if row['p'] < 0.05 else '#AAAAAA'
    ax_d.scatter(row['exp(coef)'], i, color=c, s=65, zorder=4)

ax_d.axvline(1.0, color='black', ls='--', lw=1)
ax_d.set_yticks(list(y_pos))
ax_d.set_yticklabels(cox_sum['label'], fontsize=9.5)
ax_d.set_xlabel('Hazard Ratio (95% CI)', fontsize=10)
ax_d.set_title('D  Multivariable Cox Regression\nIndependent Prognostic Value',
               fontsize=11, fontweight='bold', loc='left', pad=5)
ax_d.spines[['top','right']].set_visible(False)

max_upper = cox_sum['exp(coef) upper 95%'].max()
ax_d.set_xlim(0.2, max_upper * 2.8)
text_x = max_upper * 1.15

for i, (_, row) in enumerate(cox_sum.iterrows()):
    sig_str = ' *' if row['p'] < 0.05 else ''
    p_str   = 'p<0.001' if row['p'] < 0.001 else f'p={row["p"]:.3f}'
    ax_d.text(text_x, i, f"HR={row['exp(coef)']:.2f}, {p_str}{sig_str}",
              ha='left', va='center', fontsize=7.5)

patch_sig = mpatches.Patch(color=C_HIGH,    label='p < 0.05')
patch_ns  = mpatches.Patch(color='#AAAAAA', label='n.s.')
ax_d.legend(handles=[patch_sig, patch_ns], fontsize=8.5, frameon=True,
            framealpha=0.9, loc='lower right')

save_fig(fig_d, DIR_B, 'figure6D_cox_forest')

print("\nAll Figure 6 panels saved.")
print(f"Panel A -> {DIR_A}")
print(f"Panels B-D -> {DIR_B}")
