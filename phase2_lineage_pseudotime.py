#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 2 — Lineage-specific pseudotime
Run Palantir separately on Myeloid cells and T cells.
Saves results + plots in figures/phase2/{myeloid,tcell}/
"""
import scanpy as sc
import numpy as np
import pandas as pd
import palantir
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys, io
from datetime import datetime

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# ── Helpers ───────────────────────────────────────────────────────

def resample_zscore(df, genes, n_pts=100):
    """df: genes×pseudotime; return z-scored genes×n_pts array."""
    df_sel = df.loc[genes]
    idx    = np.linspace(0, df_sel.shape[1] - 1, n_pts).astype(int)
    df_rs  = df_sel.iloc[:, idx]
    mu  = df_rs.mean(axis=1)
    sig = df_rs.std(axis=1).clip(lower=1e-6)
    return df_rs.subtract(mu, axis=0).divide(sig, axis=0).clip(-3, 3)


def make_heatmap(gene_trends, terminal_ids, deg_pair, labels, out_path):
    """2×2 heatmap: lineage1-genes vs lineage2-genes × each lineage."""
    T1, T2 = terminal_ids[0], terminal_ids[1]
    df1 = gene_trends[T1]['trends']   # genes × pseudotime
    df2 = gene_trends[T2]['trends']

    trend_set = set(df1.index) & set(df2.index)
    top = []
    for deg in deg_pair:
        if deg is not None:
            g = [x for x in deg['names'].head(60) if x in trend_set][:20]
        else:
            g = []
        if len(g) < 5:
            g = df1.var(axis=1).sort_values(ascending=False).head(20).index.tolist()
        top.append(g)
    top_g1, top_g2 = top

    N = 100
    panels = [
        (resample_zscore(df1, top_g1), top_g1, f'{labels[0]} genes | {labels[0]} lineage', '#E64B35'),
        (resample_zscore(df2, top_g2), top_g2, f'{labels[1]} genes | {labels[1]} lineage', '#00A087'),
        (resample_zscore(df1, top_g2), top_g2, f'{labels[1]} genes | {labels[0]} lineage', '#4DBBD5'),
        (resample_zscore(df2, top_g1), top_g1, f'{labels[0]} genes | {labels[1]} lineage', '#F39B7F'),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10),
                             gridspec_kw={'hspace': 0.35, 'wspace': 0.4})
    fig.suptitle(f'Gene Trends: {labels[0]} vs {labels[1]} lineages', fontsize=12, fontweight='bold')
    for ax, (mat, genes, title, col) in zip(axes.flat, panels):
        im = ax.imshow(mat.values, aspect='auto', cmap='RdBu_r',
                       vmin=-2.5, vmax=2.5, interpolation='nearest')
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels(genes, fontsize=7)
        ax.set_xticks([0, N // 2, N - 1])
        ax.set_xticklabels(['Early\n(UC)', 'Mid', 'Late\n(CRC)'], fontsize=8)
        ax.set_title(title, fontsize=9, fontweight='bold', color=col, pad=4)
        plt.colorbar(im, ax=ax, shrink=0.8, label='Z-score', pad=0.01)
    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"   Saved: {out_path}")


def run_lineage(adata_all, cell_type, tag):
    """Full pseudotime + DEG + plots for one cell type."""
    print(f"\n{'='*65}")
    print(f"  {tag.upper()} lineage  |  {datetime.now().strftime('%H:%M:%S')}")
    print(f"{'='*65}")

    fig_dir = f'figures/phase2/{tag}'
    res_dir = f'phase2_results/{tag}'
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(res_dir, exist_ok=True)

    # ── 1. Subset ──────────────────────────────────────────────────
    types = [cell_type] if isinstance(cell_type, str) else cell_type
    mask = (adata_all.obs['cell_type_grouped'].isin(types)) & \
           (adata_all.obs['batch'].isin(['GSE125527', 'GSE132465', 'GSE132257', 'GSE150115']))
    adata = adata_all[mask].copy()
    print(f"[1] {adata.n_obs:,} cells")
    print(pd.crosstab(adata.obs['batch'], adata.obs['cell_type_grouped']).to_string())

    # ── 2. KNN + UMAP ──────────────────────────────────────────────
    print("[2] KNN + UMAP ...")
    sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=30, random_state=0)
    sc.tl.umap(adata, min_dist=0.3, random_state=0)

    # ── 3. Diffusion map ───────────────────────────────────────────
    print("[3] Diffusion map ...")
    dm_res  = palantir.utils.run_diffusion_maps(
        pd.DataFrame(adata.obsm['X_scvi'], index=adata.obs_names),
        n_components=10, knn=30,
    )
    adata.obsm['X_diffmap'] = dm_res['EigenVectors'].values
    ms_data = palantir.utils.determine_multiscale_space(dm_res)

    # ── 4. Root = UC cell at extreme DC1 ──────────────────────────
    print("[4] Selecting root ...")
    uc_mask = (adata.obs['disease'] == 'UC') & \
              (adata.obs['cell_type_grouped'].isin(types))
    dc1       = dm_res['EigenVectors'].iloc[:, 0]
    root_cell = dc1[uc_mask].idxmin()
    print(f"   Root: {root_cell}  batch={adata.obs.loc[root_cell,'batch']}")

    # ── 5. Palantir ────────────────────────────────────────────────
    print("[5] Palantir pseudotime ...")
    pr_res = palantir.core.run_palantir(
        ms_data, early_cell=root_cell,
        num_waypoints=1200, use_early_cell_as_start=True, n_jobs=4,
    )
    adata.obs['palantir_pseudotime'] = pr_res.pseudotime
    adata.obs['entropy']             = pr_res.entropy
    terminal_ids = list(pr_res.branch_probs.columns)
    n_term       = len(terminal_ids)
    print(f"   Pseudotime: {pr_res.pseudotime.min():.3f} to {pr_res.pseudotime.max():.3f}")
    print(f"   Terminals ({n_term}):")
    for t in terminal_ids:
        print(f"     {t}  {adata.obs.loc[t,'cell_type_grouped']}  {adata.obs.loc[t,'batch']}")

    # Save core results
    pr_res.pseudotime.to_csv(f'{res_dir}/pseudotime_values.csv')
    pr_res.branch_probs.to_csv(f'{res_dir}/branch_probs.csv')
    adata.obs[['batch','cell_type_grouped','palantir_pseudotime','entropy']].to_csv(
        f'{res_dir}/cell_metadata.csv')

    # ── 6. Gene trends ─────────────────────────────────────────────
    print("[6] Gene trends ...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='batch')
    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    adata_hvg.obs['palantir_pseudotime'] = adata.obs['palantir_pseudotime']
    adata_hvg.obsm['palantir_fate_probabilities'] = \
        pr_res.branch_probs.reindex(adata.obs_names)

    palantir.presults.select_branch_cells(adata_hvg, q=0.01, eps=0.01)
    gene_trends = palantir.presults.compute_gene_trends(
        adata_hvg, lineages=terminal_ids, pseudo_time_key='palantir_pseudotime',
    )
    print(f"   gene_trends keys: {list(gene_trends.keys())}")

    # ── 7. Branch DEG (UC fate vs CRC fate, disease-aware) ────────
    deg1 = deg2 = None
    BATCH_UC  = {'GSE125527', 'GSE150115'}
    BATCH_CRC = {'GSE132465', 'GSE132257'}
    t_labels  = ['UC_fate', 'CRC_fate']

    if n_term >= 2:
        print("[7] Branch DEG ...")
        bp = pr_res.branch_probs.reindex(adata.obs_names)

        # Sum probabilities across all UC terminals and all CRC terminals
        uc_term_ids  = [t for t in terminal_ids
                        if adata.obs.loc[t, 'batch'] in BATCH_UC]
        crc_term_ids = [t for t in terminal_ids
                        if adata.obs.loc[t, 'batch'] in BATCH_CRC]
        print(f"   UC  terminals ({len(uc_term_ids)}): {uc_term_ids}")
        print(f"   CRC terminals ({len(crc_term_ids)}): {crc_term_ids}")

        p_uc  = bp[uc_term_ids].sum(axis=1).values  if uc_term_ids  else np.zeros(adata.n_obs)
        p_crc = bp[crc_term_ids].sum(axis=1).values if crc_term_ids else np.zeros(adata.n_obs)

        # Fallback: if one disease group has no terminals, compare T1 vs T2 directly
        if not uc_term_ids or not crc_term_ids:
            print("   [fallback] All terminals from one disease — comparing T1 vs T2 directly")
            t_labels = [
                f"{adata.obs.loc[terminal_ids[0],'batch']}_term1",
                f"{adata.obs.loc[terminal_ids[1],'batch']}_term2",
            ]
            p_uc  = bp[terminal_ids[0]].values
            p_crc = bp[terminal_ids[1]].values
            t_labels[0] = t_labels[0].replace('UC_fate', 'term1')

        thr_uc  = np.percentile(p_uc,  70)
        thr_crc = np.percentile(p_crc, 50)
        keep = (p_uc > thr_uc) | (p_crc > thr_crc)
        adata_bp = adata[keep].copy()
        adata_bp.obs['fate_label'] = np.where(
            p_uc[keep] > p_crc[keep], 'UC_fate', 'CRC_fate'
        )
        n1 = (adata_bp.obs['fate_label'] == 'UC_fate').sum()
        n2 = (adata_bp.obs['fate_label'] == 'CRC_fate').sum()
        print(f"   UC_fate={n1:,}  CRC_fate={n2:,}")

        for k in ['neighbors', 'diffmap_evals']:
            adata_bp.uns.pop(k, None)
        sc.tl.rank_genes_groups(
            adata_bp, groupby='fate_label',
            groups=['UC_fate', 'CRC_fate'], reference='rest',
            method='wilcoxon', key_added='branch_deg',
        )
        deg1 = sc.get.rank_genes_groups_df(adata_bp, group='UC_fate',  key='branch_deg')
        deg2 = sc.get.rank_genes_groups_df(adata_bp, group='CRC_fate', key='branch_deg')
        deg1 = deg1[deg1['logfoldchanges'] > 0].head(30)
        deg2 = deg2[deg2['logfoldchanges'] > 0].head(30)

        deg1.to_csv(f'{res_dir}/branch_DEG_UC_fate.csv',  index=False)
        deg2.to_csv(f'{res_dir}/branch_DEG_CRC_fate.csv', index=False)
        print(f"   Top UC_fate:  {', '.join(deg1['names'].head(8))}")
        print(f"   Top CRC_fate: {', '.join(deg2['names'].head(8))}")
    else:
        print(f"[7] Only {n_term} terminal(s) — skipping branch DEG")

    # ── 8. Plots ───────────────────────────────────────────────────
    print("[8] Plotting ...")
    umap  = adata.obsm['X_umap']
    pt    = adata.obs['palantir_pseudotime'].values
    ent   = adata.obs['entropy'].values
    rng   = np.random.default_rng(0)
    shuf  = rng.permutation(adata.n_obs)
    BATCH = {'GSE125527': '#00A087', 'GSE132465': '#4DBBD5'}
    S, A  = 1.5, 0.5

    # Overview: batch | pseudotime | entropy | branch_prob
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    ct_str = cell_type if isinstance(cell_type, str) else '+'.join(cell_type)
    fig.suptitle(f'Pseudotime — {ct_str}  (UC root → CRC terminal)',
                 fontsize=12, fontweight='bold')

    ax = axes[0, 0]
    for bat, col in BATCH.items():
        m = adata.obs['batch'].values == bat
        ax.scatter(umap[m, 0], umap[m, 1], c=col, s=S, alpha=A, linewidths=0, label=bat)
    ax.set_title('By dataset', fontsize=10, fontweight='bold')
    ax.axis('off')
    ax.legend(markerscale=4, fontsize=8)

    ax = axes[0, 1]
    sc_pt = ax.scatter(umap[shuf, 0], umap[shuf, 1], c=pt[shuf],
                       s=S, alpha=0.6, cmap='viridis', linewidths=0)
    plt.colorbar(sc_pt, ax=ax, shrink=0.8, label='Pseudotime')
    root_idx = adata.obs_names.get_loc(root_cell)
    ax.scatter(umap[root_idx, 0], umap[root_idx, 1],
               s=180, c='red', marker='*', zorder=10, label='Root (UC)')
    for t in terminal_ids:
        tidx = adata.obs_names.get_loc(t)
        ax.scatter(umap[tidx, 0], umap[tidx, 1],
                   s=120, c='orange', marker='D', zorder=10)
    ax.set_title('Pseudotime  (* root  D terminals)', fontsize=10, fontweight='bold')
    ax.axis('off')
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    sc_ent = ax.scatter(umap[shuf, 0], umap[shuf, 1], c=ent[shuf],
                        s=S, alpha=0.6, cmap='magma', linewidths=0)
    plt.colorbar(sc_ent, ax=ax, shrink=0.8, label='Entropy')
    ax.set_title('Branch entropy', fontsize=10, fontweight='bold')
    ax.axis('off')

    ax = axes[1, 1]
    bp_vals = pr_res.branch_probs.iloc[:, 0].reindex(adata.obs_names).values
    sc_bp = ax.scatter(umap[shuf, 0], umap[shuf, 1], c=bp_vals[shuf],
                       s=S, alpha=0.6, cmap='YlOrRd', vmin=0, vmax=1, linewidths=0)
    plt.colorbar(sc_bp, ax=ax, shrink=0.8)
    tlabel = f"{adata.obs.loc[terminal_ids[0],'cell_type_grouped']} ({adata.obs.loc[terminal_ids[0],'batch']})"
    ax.set_title(f'Branch prob → {tlabel}', fontsize=9, fontweight='bold')
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(f'{fig_dir}/pseudotime_overview.png', dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"   Saved: {fig_dir}/pseudotime_overview.png")

    # Pseudotime distribution violin (UC vs CRC)
    fig, ax = plt.subplots(figsize=(6, 4))
    data_v = [pt[adata.obs['batch'].values == b] for b in ['GSE125527', 'GSE132465']]
    vp = ax.violinplot(data_v, positions=[0, 1], showmedians=True, showextrema=False)
    for body, col in zip(vp['bodies'], ['#00A087', '#4DBBD5']):
        body.set_facecolor(col); body.set_alpha(0.75)
    vp['cmedians'].set_color('black'); vp['cmedians'].set_linewidth(1.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['UC (GSE125527)', 'CRC (GSE132465)'], fontsize=10)
    ax.set_ylabel('Pseudotime', fontsize=10)
    ax.set_title(f'{ct_str}: Pseudotime UC vs CRC', fontsize=11, fontweight='bold')
    uc_med = np.median(data_v[0]); crc_med = np.median(data_v[1])
    ax.text(0, uc_med + 0.02, f'med={uc_med:.3f}', ha='center', fontsize=8)
    ax.text(1, crc_med + 0.02, f'med={crc_med:.3f}', ha='center', fontsize=8)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(f'{fig_dir}/pseudotime_distribution.png', dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"   Saved: {fig_dir}/pseudotime_distribution.png")

    # Volcano (if DEG available)
    if n_term >= 2 and deg1 is not None:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
        fig.suptitle(f'Branch-point DEG: {t_labels[0]} vs {t_labels[1]}',
                     fontsize=11, fontweight='bold')
        for ax, grp, color, lbl in [
            (ax1, 'UC_fate',  '#E64B35', t_labels[0]),
            (ax2, 'CRC_fate', '#00A087', t_labels[1]),
        ]:
            df_all = sc.get.rank_genes_groups_df(adata_bp, group=grp, key='branch_deg')
            x = df_all['logfoldchanges'].values
            y = -np.log10(df_all['pvals_adj'].values + 1e-300)
            sig = (df_all['pvals_adj'] < 0.05) & (df_all['logfoldchanges'] > 0.5)
            ax.scatter(x[~sig], y[~sig], s=3, alpha=0.3, c='#cccccc', linewidths=0)
            ax.scatter(x[sig],  y[sig],  s=5, alpha=0.7, c=color, linewidths=0)
            top10 = df_all[sig].nlargest(10, 'logfoldchanges')
            for _, row in top10.iterrows():
                ax.annotate(row['names'],
                            (row['logfoldchanges'], -np.log10(row['pvals_adj'] + 1e-300)),
                            fontsize=7, xytext=(3, 2), textcoords='offset points')
            ax.axvline(0.5, color='#888', lw=0.8, ls='--')
            ax.axhline(2,   color='#888', lw=0.8, ls='--')
            ax.set_xlabel('log2 Fold Change', fontsize=9)
            ax.set_ylabel('-log10(adj p-value)', fontsize=9)
            ax.set_title(lbl, fontsize=9, fontweight='bold', color=color)
            ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
        plt.tight_layout()
        plt.savefig(f'{fig_dir}/branch_DEG_volcano.png', dpi=200, bbox_inches='tight', facecolor='white')
        plt.close()
        print(f"   Saved: {fig_dir}/branch_DEG_volcano.png")

    # Gene trend heatmap
    if n_term >= 2 and gene_trends:
        make_heatmap(gene_trends, terminal_ids, [deg1, deg2], t_labels,
                     f'{fig_dir}/gene_trend_heatmaps.png')

    print(f"[DONE] {tag}  |  {datetime.now().strftime('%H:%M:%S')}")


# ── Main ──────────────────────────────────────────────────────────
print("=" * 65)
print(f"Phase 2 — Lineage-specific pseudotime  |  {datetime.now().strftime('%H:%M:%S')}")
print("=" * 65)

print("\nLoading merged data ...")
adata_all = sc.read_h5ad('processed_data/merged_annotated.h5ad')
print(f"  Total: {adata_all.n_obs:,} cells")

run_lineage(adata_all, 'Myeloids', 'myeloid_v2')
run_lineage(adata_all, 'T cells',  'tcell_v2')

print("\n" + "=" * 65)
print(f"All done  |  {datetime.now().strftime('%H:%M:%S')}")
print("=" * 65)
