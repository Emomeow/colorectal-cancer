# Colorectal Cancer Immune Microenvironment — Single-Cell RNA-seq Analysis

**Project**: Three-phase scRNA-seq study of immune resistance mechanisms in microsatellite-stable (MMRp) colorectal cancer (CRC) versus ulcerative colitis (UC)
**Date**: 2026-03-15

---

## Overview

This project integrates five multi-institutional scRNA-seq cohorts to characterize the tumor microenvironment (TME) in CRC and UC, identify suppressive "pre-driver" immune populations driving MMRp immune resistance, and build a patient-level immune classifier using attention-based multiple instance learning (MIL).

| Phase | Focus | Key Method |
|-------|-------|------------|
| Phase 1 | Multi-cohort integration & cell type annotation | scVI / scANVI / CellAssign |
| Phase 2 | TME characterization & differentiation trajectory | Differential abundance / Palantir / decoupleR |
| Phase 3 | Patient-level MMRp/d stratification | Gated Attention-MIL (ABMIL) |

---

## Datasets

| Dataset | Disease | Cells (raw) | Cells (after QC) | Source |
|---------|---------|-------------|------------------|--------|
| GSE132465 | CRC | 63,689 | 54,487 | GEO |
| GSE132257 | CRC | 18,409 | 15,642 | GEO |
| cellxgene | CRC | 47,107 | 9,003 | CellxGene |
| GSE125527 | UC | 84,129 | 71,866 | GEO |
| GSE150115 | UC | 32,000 | 2,138 | GEO |
| **Total** | CRC + UC | **245,334** | **153,136** | |
| GSE178341 | CRC (MIL training) | 94,196 immune cells | — | GEO |
| GSE132465 (SMC) | CRC (MIL validation) | 63,689 | — | GEO |

QC filters: mitochondrial gene fraction ≤ 10%, Scrublet doublet removal. GSE150115 additionally required ≥ 500 genes and ≥ 1,000 UMIs.

---

## Phase 1 — Data Integration & Cell Type Annotation

**Goal**: Integrate five heterogeneous cohorts into a harmonized single-cell atlas while preserving biological variation.

**Methods**:
- Batch-aware HVG selection (3,000 genes)
- **scVI** (n_latent=30, n_layers=2, batch_key=`batch`, disease covariate) trained with early stopping (patience=50, ~96 epochs)
- **scANVI** semi-supervised refinement (labels_key=`cellass_type`, max_epochs=50)
- Integration quality assessed via **iLISI / cLISI** (k=90): scANVI achieved cLISI=1.21 vs scVI 1.33, confirming superior cell-type separation
- **CellAssign** probabilistic marker-driven annotation of 10 immune and stromal lineages

**Cell type composition** (153,136 cells):

| Cell Type | Count | % |
|-----------|-------|---|
| T cells | 77,927 | 50.9% |
| B cells | 36,987 | 24.2% |
| Epithelial | 9,003 | 5.9% |
| Myeloid | 15,972 | 10.4% |
| NK cells | 7,224 | 4.7% |
| Others | 6,023 | 3.9% |

---

## Phase 2 — TME Characterization & Trajectory Analysis

**Goal**: Identify CRC-specific suppressive immune subpopulations and reconstruct their differentiation trajectory.

**Methods**:
- Fine annotation of T cells and myeloid cells using Leiden clustering, marker gene scoring, and KNN label transfer
- **Differential abundance**: per-patient subtype proportions, Mann-Whitney U test + BH FDR correction (CRC n=25 vs UC n=19)
- **Batch-corrected DEG**: scVI counterfactual batch marginalization (transform across all 5 batches) + Wilcoxon / BH FDR
- **Palantir** pseudotime: Classical Monocyte → TAM differentiation trajectory
- **TF regulon analysis**: decoupleR activity scores along pseudotime; identified core TFs driving suppressive TAM fate

**Key findings**:

| Subtype | CRC enrichment | FDR |
|---------|---------------|-----|
| TAM | **75.5×** (log2FC = +5.68) | 4.1×10⁻⁷ |
| CD8_Tex | **13.8×** (log2FC = +3.87) | 9.4×10⁻⁶ |
| CD8_Tpex | 0.05× — UC-enriched (log2FC = −4.16) | 0.022 |

TAM and CD8_Tex are established as the dominant **Pre-driver suppressive populations** in CRC. The CRC Tex/Tpex ratio = 12.8 vs UC = 0.04, indicating a deeply exhausted rather than stem-like T cell state in CRC.

---

## Phase 3 — Attention-MIL for Patient Stratification

**Goal**: Classify MMRp vs MMRd patients from single-cell immune profiles and identify the immune populations most predictive of mismatch repair status.

**Methods**:
- **Gated Attention-MIL (ABMIL)** in PyTorch: encoder (Linear 50→256→256) + gated attention + classification head
- Training: GSE178341, 62 patients (MMRd=34, MMRp=28), 94,196 immune cells; features = PCA top 50 dims of 2,000 HVGs
- 5-fold stratified cross-validation, AdamW optimizer, BCE loss, early stopping (patience=20)
- **Cross-dataset transfer**: retrained on 766 shared HVGs between GSE178341 and GSE132465; PCA projection applied to validation cohort
- **Augur-like predictability**: per-cell-type random forest AUC scores to rank immune populations by MMRp/d discriminability
- Paired tumor vs normal tissue validation (n=10 pairs, TCD8 exhaustion signature, paired Wilcoxon test)
- **Hyperparameter robustness**: one-at-a-time sweep across 6 parameters (n_PCA, hidden dim, attention dim, max cells, dropout, LR)

**Results**:

| Metric | Value |
|--------|-------|
| Out-of-fold AUC | **0.909** |
| Out-of-fold AUPRC | **0.922** |
| Fold AUCs | 1.000 / 0.929 / 1.000 / 1.000 / 0.829 |

---

## Repository Structure

```
Project_Delivery/
├── data/
│   └── metadata/                        # Batch statistics, cell type counts (CSV)
│       ├── merged_unified_batch_stats.csv
│       ├── merged_unified_celltype_by_batch.csv
│       ├── merged_unified_summary.csv
│       └── preprocessing_summary.csv
│
├── results/
│   ├── Figure1_scVI_Atlas/
│   │   ├── scVI_Integration/            # UMAP figures (panel A/B), integration tables
│   │   ├── Cell_Composition/            # Cell type composition (panel C)
│   │   ├── Marker_Expression/           # Marker dotplot (panel D)
│   │   ├── Cluster_Stability/           # LISI / bootstrap AMI (panel E/F)
│   │   ├── figure1_landscape.png/pdf    # Full Figure 1 assembly
│   │   └── script/                      # Phase 1 analysis scripts
│   │
│   ├── Figure2_Pre_driver_Trajectory/
│   │   ├── Trajectory_Analysis/         # UMAP subtypes, pseudotime, gene expression
│   │   ├── Spatial_Distribution/        # Abundance strip plots, dotplots
│   │   ├── SCENIC_Regulon/              # TF heatmaps, pseudotime correlations
│   │   └── script/                      # Phase 2 analysis scripts
│   │
│   ├── Figure3_MIL_MMR_Prediction/
│   │   ├── Model_Performance/           # ROC / PR curves, fold AUCs
│   │   ├── Attention_Analysis/          # Attention score UMAP, per-subtype scores
│   │   ├── Cross_Dataset_Validation/    # GSE132465 transfer results
│   │   └── script/                      # Phase 3 analysis scripts
│   │
│   └── Supplementary/
│       ├── Supp2_Extended_Markers/      # Extended marker dotplots
│       └── Supp3_MIL_Robustness/        # Hyperparameter sweep results
│
├── script/                              # All analysis scripts (master copy)
├── report_phase1.pdf                    # Phase 1 full report
├── report_phase2.pdf                    # Phase 2 full report
└── report_phase3.pdf                    # Phase 3 full report
```

> **Note**: Processed `.h5ad` data files (3.2 GB total) are excluded from this repository. Raw data is available from GEO (accessions listed above) and CellxGene.

---

## Environment

Python 3.12 with the following key packages:

```bash
pip install -r requirements.txt
```

| Package | Version |
|---------|---------|
| scvi-tools | 1.4.1 |
| scanpy | ≥1.9 |
| palantir | ≥1.3 |
| decoupler | ≥1.4 |
| torch | ≥2.0 |
| scikit-learn | ≥1.3 |

---

## Citation / Contact

Yue Zhan — yzhan@unmc.edu
University of Nebraska Medical Center, Department of Biostatistics
