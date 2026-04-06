# Colorectal Cancer Immune Microenvironment — Single-Cell RNA-seq Analysis

**Project**: Five-phase scRNA-seq study of immune resistance mechanisms in microsatellite-stable (MMRp) colorectal cancer (CRC) versus ulcerative colitis (UC)  
**Author**: Yue Zhan, University of Nebraska Medical Center  
**Contact**: yzhan@unmc.edu

---

## Overview

This project integrates five multi-institutional scRNA-seq cohorts to characterize the tumor microenvironment (TME) in CRC and UC, identify suppressive pre-driver immune populations driving MMRp immune evasion, dissect their intercellular communication with stromal cells, and validate clinical significance in independent bulk-RNA cohorts.

**Core story**: Chronic inflammation (UC) → CRC evolution drives **OLR1⁺ lipid-associated TAM** polarization. These TAMs crosstalk with **FAP⁺ myofibroblasts** at the tumor-stroma boundary via SPP1–CD44 and LGALS9–HAVCR2 signaling, creating a dual physical (ECM remodeling) and immune (T cell exclusion) barrier — the primary mechanism of MMRp CRC immune exclusion.

| Phase | Focus | Key Method | Key Result |
|-------|-------|------------|------------|
| 1 | Multi-cohort integration & cell annotation | scVI / scANVI / CellAssign | 153,136 cells, 10 lineages |
| 2 | TME characterization & trajectory | Palantir / decoupleR / DEG | OLR1⁺ TAMs as pre-driver population |
| 3 | Patient-level MMRp/d stratification | Gated Attention-MIL (ABMIL) | AUC = 0.909 |
| 4 | Cell communication & metabolic mechanism | LIANA / NicheNet / ssGSEA | SPP1–CD44 & LGALS9–HAVCR2 top axes |
| 5 | External clinical validation | GSE188711 / TCGA-COAD Cox | Right-sided CRC HR = 1.74, p = 0.001 |

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
| GSE144735 | CRC (border cells, Phase 4) | 9,424 | — | GEO |
| GSE188711 | CRC left/right (Phase 5) | 35,693 | 26,218 | GEO |
| TCGA-COAD | CRC bulk (Phase 5) | 514 tumor samples | 747 with OS | UCSC Xena |

QC filters: mitochondrial gene fraction ≤ 10%, Scrublet doublet removal.

---

## Phase 1 — Data Integration & Cell Type Annotation

**Goal**: Integrate five heterogeneous cohorts into a harmonized single-cell atlas while preserving biological variation.

**Methods**:
- Batch-aware HVG selection (3,000 genes)
- **scVI** (n_latent=30, n_layers=2, batch_key=`batch`, disease covariate) trained with early stopping (~96 epochs)
- **scANVI** semi-supervised refinement (labels_key=`cellass_type`, max_epochs=50)
- Integration quality: scANVI cLISI=1.21 vs scVI 1.33 (k=90), confirming superior cell-type separation
- **CellAssign** probabilistic marker-driven annotation of 10 immune and stromal lineages

**Cell type composition** (153,136 cells):

| Cell Type | Count | % |
|-----------|-------|---|
| T cells | 77,927 | 50.9% |
| Myeloid | 15,972 | 10.4% |
| B cells | 14,618 | 9.6% |
| Epithelial | 13,587 | 8.9% |
| Fibroblasts | 12,741 | 8.3% |
| NK cells | 8,923 | 5.8% |
| Endothelial | 5,203 | 3.4% |
| Plasma | 3,107 | 2.0% |
| Mast cells | 695 | 0.5% |
| pDC | 363 | 0.2% |

**Key figures**: `Project_Delivery/results/Figure1_Integration_Annotation/`

---

## Phase 2 — TME Characterization & Trajectory Analysis

**Goal**: Fine-annotate immune subpopulations, identify CRC-specific suppressive pre-driver populations, and reconstruct myeloid differentiation trajectories.

**Methods**:
- Fine annotation of T cells (77,927) and myeloid cells (15,972) into functional subtypes
- Differential abundance analysis (Mann-Whitney U + BH FDR) to identify CRC-specific enrichment
- scVI-based batch-corrected DEG (counterfactual batch marginalization, Wilcoxon + BH FDR)
- **Palantir** pseudotime reconstruction of Classical_Mono → TAM differentiation trajectory
- **decoupleR** + **SCENIC** TF regulon activity scoring to identify drivers of suppressive TAM fate

**Key findings**:
- **OLR1⁺ lipid-associated TAMs** (OLR1⁺ MMP9⁺ SPP1⁺) identified as the primary CRC-specific suppressive pre-driver population (75.5× enrichment vs UC)
- CD8_Tex (exhausted T cells) co-enriched in CRC (13.8×)
- Palantir trajectory confirms Classical_Mono → intermediate → OLR1⁺ TAM fate, driven by PPARG, MAF, BHLHE41 regulons

**Key figures**: `Project_Delivery/results/Figure2_Pre_driver_Trajectory/`

---

## Phase 3 — Attention-MIL Patient Stratification

**Goal**: Build a patient-level classifier to distinguish MMRp (immune-excluded) from MMRd (immune-inflamed) CRC using single-cell immune profiles.

**Methods**:
- **Gated Attention-MIL (ABMIL)** in PyTorch; each patient treated as a "bag" of immune cell embeddings
- Training: GSE178341 (62 patients, 94,196 cells); 5-fold stratified cross-validation
- Feature alignment via PCA on 766 shared HVGs for transfer to independent validation cohort
- Augur-like per-cell-type predictability scoring to rank immune populations by discriminability
- Hyperparameter robustness benchmarked across 6 parameters

**Key results**:
- OOF AUC = **0.909**, AUPRC = **0.922**
- OLR1⁺ TAMs and CD8_Tex carry highest attention weights in MMRp patients
- CD8 exhaustion signature validated in paired tumor vs. normal tissue (paired Wilcoxon, n=10, p < 0.001)

**Key figures**: `Project_Delivery/results/Figure3_MIL_MMR_Prediction/`

---

## Phase 4 — Cell Communication & Metabolic Mechanism Dissection

**Goal**: Dissect how OLR1⁺ TAMs and FAP⁺ myofibroblasts cooperate to build the dual physical–immune barrier.

**Methods**:
- **LIANA rank_aggregate** on 9,424 tumor-border cells (GSE144735) across 10 immune/stromal lineages
- **NicheNet** ligand prioritization: FAP⁺ myofibroblast → OLR1⁺ TAM signaling axis
- ssGSEA-like pathway scoring for five effector programs (matrix remodeling, T cell exclusion, lipid metabolism, angiogenesis, immune suppression)

**Key findings**:
- Top immunosuppressive axes: **SPP1–CD44** and **LGALS9–HAVCR2** (OLR1⁺ TAMs → T cells / NK cells)
- NicheNet prioritizes MMP2, TGFB1, POSTN from FAP⁺ myofibroblasts as drivers of OLR1⁺ TAM ECM remodeling
- OLR1⁺ TAM × FAP⁺ myofibroblast crosstalk creates a dual barrier:
  - **Physical**: MMP9/COL1A1/FN1-mediated ECM remodeling → T cell physical exclusion
  - **Immune**: LGALS9/CD274/VEGFA-mediated checkpoint suppression

**Key figures**: `Project_Delivery/results/Figure4_Communication_Metabolism/`

---

## Phase 5 — External Clinical Validation

**Goal**: Validate OLR1⁺ TAM findings in independent cohorts and assess clinical prognostic impact.

### GSE188711 — Anatomic Sidedness Validation (n=6 samples, 26,218 cells after QC)

- OLR1⁺ TAM signature reproducibly detected across independent left/right CRC samples
- High-OLR1⁺ cells show significantly elevated T cell exclusion (log2FC = −2.08, p < 10⁻¹⁰⁸) and matrix remodeling (log2FC = −1.77, p < 10⁻⁸⁵) pathway activity

### TCGA-COAD — Bulk Cohort Survival Analysis (n=747)

- OLR1⁺ TAM ssGSEA signature scored across 470 tumor samples (all 20 signature genes present)
- Right-sided CRC (anatomic surrogate for OLR1⁺ TAM-enriched TME) independently predicts worse overall survival:
  - Log-rank p = 0.0011
  - **Cox HR = 1.74** (95% CI: 1.25–2.43), p = 0.001, adjusting for stage, age, sex, and MSI status

| Cox Covariate | HR | 95% CI | p |
|---------------|----|--------|---|
| OLR1⁺ TAM score (bulk ssGSEA) | 1.02 | 0.88–1.18 | 0.79 |
| Tumor stage (I→IV) | **2.28** | 1.89–2.74 | <0.001 |
| **Right-sided CRC** | **1.74** | **1.25–2.43** | **0.001** |
| Age | 1.04 | 1.02–1.05 | <0.001 |
| Female sex | 1.01 | 0.73–1.39 | 0.97 |
| MSI-H | 1.20 | 0.69–2.10 | 0.52 |

**Key figures**: `Project_Delivery/results/Figure5_Anatomical_Heterogeneity/`, `Project_Delivery/results/Figure6_Bulk_Clinical_Validation/`

---

## Repository Structure

```
Project_Delivery/
├── results/
│   ├── Figure1_Integration_Annotation/       # Phase 1: UMAP, cell type composition
│   ├── Figure2_Pre_driver_Trajectory/        # Phase 2: trajectory, SCENIC, spatial distribution
│   ├── Figure3_MIL_MMR_Prediction/           # Phase 3: ABMIL attention maps, AUC curves
│   ├── Figure4_Communication_Metabolism/     # Phase 4: LIANA, NicheNet, pathway scores
│   ├── Figure5_Anatomical_Heterogeneity/     # Phase 5a: GSE188711 left/right CRC
│   ├── Figure6_Bulk_Clinical_Validation/     # Phase 5b: TCGA-COAD KM curves + Cox
│   └── Supplementary/                        # Supplementary Figures 1–5
└── script/                                   # Core analysis scripts
```

---

## Environment

```
Python 3.10+
scvi-tools==1.4.1    scanpy==1.10    anndata
liana                decoupler-py    palantir
lifelines            scipy           matplotlib
PyTorch >= 2.0
```

---

## Author

**Yue Zhan** — PhD Candidate, Department of Biostatistics, University of Nebraska Medical Center  
yzhan@unmc.edu | [GitHub: Emomeow](https://github.com/Emomeow)
