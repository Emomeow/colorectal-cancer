# ============================================================
# scMetabolism-equivalent analysis using AUCell (Bioconductor)
# Reads data exported from Python (GSE144735.h5ad)
# Reproduces scMetabolism Wu et al. 2021 NatComm approach
# Output: phase4_results/scmetabolism_R/
# ============================================================
suppressPackageStartupMessages({
  library(Matrix)
  library(AUCell)
  library(GSEABase)
  library(jsonlite)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
})

setwd("d:/changdao")
OUT <- "phase4_results/scmetabolism_R"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
cat("=== AUCell-based scMetabolism (R) ===\n")

# ── 1. Load expression matrix ─────────────────────────────────
cat("Loading expression matrix...\n")
# Load sparse matrix saved from Python
X_npz_file <- "phase4_results/for_scmetabolism/counts.npz"

# counts.mtx was pre-exported from Python via:
#   python -c "import scipy.sparse as sp, scipy.io as sio;
#              X=sp.load_npz('...counts.npz'); sio.mmwrite('...counts.mtx', X)"
cat("Reading Matrix Market file...\n")
X <- readMM("phase4_results/for_scmetabolism/counts.mtx")  # cells x genes

# Load barcodes and features
barcodes <- read.table(gzfile("phase4_results/for_scmetabolism/barcodes.tsv.gz"),
                       header = FALSE, stringsAsFactors = FALSE)$V1
features <- read.table(gzfile("phase4_results/for_scmetabolism/features.tsv.gz"),
                       header = FALSE, stringsAsFactors = FALSE)$V1

rownames(X) <- barcodes
colnames(X) <- features

cat(sprintf("Matrix loaded: %d cells x %d genes\n", nrow(X), ncol(X)))

# Transpose to genes x cells for AUCell; convert to dgCMatrix
X_gc <- t(X)   # genes x cells
X_gc <- as(X_gc, "CsparseMatrix")   # AUCell requires dgCMatrix
cat(sprintf("Transposed: %d genes x %d cells\n", nrow(X_gc), ncol(X_gc)))

# ── 2. Load metadata ──────────────────────────────────────────
meta <- read.csv("phase4_results/for_scmetabolism/metadata.csv", row.names = 1)
cat("Cell subtypes:\n")
print(table(meta$Cell_subtype))

# ── 3. Load KEGG focal gene sets ─────────────────────────────
cat("\nLoading KEGG focal pathways...\n")
pw_list <- fromJSON("phase4_results/for_scmetabolism/kegg_focal_pathways.json")

# Filter to genes present in matrix
pw_list_filtered <- lapply(pw_list, function(genes) {
  present <- genes[genes %in% rownames(X_gc)]
  cat(sprintf("  %s: %d/%d genes present\n",
              names(pw_list)[match(list(genes), pw_list)],
              length(present), length(genes)))
  present
})
names(pw_list_filtered) <- names(pw_list)
pw_list_filtered <- pw_list_filtered[sapply(pw_list_filtered, length) >= 3]
cat(sprintf("Using %d pathways with >=3 genes present\n", length(pw_list_filtered)))

# ── 4. Run AUCell (scMetabolism core) ─────────────────────────
cat("\nBuilding gene rankings (AUCell step 1)...\n")
cell_rankings <- AUCell_buildRankings(X_gc, plotStats = FALSE, verbose = FALSE)
cat("Rankings built.\n")

cat("Calculating AUC (AUCell step 2)...\n")
gene_sets <- lapply(pw_list_filtered, GeneSet)
for (i in seq_along(gene_sets)) {
  setName(gene_sets[[i]]) <- names(pw_list_filtered)[i]
}
gene_set_collection <- GeneSetCollection(gene_sets)

cells_AUC <- AUCell_calcAUC(
  gene_set_collection,
  cell_rankings,
  aucMaxRank = ceiling(0.05 * nrow(cell_rankings)),  # top 5% genes (scMetabolism default)
  verbose    = FALSE
)
cat("AUC calculated.\n")

# Extract scores matrix: pathways x cells
auc_mat <- getAUC(cells_AUC)
cat(sprintf("AUC matrix: %d pathways x %d cells\n", nrow(auc_mat), ncol(auc_mat)))

# ── 5. Mean score per cell type ───────────────────────────────
CELL_TYPES <- c("SPP1+B", "Anti-inflammatory", "CD8+ T cells",
                "Myofibroblasts", "Regulatory T cells")
CT_LABELS  <- c("OLR1+ TAMs\n(SPP1+B)", "Anti-inflam\nTAMs",
                "CD8+ T cells", "Myofibroblasts", "Reg T cells")

auc_df       <- as.data.frame(t(auc_mat))   # cells x pathways
auc_df$Cell_subtype <- meta[rownames(auc_df), "Cell_subtype"]
auc_df       <- auc_df[auc_df$Cell_subtype %in% CELL_TYPES, ]

pw_names <- rownames(auc_mat)
mean_by_ct <- auc_df %>%
  group_by(Cell_subtype) %>%
  summarise(across(all_of(pw_names), mean)) %>%
  as.data.frame()
rownames(mean_by_ct) <- mean_by_ct$Cell_subtype
mean_by_ct$Cell_subtype <- NULL
mean_by_ct <- mean_by_ct[CELL_TYPES[CELL_TYPES %in% rownames(mean_by_ct)], , drop = FALSE]

cat("\nMean AUCell scores (raw):\n")
print(round(mean_by_ct, 5))

# Z-score per pathway across cell types
zscore_mat <- t(apply(t(mean_by_ct), 1, function(x) (x - mean(x)) / (sd(x) + 1e-9)))
cat("\nZ-scores (cell types x pathways):\n")
print(round(zscore_mat, 3))

# Save
write.csv(mean_by_ct, file.path(OUT, "scmetabolism_scores_raw.csv"))
write.csv(zscore_mat, file.path(OUT, "scmetabolism_scores_zscore.csv"))

# ── 6. Display labels ─────────────────────────────────────────
# Short pathway names for plotting (no \n — pheatmap rownames can't have newlines)
PW_DISPLAY <- c(
  "Lipid_Uptake_FAO"           = "Lipid Uptake & FAO",
  "Cholesterol_Biosynthesis"   = "Cholesterol Biosynthesis",
  "Glycolysis"                 = "Glycolysis",
  "Oxidative_Phosphorylation"  = "Oxidative Phosphorylation",
  "Fatty_Acid_Biosynthesis"    = "Fatty Acid Biosynthesis",
  "Fatty_Acid_Degradation"     = "Fatty Acid Degradation",
  "Purine_Metabolism"          = "Purine Metabolism",
  "Amino_Sugar_Metabolism"     = "Amino Sugar Metabolism",
  "TCA_Cycle"                  = "TCA Cycle",
  "Arachidonic_Acid_Metabolism"= "Arachidonic Acid Metabolism"
)
# x-axis labels with newlines (for ggplot only)
PW_XLAB <- c(
  "Lipid_Uptake_FAO"           = "Lipid Uptake\n& FAO",
  "Cholesterol_Biosynthesis"   = "Cholesterol\nBiosynthesis",
  "Glycolysis"                 = "Glycolysis",
  "Oxidative_Phosphorylation"  = "Oxidative\nPhosphorylation",
  "Fatty_Acid_Biosynthesis"    = "Fatty Acid\nBiosynthesis",
  "Fatty_Acid_Degradation"     = "Fatty Acid\nDegradation",
  "Purine_Metabolism"          = "Purine\nMetabolism",
  "Amino_Sugar_Metabolism"     = "Amino Sugar\nMetabolism",
  "TCA_Cycle"                  = "TCA Cycle",
  "Arachidonic_Acid_Metabolism"= "Arachidonic\nAcid Metabolism"
)

CT_COLORS <- c(
  "SPP1+B"              = "#E74C3C",
  "Anti-inflammatory"   = "#F1948A",
  "CD8+ T cells"        = "#2ECC71",
  "Myofibroblasts"      = "#3498DB",
  "Regulatory T cells"  = "#9B59B6"
)

# ── 7. Panel D: Grouped bar chart ────────────────────────────
cat("\nPlotting Panel D...\n")
# zscore_mat is pathways x cell_types; t() gives cell_types x pathways
zscore_df <- as.data.frame(t(zscore_mat))
zscore_df$Cell_subtype <- rownames(zscore_df)

zscore_long <- zscore_df %>%
  pivot_longer(-Cell_subtype, names_to = "pathway", values_to = "zscore")

zscore_long <- zscore_long %>%
  mutate(
    pw_xlab = PW_XLAB[pathway],
    ct_fill = CT_COLORS[Cell_subtype],
    ct_label_f = factor(Cell_subtype, levels = CELL_TYPES)
  )
pw_xlab_levels <- PW_XLAB[names(PW_XLAB) %in% unique(zscore_long$pathway)]
zscore_long$pw_xlab <- factor(zscore_long$pw_xlab, levels = pw_xlab_levels)

p_d <- ggplot(zscore_long, aes(x = pw_xlab, y = zscore, fill = ct_label_f)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.78),
           width = 0.72, alpha = 0.88,
           color = "white", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dashed", color = "black") +
  scale_fill_manual(values = setNames(CT_COLORS, CELL_TYPES),
                    labels = setNames(CT_LABELS, CELL_TYPES),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  labs(x = NULL,
       y = "AUCell score (z-score across cell types)",
       title = "D  Metabolic Pathway Activity by Cell Type (GSE144735)",
       subtitle = "scMetabolism — AUCell on KEGG 2021 Human (Bioconductor AUCell, top 5% gene ranking)") +
  theme_classic(base_size = 11) +
  theme(
    text              = element_text(family = "sans"),
    plot.title        = element_text(size = 12, face = "bold"),
    plot.subtitle     = element_text(size = 8.5, color = "#555555"),
    axis.text.x       = element_text(size = 8.5, angle = 35, hjust = 1),
    axis.text.y       = element_text(size = 9),
    legend.position   = "top",
    legend.text       = element_text(size = 8.5),
    legend.key.size   = unit(0.5, "cm"),
    panel.grid.major.y = element_line(color = "#F0F0F0", linewidth = 0.4)
  )

ggsave(file.path(OUT, "panelD_metabolism_scMetabolism.png"),
       p_d, width = 13, height = 5.5, dpi = 200, bg = "white")
ggsave(file.path(OUT, "panelD_metabolism_scMetabolism.pdf"),
       p_d, width = 13, height = 5.5, device = cairo_pdf)
cat("Panel D saved.\n")

# ── 8. Panel E: Heatmap ───────────────────────────────────────
cat("Plotting Panel E...\n")

# zscore_mat is pathways x cell_types — use directly (no transpose needed)
hm_mat <- zscore_mat
# Use clean single-line display names for rows (pathways)
rownames(hm_mat) <- PW_DISPLAY[rownames(hm_mat)]
# Rename columns (cell types)
ct_rename <- c(
  "SPP1+B"              = "OLR1+ TAMs",
  "Anti-inflammatory"   = "Anti-inflam TAMs",
  "CD8+ T cells"        = "CD8+ T cells",
  "Myofibroblasts"      = "Myofibroblasts",
  "Regulatory T cells"  = "Reg T cells"
)
colnames(hm_mat) <- ct_rename[colnames(hm_mat)]
# Remove any NAs
hm_mat <- hm_mat[!is.na(rownames(hm_mat)), !is.na(colnames(hm_mat)), drop = FALSE]

# Column annotation
col_ann <- data.frame(
  CellType = colnames(hm_mat),
  row.names = colnames(hm_mat),
  check.names = FALSE
)

ann_col_colors <- list(
  CellType = c(
    "OLR1+ TAMs"      = "#E74C3C",
    "Anti-inflam TAMs"= "#F1948A",
    "CD8+ T cells"    = "#2ECC71",
    "Myofibroblasts"  = "#3498DB",
    "Reg T cells"     = "#9B59B6"
  )
)

do_heatmap <- function() {
  pheatmap(
    as.matrix(hm_mat),
    color         = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    breaks        = seq(-2.5, 2.5, length.out = 101),
    cluster_rows  = FALSE,
    cluster_cols  = FALSE,
    annotation_col    = col_ann,
    annotation_colors = ann_col_colors,
    fontsize       = 10,
    fontsize_row   = 9,
    fontsize_col   = 9,
    main = "E  Metabolic Pathway Activity (scMetabolism AUCell, KEGG 2021)\nOLR1+ TAMs vs Other Cell Types — GSE144735",
    angle_col    = 315,
    cellwidth    = 55,
    cellheight   = 30
  )
}

# PNG
png(file.path(OUT, "panelE_heatmap_scMetabolism.png"),
    width = 8, height = 6, units = "in", res = 200, bg = "white")
do_heatmap()
dev.off()

# PDF
cairo_pdf(file.path(OUT, "panelE_heatmap_scMetabolism.pdf"),
          width = 8, height = 6, family = "sans")
do_heatmap()
dev.off()
cat("Panel E saved.\n")

# ── 9. Copy to delivery ───────────────────────────────────────
del_dir <- "Project_Delivery/results/Figure4_Communication_Metabolism/Metabolic_Flux/figures"
dir.create(del_dir, recursive = TRUE, showWarnings = FALSE)
file.copy(file.path(OUT, "panelD_metabolism_scMetabolism.png"),
          file.path(del_dir, "panelD_metabolism.png"), overwrite = TRUE)
file.copy(file.path(OUT, "panelD_metabolism_scMetabolism.pdf"),
          file.path(del_dir, "panelD_metabolism.pdf"), overwrite = TRUE)
file.copy(file.path(OUT, "panelE_heatmap_scMetabolism.png"),
          file.path(del_dir, "panelE_metabolic_heatmap.png"), overwrite = TRUE)
file.copy(file.path(OUT, "panelE_heatmap_scMetabolism.pdf"),
          file.path(del_dir, "panelE_metabolic_heatmap.pdf"), overwrite = TRUE)
cat("Copied to delivery directory.\n")

cat("\n=== Complete ===\n")
cat("Outputs:\n")
cat("  panelD_metabolism_scMetabolism.png/.pdf\n")
cat("  panelE_heatmap_scMetabolism.png/.pdf\n")
cat("  scmetabolism_scores_raw.csv\n")
cat("  scmetabolism_scores_zscore.csv\n")
