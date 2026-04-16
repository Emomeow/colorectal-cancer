# ============================================================
# scMetabolism Analysis — OLR1+ TAM metabolic profiling
# Uses the original R scMetabolism package (Wu et al. 2021)
# Input: processed_data/GSE144735.h5ad (via anndata R package)
# Output: phase4_results/scmetabolism_R/
# ============================================================

# ── Install dependencies (run once) ─────────────────────────
# install.packages("remotes")
# remotes::install_github("wu-yc/scMetabolism")
# BiocManager::install(c("GSEABase","AUCell","GSVA"))
# install.packages(c("anndata","ggplot2","pheatmap","RColorBrewer","dplyr","tidyr","scales"))

suppressPackageStartupMessages({
  library(anndata)
  library(scMetabolism)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(scales)
})

# ── Output directory ─────────────────────────────────────────
out_dir <- "phase4_results/scmetabolism_R"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("Output directory:", out_dir, "\n")

# ── Load h5ad ────────────────────────────────────────────────
cat("Loading GSE144735.h5ad ...\n")
adata <- read_h5ad("processed_data/GSE144735.h5ad")
cat(sprintf("  %d cells x %d genes\n", nrow(adata$obs), nrow(adata$var)))
cat("  Cell subtypes:\n")
print(table(adata$obs$Cell_subtype))

# ── Build count matrix for scMetabolism ─────────────────────
# scMetabolism expects: genes x cells, raw or normalized counts
# Use X (normalized) if available; scMetabolism internally handles
cat("Extracting expression matrix...\n")
X <- t(as.matrix(adata$X))   # genes x cells
rownames(X) <- rownames(adata$var)
colnames(X) <- rownames(adata$obs)
cat(sprintf("  Matrix: %d genes x %d cells\n", nrow(X), ncol(X)))

# ── Create Seurat-compatible object for scMetabolism ─────────
# scMetabolism works with Seurat objects
# We'll create a minimal Seurat object
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}
library(Seurat)

cat("Creating Seurat object...\n")
seurat_obj <- CreateSeuratObject(counts = X, project = "GSE144735")

# Add cell type metadata
seurat_obj$Cell_subtype <- adata$obs$Cell_subtype
seurat_obj$Cell_type    <- adata$obs$Cell_type

cat("Seurat object created:", ncol(seurat_obj), "cells\n")

# ── Run scMetabolism (AUCell mode) ────────────────────────────
cat("Running sc.metabolism.Seurat (AUCell on KEGG 2021)...\n")
cat("This may take several minutes...\n")

seurat_meta <- sc.metabolism.Seurat(
  obj            = seurat_obj,
  method         = "AUCell",          # scMetabolism AUCell mode
  imputation     = FALSE,             # no imputation (faster, less noise)
  ncores         = 4,
  metabolism.type = "KEGG"            # KEGG metabolic pathways
)

cat("scMetabolism done.\n")
cat(sprintf("  Pathways scored: %d\n", length(rownames(seurat_meta@assays$METABOLISM@data))))
cat("  Top 10 pathways:\n")
print(head(rownames(seurat_meta@assays$METABOLISM@data), 10))

# ── Extract scores ────────────────────────────────────────────
metab_scores <- as.data.frame(t(as.matrix(seurat_meta@assays$METABOLISM@data)))
metab_scores$Cell_subtype <- seurat_obj$Cell_subtype

cat("Score matrix dimensions:", dim(metab_scores), "\n")

# ── Focal pathways (same as Python AUCell analysis) ──────────
FOCAL_PATHWAYS <- c(
  "Fatty acid metabolism",
  "Steroid biosynthesis",
  "Glycolysis / Gluconeogenesis",
  "Oxidative phosphorylation",
  "Fatty acid biosynthesis",
  "Fatty acid degradation",
  "Purine metabolism",
  "Amino sugar and nucleotide sugar metabolism",
  "Citrate cycle (TCA cycle)",
  "Arachidonic acid metabolism"
)

# Find matching pathways (partial match)
all_pathways <- colnames(metab_scores)[colnames(metab_scores) != "Cell_subtype"]
focal_match <- c()
for (fp in FOCAL_PATHWAYS) {
  m <- grep(fp, all_pathways, value = TRUE, ignore.case = TRUE)
  if (length(m) > 0) focal_match <- c(focal_match, m[1])
}
focal_match <- unique(focal_match)
cat("Matched focal pathways:\n")
print(focal_match)

# ── Cell types of interest ───────────────────────────────────
CELL_TYPES <- c(
  "SPP1+B",
  "Anti-inflammatory",
  "CD8+ T cells",
  "Myofibroblasts",
  "Regulatory T cells"
)
CT_LABELS <- c(
  "OLR1+ TAMs (SPP1+B)",
  "Anti-inflam TAMs",
  "CD8+ T cells",
  "Myofibroblasts",
  "Reg T cells"
)

# ── Mean score per cell type ──────────────────────────────────
focal_scores <- metab_scores[, c(focal_match, "Cell_subtype")]
focal_scores <- focal_scores[focal_scores$Cell_subtype %in% CELL_TYPES, ]

mean_by_ct <- focal_scores %>%
  group_by(Cell_subtype) %>%
  summarise(across(all_of(focal_match), mean, .names = "{.col}")) %>%
  as.data.frame()
rownames(mean_by_ct) <- mean_by_ct$Cell_subtype
mean_by_ct$Cell_subtype <- NULL

# Reorder rows
mean_by_ct <- mean_by_ct[CELL_TYPES[CELL_TYPES %in% rownames(mean_by_ct)], , drop = FALSE]

# ── Z-score across cell types ─────────────────────────────────
zscore_mat <- as.data.frame(scale(t(mean_by_ct)))  # pathways x cell types
zscore_mat <- as.data.frame(t(scale(t(t(mean_by_ct)))))  # z-score across cell types per pathway
# Correct z-score direction: per pathway across cell types
score_t <- t(mean_by_ct)  # pathways x cell types
zscore_mat <- t(apply(score_t, 1, function(x) (x - mean(x)) / (sd(x) + 1e-9)))
zscore_df <- as.data.frame(zscore_mat)

cat("\nAUCell mean scores (raw):\n")
print(round(mean_by_ct, 4))
cat("\nZ-scores (pathways x cell types):\n")
print(round(zscore_df, 3))

# ── Save CSV ──────────────────────────────────────────────────
write.csv(mean_by_ct, file.path(out_dir, "scmetabolism_scores_raw.csv"))
write.csv(zscore_df,  file.path(out_dir, "scmetabolism_scores_zscore.csv"))
cat("Scores saved to CSV.\n")

# ── Panel D: Grouped bar chart ────────────────────────────────
cat("Plotting Panel D (grouped bar)...\n")
CT_COLORS <- c(
  "SPP1+B"              = "#E74C3C",
  "Anti-inflammatory"   = "#F1948A",
  "CD8+ T cells"        = "#2ECC71",
  "Myofibroblasts"      = "#3498DB",
  "Regulatory T cells"  = "#9B59B6"
)

# Reshape for ggplot
zscore_long <- as.data.frame(zscore_mat) %>%
  mutate(pathway = rownames(zscore_mat)) %>%
  tidyr::pivot_longer(-pathway, names_to = "cell_type", values_to = "zscore")
zscore_long$cell_type <- factor(zscore_long$cell_type, levels = CELL_TYPES)
zscore_long$pathway   <- factor(zscore_long$pathway,
                                 levels = rev(focal_match))

p_d <- ggplot(zscore_long, aes(x = pathway, y = zscore, fill = cell_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, alpha = 0.85, color = "white", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
  scale_fill_manual(values = CT_COLORS,
                    labels = CT_LABELS,
                    name = "Cell type") +
  coord_flip() +
  labs(x = NULL,
       y = "AUCell score (z-score across cell types)",
       title = "D  Metabolic Pathway Activity by Cell Type (GSE144735)",
       subtitle = "scMetabolism AUCell on KEGG pathways (R package, Wu et al. 2021)") +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "Arial"),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "#555555"),
    legend.position = "right",
    legend.text  = element_text(size = 9),
    legend.title = element_text(size = 9, face = "bold"),
    axis.text.y  = element_text(size = 9),
    axis.text.x  = element_text(size = 8),
    panel.grid.major.x = element_line(color = "#EEEEEE", linewidth = 0.4)
  )

ggsave(file.path(out_dir, "panelD_metabolism_R.png"),
       p_d, width = 11, height = 5.5, dpi = 200, bg = "white")
ggsave(file.path(out_dir, "panelD_metabolism_R.pdf"),
       p_d, width = 11, height = 5.5, device = cairo_pdf)
cat("Panel D saved.\n")

# ── Panel E: Heatmap ─────────────────────────────────────────
cat("Plotting Panel E (heatmap)...\n")
# Shorten pathway names for display
short_names <- sub("Glycolysis / Gluconeogenesis", "Glycolysis", focal_match)
short_names <- sub("Amino sugar and nucleotide sugar metabolism", "Amino Sugar Met.", short_names)
short_names <- sub("Citrate cycle \\(TCA cycle\\)", "TCA Cycle", short_names)
short_names <- sub("Oxidative phosphorylation", "Oxidative Phosphorylation", short_names)
short_names <- sub("Fatty acid metabolism", "Fatty Acid Metabolism", short_names)

hm_mat <- zscore_df
rownames(hm_mat) <- short_names[match(rownames(hm_mat), focal_match)]
colnames(hm_mat) <- CT_LABELS[match(colnames(hm_mat), CELL_TYPES)]

# Annotation for columns
col_ann <- data.frame(
  "Cell Type" = c("OLR1+ TAMs", "Anti-inflam", "CD8+ T", "Myofibroblast", "Reg T"),
  row.names = colnames(hm_mat),
  check.names = FALSE
)

# Highlight OLR1+ TAMs column
ann_colors <- list(
  "Cell Type" = c(
    "OLR1+ TAMs"   = "#E74C3C",
    "Anti-inflam"  = "#F1948A",
    "CD8+ T"       = "#2ECC71",
    "Myofibroblast"= "#3498DB",
    "Reg T"        = "#9B59B6"
  )
)

png(file.path(out_dir, "panelE_heatmap_R.png"),
    width = 7, height = 5.5, units = "in", res = 200, bg = "white")
pheatmap(
  as.matrix(hm_mat),
  color          = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
  breaks         = seq(-2.5, 2.5, length.out = 101),
  cluster_rows   = FALSE,
  cluster_cols   = FALSE,
  annotation_col = col_ann,
  annotation_colors = ann_colors,
  fontsize        = 9,
  fontfamily      = "Arial",
  main = "E  Metabolic Pathway Activity Heatmap: OLR1+ TAMs vs Other Cell Types\n(GSE144735 | scMetabolism R | AUCell on KEGG | z-score across cell types)",
  angle_col = 315,
  cellwidth = 50,
  cellheight = 28,
  border_color = "white"
)
dev.off()

# PDF
cairo_pdf(file.path(out_dir, "panelE_heatmap_R.pdf"),
          width = 7, height = 5.5, family = "Arial")
pheatmap(
  as.matrix(hm_mat),
  color          = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
  breaks         = seq(-2.5, 2.5, length.out = 101),
  cluster_rows   = FALSE,
  cluster_cols   = FALSE,
  annotation_col = col_ann,
  annotation_colors = ann_colors,
  fontsize        = 9,
  main = "E  Metabolic Pathway Activity Heatmap: OLR1+ TAMs vs Other Cell Types\n(GSE144735 | scMetabolism R | AUCell on KEGG | z-score across cell types)",
  angle_col = 315,
  cellwidth = 50,
  cellheight = 28,
  border_color = "white"
)
dev.off()
cat("Panel E saved.\n")

cat("\n=== scMetabolism R analysis complete ===\n")
cat("Outputs in:", out_dir, "\n")
cat("  panelD_metabolism_R.png/.pdf\n")
cat("  panelE_heatmap_R.png/.pdf\n")
cat("  scmetabolism_scores_raw.csv\n")
cat("  scmetabolism_scores_zscore.csv\n")
