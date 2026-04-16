# ============================================================
# Bulk RNA-seq Deconvolution — dtangle
# Reference: scRNA-seq mean profiles (GSE144735, 40 cell types)
# Target: TCGA-COAD bulk RNA-seq tumor samples
# Output: phase5_results/dtangle_proportions.csv
#         Figure 6 Panel A figures
# ============================================================
suppressPackageStartupMessages({
  library(dtangle)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
})

setwd("d:/changdao")
OUT_DATA <- "phase5_results"
OUT_FIG  <- "Project_Delivery/results/Figure6_Bulk_Clinical_Validation/Deconvolution/figures"
dir.create(OUT_FIG, recursive=TRUE, showWarnings=FALSE)

# ── 1. Load reference matrix (genes x cell types) ────────────
cat("Loading scRNA-seq reference matrix...\n")
ref <- read.csv(gzfile("phase5_results/for_dtangle/reference_matrix.csv.gz"),
                row.names=1, check.names=FALSE)
ref <- as.matrix(ref)
cat(sprintf("Reference: %d genes x %d cell types\n", nrow(ref), ncol(ref)))

# ── 2. Load TCGA-COAD bulk expression ────────────────────────
cat("Loading TCGA-COAD bulk expression...\n")
# Read compressed TSV
con <- gzfile("data/TCGA-COAD.star_counts.tsv.gz", "rt")
header <- readLines(con, n=1)
samples <- strsplit(header, "\t")[[1]][-1]
close(con)

bulk <- read.table(gzfile("data/TCGA-COAD.star_counts.tsv.gz"),
                   header=TRUE, row.names=1, sep="\t",
                   check.names=FALSE)
cat(sprintf("Bulk: %d genes x %d samples\n", nrow(bulk), ncol(bulk)))

# Strip Ensembl version numbers
rownames(bulk) <- sub("\\.[0-9]+$", "", rownames(bulk))

# ── 3. Map TCGA ENSG IDs → gene symbols ──────────────────────
cat("Mapping ENSG IDs to gene symbols...\n")
gene_map <- read.csv("phase5_results/for_dtangle/ensg_to_symbol.csv",
                     row.names=1, stringsAsFactors=FALSE)
ensg_to_sym <- setNames(gene_map$symbol, rownames(gene_map))

# Strip version from TCGA gene IDs and map to symbols
tcga_ensg_clean <- sub("\\.[0-9]+$", "", rownames(bulk))
tcga_sym <- ifelse(tcga_ensg_clean %in% names(ensg_to_sym),
                   ensg_to_sym[tcga_ensg_clean],
                   tcga_ensg_clean)
rownames(bulk) <- make.unique(tcga_sym)
# Remove duplicated symbols (keep first occurrence)
dup_mask <- duplicated(tcga_sym)
bulk <- bulk[!dup_mask, ]
rownames(bulk) <- tcga_sym[!dup_mask]
cat(sprintf("Bulk after symbol mapping: %d genes\n", nrow(bulk)))

common_genes <- intersect(rownames(bulk), rownames(ref))
cat(sprintf("Common genes with reference: %d\n", length(common_genes)))
if (length(common_genes) < 100) stop("Too few common genes.")

# ── 4. Filter to tumor samples (01A/01B) ─────────────────────
all_samples <- colnames(bulk)
tumor_samples <- all_samples[grepl("-01[AB]$", all_samples)]
cat(sprintf("Tumor samples: %d\n", length(tumor_samples)))

bulk_tumor <- bulk[common_genes, tumor_samples]
ref_common <- ref[common_genes, ]

cat(sprintf("Working matrix: %d genes x %d samples\n",
            nrow(bulk_tumor), ncol(bulk_tumor)))

# ── 5. Run dtangle deconvolution ─────────────────────────────
cat("\nRunning dtangle deconvolution...\n")
cat("(This estimates cell-type proportions from bulk RNA-seq)\n")

# dtangle input: rows = samples, cols = genes (for both mixture and reference)
Y <- t(as.matrix(bulk_tumor))    # samples x genes (mixture)
ref_t <- t(ref_common)           # cell_types x genes -> need genes x celltypes
# dtangle expects: Y = samples x genes, pure_samples = list of cell-type indices

# Use dtangle in marker-based mode with pure reference
# Combine: rows = [samples, cell_types], cols = genes
n_ct <- nrow(ref_t)
n_samp <- nrow(Y)
combined <- rbind(Y, ref_t)   # (n_samp + n_ct) x genes

# Cell type indices in combined matrix
pure_samples_idx <- lapply(seq(n_samp + 1, n_samp + n_ct), function(i) i)
names(pure_samples_idx) <- rownames(ref_t)

cat("Finding marker genes...\n")
markers <- find_markers(
  Y            = combined,
  pure_samples = pure_samples_idx,
  data_type    = "rna-seq",
  marker_method = "ratio"
)

cat("Running dtangle estimation...\n")
result <- dtangle(
  Y            = combined,
  pure_samples = pure_samples_idx,
  data_type    = "rna-seq",
  markers      = markers
)

# Extract proportions for bulk tumor samples only
props <- result$estimates[1:n_samp, ]
rownames(props) <- rownames(Y)
colnames(props) <- rownames(ref_t)

cat(sprintf("Proportions estimated: %d samples x %d cell types\n",
            nrow(props), ncol(props)))

# Save
write.csv(props, file.path(OUT_DATA, "dtangle_proportions.csv"))
cat("Proportions saved.\n")

# Summary statistics
cat("\nMean proportions across samples (top cell types):\n")
mean_props <- sort(colMeans(props), decreasing=TRUE)
print(round(head(mean_props, 12), 4))

# ── 6. OLR1+ TAM proportion (SPP1+B) ─────────────────────────
olr1_prop <- props[, "SPP1+B"]
cat(sprintf("\nOLR1+ TAM (SPP1+B) proportion: mean=%.4f, range=[%.4f, %.4f]\n",
            mean(olr1_prop), min(olr1_prop), max(olr1_prop)))

# ── 7. Load clinical data for Panel A ─────────────────────────
cat("\nLoading clinical/survival data...\n")
clin <- read.table("data/TCGA-COAD_clinicalMatrix.tsv", sep="\t",
                   header=TRUE, stringsAsFactors=FALSE, quote="")
surv_all <- read.table("data/TCGA-COAD_survival_pancan.tsv", sep="\t",
                       header=TRUE, stringsAsFactors=FALSE)
coad_surv <- surv_all[surv_all$cancer.type.abbreviation == "COAD", ]

# Patient barcode (first 12 chars) from sample names
sample_patient <- sub("(-[0-9]+[A-Z])+$", "", rownames(props))
# Standardise to TCGA-XX-XXXX
sample_patient <- substr(rownames(props), 1, 12)

prop_df <- as.data.frame(props)
prop_df$sample   <- rownames(props)
prop_df$patient  <- sample_patient

# Merge survival
coad_surv$patient <- substr(coad_surv$sample, 1, 12)
merged <- merge(prop_df, coad_surv[, c("patient","OS","OS.time")], by="patient", all.x=TRUE)

# Sidedness from clinical matrix
clin$patient <- substr(clin$sampleID, 1, 12)
right_locs <- c("Cecum","Ascending Colon","Transverse Colon","Hepatic Flexure")
left_locs  <- c("Sigmoid Colon","Descending Colon","Splenic Flexure","Rectosigmoid Junction")
clin$side <- ifelse(
  sapply(clin$anatomic_neoplasm_subdivision, function(x)
    any(sapply(right_locs, function(r) grepl(r, x, ignore.case=TRUE)))),
  "Right",
  ifelse(
    sapply(clin$anatomic_neoplasm_subdivision, function(x)
      any(sapply(left_locs,  function(l) grepl(l, x, ignore.case=TRUE)))),
    "Left", NA
  )
)

# Select available columns
clin_cols <- c("patient","side")
for (cn in c("MSI_updated_Oct62011","ajcc_pathologic_tumor_stage","age_at_initial_pathologic_diagnosis")) {
  if (cn %in% names(clin)) clin_cols <- c(clin_cols, cn)
}
merged <- merge(merged, clin[, clin_cols],
                by="patient", all.x=TRUE)
cat(sprintf("Merged samples with clinical: %d\n", nrow(merged)))
write.csv(merged, file.path(OUT_DATA, "dtangle_with_clinical.csv"), row.names=FALSE)

# ── 8. Panel A: Stacked bar — cell-type proportions sorted by OLR1+ TAM ──
cat("\nPlotting Panel A...\n")

# Top cell types to display (collapse rest into "Other")
top_cts <- names(sort(colMeans(props), decreasing=TRUE))[1:10]
if (!"SPP1+B" %in% top_cts) top_cts <- c("SPP1+B", top_cts[1:9])

plot_df <- as.data.frame(props[, top_cts])
plot_df$Other   <- 1 - rowSums(plot_df)
plot_df$Other   <- pmax(plot_df$Other, 0)
plot_df$sample  <- rownames(props)
plot_df$OLR1_prop <- props[, "SPP1+B"]

# Sort by OLR1+ TAM proportion
plot_df <- plot_df[order(plot_df$OLR1_prop), ]
plot_df$x_idx <- seq_len(nrow(plot_df))

# Sidedness annotation
side_vec <- merged$side[match(substr(plot_df$sample, 1, 12),
                               substr(merged$sample, 1, 12))]
plot_df$side <- side_vec

# Long format for ggplot
long_df <- plot_df %>%
  pivot_longer(cols = c(all_of(top_cts), "Other"),
               names_to = "cell_type", values_to = "proportion") %>%
  mutate(cell_type = factor(cell_type,
                             levels = c(top_cts, "Other")))

# Colors
ct_colors <- c(
  "SPP1+B"             = "#E74C3C",
  "Anti-inflammatory"  = "#F1948A",
  "Pro-inflammatory"   = "#E59866",
  "Epithelial cells"   = "#9B59B6",
  "Stromal 3"          = "#3498DB",
  "Myofibroblasts"     = "#2980B9",
  "CD8+ T cells"       = "#2ECC71",
  "CD4+ T cells"       = "#27AE60",
  "Regulatory T cells" = "#1ABC9C",
  "IgA+ Plasma"        = "#F39C12",
  "Other"              = "#D5D8DC"
)
used_colors <- ct_colors[levels(long_df$cell_type)]
used_colors[is.na(used_colors)] <- "#AAAAAA"

p_a <- ggplot(long_df, aes(x = x_idx, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 1, position = "stack") +
  scale_fill_manual(values = used_colors, name = "Cell type") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0,0)) +
  labs(
    x = paste0("TCGA-COAD tumor samples (n=", nrow(plot_df),
                ")\nsorted by OLR1+ TAM proportion"),
    y = "Estimated cell-type proportion",
    title = "A  Cellular Composition of TCGA-COAD Tumors (dtangle Deconvolution)",
    subtitle = paste0("Reference: GSE144735 scRNA-seq (40 cell types) | ",
                      "Bulk: TCGA-COAD star_counts | ",
                      "Red = OLR1+ TAMs (SPP1+B)")
  ) +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(family = "sans"),
    plot.title    = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 8, color = "#555555"),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    legend.position = "right",
    legend.text   = element_text(size = 8),
    legend.key.size = unit(0.45, "cm")
  )

ggsave(file.path(OUT_FIG, "panelA_deconvolution.png"),
       p_a, width=12, height=5.5, dpi=200, bg="white")
ggsave(file.path(OUT_FIG, "panelA_deconvolution.pdf"),
       p_a, width=12, height=5.5, device=cairo_pdf)
cat("Panel A saved.\n")

# ── 9. Supplementary: OLR1+ TAM proportion by sidedness ───────
if (any(!is.na(merged$side))) {
  side_df <- merged[!is.na(merged$side), ]
  cat(sprintf("\nOLR1+ TAM by side: Right n=%d, Left n=%d\n",
              sum(side_df$side=="Right"), sum(side_df$side=="Left")))

  r_v <- side_df$`SPP1+B`[side_df$side=="Right"]
  l_v <- side_df$`SPP1+B`[side_df$side=="Left"]
  wt  <- wilcox.test(r_v, l_v)
  cat(sprintf("Right mean=%.4f, Left mean=%.4f, Wilcoxon p=%.4f\n",
              mean(r_v, na.rm=TRUE), mean(l_v, na.rm=TRUE), wt$p.value))

  side_box <- ggplot(side_df[!is.na(side_df$`SPP1+B`),],
                     aes(x=side, y=`SPP1+B`, fill=side)) +
    geom_boxplot(alpha=0.8, outlier.size=0.8) +
    geom_jitter(width=0.2, size=1, alpha=0.5) +
    scale_fill_manual(values=c("Right"="#E74C3C","Left"="#3498DB"), guide="none") +
    annotate("text", x=1.5, y=max(side_df$`SPP1+B`, na.rm=TRUE)*1.05,
             label=sprintf("p=%.3f", wt$p.value), size=3.5) +
    labs(x="Anatomic side", y="OLR1+ TAM proportion (dtangle)",
         title="OLR1+ TAM Proportion by Sidedness\n(TCGA-COAD, dtangle deconvolution)") +
    theme_classic(base_size=10) +
    theme(text=element_text(family="sans"),
          plot.title=element_text(size=10, face="bold"))

  ggsave(file.path(OUT_FIG, "panelA_olr1_by_side.png"),
         side_box, width=5, height=5, dpi=200, bg="white")
  ggsave(file.path(OUT_FIG, "panelA_olr1_by_side.pdf"),
         side_box, width=5, height=5, device=cairo_pdf)
  cat("Sidedness boxplot saved.\n")
}

cat("\n=== dtangle deconvolution complete ===\n")
