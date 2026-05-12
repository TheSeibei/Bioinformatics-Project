# scripts/04_differential_expression.R

message("Starting 04_differential_expression.R")

# ---- libraries ----
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

required_pkgs <- c("limma", "dplyr", "tibble", "rio", "ggplot2")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("limma")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

# ---- directories ----
resultsdir <- here("results")
plotsdir <- here("plots")

dir.create(resultsdir, showWarnings = FALSE, recursive = TRUE)
dir.create(plotsdir, showWarnings = FALSE, recursive = TRUE)

# ---- helper functions ----

clean_group_factor <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA
  factor(x, levels = c("CONTROL", "MCI", "AD"))
}

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_\\-]", "_", x)
}

plot_volcano <- function(tt, outfile, title_txt, lfc_cutoff = 1, p_cutoff = 0.05) {
  df <- tt
  
  if (!"SYMBOL" %in% colnames(df)) {
    df$SYMBOL <- rownames(df)
  }
  
  df <- df %>%
    mutate(
      negLog10P = -log10(P.Value),
      status = case_when(
        adj.P.Val < p_cutoff & logFC >= lfc_cutoff ~ "Up",
        adj.P.Val < p_cutoff & logFC <= -lfc_cutoff ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  p <- ggplot(df, aes(x = logFC, y = negLog10P, color = status)) +
    geom_point(alpha = 0.7, size = 1.6) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
    scale_color_manual(values = c("Down" = "#2C7BB6", "NS" = "grey70", "Up" = "#D7191C")) +
    labs(
      title = title_txt,
      x = "log2 fold change",
      y = "-log10(p-value)"
    ) +
    theme_bw(base_size = 12)
  
  ggplot2::ggsave(outfile, p, width = 8, height = 6, dpi = 180)
}

plot_heatmap_base <- function(mat, pheno, genes, outfile, title_txt) {
  genes <- unique(genes)
  genes <- genes[genes %in% rownames(mat)]
  
  if (length(genes) < 2) {
    message("Skipping heatmap for ", title_txt, " because fewer than 2 genes are available.")
    return(invisible(NULL))
  }
  
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    install.packages("pheatmap")
  }
  
  submat <- mat[genes, , drop = FALSE]
  submat <- t(scale(t(submat)))
  submat[is.na(submat)] <- 0
  
  ord <- order(pheno$GROUP)
  submat <- submat[, ord, drop = FALSE]
  pheno <- pheno[ord, , drop = FALSE]
  
  ann_col <- data.frame(GROUP = pheno$GROUP)
  rownames(ann_col) <- rownames(pheno)
  
  pheatmap::pheatmap(
    mat = submat,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    annotation_col = ann_col,
    main = title_txt,
    filename = outfile,
    width = 8,
    height = 10
  )
}


run_one_contrast <- function(fit2, coef_name, expr_mat, pheno, out_prefix,
                             p_cutoff = 0.05, lfc_cutoff = 1) {
  message("Running contrast: ", coef_name)
  
  tt <- limma::topTable(
    fit2,
    coef = coef_name,
    number = Inf,
    sort.by = "P"
  )
  
  if (!"SYMBOL" %in% colnames(tt)) {
    tt$SYMBOL <- rownames(tt)
  }
  
  tt <- tt %>%
    dplyr::select(
      SYMBOL, logFC, AveExpr, t, P.Value, adj.P.Val, B,
      everything()
    )
  
  out_csv_all <- here(resultsdir, paste0("DEG_", out_prefix, "_all.csv"))
  out_csv_sig <- here(resultsdir, paste0("DEG_", out_prefix, "_sig_FDR05.csv"))
  out_rds_all <- here(resultsdir, paste0("DEG_", out_prefix, "_all.rds"))
  
  saveRDS(tt, out_rds_all)
  rio::export(tt, out_csv_all)
  
  tt_sig <- tt %>%
    dplyr::filter(adj.P.Val < p_cutoff)
  
  rio::export(tt_sig, out_csv_sig)
  
  message(out_prefix, ": significant genes at FDR < ", p_cutoff, " = ", nrow(tt_sig))
  
  plot_volcano(
    tt = tt,
    outfile = here(plotsdir, paste0("volcano_", out_prefix, ".png")),
    title_txt = paste0("Differential expression: ", coef_name),
    lfc_cutoff = lfc_cutoff,
    p_cutoff = p_cutoff
  )
  
  top_heatmap_genes <- tt %>%
    dplyr::filter(adj.P.Val < p_cutoff) %>%
    dplyr::arrange(adj.P.Val) %>%
    dplyr::slice_head(n = 50) %>%
    dplyr::pull(SYMBOL)
  
  plot_heatmap_base(
    mat = expr_mat,
    pheno = pheno,
    genes = top_heatmap_genes,
    outfile = here(plotsdir, paste0("heatmap_top50_", out_prefix, ".png")),
    title_txt = paste0("Top DE genes: ", coef_name)
  )
  
  invisible(tt)
}

# ---- load data ----
message("Loading merged expression and phenotype data...")

expr_merged_combat <- readRDS(here(resultsdir, "expr_merged_combat.rds"))
pheno_merged <- readRDS(here(resultsdir, "pheno_merged.rds"))

# optional, only if present
feature_file <- here(resultsdir, "featureData_merged.rds")
featureData_merged <- NULL
if (file.exists(feature_file)) {
  featureData_merged <- readRDS(feature_file)
}

# ---- sanity checks ----
expr_merged_combat <- as.matrix(expr_merged_combat)
mode(expr_merged_combat) <- "numeric"

stopifnot(identical(colnames(expr_merged_combat), rownames(pheno_merged)))

if (!"GROUP" %in% colnames(pheno_merged)) {
  stop("GROUP column not found in pheno_merged.")
}

pheno_merged$GROUP <- clean_group_factor(pheno_merged$GROUP)

keep <- !is.na(pheno_merged$GROUP)
if (sum(!keep) > 0) {
  message("Removing ", sum(!keep), " samples with missing GROUP.")
}

expr_merged_combat <- expr_merged_combat[, keep, drop = FALSE]
pheno_merged <- pheno_merged[keep, , drop = FALSE]

stopifnot(identical(colnames(expr_merged_combat), rownames(pheno_merged)))

message("Samples per group:")
print(table(pheno_merged$GROUP, useNA = "ifany"))

if (any(table(pheno_merged$GROUP) < 2)) {
  stop("At least one group has fewer than 2 samples. Cannot run reliable DE.")
}

# ---- design matrix ----
design <- model.matrix(~ 0 + GROUP, data = pheno_merged)
colnames(design) <- sub("^GROUP", "", colnames(design))

message("Design matrix columns:")
print(colnames(design))

print(table(pheno_merged$GROUP))

# ---- limma fit ----
fit <- limma::lmFit(expr_merged_combat, design)

contrast_matrix <- limma::makeContrasts(
  AD_vs_CONTROL = AD - CONTROL,
  MCI_vs_CONTROL = MCI - CONTROL,
  AD_vs_MCI = AD - MCI,
  levels = design
)

fit2 <- limma::contrasts.fit(fit, contrast_matrix)
fit2 <- limma::eBayes(fit2)

# ---- save fit objects ----
saveRDS(design, here(resultsdir, "DE_design_matrix.rds"))
saveRDS(contrast_matrix, here(resultsdir, "DE_contrast_matrix.rds"))
saveRDS(fit2, here(resultsdir, "DE_limma_fit.rds"))

# ---- run contrasts ----
res_ad_vs_control <- run_one_contrast(
  fit2 = fit2,
  coef_name = "AD_vs_CONTROL",
  expr_mat = expr_merged_combat,
  pheno = pheno_merged,
  out_prefix = "AD_vs_CONTROL"
)

res_mci_vs_control <- run_one_contrast(
  fit2 = fit2,
  coef_name = "MCI_vs_CONTROL",
  expr_mat = expr_merged_combat,
  pheno = pheno_merged,
  out_prefix = "MCI_vs_CONTROL"
)

res_ad_vs_mci <- run_one_contrast(
  fit2 = fit2,
  coef_name = "AD_vs_MCI",
  expr_mat = expr_merged_combat,
  pheno = pheno_merged,
  out_prefix = "AD_vs_MCI"
)

# ---- summary table ----
make_summary_row <- function(tt, contrast_name, p_cutoff = 0.05, lfc_cutoff = 1) {
  data.frame(
    contrast = contrast_name,
    n_total = nrow(tt),
    n_fdr_0_05 = sum(tt$adj.P.Val < p_cutoff, na.rm = TRUE),
    n_fdr_0_05_abs_logFC_1 = sum(tt$adj.P.Val < p_cutoff & abs(tt$logFC) >= lfc_cutoff, na.rm = TRUE),
    top_gene = if (nrow(tt) > 0) tt$SYMBOL[1] else NA_character_,
    top_gene_adjP = if (nrow(tt) > 0) tt$adj.P.Val[1] else NA_real_,
    stringsAsFactors = FALSE
  )
}

summary_tbl <- dplyr::bind_rows(
  make_summary_row(res_ad_vs_control, "AD_vs_CONTROL"),
  make_summary_row(res_mci_vs_control, "MCI_vs_CONTROL"),
  make_summary_row(res_ad_vs_mci, "AD_vs_MCI")
)

print(summary_tbl)

saveRDS(summary_tbl, here(resultsdir, "DE_summary.rds"))
rio::export(summary_tbl, here(resultsdir, "DE_summary.csv"))

# ---- optional: ranked gene lists for GSEA ----
write_ranked_list <- function(tt, outfile) {
  ranks <- tt[, c("SYMBOL", "t")]
  ranks <- ranks[!is.na(ranks$SYMBOL) & !duplicated(ranks$SYMBOL), , drop = FALSE]
  ranks <- ranks[order(ranks$t, decreasing = TRUE), , drop = FALSE]
  rio::export(ranks, outfile)
}

write_ranked_list(res_ad_vs_control, here(resultsdir, "ranked_AD_vs_CONTROL_tstat.csv"))
write_ranked_list(res_mci_vs_control, here(resultsdir, "ranked_MCI_vs_CONTROL_tstat.csv"))
write_ranked_list(res_ad_vs_mci, here(resultsdir, "ranked_AD_vs_MCI_tstat.csv"))

message("Finished 04_differential_expression.R")
message("Outputs saved in results/ and plots/")
