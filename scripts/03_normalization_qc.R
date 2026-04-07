message("Starting 03_normalization.R")

# ---- libraries ----
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

required_pkgs <- c("dplyr", "tibble", "rio", "limma", "sva", "ggplot2", "stringr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("limma", "sva")) {
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
datadir <- here("data")
resultsdir <- here("results")
plotsdir <- here("plots")

dir.create(resultsdir, showWarnings = FALSE, recursive = TRUE)
dir.create(plotsdir, showWarnings = FALSE, recursive = TRUE)

# ---- helper functions ----

is_log_scale <- function(mat) {
  x <- as.numeric(mat)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(TRUE)
  
  qs <- stats::quantile(x, probs = c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)
  
  # typical heuristic used for microarray expression objects
  loglike <- (qs[5] < 100 && qs[6] < 100)
  return(loglike)
}

safe_log2_transform <- function(mat) {
  if (is_log_scale(mat)) {
    message("Matrix already appears to be on log scale. Skipping log2 transform.")
    return(mat)
  }
  
  message("Matrix appears NOT to be on log scale. Applying log2(x + 1).")
  mat[mat < 0] <- NA
  mat <- log2(mat + 1)
  return(mat)
}

plot_boxplot <- function(mat, outfile, title_txt) {
  grDevices::png(outfile, width = 2200, height = 1200, res = 180)
  boxplot(
    mat,
    las = 2,
    outline = FALSE,
    cex.axis = 0.45,
    main = title_txt,
    ylab = "Expression"
  )
  grDevices::dev.off()
}

plot_density <- function(mat, outfile, title_txt) {
  grDevices::png(outfile, width = 1800, height = 1200, res = 180)
  dens_list <- apply(mat, 2, function(x) density(x[is.finite(x)], na.rm = TRUE))
  xlim <- range(vapply(dens_list, function(d) range(d$x), numeric(2)), na.rm = TRUE)
  ylim <- range(vapply(dens_list, function(d) range(d$y), numeric(2)), na.rm = TRUE)
  
  plot(
    dens_list[[1]],
    main = title_txt,
    xlab = "Expression",
    ylab = "Density",
    xlim = xlim,
    ylim = ylim,
    lwd = 1
  )
  if (length(dens_list) > 1) {
    for (i in 2:length(dens_list)) {
      lines(dens_list[[i]], lwd = 1)
    }
  }
  grDevices::dev.off()
}

make_pca_df <- function(mat, pheno, color_var = "BATCH") {
  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
  
  df <- data.frame(
    SAMPLEID = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )
  
  pheno2 <- pheno %>%
    tibble::rownames_to_column("SAMPLEID")
  
  df <- df %>%
    dplyr::left_join(pheno2, by = "SAMPLEID")
  
  attr(df, "pc1_var") <- round(100 * var_expl[1], 2)
  attr(df, "pc2_var") <- round(100 * var_expl[2], 2)
  
  return(df)
}

plot_pca <- function(mat, pheno, outfile, title_txt, color_var = "BATCH", shape_var = "GROUP") {
  df <- make_pca_df(mat, pheno, color_var = color_var)
  pc1_var <- attr(df, "pc1_var")
  pc2_var <- attr(df, "pc2_var")
  
  keep_cols <- c("PC1", "PC2", color_var, shape_var)
  df <- df[complete.cases(df[, keep_cols, drop = FALSE]), , drop = FALSE]
  
  p <- ggplot(
    df,
    aes(x = PC1, y = PC2, color = .data[[color_var]], shape = .data[[shape_var]])
  ) +
    geom_point(size = 3, alpha = 0.9) +
    labs(
      title = title_txt,
      x = paste0("PC1 (", pc1_var, "%)"),
      y = paste0("PC2 (", pc2_var, "%)")
    ) +
    theme_bw(base_size = 12)
  
  ggplot2::ggsave(outfile, p, width = 9, height = 6, dpi = 180)
}

collapse_to_symbol <- function(expr_mat, feature_df, dataset_name = "dataset") {
  stopifnot(nrow(expr_mat) == nrow(feature_df))
  stopifnot(identical(rownames(expr_mat), rownames(feature_df)))
  
  tmp <- data.frame(
    FEATUREID = rownames(feature_df),
    SYMBOL = feature_df$SYMBOL,
    stringsAsFactors = FALSE
  )
  
  keep <- !is.na(tmp$SYMBOL) & trimws(tmp$SYMBOL) != ""
  tmp <- tmp[keep, , drop = FALSE]
  expr_mat <- expr_mat[keep, , drop = FALSE]
  
  message(dataset_name, ": probes with non-missing gene symbol = ", nrow(expr_mat))
  
  # limma::avereps averages probes that map to the same ID
  collapsed <- limma::avereps(expr_mat, ID = tmp$SYMBOL)
  
  message(dataset_name, ": unique gene symbols after collapsing = ", nrow(collapsed))
  return(collapsed)
}

# ---- load data ----
message("Loading annotated assay, feature, and pheno data...")

assayData_63060 <- readRDS(here("results", "assayData_63060_annotated.rds"))
assayData_63061 <- readRDS(here("results", "assayData_63061_annotated.rds"))

featureData_63060 <- readRDS(here("results", "featureData_63060.rds"))
featureData_63061 <- readRDS(here("results", "featureData_63061.rds"))

phenoData_63060 <- readRDS(here("results", "phenoData_63060_raw.rds"))
phenoData_63061 <- readRDS(here("results", "phenoData_63061_raw.rds"))

# ---- sanity checks ----
stopifnot(identical(rownames(assayData_63060), rownames(featureData_63060)))
stopifnot(identical(rownames(assayData_63061), rownames(featureData_63061)))

stopifnot(identical(colnames(assayData_63060), rownames(phenoData_63060)))
stopifnot(identical(colnames(assayData_63061), rownames(phenoData_63061)))

message("Dimensions before normalization:")
message("GSE63060 assay: ", paste(dim(assayData_63060), collapse = " x "))
message("GSE63061 assay: ", paste(dim(assayData_63061), collapse = " x "))

# ---- QC before normalization ----
plot_boxplot(
  assayData_63060,
  here("plots", "boxplot_GSE63060_before_normalization.png"),
  "GSE63060 before normalization"
)

plot_boxplot(
  assayData_63061,
  here("plots", "boxplot_GSE63061_before_normalization.png"),
  "GSE63061 before normalization"
)

plot_density(
  assayData_63060,
  here("plots", "density_GSE63060_before_normalization.png"),
  "GSE63060 before normalization"
)

plot_density(
  assayData_63061,
  here("plots", "density_GSE63061_before_normalization.png"),
  "GSE63061 before normalization"
)

# ---- log transform if needed ----
assayData_63060 <- safe_log2_transform(assayData_63060)
assayData_63061 <- safe_log2_transform(assayData_63061)

# ---- within-dataset normalization ----
message("Applying quantile normalization within each dataset...")

assayNorm_63060 <- limma::normalizeBetweenArrays(assayData_63060, method = "quantile")
assayNorm_63061 <- limma::normalizeBetweenArrays(assayData_63061, method = "quantile")

# ---- QC after within-dataset normalization ----
plot_boxplot(
  assayNorm_63060,
  here("plots", "boxplot_GSE63060_after_quantile_normalization.png"),
  "GSE63060 after quantile normalization"
)

plot_boxplot(
  assayNorm_63061,
  here("plots", "boxplot_GSE63061_after_quantile_normalization.png"),
  "GSE63061 after quantile normalization"
)

plot_density(
  assayNorm_63060,
  here("plots", "density_GSE63060_after_quantile_normalization.png"),
  "GSE63060 after quantile normalization"
)

plot_density(
  assayNorm_63061,
  here("plots", "density_GSE63061_after_quantile_normalization.png"),
  "GSE63061 after quantile normalization"
)

# ---- collapse probes to gene symbols ----
message("Collapsing probes to gene symbols...")

exprGene_63060 <- collapse_to_symbol(
  expr_mat = assayNorm_63060,
  feature_df = featureData_63060,
  dataset_name = "GSE63060"
)

exprGene_63061 <- collapse_to_symbol(
  expr_mat = assayNorm_63061,
  feature_df = featureData_63061,
  dataset_name = "GSE63061"
)

# ---- keep common genes ----
common_genes <- intersect(rownames(exprGene_63060), rownames(exprGene_63061))
message("Common genes between GSE63060 and GSE63061: ", length(common_genes))

if (length(common_genes) < 1000) {
  warning("Low number of common genes detected: ", length(common_genes))
}

exprGene_63060 <- exprGene_63060[common_genes, , drop = FALSE]
exprGene_63061 <- exprGene_63061[common_genes, , drop = FALSE]

# ---- merge datasets ----
message("Merging datasets by common gene symbols...")

expr_merged_precombat <- cbind(exprGene_63060, exprGene_63061)

pheno_merged <- dplyr::bind_rows(
  phenoData_63060 %>% tibble::rownames_to_column("SAMPLEID"),
  phenoData_63061 %>% tibble::rownames_to_column("SAMPLEID")
) %>%
  tibble::column_to_rownames("SAMPLEID")

pheno_merged <- pheno_merged[colnames(expr_merged_precombat), , drop = FALSE]

stopifnot(identical(colnames(expr_merged_precombat), rownames(pheno_merged)))

# ---- create merged feature data ----
featureData_merged <- data.frame(
  SYMBOL = rownames(expr_merged_precombat),
  stringsAsFactors = FALSE
)
rownames(featureData_merged) <- featureData_merged$SYMBOL

# ---- PCA before batch correction ----
plot_pca(
  mat = expr_merged_precombat,
  pheno = pheno_merged,
  outfile = here("plots", "PCA_merged_before_ComBat_color_BATCH.png"),
  title_txt = "Merged data before ComBat",
  color_var = "BATCH",
  shape_var = "GROUP"
)

plot_pca(
  mat = expr_merged_precombat,
  pheno = pheno_merged,
  outfile = here("plots", "PCA_merged_before_ComBat_color_GROUP.png"),
  title_txt = "Merged data before ComBat",
  color_var = "GROUP",
  shape_var = "BATCH"
)

# ---- batch correction with ComBat ----
message("Running ComBat batch correction...")

batch <- pheno_merged$BATCH
if (any(is.na(batch))) {
  stop("BATCH contains NA values in merged phenotype data.")
}

if ("GROUP" %in% colnames(pheno_merged)) {
  keep_combat <- !is.na(pheno_merged$GROUP)
  
  if (sum(!keep_combat) > 0) {
    message("Excluding ", sum(!keep_combat), " samples from ComBat due to missing GROUP.")
  }
  
  expr_merged_precombat_combat <- expr_merged_precombat[, keep_combat, drop = FALSE]
  pheno_merged_combat <- pheno_merged[keep_combat, , drop = FALSE]
  batch_combat <- pheno_merged_combat$BATCH
  mod <- model.matrix(~ GROUP, data = pheno_merged_combat)
  
  stopifnot(
    ncol(expr_merged_precombat_combat) == nrow(pheno_merged_combat),
    length(batch_combat) == nrow(pheno_merged_combat),
    nrow(mod) == nrow(pheno_merged_combat)
  )
  
  expr_merged_combat <- sva::ComBat(
    dat = as.matrix(expr_merged_precombat_combat),
    batch = batch_combat,
    mod = mod,
    par.prior = TRUE,
    prior.plots = FALSE
  )
  
  pheno_merged_final <- pheno_merged_combat
} else {
  expr_merged_combat <- sva::ComBat(
    dat = as.matrix(expr_merged_precombat),
    batch = batch,
    mod = NULL,
    par.prior = TRUE,
    prior.plots = FALSE
  )
  
  pheno_merged_final <- pheno_merged
}

# ---- PCA after batch correction ----
plot_pca(
  mat = expr_merged_combat,
  pheno = pheno_merged_final,
  outfile = here("plots", "PCA_merged_after_ComBat_color_BATCH.png"),
  title_txt = "Merged data after ComBat",
  color_var = "BATCH",
  shape_var = "GROUP"
)

plot_pca(
  mat = expr_merged_combat,
  pheno = pheno_merged_final,
  outfile = here("plots", "PCA_merged_after_ComBat_color_GROUP.png"),
  title_txt = "Merged data after ComBat",
  color_var = "GROUP",
  shape_var = "BATCH"
)

plot_boxplot(
  expr_merged_precombat,
  here("plots", "boxplot_merged_before_ComBat.png"),
  "Merged data before ComBat"
)

plot_boxplot(
  expr_merged_combat,
  here("plots", "boxplot_merged_after_ComBat.png"),
  "Merged data after ComBat"
)

plot_density(
  expr_merged_precombat,
  here("plots", "density_merged_before_ComBat.png"),
  "Merged data before ComBat"
)

plot_density(
  expr_merged_combat,
  here("plots", "density_merged_after_ComBat.png"),
  "Merged data after ComBat"
)

# ---- save outputs ----
message("Saving normalized and batch-corrected outputs...")

saveRDS(assayNorm_63060, here("results", "assayData_63060_normalized.rds"))
saveRDS(assayNorm_63061, here("results", "assayData_63061_normalized.rds"))

saveRDS(exprGene_63060, here("results", "exprGene_63060_normalized.rds"))
saveRDS(exprGene_63061, here("results", "exprGene_63061_normalized.rds"))

saveRDS(expr_merged_precombat, here("results", "expr_merged_precombat.rds"))
saveRDS(expr_merged_combat, here("results", "expr_merged_combat.rds"))

saveRDS(pheno_merged_final, here("results", "pheno_merged.rds"))
saveRDS(featureData_merged, here("results", "featureData_merged.rds"))

rio::export(
  tibble::rownames_to_column(as.data.frame(exprGene_63060), "SYMBOL"),
  here("results", "exprGene_63060_normalized.csv")
)

rio::export(
  tibble::rownames_to_column(as.data.frame(exprGene_63061), "SYMBOL"),
  here("results", "exprGene_63061_normalized.csv")
)

rio::export(
  tibble::rownames_to_column(as.data.frame(expr_merged_precombat), "SYMBOL"),
  here("results", "expr_merged_precombat.csv")
)

rio::export(
  tibble::rownames_to_column(as.data.frame(expr_merged_combat), "SYMBOL"),
  here("results", "expr_merged_combat.csv")
)

rio::export(
  tibble::rownames_to_column(pheno_merged_final, "SAMPLEID"),
  here("results", "pheno_merged.csv")
)

rio::export(
  featureData_merged,
  here("results", "featureData_merged.csv")
)

# ---- summary table ----
summary_tbl <- data.frame(
  object = c(
    "assayData_63060_normalized",
    "assayData_63061_normalized",
    "exprGene_63060_normalized",
    "exprGene_63061_normalized",
    "expr_merged_precombat",
    "expr_merged_combat"
  ),
  n_features = c(
    nrow(assayNorm_63060),
    nrow(assayNorm_63061),
    nrow(exprGene_63060),
    nrow(exprGene_63061),
    nrow(expr_merged_precombat),
    nrow(expr_merged_combat)
  ),
  n_samples = c(
    ncol(assayNorm_63060),
    ncol(assayNorm_63061),
    ncol(exprGene_63060),
    ncol(exprGene_63061),
    ncol(expr_merged_precombat),
    ncol(expr_merged_combat)
  ),
  stringsAsFactors = FALSE
)

rio::export(summary_tbl, here("results", "normalization_summary.csv"))

message("Finished 03_normalization.R")
message("Main final object for downstream analysis: results/expr_merged_combat.rds")
message("Phenotype table: results/pheno_merged.rds")
message("Feature table: results/featureData_merged.rds")
