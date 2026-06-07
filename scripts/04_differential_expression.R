# scripts/04_differential_expression.R

message("Starting 04_differential_expression.R")

# ---- libraries ----
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

required_pkgs <- c("limma", "dplyr", "tibble", "rio", "ggplot2", "ggrepel")
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

plot_volcano <- function(tt, outfile, title_txt, p_cutoff = 0.05,
                         xlim = NULL, ylim = NULL,
                         label_genes = NULL, label_n = 2) {
  df <- tt
  
  if (!"SYMBOL" %in% colnames(df)) {
    df$SYMBOL <- rownames(df)
  }
  
  # Keine logFC-Schranke mehr: Signifikanz allein bestimmt die Markierung,
  # die Farbe ergibt sich aus der Richtung der Regulation.
  df <- df %>%
    mutate(
      negLog10P = -log10(P.Value),
      status = case_when(
        adj.P.Val < p_cutoff & logFC > 0 ~ "Up",
        adj.P.Val < p_cutoff & logFC < 0 ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  # NS zuerst zeichnen, damit die farbigen (signifikanten) Punkte oben liegen
  df$status <- factor(df$status, levels = c("NS", "Down", "Up"))
  df <- df[order(df$status), , drop = FALSE]
  
  # Schwellenlinie auf der rohen-p-Achse setzen, wo die FDR-Grenze tatsaechlich
  # liegt (groesster roher p-Wert unter den FDR-signifikanten Genen). Da die
  # BH-Korrektur monoton im Rang ist, ist alles oberhalb dieser Linie auch
  # FDR-signifikant und damit gefaerbt -- Linie und Einfaerbung passen zusammen.
  sig <- df$adj.P.Val < p_cutoff
  y_line <- if (any(sig, na.rm = TRUE)) {
    -log10(max(df$P.Value[sig], na.rm = TRUE))
  } else {
    -log10(p_cutoff)
  }
  
  # Zu beschriftende Gene: entweder explizit per label_genes vorgegeben,
  # sonst die label_n signifikantesten (kleinstes adj.P.Val) dieses Kontrasts.
  if (!is.null(label_genes)) {
    lab_df <- df[df$SYMBOL %in% label_genes, , drop = FALSE]
  } else {
    lab_df <- df[order(df$adj.P.Val), , drop = FALSE]
    lab_df <- head(lab_df, label_n)
  }
  
  p <- ggplot(df, aes(x = logFC, y = negLog10P, color = status)) +
    geom_point(alpha = 0.7, size = 1.6) +
    geom_hline(yintercept = y_line, linetype = "dashed") +
    scale_color_manual(
      values = c("Down" = "#2C7BB6", "NS" = "grey70", "Up" = "#D7191C"),
      breaks = c("Up", "Down", "NS"),
      labels = c("Up" = "hochreguliert", "Down" = "runterreguliert", "NS" = "n.s.")
    ) +
    labs(
      title = title_txt,
      x = "log2 fold change",
      y = "-log10(p-value)",
      color = "status",
      caption = paste0("gestrichelte Linie: FDR < ", p_cutoff)
    ) +
    theme_bw(base_size = 12)
  
  # Gen-Labels (mit Verbindungslinien, schwarz, ohne Legendeneintrag)
  if (nrow(lab_df) > 0) {
    p <- p +
      ggrepel::geom_text_repel(
        data = lab_df,
        aes(x = logFC, y = negLog10P, label = SYMBOL),
        inherit.aes = FALSE,
        color = "black",
        size = 3.5,
        fontface = "bold",
        box.padding = 0.6,
        point.padding = 0.3,
        min.segment.length = 0,
        max.overlaps = Inf,
        seed = 42
      )
  }
  
  # Synchrone Achsen ueber alle Kontraste hinweg (coord_cartesian schneidet
  # nur die Ansicht zu, verwirft aber keine Punkte).
  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  }
  
  ggplot2::ggsave(outfile, p, width = 8, height = 6, dpi = 180)
}

plot_heatmap_base <- function(mat, pheno, genes, outfile, title_txt, z_cap = 2.5) {
  genes <- unique(genes)
  genes <- genes[genes %in% rownames(mat)]
  
  if (length(genes) < 2) {
    message("Skipping heatmap for ", title_txt, " because fewer than 2 genes are available.")
    return(invisible(NULL))
  }
  
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    install.packages("pheatmap")
  }
  
  # z-Transformation pro Gen (zeilenweise) ueber die Proben.
  submat <- mat[genes, , drop = FALSE]
  submat <- t(scale(t(submat)))
  submat[is.na(submat)] <- 0
  
  # Extreme z-Werte symmetrisch kappen, damit einzelne Ausreisser die
  # Farbskala nicht auseinanderziehen und der Grossteil nicht im neutralen
  # Bereich verschwindet.
  submat[submat >  z_cap] <-  z_cap
  submat[submat < -z_cap] <- -z_cap
  
  ord <- order(pheno$GROUP)
  submat <- submat[, ord, drop = FALSE]
  pheno <- pheno[ord, , drop = FALSE]
  
  ann_col <- data.frame(GROUP = pheno$GROUP)
  rownames(ann_col) <- rownames(pheno)
  
  # Divergierende Skala, symmetrisch um 0 zentriert: Weiss = neutral
  # (statt Gelb), Blau = runterreguliert, Rot = hochreguliert.
  heat_breaks <- seq(-z_cap, z_cap, length.out = 101)
  heat_colors <- colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(length(heat_breaks) - 1)
  
  pheatmap::pheatmap(
    mat = submat,
    color = heat_colors,
    breaks = heat_breaks,
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
                             p_cutoff = 0.05, lfc_cutoff = 1,
                             volcano_xlim = NULL, volcano_ylim = NULL,
                             volcano_label_genes = NULL, volcano_label_n = 2) {
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
    p_cutoff = p_cutoff,
    xlim = volcano_xlim,
    ylim = volcano_ylim,
    label_genes = volcano_label_genes,
    label_n = volcano_label_n
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

# Gemeinsame Volcano-Achsen ueber alle Kontraste, damit die Signalstaerke
# direkt vergleichbar ist (sonst hat jeder Plot seine eigene y-Skala).
volcano_coefs <- c("AD_vs_CONTROL", "MCI_vs_CONTROL", "AD_vs_MCI")

global_logFC <- numeric(0)
global_negLog10P <- numeric(0)
for (cf in volcano_coefs) {
  tt_tmp <- limma::topTable(fit2, coef = cf, number = Inf, sort.by = "P")
  global_logFC <- c(global_logFC, tt_tmp$logFC)
  global_negLog10P <- c(global_negLog10P, -log10(tt_tmp$P.Value))
}

x_abs <- max(abs(global_logFC), na.rm = TRUE)
volcano_xlim <- c(-x_abs, x_abs)                       # symmetrisch um 0
volcano_ylim <- c(0, max(global_negLog10P, na.rm = TRUE))

message("Shared volcano x-limits: ", paste(round(volcano_xlim, 3), collapse = " ... "))
message("Shared volcano y-limits: ", paste(round(volcano_ylim, 3), collapse = " ... "))

# Feste Gene, die in ALLEN Volcano-Plots beschriftet werden, damit man
# dieselben Gene ueber die Kontraste hinweg wiederfindet und ihre
# Verschiebung vergleichen kann. Ein Gen pro Pathview-Karte, eindeutig
# zugeordnet:
#   NDUFA1 = Oxidative phosphorylation (hsa00190), Atmungskette Komplex I
#   PSMA3  = Alzheimer disease (hsa05010), Proteasom-Arm -- liegt NICHT in
#            hsa00190, daher eindeutig der Alzheimer-Karte zuzuordnen.
volcano_label_genes <- c("NDUFA1", "PSMA3")

res_ad_vs_control <- run_one_contrast(
  fit2 = fit2,
  coef_name = "AD_vs_CONTROL",
  expr_mat = expr_merged_combat,
  pheno = pheno_merged,
  out_prefix = "AD_vs_CONTROL",
  volcano_xlim = volcano_xlim,
  volcano_ylim = volcano_ylim,
  volcano_label_genes = volcano_label_genes
)

res_mci_vs_control <- run_one_contrast(
  fit2 = fit2,
  coef_name = "MCI_vs_CONTROL",
  expr_mat = expr_merged_combat,
  pheno = pheno_merged,
  out_prefix = "MCI_vs_CONTROL",
  volcano_xlim = volcano_xlim,
  volcano_ylim = volcano_ylim,
  volcano_label_genes = volcano_label_genes
)

res_ad_vs_mci <- run_one_contrast(
  fit2 = fit2,
  coef_name = "AD_vs_MCI",
  expr_mat = expr_merged_combat,
  pheno = pheno_merged,
  out_prefix = "AD_vs_MCI",
  volcano_xlim = volcano_xlim,
  volcano_ylim = volcano_ylim,
  volcano_label_genes = volcano_label_genes
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