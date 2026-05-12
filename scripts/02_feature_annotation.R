# scripts/02_feature_annotation.R

message("Starting 02_feature_annotation.R")

# ---- libraries ----
library(here)
library(GEOquery)
library(dplyr)
library(stringr)
library(tibble)
library(rio)

# ---- directories ----
datadir <- here("data")
resultsdir <- here("results")

# ---- load assay data ----
assayData_63060 <- readRDS(here(resultsdir, "assayData_63060_raw.rds"))
assayData_63061 <- readRDS(here(resultsdir, "assayData_63061_raw.rds"))

# ---- download/load GPL annotation ----
message("Downloading GPL annotation tables...")

gpl6947 <- getGEO("GPL6947", AnnotGPL = TRUE)
gpl10558 <- getGEO("GPL10558", AnnotGPL = TRUE)

table6947 <- Table(gpl6947)
table10558 <- Table(gpl10558)

# ---- helper: clean annotation ----
clean_annotation <- function(tbl) {
  
  colnames(tbl) <- make.names(colnames(tbl))
  
  # helper to safely pick a column
  pick_col <- function(patterns, cols) {
    for (p in patterns) {
      hit <- grep(p, cols, ignore.case = TRUE, value = TRUE)
      if (length(hit) > 0) return(hit[1])
    }
    return(NA_character_)
  }
  
  cols <- colnames(tbl)
  
  id_col <- pick_col(c("^ID$", "ID_REF", "Probe"), cols)
  symbol_col <- pick_col(c("Gene.Symbol", "Symbol", "GENE_SYMBOL"), cols)
  entrez_col <- pick_col(c("Entrez", "ENTREZ"), cols)
  gene_col <- pick_col(c("Gene.Title", "Gene.Name", "GENE_NAME"), cols)
  
  message("Detected columns:")
  message("ID: ", id_col)
  message("SYMBOL: ", symbol_col)
  message("ENTREZ: ", entrez_col)
  message("GENENAME: ", gene_col)
  
  # build dataframe safely
  out <- data.frame(
    FEATUREID = tbl[[id_col]],
    SYMBOL = if (!is.na(symbol_col)) tbl[[symbol_col]] else NA,
    ENTREZ = if (!is.na(entrez_col)) tbl[[entrez_col]] else NA,
    GENENAME = if (!is.na(gene_col)) tbl[[gene_col]] else NA,
    stringsAsFactors = FALSE
  )
  
  # clean multi-mapping
  out <- out %>%
    mutate(
      SYMBOL = ifelse(is.na(SYMBOL), NA,
                      str_split(SYMBOL, " // ") %>% sapply(`[`, 1)),
      ENTREZ = ifelse(is.na(ENTREZ), NA,
                      str_split(ENTREZ, " // ") %>% sapply(`[`, 1))
    )
  
  return(out)
}

annot_6947 <- clean_annotation(table6947)
annot_10558 <- clean_annotation(table10558)

# ---- match annotation to assay ----
featureData_63060 <- annot_6947 %>%
  filter(FEATUREID %in% rownames(assayData_63060)) %>%
  distinct(FEATUREID, .keep_all = TRUE) %>%
  tibble::column_to_rownames("FEATUREID")

featureData_63061 <- annot_10558 %>%
  filter(FEATUREID %in% rownames(assayData_63061)) %>%
  distinct(FEATUREID, .keep_all = TRUE) %>%
  tibble::column_to_rownames("FEATUREID")

# ---- align order ----
featureData_63060 <- featureData_63060[rownames(assayData_63060), , drop = FALSE]
featureData_63061 <- featureData_63061[rownames(assayData_63061), , drop = FALSE]

# ---- QC ----
message("Checking annotation coverage...")

message("GSE63060 annotated probes: ",
        sum(!is.na(featureData_63060$SYMBOL)), "/",
        nrow(featureData_63060))

message("GSE63061 annotated probes: ",
        sum(!is.na(featureData_63061$SYMBOL)), "/",
        nrow(featureData_63061))

# ---- optional: filter unannotated probes ----
keep_63060 <- !is.na(featureData_63060$SYMBOL) & featureData_63060$SYMBOL != ""
keep_63061 <- !is.na(featureData_63061$SYMBOL) & featureData_63061$SYMBOL != ""

assayData_63060 <- assayData_63060[keep_63060, ]
assayData_63061 <- assayData_63061[keep_63061, ]

featureData_63060 <- featureData_63060[keep_63060, ]
featureData_63061 <- featureData_63061[keep_63061, ]

# ---- save outputs ----
saveRDS(featureData_63060, here(resultsdir, "featureData_63060.rds"))
saveRDS(featureData_63061, here(resultsdir, "featureData_63061.rds"))

saveRDS(assayData_63060, here(resultsdir, "assayData_63060_annotated.rds"))
saveRDS(assayData_63061, here(resultsdir, "assayData_63061_annotated.rds"))

rio::export(
  tibble::rownames_to_column(featureData_63060, "FEATUREID"),
  here(resultsdir, "featureData_63060.csv")
)

rio::export(
  tibble::rownames_to_column(featureData_63061, "FEATUREID"),
  here(resultsdir, "featureData_63061.csv")
)

message("Finished 02_feature_annotation.R")
message("Next step: normalization + batch correction (03_normalization.R)")