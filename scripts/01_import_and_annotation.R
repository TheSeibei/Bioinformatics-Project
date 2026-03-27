# scripts/01_import_and_annotation.R

message("Starting 01_import_and_annotation.R")

# ---- libraries ----
libraries <- unlist(read.table("libraries.txt", stringsAsFactors = FALSE))

for (pkg in libraries) {
  if (!require(pkg, character.only = TRUE)) {
    message(paste("Installing missing package:", pkg))
    
    if (pkg %in% c("GEOquery", "limma", "sva", "AnnotationDbi", "org.Hs.eg.db")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    
    library(pkg, character.only = TRUE)
  }
}

# ---- directories ----
dir.create(here("data"), showWarnings = FALSE)
dir.create(here("results"), showWarnings = FALSE)
dir.create(here("plots"), showWarnings = FALSE)
dir.create(here("scripts"), showWarnings = FALSE)

message("Project root detected by here(): ", here())

# ---- helper function to read GEO-like txt files ----
read_geo_table <- function(path) {
  # detect header row beginning with ID_REF if GEO metadata lines are present
  lines <- readLines(path, n = 200)
  id_ref_line <- grep("^ID_REF\\b", lines)
  
  if (length(id_ref_line) > 0) {
    df <- read.delim(
      path,
      skip = id_ref_line[1] - 1,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    df <- read.delim(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  
  df
}

# ---- input files ----
# rename these filenames if yours differ
file_63060 <- here("data", "GSE63060_non-normalized.txt")
file_63061 <- here("data", "GSE63061_non-normalized.txt")

stopifnot(file.exists(file_63060))
stopifnot(file.exists(file_63061))

# ---- import raw expression tables ----
raw_63060 <- read_geo_table(file_63060)
raw_63061 <- read_geo_table(file_63061)

message("Imported raw tables.")
message("GSE63060 dim: ", paste(dim(raw_63060), collapse = " x "))
message("GSE63061 dim: ", paste(dim(raw_63061), collapse = " x "))

# inspect first columns
print(colnames(raw_63060)[1:min(10, ncol(raw_63060))])
print(colnames(raw_63061)[1:min(10, ncol(raw_63061))])

# ---- build assayData ----
# assumes first column is probe ID, often ID_REF
make_assay_matrix <- function(df) {
  probe_col <- colnames(df)[1]
  
  assay <- df %>%
    rename(FEATUREID = all_of(probe_col)) %>%
    group_by(FEATUREID) %>%
    summarise(across(everything(), ~ suppressWarnings(max(as.numeric(.x), na.rm = TRUE)))) %>%
    ungroup() %>%
    column_to_rownames("FEATUREID") %>%
    as.matrix()
  
  mode(assay) <- "numeric"
  assay
}

assayData_63060 <- make_assay_matrix(raw_63060)
assayData_63061 <- make_assay_matrix(raw_63061)

message("Constructed assay matrices.")
message("assayData_63060 dim: ", paste(dim(assayData_63060), collapse = " x "))
message("assayData_63061 dim: ", paste(dim(assayData_63061), collapse = " x "))

# ---- get sample annotation from GEO ----
# this downloads/reads the series matrix metadata
gse63060 <- getGEO("GSE63060", GSEMatrix = TRUE)
gse63061 <- getGEO("GSE63061", GSEMatrix = TRUE)

pheno_raw_63060 <- pData(gse63060[[1]])
pheno_raw_63061 <- pData(gse63061[[1]])

# save raw phenotype tables for manual inspection
rio::export(
  tibble::rownames_to_column(pheno_raw_63060, "SAMPLEID"),
  here("results", "pheno_raw_GSE63060.csv")
)

rio::export(
  tibble::rownames_to_column(pheno_raw_63061, "SAMPLEID"),
  here("results", "pheno_raw_GSE63061.csv")
)

message("Saved raw phenotype tables to results/")

# ---- create preliminary phenoData ----
# IMPORTANT:
# inspect the exported raw phenotype csv files and adjust GROUP extraction if needed

extract_group <- function(x) {
  case_when(
    grepl("alzheimer", x, ignore.case = TRUE) ~ "AD",
    grepl("\\bmci\\b", x, ignore.case = TRUE) ~ "MCI",
    grepl("control", x, ignore.case = TRUE) ~ "CONTROL",
    grepl("\\bctl\\b", x, ignore.case = TRUE) ~ "CONTROL",
    grepl("normal", x, ignore.case = TRUE) ~ "CONTROL",
    TRUE ~ NA_character_
  )
}

collapse_characteristics <- function(df) {
  char_cols <- grep("^characteristics_ch1", colnames(df), value = TRUE)
  
  if (length(char_cols) == 0) {
    return(rep(NA_character_, nrow(df)))
  }
  
  apply(df[, char_cols, drop = FALSE], 1, function(z) paste(z, collapse = " | "))
}

phenoData_63060 <- pheno_raw_63060 %>%
  tibble::rownames_to_column("SAMPLEID") %>%
  mutate(
    CHARACTERISTICS_ALL = collapse_characteristics(pheno_raw_63060),
    GROUP = extract_group(CHARACTERISTICS_ALL),
    BATCH = "GSE63060",
    STUDY = "GSE63060"
  ) %>%
  select(SAMPLEID, GROUP, BATCH, STUDY, everything())

phenoData_63061 <- pheno_raw_63061 %>%
  tibble::rownames_to_column("SAMPLEID") %>%
  mutate(
    CHARACTERISTICS_ALL = collapse_characteristics(pheno_raw_63061),
    GROUP = extract_group(CHARACTERISTICS_ALL),
    BATCH = "GSE63061",
    STUDY = "GSE63061"
  ) %>%
  select(SAMPLEID, GROUP, BATCH, STUDY, everything())

# ---- align sample names between assayData and phenoData ----
common_63060 <- intersect(colnames(assayData_63060), phenoData_63060$SAMPLEID)
common_63061 <- intersect(colnames(assayData_63061), phenoData_63061$SAMPLEID)

assayData_63060 <- assayData_63060[, common_63060, drop = FALSE]
assayData_63061 <- assayData_63061[, common_63061, drop = FALSE]

phenoData_63060 <- phenoData_63060 %>%
  filter(SAMPLEID %in% common_63060) %>%
  mutate(rowname = SAMPLEID) %>%
  column_to_rownames("rowname")

phenoData_63061 <- phenoData_63061 %>%
  filter(SAMPLEID %in% common_63061) %>%
  mutate(rowname = SAMPLEID) %>%
  column_to_rownames("rowname")

phenoData_63060 <- phenoData_63060[colnames(assayData_63060), , drop = FALSE]
phenoData_63061 <- phenoData_63061[colnames(assayData_63061), , drop = FALSE]

stopifnot(identical(colnames(assayData_63060), rownames(phenoData_63060)))
stopifnot(identical(colnames(assayData_63061), rownames(phenoData_63061)))

# ---- save intermediate objects ----
saveRDS(raw_63060, here("results", "raw_63060.rds"))
saveRDS(raw_63061, here("results", "raw_63061.rds"))

saveRDS(assayData_63060, here("results", "assayData_63060_raw.rds"))
saveRDS(assayData_63061, here("results", "assayData_63061_raw.rds"))

saveRDS(phenoData_63060, here("results", "phenoData_63060_raw.rds"))
saveRDS(phenoData_63061, here("results", "phenoData_63061_raw.rds"))

rio::export(
  tibble::rownames_to_column(phenoData_63060, "ROWNAME"),
  here("results", "phenoData_63060_raw.csv")
)

rio::export(
  tibble::rownames_to_column(phenoData_63061, "ROWNAME"),
  here("results", "phenoData_63061_raw.csv")
)

message("Finished import and preliminary annotation.")
message("Next: inspect phenoData csv files and then build featureData from GPL annotation.")