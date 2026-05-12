# scripts/01_import_and_annotation.R

message("Starting 01_import_and_annotation.R")

# ---- load packages from libraries.txt ----
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

libraries_file <- here("libraries.txt")

if (!file.exists(libraries_file)) {
  stop("libraries.txt not found in project root.")
}

libraries <- unlist(read.table(libraries_file, stringsAsFactors = FALSE))

bioc_packages <- c(
  "GEOquery",
  "limma",
  "sva",
  "AnnotationDbi",
  "org.Hs.eg.db"
)

for (pkg in libraries) {
  if (!require(pkg, character.only = TRUE)) {
    message(paste("Installing missing package:", pkg))
    
    if (pkg %in% bioc_packages) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg)
    }
    
    library(pkg, character.only = TRUE)
  }
}

# ---- directories ----
datadir <- here("data")
resultsdir <- here("results")
plotsdir <- here("plots")
scriptsdir <- here("scripts")

dir.create(datadir, showWarnings = FALSE)
dir.create(resultsdir, showWarnings = FALSE)
dir.create(plotsdir, showWarnings = FALSE)
dir.create(scriptsdir, showWarnings = FALSE)

message("Project root detected by here(): ", here())

# ---- input files ----
raw_file_63060 <- here(datadir, "GSE63060_non-normalized.txt")
raw_file_63061 <- here(datadir, "GSE63061_non-normalized.txt")

series_file_63060 <- here(datadir, "GSE63060_series_matrix.txt")
series_file_63061 <- here(datadir, "GSE63061_series_matrix.txt")

stopifnot(file.exists(raw_file_63060))
stopifnot(file.exists(raw_file_63061))
stopifnot(file.exists(series_file_63060))
stopifnot(file.exists(series_file_63061))

# ---- helper functions ----

read_non_normalized_assay <- function(path) {
  message("Reading assay file: ", path)
  
  df <- read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fill = TRUE,
    quote = ""
  )
  
  if (ncol(df) < 2) {
    stop("Assay file appears to have fewer than 2 columns: ", path)
  }
  
  # remove unnamed columns
  bad_cols <- which(is.na(colnames(df)) | trimws(colnames(df)) == "")
  if (length(bad_cols) > 0) {
    message("Removing unnamed columns at positions: ", paste(bad_cols, collapse = ", "))
    df <- df[, -bad_cols, drop = FALSE]
  }
  
  # first column contains probe IDs, e.g. ID_REF / IND_ID
  colnames(df)[1] <- "FEATUREID"
  
  # remove empty feature rows
  df <- df[!is.na(df$FEATUREID) & trimws(df$FEATUREID) != "", , drop = FALSE]
  
  # convert expression columns to numeric
  for (j in 2:ncol(df)) {
    df[[j]] <- suppressWarnings(as.numeric(df[[j]]))
  }
  
  # collapse duplicate probes if present
  df <- tibble::as_tibble(df, .name_repair = "unique") %>%
    dplyr::group_by(FEATUREID) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::everything(),
        ~ {
          x <- .x[!is.na(.x)]
          if (length(x) == 0) NA_real_ else max(x)
        }
      ),
      .groups = "drop"
    )
  
  assay <- df %>%
    tibble::column_to_rownames("FEATUREID") %>%
    as.matrix()
  
  mode(assay) <- "numeric"
  
  return(assay)
}

clean_raw_colnames <- function(x) {
  x <- gsub("\\.AVG_Signal$", "", x)
  x <- gsub("\\.Detection Pval$", "", x)
  x <- gsub("\\.Pval$", "", x)
  x
}

read_series_metadata <- function(path) {
  message("Reading series matrix metadata: ", path)
  
  lines <- readLines(path, warn = FALSE)
  
  get_first_matching_line <- function(pattern) {
    hit <- grep(pattern, lines, value = TRUE)
    if (length(hit) == 0) {
      return(NULL)
    }
    hit[1]
  }
  
  split_tab <- function(x) {
    strsplit(x, "\t", fixed = TRUE)[[1]]
  }
  
  sample_title_line <- get_first_matching_line("^!Sample_title")
  sample_geo_line <- get_first_matching_line("^!Sample_geo_accession")
  sample_platform_line <- get_first_matching_line("^!Sample_platform_id")
  
  if (is.null(sample_title_line)) {
    stop("Could not find !Sample_title in series matrix file: ", path)
  }
  if (is.null(sample_geo_line)) {
    stop("Could not find !Sample_geo_accession in series matrix file: ", path)
  }
  if (is.null(sample_platform_line)) {
    stop("Could not find !Sample_platform_id in series matrix file: ", path)
  }
  
  sample_titles <- split_tab(sample_title_line)[-1]
  sample_titles <- gsub('"', '', sample_titles)
  sample_titles <- trimws(sample_titles)
  sample_geo <- split_tab(sample_geo_line)[-1]
  sample_platform <- split_tab(sample_platform_line)[-1]
  
  char_lines <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
  
  if (length(char_lines) == 0) {
    characteristics_all <- rep(NA_character_, length(sample_titles))
  } else {
    char_list <- lapply(char_lines, function(x) split_tab(x)[-1])
    
    # pad rows if needed
    max_len <- max(vapply(char_list, length, integer(1)))
    char_list <- lapply(char_list, function(x) {
      length(x) <- max_len
      x
    })
    
    char_mat <- do.call(rbind, char_list)
    char_df <- as.data.frame(char_mat, stringsAsFactors = FALSE)
    colnames(char_df) <- sample_titles[seq_len(ncol(char_df))]
    
    characteristics_all <- apply(char_df, 2, function(z) {
      paste(z[!is.na(z) & z != ""], collapse = " | ")
    })
  }
  
  pheno <- data.frame(
    SAMPLEID = sample_titles,
    GSM = sample_geo,
    PLATFORM = sample_platform,
    CHARACTERISTICS_ALL = characteristics_all,
    stringsAsFactors = FALSE
  )
  
  pheno
}

extract_group <- function(x) {
  dplyr::case_when(
    grepl("status:\\s*ad\\b", x, ignore.case = TRUE) ~ "AD",
    grepl("status:\\s*alzheimer", x, ignore.case = TRUE) ~ "AD",
    grepl("status:\\s*mci\\b", x, ignore.case = TRUE) ~ "MCI",
    grepl("status:\\s*control", x, ignore.case = TRUE) ~ "CONTROL",
    grepl("status:\\s*ctl\\b", x, ignore.case = TRUE) ~ "CONTROL",
    TRUE ~ NA_character_
  )
}

extract_characteristic_value <- function(x, key) {
  pattern <- paste0(key, ":\\s*([^|]+)")
  out <- stringr::str_match(x, pattern)[, 2]
  out <- trimws(out)
  out[out == ""] <- NA_character_
  out
}

# ---- import assay data from non-normalized files ----
assayData_63060 <- read_non_normalized_assay(raw_file_63060)
assayData_63061 <- read_non_normalized_assay(raw_file_63061)

colnames(assayData_63060) <- clean_raw_colnames(colnames(assayData_63060))
colnames(assayData_63061) <- clean_raw_colnames(colnames(assayData_63061))

message("Constructed assay matrices.")
message("assayData_63060 dim: ", paste(dim(assayData_63060), collapse = " x "))
message("assayData_63061 dim: ", paste(dim(assayData_63061), collapse = " x "))

message("First assay columns GSE63060:")
print(colnames(assayData_63060)[1:min(10, ncol(assayData_63060))])

message("First assay columns GSE63061:")
print(colnames(assayData_63061)[1:min(10, ncol(assayData_63061))])

# ---- import pheno data from Series Matrix files ----
phenoData_63060 <- read_series_metadata(series_file_63060) %>%
  dplyr::mutate(
    GROUP = extract_group(CHARACTERISTICS_ALL),
    STATUS = extract_characteristic_value(CHARACTERISTICS_ALL, "status"),
    ETHNICITY = extract_characteristic_value(CHARACTERISTICS_ALL, "ethnicity"),
    AGE = suppressWarnings(as.numeric(extract_characteristic_value(CHARACTERISTICS_ALL, "age"))),
    GENDER = extract_characteristic_value(CHARACTERISTICS_ALL, "gender"),
    TISSUE = extract_characteristic_value(CHARACTERISTICS_ALL, "tissue"),
    BATCH = "GSE63060",
    STUDY = "GSE63060"
  ) %>%
  dplyr::select(
    SAMPLEID, GSM, GROUP, STATUS, AGE, GENDER, ETHNICITY, TISSUE,
    BATCH, STUDY, PLATFORM, CHARACTERISTICS_ALL
  )

phenoData_63061 <- read_series_metadata(series_file_63061) %>%
  dplyr::mutate(
    GROUP = extract_group(CHARACTERISTICS_ALL),
    STATUS = extract_characteristic_value(CHARACTERISTICS_ALL, "status"),
    ETHNICITY = extract_characteristic_value(CHARACTERISTICS_ALL, "ethnicity"),
    AGE = suppressWarnings(as.numeric(extract_characteristic_value(CHARACTERISTICS_ALL, "age"))),
    GENDER = extract_characteristic_value(CHARACTERISTICS_ALL, "gender"),
    TISSUE = extract_characteristic_value(CHARACTERISTICS_ALL, "tissue"),
    BATCH = "GSE63061",
    STUDY = "GSE63061"
  ) %>%
  dplyr::select(
    SAMPLEID, GSM, GROUP, STATUS, AGE, GENDER, ETHNICITY, TISSUE,
    BATCH, STUDY, PLATFORM, CHARACTERISTICS_ALL
  )

message("Constructed preliminary phenoData.")
message("phenoData_63060 dim: ", paste(dim(phenoData_63060), collapse = " x "))
message("phenoData_63061 dim: ", paste(dim(phenoData_63061), collapse = " x "))

message("Group counts GSE63060:")
print(table(phenoData_63060$GROUP, useNA = "ifany"))

message("Group counts GSE63061:")
print(table(phenoData_63061$GROUP, useNA = "ifany"))


# ---- align assay and pheno by SAMPLEID ----
common_63060 <- intersect(colnames(assayData_63060), phenoData_63060$SAMPLEID)
common_63061 <- intersect(colnames(assayData_63061), phenoData_63061$SAMPLEID)

message("Matched samples GSE63060: ", length(common_63060))
message("Matched samples GSE63061: ", length(common_63061))

if (length(common_63060) == 0) {
  stop("No matching sample IDs found for GSE63060 between assayData and phenoData.")
}

if (length(common_63061) == 0) {
  stop("No matching sample IDs found for GSE63061 between assayData and phenoData.")
}

assayData_63060 <- assayData_63060[, common_63060, drop = FALSE]
assayData_63061 <- assayData_63061[, common_63061, drop = FALSE]

phenoData_63060 <- as.data.frame(phenoData_63060)
rownames(phenoData_63060) <- NULL
phenoData_63060 <- phenoData_63060 %>%
  dplyr::filter(SAMPLEID %in% common_63060) %>%
  tibble::column_to_rownames("SAMPLEID")

phenoData_63061 <- as.data.frame(phenoData_63061)
rownames(phenoData_63061) <- NULL
phenoData_63061 <- phenoData_63061 %>%
  dplyr::filter(SAMPLEID %in% common_63061) %>%
  tibble::column_to_rownames("SAMPLEID")

phenoData_63060 <- phenoData_63060[colnames(assayData_63060), , drop = FALSE]
phenoData_63061 <- phenoData_63061[colnames(assayData_63061), , drop = FALSE]

stopifnot(identical(colnames(assayData_63060), rownames(phenoData_63060)))
stopifnot(identical(colnames(assayData_63061), rownames(phenoData_63061)))

# ---- optional: save sample mapping tables ----
sample_map_63060 <- phenoData_63060 %>%
  tibble::rownames_to_column("SAMPLEID") %>%
  dplyr::select(SAMPLEID, GSM, GROUP, BATCH, STUDY, PLATFORM)

sample_map_63061 <- phenoData_63061 %>%
  tibble::rownames_to_column("SAMPLEID") %>%
  dplyr::select(SAMPLEID, GSM, GROUP, BATCH, STUDY, PLATFORM)

rio::export(sample_map_63060, here(resultsdir, "sample_map_63060.csv"))
rio::export(sample_map_63061, here(resultsdir, "sample_map_63061.csv"))

# ---- save intermediate objects ----
saveRDS(assayData_63060, here(resultsdir, "assayData_63060_raw.rds"))
saveRDS(assayData_63061, here(resultsdir, "assayData_63061_raw.rds"))

saveRDS(phenoData_63060, here(resultsdir, "phenoData_63060_raw.rds"))
saveRDS(phenoData_63061, here(resultsdir, "phenoData_63061_raw.rds"))

rio::export(
  tibble::rownames_to_column(phenoData_63060, "ROWNAME"),
  here(resultsdir, "phenoData_63060_raw.csv")
)

rio::export(
  tibble::rownames_to_column(phenoData_63061, "ROWNAME"),
  here(resultsdir, "phenoData_63061_raw.csv")
)

message("Finished 01_import_and_annotation.R")
message("Next steps:")
message("1) inspect phenoData csv files in results/")
message("2) download/load GPL6947 and GPL10558 annotation tables")
message("3) build featureData in 02_feature_annotation.R")