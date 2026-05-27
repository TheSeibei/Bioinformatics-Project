# scripts/05_functional_enrichment.R

message("Starting 05_functional_enrichment.R")

# ---- libraries ----
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

required_pkgs <- c(
  "clusterProfiler",
  "enrichplot",
  "org.Hs.eg.db",
  "DOSE",
  "dplyr",
  "ggplot2",
  "rio",
  "tibble"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    if (pkg %in% c(
      "clusterProfiler",
      "enrichplot",
      "org.Hs.eg.db",
      "DOSE"
    )) {
      
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

dir.create(resultsdir, showWarnings = FALSE)
dir.create(plotsdir, showWarnings = FALSE)

# ---- helper functions ----

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_\\-]", "_", x)
}

convert_symbols_to_entrez <- function(symbols) {
  
  conv <- bitr(
    symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  conv <- conv[!duplicated(conv$SYMBOL), ]
  
  return(conv)
}

run_go_enrichment <- function(entrez_ids, comparison_name) {
  
  ego <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message("No GO enrichment found for ", comparison_name)
    return(NULL)
  }
  
  # save table
  rio::export(
    as.data.frame(ego),
    here(resultsdir,
         paste0("GO_BP_", safe_filename(comparison_name), ".csv"))
  )
  
  # dotplot
  p <- dotplot(
    ego,
    showCategory = 20,
    font.size = 10,
    title = paste0("GO Biological Process: ", comparison_name)
  )
  
  ggplot2::ggsave(
    filename = here(
      plotsdir,
      paste0("DOTPLOT_GO_BP_", safe_filename(comparison_name), ".png")
    ),
    plot = p,
    width = 10,
    height = 8,
    dpi = 200
  )
  
  return(ego)
}

run_kegg_enrichment <- function(entrez_ids, comparison_name) {
  
  ekegg <- enrichKEGG(
    gene = entrez_ids,
    organism = "hsa",
    pvalueCutoff = 0.05
  )
  
  if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
    message("No KEGG enrichment found for ", comparison_name)
    return(NULL)
  }
  
  # convert IDs to readable gene symbols
  ekegg <- setReadable(
    ekegg,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID"
  )
  
  # save results
  rio::export(
    as.data.frame(ekegg),
    here(resultsdir,
         paste0("KEGG_", safe_filename(comparison_name), ".csv"))
  )
  
  # dotplot
  p <- dotplot(
    ekegg,
    showCategory = 20,
    font.size = 10,
    title = paste0("KEGG Pathways: ", comparison_name)
  )
  
  ggplot2::ggsave(
    filename = here(
      plotsdir,
      paste0("DOTPLOT_KEGG_", safe_filename(comparison_name), ".png")
    ),
    plot = p,
    width = 10,
    height = 8,
    dpi = 200
  )
  
  return(ekegg)
}

# ---- find DEG result files ----

deg_files <- list.files(
  resultsdir,
  pattern = "^DEG_.*_sig_FDR05\\.csv$",
  full.names = TRUE
)

message("Found DEG files:")
print(basename(deg_files))

# ---- run enrichment for each contrast ----

for (f in deg_files) {
  
  comparison_name <- tools::file_path_sans_ext(basename(f))
  comparison_name <- gsub("^DEG_", "", comparison_name)
  comparison_name <- gsub("_sig_FDR05$", "", comparison_name)
  
  message("Processing: ", comparison_name)
  
  deg <- rio::import(f)
  
  if (!"SYMBOL" %in% colnames(deg)) {
    message("Skipping ", comparison_name, " because SYMBOL column missing.")
    next
  }
  
  symbols <- unique(deg$SYMBOL)
  symbols <- symbols[!is.na(symbols)]
  symbols <- symbols[symbols != ""]
  
  message("Number of significant symbols: ", length(symbols))
  
  # convert to ENTREZ
  conv <- convert_symbols_to_entrez(symbols)
  
  if (nrow(conv) == 0) {
    message("No ENTREZ IDs found for ", comparison_name)
    next
  }
  
  entrez_ids <- unique(conv$ENTREZID)
  
  message("Mapped ENTREZ IDs: ", length(entrez_ids))
  
  # GO enrichment
  run_go_enrichment(entrez_ids, comparison_name)
  
  # KEGG enrichment
  run_kegg_enrichment(entrez_ids, comparison_name)
}

message("Finished 05_functional_enrichment.R")