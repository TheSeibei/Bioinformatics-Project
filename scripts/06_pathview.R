# scripts/06_pathview.R
#
# Visualisiert zwei KEGG-Pathways mit pathview, eingefaerbt nach den logFC
# aus dem Kontrast AD_vs_CONTROL:
#   - hsa00190  Oxidative phosphorylation
#   - hsa05010  Alzheimer disease
#
# Hinweis: pathview laedt die Pathway-Karten zur Laufzeit von KEGG
# (www.kegg.jp) herunter -- dieses Skript braucht also eine Internetverbindung.

message("Starting 06_pathview.R")

# ---- libraries ----
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

required_pkgs <- c("pathview", "org.Hs.eg.db", "clusterProfiler", "dplyr", "rio")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("pathview", "org.Hs.eg.db", "clusterProfiler")) {
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
plotsdir   <- here("plots")
kegg_cache <- here("results", "kegg_cache")   # Zwischenspeicher fuer KGML/PNG

dir.create(plotsdir, showWarnings = FALSE, recursive = TRUE)
dir.create(kegg_cache, showWarnings = FALSE, recursive = TRUE)

# ---- settings ----
contrast_prefix <- "AD_vs_CONTROL"
pathways <- c(
  "00190" = "Oxidative phosphorylation",
  "05010" = "Alzheimer disease"
)

# ---- load DE results ----
deg_rds <- here(resultsdir, paste0("DEG_", contrast_prefix, "_all.rds"))
deg_csv <- here(resultsdir, paste0("DEG_", contrast_prefix, "_all.csv"))

if (file.exists(deg_rds)) {
  deg <- readRDS(deg_rds)
} else if (file.exists(deg_csv)) {
  deg <- rio::import(deg_csv)
} else {
  stop("DE-Ergebnisse fuer ", contrast_prefix, " nicht gefunden. ",
       "Bitte zuerst 04_differential_expression.R laufen lassen.")
}

if (!all(c("SYMBOL", "logFC") %in% colnames(deg))) {
  stop("Spalten SYMBOL und/oder logFC fehlen in den DE-Ergebnissen.")
}

# ---- map gene symbols to ENTREZ (KEGG-Standard-ID fuer hsa) ----
symbols <- unique(deg$SYMBOL[!is.na(deg$SYMBOL) & deg$SYMBOL != ""])

conv <- clusterProfiler::bitr(
  symbols,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
conv <- conv[!duplicated(conv$SYMBOL), ]

deg_mapped <- merge(deg[, c("SYMBOL", "logFC")], conv, by = "SYMBOL")

# Falls mehrere Symbole auf dieselbe ENTREZ-ID fallen: das mit dem
# staerksten Effekt (groesstes |logFC|) behalten.
deg_mapped <- deg_mapped %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::slice_max(abs(logFC), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

# Benannter logFC-Vektor: Namen = ENTREZ-IDs, Werte = logFC.
gene_vec <- deg_mapped$logFC
names(gene_vec) <- deg_mapped$ENTREZID

message("Genes mapped to ENTREZ for pathview: ", length(gene_vec))

# Symmetrische Farbskala um 0 (gleicher Betrag fuer hoch- und runterreguliert).
fc_limit <- ceiling(max(abs(gene_vec), na.rm = TRUE) * 10) / 10
if (!is.finite(fc_limit) || fc_limit <= 0) fc_limit <- 1

# ---- pathview helper ----
# pathview schreibt die fertige PNG ins aktuelle Arbeitsverzeichnis, daher
# wechseln wir kurz nach plots/ und stellen das wd anschliessend wieder her.
run_pathview <- function(gene_vec, pid, pname, cache, outdir, suffix, limit) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(outdir)
  
  message("Pathview: hsa", pid, " (", pname, ")")
  
  tryCatch(
    {
      pathview::pathview(
        gene.data  = gene_vec,
        pathway.id = pid,
        species    = "hsa",
        out.suffix = suffix,
        kegg.dir   = cache,
        limit      = list(gene = limit),
        bins       = list(gene = 20),
        low        = list(gene = "#2C7BB6"),  # blau  = runterreguliert
        mid        = list(gene = "#F7F7F7"),  # weiss = neutral
        high       = list(gene = "#D7191C")   # rot   = hochreguliert
      )
    },
    error = function(e) {
      message("  -> fehlgeschlagen fuer hsa", pid, ": ", conditionMessage(e))
    }
  )
}

# ---- run pathview for both pathways ----
for (pid in names(pathways)) {
  run_pathview(
    gene_vec = gene_vec,
    pid      = pid,
    pname    = pathways[[pid]],
    cache    = kegg_cache,
    outdir   = plotsdir,
    suffix   = contrast_prefix,
    limit    = fc_limit
  )
}

message("Finished 06_pathview.R")
message("Output PNGs (z.B. hsa05010.", contrast_prefix, ".png) liegen in: ", plotsdir)