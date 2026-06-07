# run_all.R
library(here)

scripts <- list.files(
  here("scripts"),
  pattern = "^[0-9]{2}_.*\\.R$",   # 01_..., 02_..., usw.
  full.names = TRUE
)
scripts <- sort(scripts)            # garantiert die Reihenfolge 01 -> 05

for (s in scripts) {
  message("\n========================================")
  message("Running: ", basename(s))
  message("========================================\n")
  t0 <- Sys.time()
  
  tryCatch(
    source(s, echo = FALSE),
    error = function(e) {
      stop("Abbruch bei ", basename(s), ": ", conditionMessage(e), call. = FALSE)
    }
  )
  
  message(basename(s), " fertig in ",
          round(difftime(Sys.time(), t0, units = "mins"), 2), " min")
}

message("\nAlle Skripte erfolgreich durchgelaufen.")