# =============================================================================
# 00_run_all.R  —  master replication script
#
# Which margin absorbs the shock? Commodity terms of trade, external adjustment
# and inflation in resource-rich economies.
#
# USAGE
#   setwd("path/to/Which-margin-absorbs-the-shock-")
#   source("code/00_run_all.R")
#
# WHAT IT DOES
#   1. checks that the working directory is the repository root
#   2. checks packages and reports versions against those used in the paper
#   3. runs the main estimation pipeline   (code/02_replicate_empirical.R)
#   4. runs the robustness exercise        (code/03_robustness_pricesetters.R)
#   5. writes a session log to output/logs/
#
# RUNTIME
#   Cold run: several hours on one core. Every stage caches to
#   output/tables/raw_*.csv, so a second run rebuilds all tables and figures in
#   seconds. Delete a raw_*.csv, or set FORCE <- TRUE in a script, to re-estimate.
# =============================================================================

t_start <- Sys.time()

# ---------------------------------------------------------------- 1. location
need <- c("code/02_replicate_empirical.R", "data/derived/view_dml_lp.csv")
if (!all(file.exists(need))) {
  stop("Run this from the repository root. Missing: ",
       paste(need[!file.exists(need)], collapse = ", "),
       "\nCurrent working directory: ", getwd(), call. = FALSE)
}
for (d in c("output/tables", "output/figures", "output/logs"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------- 2. packages
REQUIRED <- c("data.table", "DoubleML", "mlr3", "mlr3learners",
              "ranger", "xgboost", "fixest", "ggplot2")
missing <- REQUIRED[!vapply(REQUIRED, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       "\nInstall with:\n  install.packages(c(",
       paste0('"', missing, '"', collapse = ", "), "))", call. = FALSE)
}

# versions used to produce the published results
PAPER <- c(R = "4.5.0", data.table = "1.17", DoubleML = "1.0", mlr3 = "1.0",
           mlr3learners = "0.9", ranger = "0.17", xgboost = "1.7",
           fixest = "0.12", ggplot2 = "3.5")

message("\n", strrep("=", 74))
message("REPLICATION: Which margin absorbs the shock?")
message(strrep("=", 74))
message(sprintf("  %-16s %-12s %s", "component", "installed", "used in paper"))
message(sprintf("  %-16s %-12s %s", "R",
                paste(R.version$major, R.version$minor, sep = "."), PAPER[["R"]]))
for (p in REQUIRED)
  message(sprintf("  %-16s %-12s %s", p,
                  as.character(utils::packageVersion(p)),
                  if (p %in% names(PAPER)) PAPER[[p]] else "-"))
message("")
message("  Results are seed-controlled and reproduce exactly on these versions.")
message("  Minor differences are possible across ranger releases, whose internal")
message("  RNG has changed; coefficients are stable to three decimal places")
message("  across the versions we have tested.")
message(strrep("=", 74), "\n")

# ------------------------------------------------------------------ 3. logging
log_path <- file.path("output", "logs",
                      format(Sys.time(), "run_%Y%m%d_%H%M%S.log"))
con <- file(log_path, open = "wt")
sink(con, split = TRUE); sink(con, type = "message")
on.exit({ sink(type = "message"); sink(); close(con) }, add = TRUE)

cat("Replication run started", format(t_start), "\n")
cat(R.version.string, "|", Sys.info()[["sysname"]], "\n\n")

# --------------------------------------------------------------- 4. pipeline
run_stage <- function(script, label) {
  cat("\n", strrep("-", 74), "\n", label, "\n", strrep("-", 74), "\n", sep = "")
  t0 <- Sys.time()
  ok <- tryCatch({ source(script, echo = FALSE, local = new.env()); TRUE },
                 error = function(e) { message("FAILED: ", conditionMessage(e)); FALSE })
  cat(sprintf("\n%s: %s in %.1f min\n", label, if (ok) "completed" else "FAILED",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  ok
}

ok1 <- run_stage("code/02_replicate_empirical.R",
                 "STAGE 1  Main estimation: Tables 1-6, Appendix A1-A2, Figures 1-6")
ok2 <- run_stage("code/03_robustness_pricesetters.R",
                 "STAGE 2  Robustness: dominant commodity price-setters (Table 7)")
ok3 <- run_stage("code/04_robustness_learner_composition.R",
                 "STAGE 3  Robustness: nuisance learner and composition controls (Table 8)")
ok4 <- run_stage("code/05_fix_twfe_band.R",
                 "STAGE 4  Corrected two-way FE band average (supersedes the FE row of Stage 3)")

# ------------------------------------------------------------------ 5. summary
cat("\n", strrep("=", 74), "\n", sep = "")
cat("DONE in ", sprintf("%.1f", as.numeric(difftime(Sys.time(), t_start,
                                                    units = "mins"))), " minutes\n", sep = "")
cat("  main estimation : ", if (ok1) "ok" else "FAILED", "\n", sep = "")
cat("  price-setters   : ", if (ok2) "ok" else "FAILED", "\n", sep = "")
cat("  learner/composition: ", if (ok3) "ok" else "FAILED", "\n", sep = "")
cat("  corrected FE band  : ", if (ok4) "ok" else "FAILED", "\n\n", sep = "")
cat("  tables   -> output/tables/    (ALL_TABLES.txt, tab*.csv)\n")
cat("  figures  -> output/figures/   (Fig1-Fig6, .tiff at 600 dpi and .png)\n")
cat("  log      -> ", log_path, "\n\n", sep = "")
cat("Check your numbers against docs/results_summary.md.\n")
cat("If they differ by more than rounding, please open an issue.\n")
