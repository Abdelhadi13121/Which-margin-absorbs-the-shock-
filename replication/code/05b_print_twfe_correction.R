# =============================================================================
# 05b_print_twfe_correction.R
#
# The estimation in 05_fix_twfe_band.R completed and wrote
# output/tables/rob2_A_twfe_corrected.csv before the final print block failed on
# a naming bug (TH was labelled with INTER, "shk_fx", while the comparison table
# indexed tb[["fx"]]). 05 has since been fixed, but there is no need to re-run
# it: this reads the saved results and produces the report.
#
#   source("code/05b_print_twfe_correction.R")
#
# Takes seconds. Requires only the CSV that the earlier run already wrote.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

f <- "output/tables/rob2_A_twfe_corrected.csv"
if (!file.exists(f)) stop("Not found: ", f,
  "\nRun code/05_fix_twfe_band.R first (it is now fixed).", call. = FALSE)
out <- fread(f)

stars <- function(p) ifelse(is.na(p), "", ifelse(p < .01, "***",
                     ifelse(p < .05, "**", ifelse(p < .10, "*", ""))))
fx <- out[kind == "estimate" & term == "fx"]

sink("output/tables/CORRECTION_twfe_band.txt")
cat("CORRECTED two-way fixed-effects band average, h = 12-24\n")
cat("generated", format(Sys.time()), "\n")
cat(strrep("=", 74), "\n\n")
cat("The earlier version of this column averaged per-horizon covariance matrices\n")
cat("as V = sum(V_h)/H^2, which assumes the horizons are independent. They are not:\n")
cat("they share regressors and overlapping 12-month windows. That expression drops\n")
cat("every cross-horizon covariance term and understates the standard errors by\n")
cat("roughly a factor of H. The reported t of -11.79 on fx was an artefact of it.\n\n")
cat("This version stacks the per-horizon influence functions, averages them, and\n")
cat("clusters at country level, so the cross-horizon covariance enters the variance\n")
cat("of the average. It is the same construction used for the DML band estimates.\n\n")
cat(sprintf("Sample: %d observations, %d clusters\n\n", out$n[1], out$G[1]))

cat("Estimates\n")
print(out[kind == "estimate", .(term, coef = round(coef, 4), se = round(se, 4),
      t = round(t, 2), p = round(p, 4), sig = stars(p),
      CI = sprintf("[%.4f, %.4f]", lo, hi))])
cat("\nContrasts\n")
print(out[kind == "contrast", .(term, diff = round(coef, 4), se = round(se, 4),
      t = round(t, 2), p = round(p, 4), sig = stars(p))])
cat("\nJoint equality\n")
print(out[kind == "joint", .(W = round(coef, 3), p = round(p, 4), sig = stars(p))])

cat("\n\nCOMPARISON, fx coefficient, long band\n\n")
print(data.table(
  Specification = c("DML random forest (paper, Table 5)",
                    "DML gradient boosting",
                    "Two-way FE, as first reported (INVALID)",
                    "Two-way FE, corrected"),
  coef = c(-0.0417, -0.0275, -0.0353, round(fx$coef, 4)),
  se   = c( 0.0193,  0.0269,  0.0030, round(fx$se,   4)),
  t    = c(-2.17,   -1.02,   -11.79,  round(fx$t,    2)),
  p    = c( 0.036,   0.313,   0.000,  round(fx$p,    4))))
cat("\nReport the corrected FE column; do not report the original one.\n")
sink()

cat(readLines("output/tables/CORRECTION_twfe_band.txt"), sep = "\n")
cat("\n\nWritten: output/tables/CORRECTION_twfe_band.txt\n")
