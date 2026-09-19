# =============================================================================
# 05_fix_twfe_band.R
#
# CORRECTION to Panel A of 04_robustness_learner_composition.R.
#
# WHAT WAS WRONG
#   band_fe() computed the variance of the band average as
#       V = Reduce(`+`, V_h) / H^2
#   which is the formula for averaging INDEPENDENT estimates. The thirteen
#   horizons in the long band are heavily correlated — they share regressors and
#   overlapping 12-month windows — so that expression omits every cross-horizon
#   covariance term and understates the variance by roughly a factor of H.
#   In the run of 2026-09-19 it produced t = -11.79 on the fx coefficient, which
#   is not a credible standard error for a macro panel of 46 clusters.
#
# WHAT THIS DOES INSTEAD
#   The same construction already used for the DML band averages: stack the
#   per-horizon influence functions, average them, and cluster the resulting
#   n x M matrix at country level. For OLS with two-way fixed effects the
#   influence function of coefficient k in horizon h is
#
#       IF_h,k = (X'X/n)^{-1}[k, ] %*% (x_i * e_i,h)
#
#   evaluated on the double-demeaned design. Averaging across horizons and then
#   clustering gives a variance that carries the cross-horizon covariance, which
#   is what the band average needs.
#
#   Because the horizon samples differ slightly (each h loses observations at the
#   end of each country's series), the band is estimated on a FIXED sample: the
#   rows for which every horizon in the band is observed. This matches the DML
#   band procedure and is what makes the influence functions conformable.
#
# RUNTIME  A few minutes. No machine learning, nothing cached from the DML runs
#          is touched or needed.
#
# USAGE (from the repository root)
#   setwd("path/to/Which-margin-absorbs-the-shock-")
#   source("code/05_fix_twfe_band.R")
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

PATH      <- "data/derived/view_dml_lp.csv"
N_LAGS    <- 6L
BUST_Q    <- 0.25
H_COND    <- 12L
LONG_BAND <- 12:24

CONTROLS_LEAN <- c("d_log_fx_12m", "fx_vol_12m",
                   "d_log_exports_usd_12m", "d_log_imports_usd_12m")
MARGINS <- c(imports  = "d_log_imports_usd_12m",
             reserves = "d_log_reserves_12m",
             fx       = "d_log_fx_12m")
INTER   <- paste0("shk_", names(MARGINS))
own_lags <- function(v) paste0(v, "_L", seq_len(N_LAGS))

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------- data ------
# Identical preparation to 02 and 04, so the sample matches exactly.
d <- fread(PATH)
d[, date := as.IDate(date)]
d[, mo := year(date) * 12L + month(date)]
d[, cl := .GRP, by = iso3]
setorder(d, iso3, date)

d[, res_pos := fifelse(reserves_exgold_usd > 0, reserves_exgold_usd, NA_real_)]
d[, d_log_reserves_12m := 100 * (log(res_pos) - log(shift(res_pos, 12))), by = iso3]
d[, mo_12 := shift(mo, 12), by = iso3]
d[is.na(mo_12) | (mo - mo_12) != 12L, d_log_reserves_12m := NA_real_]
d[, c("res_pos", "mo_12") := NULL]
d[, infl_z := (infl_yoy_w - mean(infl_yoy_w, na.rm = TRUE)) /
              sd(infl_yoy_w, na.rm = TRUE), by = iso3]

bust_cut <- quantile(d$d_shock_12m, BUST_Q, na.rm = TRUE)
d[, bust := as.integer(d_shock_12m <= bust_cut)]
for (m in names(MARGINS)) {
  yc <- MARGINS[[m]]
  d[, (paste0("fwd_", m)) := shift(get(yc), n = -H_COND, type = "lag"), by = iso3]
  d[, mo_f := shift(mo, n = -H_COND, type = "lag"), by = iso3]
  d[is.na(mo_f) | (mo_f - mo) != H_COND, (paste0("fwd_", m)) := NA_real_]
}
d[, mo_f := NULL]
ABS <- c(imports = -1, reserves = -1, fx = +1)
ep <- d[bust == 1 & !is.na(fwd_imports) & !is.na(fwd_reserves) & !is.na(fwd_fx)]
for (m in names(MARGINS)) {
  v <- ep[[paste0("fwd_", m)]] * ABS[[m]]
  ep[, (paste0("z_", m)) := (v - mean(v)) / sd(v)]
}
for (m in names(MARGINS)) {
  zz <- paste0("z_", m)
  ep[, (paste0("hist_", m)) := (sum(get(zz)) - get(zz)) / pmax(.N - 1L, 1L), by = iso3]
}
hc <- paste0("hist_", names(MARGINS))
ep[, margin_class_used := names(MARGINS)[max.col(as.matrix(.SD), ties.method = "first")],
   .SDcols = hc]
ep[, n_ep := .N, by = iso3]
ep[n_ep < 12L, margin_class_used := NA_character_]
d <- merge(d, ep[, .(iso3, date, margin_class_used)], by = c("iso3","date"), all.x = TRUE)
for (m in names(MARGINS))
  d[, (paste0("shk_", m)) := shock_neg *
      as.integer(!is.na(margin_class_used) & margin_class_used == m)]
for (v in c("infl_z", INTER)) for (L in seq_len(N_LAGS))
  d[, paste0(v, "_L", L) := shift(get(v), L), by = iso3]

make_y <- function(w, ycol, h, nm) {
  w[, (nm) := shift(get(ycol), n = -h, type = "lag"), by = iso3]
  w[, mo_h := shift(mo, n = -h, type = "lag"), by = iso3]
  w[is.na(mo_h) | (mo_h - mo) != h, (nm) := NA_real_]
  w[, mo_h := NULL]; w
}
demean <- function(dt, cols) {
  for (v in cols) {
    dt[, (v) := get(v) - mean(get(v), na.rm = TRUE), by = iso3]
    dt[, (v) := get(v) - mean(get(v), na.rm = TRUE), by = date]
  }
  dt
}

# ------------------------------------------------- fixed-sample band design --
hs <- LONG_BAND
w  <- copy(d); ycols <- character(0)
for (h in hs) { nm <- sprintf("y_%02d", h); w <- make_y(w, "infl_z", h, nm)
                ycols <- c(ycols, nm) }
xall <- c(CONTROLS_LEAN, "shock_pos", INTER,
          own_lags("infl_z"), unlist(lapply(INTER, own_lags)))
w <- na.omit(w[, c("iso3","date","cl", ycols, xall), with = FALSE])
n <- nrow(w); G <- uniqueN(w$cl); cl <- w$cl
cat(sprintf("fixed-sample band: n = %d, clusters = %d\n", n, G))

# Two-way fixed effects are removed by double demeaning, exactly as in the DML
# pipeline, so the OLS below is the within estimator.
w <- demean(w, c(ycols, xall))
X <- as.matrix(w[, xall, with = FALSE])
k_idx <- match(INTER, xall)
XtX_inv <- solve(crossprod(X) / n)
M <- length(INTER); H <- length(hs)

TH  <- matrix(NA_real_, H, M, dimnames = list(NULL, names(MARGINS)))
IFs <- matrix(0, n, M)

for (i in seq_along(hs)) {
  y  <- w[[ycols[i]]]
  b  <- qr.solve(X, y)
  e  <- as.numeric(y - X %*% b)
  TH[i, ] <- b[k_idx]
  # influence function for each target coefficient at this horizon
  A <- XtX_inv[k_idx, , drop = FALSE]        # M x p
  IFs <- IFs + (X * e) %*% t(A)              # n x M, accumulated across horizons
}
IFs <- IFs / H
tb  <- colMeans(TH)

# cluster-robust covariance of the band average, with finite-cluster correction
Vb <- crossprod(rowsum(IFs, cl)) / n^2 * (G / (G - 1))
se <- sqrt(diag(Vb)); tc <- qt(0.975, G - 1)
# name both by margin so the comparison table below can index them by name
names(tb) <- names(se) <- names(MARGINS)

est <- data.table(spec = "Two-way FE (no ML), corrected",
                  term = names(MARGINS), coef = tb, se = se, t = tb / se,
                  p = 2 * pt(-abs(tb / se), df = G - 1),
                  lo = tb - tc * se, hi = tb + tc * se,
                  n = n, G = G, kind = "estimate")

ct <- rbindlist(lapply(combn(M, 2, simplify = FALSE), function(pp) {
  i <- pp[1]; j <- pp[2]; dif <- tb[i] - tb[j]
  s2 <- sqrt(Vb[i,i] + Vb[j,j] - 2*Vb[i,j])
  data.table(spec = "Two-way FE (no ML), corrected",
             term = paste(names(MARGINS)[i], "-", names(MARGINS)[j]),
             coef = dif, se = s2, t = dif/s2, p = 2*pt(-abs(dif/s2), df = G-1),
             lo = dif - tc*s2, hi = dif + tc*s2, n = n, G = G, kind = "contrast")
}))

R  <- rbind(c(1,-1,0), c(0,1,-1)); Rb <- R %*% tb
W  <- as.numeric(t(Rb) %*% solve(R %*% Vb %*% t(R)) %*% Rb)
jt <- data.table(spec = "Two-way FE (no ML), corrected", term = "JOINT equality",
                 coef = W, se = NA_real_, t = NA_real_, p = 1 - pchisq(W, 2),
                 lo = NA_real_, hi = NA_real_, n = n, G = G, kind = "joint")

out <- rbind(est, ct, jt)
fwrite(out, "output/tables/rob2_A_twfe_corrected.csv")

# ---------------------------------------------------------------- report ----
stars <- function(p) ifelse(is.na(p), "", ifelse(p < .01, "***",
                     ifelse(p < .05, "**", ifelse(p < .10, "*", ""))))

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
cat(sprintf("Sample: %d observations, %d clusters, critical value %.3f\n\n", n, G, tc))

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
cmp <- data.table(
  Specification = c("DML random forest (paper, Table 5)",
                    "DML gradient boosting",
                    "Two-way FE, as first reported (INVALID)",
                    "Two-way FE, corrected"),
  coef = c(-0.0417, -0.0275, -0.0353, round(tb[["fx"]], 4)),
  se   = c(0.0193, 0.0269, 0.0030, round(se[["fx"]], 4)),
  t    = c(-2.17, -1.02, -11.79, round((tb/se)[["fx"]], 2)))
print(cmp)
cat("\nThe point estimate is stable across all three valid specifications. What\n")
cat("differs is precision. Report the corrected FE column; do not report the\n")
cat("original one.\n")
sink()

cat(readLines("output/tables/CORRECTION_twfe_band.txt"), sep = "\n")
cat("\n\nWritten: output/tables/CORRECTION_twfe_band.txt and rob2_A_twfe_corrected.csv\n")
