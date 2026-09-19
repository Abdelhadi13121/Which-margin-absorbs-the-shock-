# =============================================================================
# robustness_pricesetters.R
#
# Referee-anticipating robustness check: are the results driven by countries
# large enough to influence the world prices that build the shock?
#
# The Gruss-Kebhaj index fixes trade weights at pre-sample shares, so for the
# shock to be endogenous the causal path has to run from a country's CURRENT
# domestic conditions to WORLD commodity prices. That is plausible for Saudi
# Arabia in oil and Chile in copper, weaker for Russia, and largely absent
# elsewhere. Dropping the large producers is therefore a partial answer, not a
# clean identification argument, and the paper should say so.
#
# TWO VARIANTS
#   conservative  drop SAU, RUS, CHL          -> the headline robustness table
#   aggressive    drop the ten largest        -> reported as a further check
#
# WHAT IT RE-ESTIMATES (not the whole pipeline)
#   1. baseline inflation response, random forest, h = -6..24
#   2. band-averaged margin-conditional response, long band h = 12..24,
#      with contrasts and the joint equality test
# Everything else in the paper is unchanged by the exclusion and is not re-run.
#
# USAGE
#   setwd("path/to/repo"); source("robustness_pricesetters.R")
#
# Results cache to out/rob_*.csv exactly as the main script does, so a second
# run rebuilds the table in seconds.
#
# NOT RUN BEFORE DELIVERY: no R in the drafting environment.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(DoubleML); library(mlr3); library(mlr3learners)
})
lgr::get_logger("mlr3")$set_threshold("warn")

SEED      <- 20260914L
PATH      <- "data/derived/view_dml_lp.csv"
FORCE     <- FALSE
N_FOLDS   <- 5L
N_REP     <- 1L          # 1 repetition: this is a robustness check, not the headline
MIN_OBS   <- 400L
MIN_TREAT <- 60L         # lowered from 80: the excluded samples are smaller
N_LAGS    <- 6L
H_MIN     <- -6L
H_MAX     <- 24L
BUST_Q    <- 0.25
H_COND    <- 12L
LONG_BAND <- 12:24

# Exclusion sets. Chosen on share of world production of the country's main
# export commodity, not on sample influence.
DROP <- list(
  full         = character(0),
  conservative = c("SAU", "RUS", "CHL"),
  aggressive   = c("SAU", "RUS", "IRN", "IRQ", "ARE", "KWT", "CHL", "PER",
                   "QAT", "DZA")
)

CONTROLS_LEAN <- c("d_log_fx_12m", "fx_vol_12m",
                   "d_log_exports_usd_12m", "d_log_imports_usd_12m")
MARGINS <- c(imports  = "d_log_imports_usd_12m",
             reserves = "d_log_reserves_12m",
             fx       = "d_log_fx_12m")

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
cached <- function(f, expr) {
  fp <- file.path("out", f)
  if (!FORCE && file.exists(fp)) { message("  [cache] ", f); return(fread(fp)) }
  r <- force(expr); fwrite(r, fp); r
}

# -----------------------------------------------------------------------------
# Data preparation, identical to the main script
# -----------------------------------------------------------------------------

prep <- function() {
  d <- fread(PATH)
  d[, date := as.IDate(date)]
  d[, mo := year(date) * 12L + month(date)]
  setorder(d, iso3, date)
  stopifnot(nrow(d) == 21560L, uniqueN(d$iso3) == 49L)

  d[, res_pos := fifelse(reserves_exgold_usd > 0, reserves_exgold_usd, NA_real_)]
  d[, d_log_reserves_12m := 100 * (log(res_pos) - log(shift(res_pos, 12))), by = iso3]
  d[, mo_12 := shift(mo, 12), by = iso3]
  d[is.na(mo_12) | (mo - mo_12) != 12L, d_log_reserves_12m := NA_real_]
  d[, c("res_pos", "mo_12") := NULL]

  d[, infl_z := (infl_yoy_w - mean(infl_yoy_w, na.rm = TRUE)) /
                sd(infl_yoy_w, na.rm = TRUE), by = iso3]

  # Classification is built on the FULL sample and then carried into the
  # restricted samples. Rebuilding it inside each subsample would change the
  # standardisation and confound the exclusion with a different classification.
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

  lagv <- unique(c("infl_yoy_w", "infl_z", "shock_neg", "shock_pos",
                   paste0("shk_", names(MARGINS))))
  for (v in lagv) for (L in seq_len(N_LAGS))
    d[, paste0(v, "_L", L) := shift(get(v), L), by = iso3]
  d[]
}

own_lags <- function(v) paste0(v, "_L", seq_len(N_LAGS))
INTER <- paste0("shk_", names(MARGINS))

demean <- function(dt, cols) {
  for (v in cols) {
    dt[, (v) := get(v) - mean(get(v), na.rm = TRUE), by = iso3]
    dt[, (v) := get(v) - mean(get(v), na.rm = TRUE), by = date]
  }
  dt
}
make_y <- function(w, ycol, h, nm = "y_h") {
  w[, (nm) := shift(get(ycol), n = -h, type = "lag"), by = iso3]
  w[, mo_h := shift(mo, n = -h, type = "lag"), by = iso3]
  w[is.na(mo_h) | (mo_h - mo) != h, (nm) := NA_real_]
  w[, mo_h := NULL]; w
}

lrn_rf <- lrn("regr.ranger", num.trees = 500, min.node.size = 10,
              max.depth = 8, num.threads = 1)
cluster_folds <- function(cl, nf, nr, seed) {
  ids <- sort(unique(cl))
  lapply(seq_len(nr), function(r) {
    set.seed(seed + r)
    a <- sample(rep_len(seq_len(nf), length(ids))); names(a) <- as.character(ids)
    f <- a[as.character(cl)]
    list(train_ids = lapply(seq_len(nf), function(k) which(f != k)),
         test_ids  = lapply(seq_len(nf), function(k) which(f == k)))
  })
}
FOLD_FALLBACK <- FALSE
fit_plr <- function(df, ycol, dcol, xcols, cl) {
  mk <- function() double_ml_data_from_data_frame(df, y_col = ycol, d_cols = dcol, x_cols = xcols)
  o <- DoubleMLPLR$new(mk(), ml_l = lrn_rf$clone(), ml_m = lrn_rf$clone(),
                       n_folds = N_FOLDS, n_rep = N_REP,
                       score = "partialling out", draw_sample_splitting = FALSE)
  ok <- tryCatch({ o$set_sample_splitting(cluster_folds(cl, N_FOLDS, N_REP, SEED)); TRUE },
                 error = function(e) FALSE)
  if (!ok) {
    if (!FOLD_FALLBACK) { warning("default folds used", immediate. = TRUE); FOLD_FALLBACK <<- TRUE }
    set.seed(SEED)
    o <- DoubleMLPLR$new(mk(), ml_l = lrn_rf$clone(), ml_m = lrn_rf$clone(),
                         n_folds = N_FOLDS, n_rep = N_REP, score = "partialling out")
  }
  o$fit(); o
}
influence <- function(o, k = 1L) {
  a <- o$psi_a[, 1, k]; b <- o$psi_b[, 1, k]
  (a * (-mean(b) / mean(a)) + b) / mean(a)
}
clustV <- function(IF, cl) {
  n <- nrow(IF); G <- length(unique(cl))
  crossprod(rowsum(IF, cl)) / n^2 * (G / (G - 1))
}

# -----------------------------------------------------------------------------
# 1. Baseline response on a restricted sample
# -----------------------------------------------------------------------------

baseline_on <- function(dat, tag) {
  rbindlist(lapply(H_MIN:H_MAX, function(h) {
    cat(sprintf("  [%s] baseline h=%3d\n", tag, h)); flush.console()
    w <- make_y(copy(dat), "infl_yoy_w", h)
    xc <- c(CONTROLS_LEAN, "shock_pos", own_lags("infl_yoy_w"),
            own_lags("shock_neg"), own_lags("shock_pos"))
    w <- na.omit(w[, c("iso3","date","cl","y_h","shock_neg", xc), with = FALSE])
    if (nrow(w) < MIN_OBS) return(NULL)
    n <- nrow(w); G <- uniqueN(w$cl)
    w <- demean(w, c("y_h","shock_neg", xc)); cl <- w$cl
    w[, c("iso3","date","cl") := NULL]
    o <- fit_plr(as.data.frame(w), "y_h", "shock_neg", xc, cl)
    th <- o$coef[[1]]; se <- sqrt(clustV(cbind(influence(o)), cl)[1,1])
    tc <- qt(0.975, G - 1)
    data.table(sample = tag, h = h, coef = th, se_cl = se, t = th/se,
               p = 2*pt(-abs(th/se), df = G-1),
               lo = th - tc*se, hi = th + tc*se,
               sig = abs(th/se) > tc, n_eff = n, G = G)
  }), fill = TRUE)
}

# -----------------------------------------------------------------------------
# 2. Long-band conditional response on a restricted sample
# -----------------------------------------------------------------------------

band_on <- function(dat, tag) {
  hs <- LONG_BAND
  w <- copy(dat); yc <- character(0)
  for (h in hs) { nm <- sprintf("y_%02d", h); w <- make_y(w, "infl_z", h, nm); yc <- c(yc, nm) }
  xall <- c(CONTROLS_LEAN, "shock_pos", INTER,
            own_lags("infl_z"), unlist(lapply(INTER, own_lags)))
  w <- na.omit(w[, c("iso3","date","cl", yc, xall), with = FALSE])
  ntr <- sapply(INTER, function(v) sum(w[[v]] < 0))
  cat(sprintf("  [%s] long band n=%d treated %s\n", tag, nrow(w), paste(ntr, collapse = "/")))
  if (nrow(w) < MIN_OBS || any(ntr < MIN_TREAT)) {
    warning(sprintf("[%s] band skipped", tag), immediate. = TRUE); return(NULL)
  }
  n <- nrow(w); G <- uniqueN(w$cl)
  w <- demean(w, c(yc, xall)); cl <- w$cl
  w[, c("iso3","date","cl") := NULL]; df <- as.data.frame(w)

  H <- length(hs); M <- length(INTER)
  TH <- matrix(NA_real_, H, M); IFs <- matrix(0, n, M)
  for (i in seq_along(hs)) { cat(sprintf("    h=%2d\n", hs[i])); flush.console()
    for (j in seq_along(INTER)) {
      v <- INTER[j]
      xc <- c(CONTROLS_LEAN, "shock_pos", setdiff(INTER, v),
              own_lags("infl_z"), unlist(lapply(INTER, own_lags)))
      o <- fit_plr(df, yc[i], v, xc, cl)
      TH[i, j] <- o$coef[[1]]; IFs[, j] <- IFs[, j] + influence(o)
    } }
  V <- clustV(IFs / H, cl); tb <- colMeans(TH); se <- sqrt(diag(V)); tc <- qt(0.975, G-1)

  est <- data.table(sample = tag, term = names(MARGINS), coef = tb, se_cl = se,
                    t = tb/se, p = 2*pt(-abs(tb/se), df = G-1),
                    lo = tb - tc*se, hi = tb + tc*se, sig = abs(tb/se) > tc,
                    n_eff = n, G = G, n_treated = as.integer(ntr), kind = "estimate")
  ct <- rbindlist(lapply(combn(M, 2, simplify = FALSE), function(pp) {
    i <- pp[1]; j <- pp[2]; dif <- tb[i] - tb[j]
    s2 <- sqrt(V[i,i] + V[j,j] - 2*V[i,j])
    data.table(sample = tag, term = paste(names(MARGINS)[i], "-", names(MARGINS)[j]),
               coef = dif, se_cl = s2, t = dif/s2, p = 2*pt(-abs(dif/s2), df = G-1),
               lo = dif - tc*s2, hi = dif + tc*s2, sig = abs(dif/s2) > tc,
               n_eff = n, G = G, n_treated = NA_integer_, kind = "contrast")
  }))
  R <- rbind(c(1,-1,0), c(0,1,-1)); Rb <- R %*% tb
  W <- as.numeric(t(Rb) %*% solve(R %*% V %*% t(R)) %*% Rb)
  jt <- data.table(sample = tag, term = "JOINT equality", coef = W, se_cl = NA_real_,
                   t = NA_real_, p = 1 - pchisq(W, 2), lo = NA_real_, hi = NA_real_,
                   sig = (1 - pchisq(W, 2)) < .05, n_eff = n, G = G,
                   n_treated = NA_integer_, kind = "joint")
  rbind(est, ct, jt)
}

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------

D0 <- prep()

base_all <- rbindlist(lapply(names(DROP), function(tag) {
  cached(sprintf("rob_baseline_%s.csv", tag), {
    dd <- D0[!iso3 %in% DROP[[tag]]]
    dd[, cl := .GRP, by = iso3]
    baseline_on(dd, tag)
  })
}), fill = TRUE)

band_all <- rbindlist(lapply(names(DROP), function(tag) {
  cached(sprintf("rob_band_%s.csv", tag), {
    dd <- D0[!iso3 %in% DROP[[tag]]]
    dd[, cl := .GRP, by = iso3]
    band_on(dd, tag)
  })
}), fill = TRUE)

# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------

stars <- function(p) ifelse(is.na(p), "", ifelse(p < .01, "***",
                     ifelse(p < .05, "**", ifelse(p < .10, "*", ""))))

sink("output/tables/ROBUSTNESS_pricesetters.txt")
cat("ROBUSTNESS: excluding dominant commodity price-setters\n")
cat("generated", format(Sys.time()), "\n")
cat(strrep("=", 76), "\n\n")
cat("Exclusion sets\n")
for (k in names(DROP))
  cat(sprintf("  %-13s %s\n", k,
      if (length(DROP[[k]])) paste(DROP[[k]], collapse = ", ") else "(none)"))

cat("\n\nPANEL A. Baseline inflation response, selected horizons\n")
pa <- dcast(base_all[h %in% c(0,3,4,5,6,7,12,24)], h ~ sample,
            value.var = c("coef","t"))
print(pa[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
cat("\nsignificant horizons (5%, h >= 0) by sample:\n")
print(base_all[h >= 0 & sig == TRUE, .(horizons = paste(h, collapse = ", ")), by = sample])
cat("\nclusters and sample size at h = 6:\n")
print(base_all[h == 6, .(sample, n_eff, G)])

cat("\n\nPANEL B. Long band (h = 12-24), margin-conditional response\n")
pb <- band_all[kind == "estimate"]
print(pb[, .(sample, term, coef = round(coef, 4), se = round(se_cl, 4),
             t = round(t, 2), p = round(p, 4), sig = stars(p),
             n_eff, G, n_treated)])

cat("\n\nPANEL C. Contrasts and joint test, long band\n")
pc <- band_all[kind %in% c("contrast","joint")]
print(pc[, .(sample, term, value = round(coef, 4), se = round(se_cl, 4),
             t = round(t, 2), p = round(p, 4), sig = stars(p))])

cat("\n\nHOW TO READ THIS\n")
cat("The conservative exclusion is the one that belongs in the paper. If the\n")
cat("h = 3-6 block and the long-band exchange-rate coefficient survive there,\n")
cat("the results are not driven by countries large enough to move world prices.\n")
cat("The aggressive exclusion drops a fifth of the panel and loses power by\n")
cat("construction; treat a loss of significance there as uninformative unless\n")
cat("the point estimates also move substantially.\n")
cat("\nNote what this check can and cannot establish. The Gruss-Kebhaj weights\n")
cat("are fixed at pre-sample shares, so endogeneity requires a path from current\n")
cat("domestic conditions to world prices. Excluding large producers weakens that\n")
cat("path but does not close it. Report the check as partial.\n")
sink()

cat(readLines("output/tables/ROBUSTNESS_pricesetters.txt"), sep = "\n")
cat("\n\nWritten: out/ROBUSTNESS_pricesetters.txt and out/rob_*.csv\n")
