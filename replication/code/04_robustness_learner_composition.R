# =============================================================================
# 04_robustness_learner_composition.R
#
# Two referee-anticipating exercises on the paper's CENTRAL result, the
# margin-conditional inflation response over the long band (h = 12-24).
#
#   PANEL A  Estimator dependence.
#            The baseline average response in Table 2 is significant under the
#            random forest and at no horizon under gradient boosting. Nobody has
#            yet checked whether the MARGIN result depends on the learner. This
#            re-estimates the long band under (i) random forest, (ii) gradient
#            boosting, and (iii) a conventional two-way fixed-effects panel LP
#            with the same controls and no machine learning at all.
#
#   PANEL B  Composition versus channel.
#            Countries classified as exchange-rate absorbers may simply be the
#            economies with weaker nominal anchors. This re-estimates the long
#            band adding, one at a time and then jointly, three predetermined
#            country characteristics: peg status, reserve adequacy, and mean
#            historical inflation. If the FX coefficient survives, the result is
#            harder to read as a nominal-regime proxy.
#
# WHAT THIS CAN AND CANNOT SHOW
#   Neither exercise makes the margin comparison causal. The classification is
#   predetermined with respect to the current shock, which rules out reverse
#   causation from the outcome, but it does not rule out that FX-absorbing
#   countries differ in ways not captured by the three controls used here.
#   Report the exercise as narrowing the composition objection, not closing it.
#
# USAGE (from the repository root)
#   setwd("path/to/Which-margin-absorbs-the-shock-")
#   source("code/04_robustness_learner_composition.R")
#
# RUNTIME  Roughly 4-6 hours cold on one core; caches to output/tables/rob2_*.csv
#          and rebuilds in seconds afterwards.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(DoubleML); library(mlr3); library(mlr3learners)
  library(fixest)
})
lgr::get_logger("mlr3")$set_threshold("warn")

SEED      <- 20260914L
PATH      <- "data/derived/view_dml_lp.csv"
FORCE     <- FALSE
N_FOLDS   <- 5L
N_REP     <- 1L
MIN_OBS   <- 400L
MIN_TREAT <- 60L
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

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
cached <- function(f, expr) {
  fp <- file.path("output/tables", f)
  if (!FORCE && file.exists(fp)) { message("  [cache] ", f); return(fread(fp)) }
  r <- force(expr); fwrite(r, fp); r
}
own_lags <- function(v) paste0(v, "_L", seq_len(N_LAGS))

# -----------------------------------------------------------------------------
# Data preparation — identical to 02_replicate_empirical.R, plus the three
# composition variables built in Panel B.
# -----------------------------------------------------------------------------

prep <- function() {
  d <- fread(PATH)
  d[, date := as.IDate(date)]
  d[, mo := year(date) * 12L + month(date)]
  d[, cl := .GRP, by = iso3]
  setorder(d, iso3, date)
  stopifnot(nrow(d) == 21560L, uniqueN(d$iso3) == 49L)

  d[, res_pos := fifelse(reserves_exgold_usd > 0, reserves_exgold_usd, NA_real_)]
  d[, d_log_reserves_12m := 100 * (log(res_pos) - log(shift(res_pos, 12))), by = iso3]
  d[, mo_12 := shift(mo, 12), by = iso3]
  d[is.na(mo_12) | (mo - mo_12) != 12L, d_log_reserves_12m := NA_real_]
  d[, c("res_pos", "mo_12") := NULL]

  d[, infl_z := (infl_yoy_w - mean(infl_yoy_w, na.rm = TRUE)) /
                sd(infl_yoy_w, na.rm = TRUE), by = iso3]

  # ---- margin classification, ex ante, leave-one-episode-out ----------------
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

  # ---- PANEL B: three predetermined country characteristics ----------------
  # All are computed on the FULL country history and are therefore fixed within
  # country. They enter as interactions with the negative partial sum, so they
  # are not absorbed by the country fixed effect.
  d[, peg_share  := mean(state_peg, na.rm = TRUE), by = iso3]          # 0-1
  d[, res_adeq   := mean(reserves_months, na.rm = TRUE), by = iso3]    # months
  d[, infl_hist  := mean(infl_yoy_w, na.rm = TRUE), by = iso3]         # pp
  for (v in c("peg_share","res_adeq","infl_hist")) {
    x <- d[[v]]
    d[, (v) := (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)]      # standardised
    d[is.na(get(v)), (v) := 0]                                         # neutral if missing
    d[, (paste0("x_", v)) := shock_neg * get(v)]                       # interaction
  }

  lagv <- unique(c("infl_z", "shock_neg", "shock_pos", INTER,
                   paste0("x_", c("peg_share","res_adeq","infl_hist"))))
  for (v in lagv) for (L in seq_len(N_LAGS))
    d[, paste0(v, "_L", L) := shift(get(v), L), by = iso3]
  d[]
}

demean <- function(dt, cols) {
  for (v in cols) {
    dt[, (v) := get(v) - mean(get(v), na.rm = TRUE), by = iso3]
    dt[, (v) := get(v) - mean(get(v), na.rm = TRUE), by = date]
  }
  dt
}
make_y <- function(w, ycol, h, nm) {
  w[, (nm) := shift(get(ycol), n = -h, type = "lag"), by = iso3]
  w[, mo_h := shift(mo, n = -h, type = "lag"), by = iso3]
  w[is.na(mo_h) | (mo_h - mo) != h, (nm) := NA_real_]
  w[, mo_h := NULL]; w
}

LRN <- list(
  rf  = lrn("regr.ranger",  num.trees = 500, min.node.size = 10,
            max.depth = 8, num.threads = 1),
  xgb = lrn("regr.xgboost", nrounds = 300, eta = 0.05, max_depth = 5,
            subsample = 0.8, colsample_bytree = 0.8, nthread = 1)
)
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
fit_plr <- function(df, ycol, dcol, xcols, cl, learner) {
  mk <- function() double_ml_data_from_data_frame(df, y_col = ycol,
                                                  d_cols = dcol, x_cols = xcols)
  o <- DoubleMLPLR$new(mk(), ml_l = learner$clone(), ml_m = learner$clone(),
                       n_folds = N_FOLDS, n_rep = N_REP,
                       score = "partialling out", draw_sample_splitting = FALSE)
  ok <- tryCatch({ o$set_sample_splitting(cluster_folds(cl, N_FOLDS, N_REP, SEED)); TRUE },
                 error = function(e) FALSE)
  if (!ok) {
    if (!FOLD_FALLBACK) { warning("default folds used", immediate. = TRUE); FOLD_FALLBACK <<- TRUE }
    set.seed(SEED)
    o <- DoubleMLPLR$new(mk(), ml_l = learner$clone(), ml_m = learner$clone(),
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
# Long-band estimator, parameterised by learner and by extra controls
# -----------------------------------------------------------------------------

band_dml <- function(dat, learner, extra = character(0), tag) {
  hs <- LONG_BAND
  w <- copy(dat); yc <- character(0)
  for (h in hs) { nm <- sprintf("y_%02d", h); w <- make_y(w, "infl_z", h, nm); yc <- c(yc, nm) }
  xall <- c(CONTROLS_LEAN, "shock_pos", INTER, extra,
            own_lags("infl_z"), unlist(lapply(INTER, own_lags)),
            unlist(lapply(extra, own_lags)))
  xall <- unique(xall)
  w <- na.omit(w[, c("iso3","date","cl", yc, xall), with = FALSE])
  ntr <- sapply(INTER, function(v) sum(w[[v]] < 0))
  cat(sprintf("  [%s] n=%d treated %s\n", tag, nrow(w), paste(ntr, collapse = "/")))
  if (nrow(w) < MIN_OBS || any(ntr < MIN_TREAT)) {
    warning(sprintf("[%s] skipped: insufficient treated mass", tag), immediate. = TRUE)
    return(NULL)
  }
  n <- nrow(w); G <- uniqueN(w$cl)
  w <- demean(w, c(yc, xall)); cl <- w$cl
  w[, c("iso3","date","cl") := NULL]; df <- as.data.frame(w)

  H <- length(hs); M <- length(INTER)
  TH <- matrix(NA_real_, H, M); IFs <- matrix(0, n, M)
  for (i in seq_along(hs)) { cat(sprintf("    h=%2d\n", hs[i])); flush.console()
    for (j in seq_along(INTER)) {
      v <- INTER[j]
      xc <- setdiff(xall, v)
      o <- fit_plr(df, yc[i], v, xc, cl, learner)
      TH[i, j] <- o$coef[[1]]; IFs[, j] <- IFs[, j] + influence(o)
    } }
  V <- clustV(IFs / H, cl); tb <- colMeans(TH); se <- sqrt(diag(V)); tc <- qt(.975, G-1)

  est <- data.table(spec = tag, term = names(MARGINS), coef = tb, se = se,
                    t = tb/se, p = 2*pt(-abs(tb/se), df = G-1),
                    lo = tb - tc*se, hi = tb + tc*se, n = n, G = G, kind = "estimate")
  ct <- rbindlist(lapply(combn(M, 2, simplify = FALSE), function(pp) {
    i <- pp[1]; j <- pp[2]; dif <- tb[i] - tb[j]
    s2 <- sqrt(V[i,i] + V[j,j] - 2*V[i,j])
    data.table(spec = tag, term = paste(names(MARGINS)[i], "-", names(MARGINS)[j]),
               coef = dif, se = s2, t = dif/s2, p = 2*pt(-abs(dif/s2), df = G-1),
               lo = dif - tc*s2, hi = dif + tc*s2, n = n, G = G, kind = "contrast")
  }))
  R <- rbind(c(1,-1,0), c(0,1,-1)); Rb <- R %*% tb
  W <- as.numeric(t(Rb) %*% solve(R %*% V %*% t(R)) %*% Rb)
  jt <- data.table(spec = tag, term = "JOINT equality", coef = W, se = NA_real_,
                   t = NA_real_, p = 1 - pchisq(W, 2), lo = NA_real_, hi = NA_real_,
                   n = n, G = G, kind = "joint")
  rbind(est, ct, jt)
}

# -----------------------------------------------------------------------------
# Conventional two-way FE benchmark on the same band — no machine learning
#
# !! SUPERSEDED. The variance below averages per-horizon covariance matrices as
# !! sum(V_h)/H^2, which assumes independence across horizons. The horizons share
# !! regressors and overlapping 12-month windows, so that expression drops every
# !! cross-horizon covariance term and understates the standard errors by roughly
# !! a factor of H. It produced t = -11.79 on fx in the run of 2026-09-19.
# !! Use code/05_fix_twfe_band.R, which stacks the per-horizon influence
# !! functions as the DML band estimator does. This function is retained only so
# !! the earlier output remains reproducible; do not report its inference.

band_fe <- function(dat, tag = "twfe") {
  hs <- LONG_BAND
  w <- copy(dat); yc <- character(0)
  for (h in hs) { nm <- sprintf("y_%02d", h); w <- make_y(w, "infl_z", h, nm); yc <- c(yc, nm) }
  rhs <- paste(c(INTER, "shock_pos", CONTROLS_LEAN,
                 own_lags("infl_z"), unlist(lapply(INTER, own_lags))), collapse = " + ")
  est <- lapply(yc, function(y) {
    f <- as.formula(paste(y, "~", rhs, "| iso3 + date"))
    m <- feols(f, data = w, cluster = ~iso3)
    list(b = coef(m)[INTER], V = vcov(m)[INTER, INTER, drop = FALSE], G = m$fixef_sizes[["iso3"]])
  })
  H <- length(est); M <- length(INTER)
  tb <- rowMeans(sapply(est, function(e) e$b))
  # INVALID: assumes independence across horizons. See the note above.
  V  <- Reduce(`+`, lapply(est, `[[`, "V")) / H^2
  G  <- est[[1]]$G; se <- sqrt(diag(V)); tc <- qt(.975, G-1)
  out <- data.table(spec = tag, term = names(MARGINS), coef = tb, se = se,
                    t = tb/se, p = 2*pt(-abs(tb/se), df = G-1),
                    lo = tb - tc*se, hi = tb + tc*se, n = NA_integer_, G = G,
                    kind = "estimate")
  ct <- rbindlist(lapply(combn(M, 2, simplify = FALSE), function(pp) {
    i <- pp[1]; j <- pp[2]; dif <- tb[i] - tb[j]
    s2 <- sqrt(V[i,i] + V[j,j] - 2*V[i,j])
    data.table(spec = tag, term = paste(names(MARGINS)[i], "-", names(MARGINS)[j]),
               coef = dif, se = s2, t = dif/s2, p = 2*pt(-abs(dif/s2), df = G-1),
               lo = dif - tc*s2, hi = dif + tc*s2, n = NA_integer_, G = G, kind = "contrast")
  }))
  rbind(out, ct)
}

# =============================================================================
# Run
# =============================================================================

D <- prep()

message("\nPANEL A: estimator dependence")
A <- rbindlist(list(
  cached("rob2_A_rf.csv",   band_dml(D, LRN$rf,  character(0), "DML random forest")),
  cached("rob2_A_xgb.csv",  band_dml(D, LRN$xgb, character(0), "DML gradient boosting")),
  cached("rob2_A_twfe.csv", band_fe(D, "Two-way FE (no ML)"))
), fill = TRUE)

message("\nPANEL B: composition controls")
CTRL <- list(
  "+ peg status"        = "x_peg_share",
  "+ reserve adequacy"  = "x_res_adeq",
  "+ inflation history" = "x_infl_hist",
  "+ all three"         = c("x_peg_share","x_res_adeq","x_infl_hist")
)
B <- rbindlist(lapply(names(CTRL), function(k)
  cached(sprintf("rob2_B_%s.csv", gsub("[^a-z]", "", tolower(k))),
         band_dml(D, LRN$rf, CTRL[[k]], k))), fill = TRUE)

# -----------------------------------------------------------------------------
# Report
# -----------------------------------------------------------------------------

stars <- function(p) ifelse(is.na(p), "", ifelse(p < .01, "***",
                     ifelse(p < .05, "**", ifelse(p < .10, "*", ""))))
show <- function(D, kinds) D[kind %in% kinds,
  .(spec, term, coef = round(coef, 4), se = round(se, 4), t = round(t, 2),
    p = round(p, 4), sig = stars(p), n, G)]

sink("output/tables/ROBUSTNESS_learner_composition.txt")
cat("ROBUSTNESS: estimator dependence and composition controls\n")
cat("Long band, h = 12-24. Outcome: within-country standardised inflation.\n")
cat("generated", format(Sys.time()), "\n"); cat(strrep("=", 78), "\n")

cat("\nBENCHMARK (paper, Table 5): fx -0.0417, se 0.0193, t -2.17, p 0.036\n")
cat("                            imports 0.0133 | reserves 0.0131\n")
cat("                            joint W = 8.45, p = 0.015\n")

cat("\n\nPANEL A. Does the result depend on the nuisance learner?\n\n")
print(show(A, "estimate"))
cat("\nContrasts:\n"); print(show(A, "contrast"))
cat("\nJoint equality tests:\n")
print(A[kind == "joint", .(spec, W = round(coef, 3), p = round(p, 4), sig = stars(p))])
cat("\nRead: if the fx coefficient keeps its sign, rough magnitude and significance\n")
cat("under gradient boosting and under two-way fixed effects, the central result is\n")
cat("not an artefact of the random forest. The two-way FE column uses no machine\n")
cat("learning at all. NOTE: the standard errors in the two-way FE row are INVALID\n")
cat("— they assume independence across horizons. See code/05_fix_twfe_band.R for\n")
cat("the corrected column; use the point estimate here, not the inference.\n")

cat("\n\nPANEL B. Composition or channel?\n\n")
print(show(B, "estimate"))
cat("\nContrasts:\n"); print(show(B, "contrast"))
cat("\nJoint equality tests:\n")
print(B[kind == "joint", .(spec, W = round(coef, 3), p = round(p, 4), sig = stars(p))])
cat("\nEach specification adds the interaction of the negative partial sum with a\n")
cat("predetermined, country-fixed characteristic: the share of the sample spent on a\n")
cat("peg, mean reserve adequacy in months of imports, and mean historical inflation.\n")
cat("All three are standardised across countries. If the fx coefficient survives, the\n")
cat("differential is harder to read as a proxy for the nominal regime.\n")
cat("\nThis narrows the composition objection; it does not close it. Countries may\n")
cat("differ on dimensions none of these three captures, and the margin comparison\n")
cat("remains a comparison across historically chosen adjustment patterns rather than\n")
cat("an experiment.\n")

cat("\n\nCOUNTRY CHARACTERISTICS BY ASSIGNED MARGIN\n\n")
cc <- unique(D[!is.na(margin_class_used),
               .(iso3, margin_class_used, peg_share, res_adeq, infl_hist)])
print(cc[, .(countries = .N,
             peg = round(mean(peg_share), 2),
             reserves = round(mean(res_adeq), 2),
             infl = round(mean(infl_hist), 2)), by = margin_class_used])
cat("\n(values are standardised; positive = above the cross-country mean)\n")
sink()

cat(readLines("output/tables/ROBUSTNESS_learner_composition.txt"), sep = "\n")
cat("\n\nWritten: output/tables/ROBUSTNESS_learner_composition.txt and rob2_*.csv\n")
