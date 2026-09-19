# =============================================================================
# replicate_empirical.R
#
# Full replication of the empirical section: every table and every figure.
#
# USAGE
#   setwd("path/to/Which-margin-absorbs-the-shock-")
#   source("replicate_empirical.R")
#
# The script is STAGED and CACHED. Each estimation stage writes its results to
# out/*.csv; on a later run it reloads that file instead of re-estimating. So
# the first full run is long (several hours on one core) and every run after
# that rebuilds all tables and figures in seconds. To force a stage to re-run,
# delete its csv from out/, or set FORCE <- TRUE below.
#
# OUTPUTS
#   out/tab1_descriptives.csv     Table 1
#   out/tab2_baseline.csv         Table 2
#   out/tab3_symmetry.csv         Table 3
#   out/tab4_margins.csv          Table 4
#   out/tab5_bands.csv            Table 5
#   out/tab6_contrasts.csv        Table 6
#   out/tabA1_countries.csv       Appendix A1
#   out/tabA2_placebo.csv         Appendix A2
#   out/ALL_TABLES.txt            every table, formatted, in one file
#   fig/fig1..fig6.tiff           600 dpi, greyscale-safe
#   out/raw_*.csv                 full horizon-by-horizon estimates behind each figure
#
# Specifications are identical to the estimation scripts v5 to v9: partially
# linear DML, ranger (and xgboost for the baseline), five-fold cross-fitting
# drawn over countries, country-level clustered inference from the orthogonal
# score, t with G-1 degrees of freedom.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(DoubleML); library(mlr3); library(mlr3learners)
  library(ggplot2); library(fixest)
})
lgr::get_logger("mlr3")$set_threshold("warn")

# ---------------------------------------------------------------- 0. Config
SEED      <- 20260914L
PATH      <- "data/derived/view_dml_lp.csv"
FORCE     <- FALSE          # TRUE re-runs every stage from scratch
N_FOLDS   <- 5L
N_REP_BASE<- 3L             # baseline uses 3 repetitions, other stages 1
N_REP     <- 1L
MIN_OBS   <- 400L
MIN_TREAT <- 80L
N_LAGS    <- 6L
H_MIN     <- -6L
H_MAX     <- 24L
BUST_Q    <- 0.25
Q_CUT     <- 0.75           # large-shock threshold, Section 5.3
H_COND    <- 12L
BANDS     <- list(short = 0:6, medium = 7:11, long = 12:24)

CONTROLS_LEAN <- c("d_log_fx_12m", "fx_vol_12m",
                   "d_log_exports_usd_12m", "d_log_imports_usd_12m")
MARGINS <- c(imports  = "d_log_imports_usd_12m",
             reserves = "d_log_reserves_12m",
             fx       = "d_log_fx_12m")

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

cached <- function(file, expr) {
  fp <- file.path("output/tables", file)
  if (!FORCE && file.exists(fp)) { message("  [cache] ", file); return(fread(fp)) }
  res <- force(expr); fwrite(res, fp); res
}

# ------------------------------------------------------------------ 1. Data
message("Stage 1: data")
d <- fread(PATH)
d[, date := as.IDate(date)]
d[, mo := year(date) * 12L + month(date)]
d[, cl := .GRP, by = iso3]
setorder(d, iso3, date)
stopifnot(nrow(d) == 21560L, uniqueN(d$iso3) == 49L)

# reserves margin
d[, res_pos := fifelse(reserves_exgold_usd > 0, reserves_exgold_usd, NA_real_)]
d[, d_log_reserves_12m := 100 * (log(res_pos) - log(shift(res_pos, 12))), by = iso3]
d[, mo_12 := shift(mo, 12), by = iso3]
d[is.na(mo_12) | (mo - mo_12) != 12L, d_log_reserves_12m := NA_real_]
d[, c("res_pos", "mo_12") := NULL]

# within-country standardised inflation
d[, infl_z := (infl_yoy_w - mean(infl_yoy_w, na.rm = TRUE)) /
              sd(infl_yoy_w, na.rm = TRUE), by = iso3]

# size split, Section 5.3
d[, thr := quantile(abs(d_shock_12m), Q_CUT, na.rm = TRUE)]
d[, big := as.integer(abs(d_shock_12m) >= thr)]
d[, shock_big := d_shock_12m * big]
d[, shock_small := d_shock_12m * (1L - big)]

LAGGABLE <- unique(c("infl_yoy_w", "infl_z", "shock_neg", "shock_pos",
                     "shock_big", "shock_small", unname(MARGINS)))
for (v in LAGGABLE) for (L in seq_len(N_LAGS))
  d[, paste0(v, "_L", L) := shift(get(v), L), by = iso3]
own_lags <- function(v) paste0(v, "_L", seq_len(N_LAGS))

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

# ------------------------------------------------- 2. Estimation machinery
lrn_rf  <- lrn("regr.ranger", num.trees = 500, min.node.size = 10,
               max.depth = 8, num.threads = 1)
lrn_xgb <- lrn("regr.xgboost", nrounds = 300, eta = 0.05, max_depth = 5,
               subsample = 0.8, colsample_bytree = 0.8, nthread = 1)

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
fit_plr <- function(df, ycol, dcol, xcols, cl, learner, nrep) {
  mk <- function() double_ml_data_from_data_frame(df, y_col = ycol,
                                                  d_cols = dcol, x_cols = xcols)
  obj <- DoubleMLPLR$new(mk(), ml_l = learner$clone(), ml_m = learner$clone(),
                         n_folds = N_FOLDS, n_rep = nrep,
                         score = "partialling out", draw_sample_splitting = FALSE)
  ok <- tryCatch({ obj$set_sample_splitting(cluster_folds(cl, N_FOLDS, nrep, SEED)); TRUE },
                 error = function(e) FALSE)
  if (!ok) {
    if (!FOLD_FALLBACK) {
      warning("custom fold assignment rejected; default folds used. Inference is ",
              "still cluster-robust.", immediate. = TRUE); FOLD_FALLBACK <<- TRUE
    }
    set.seed(SEED)
    obj <- DoubleMLPLR$new(mk(), ml_l = learner$clone(), ml_m = learner$clone(),
                           n_folds = N_FOLDS, n_rep = nrep, score = "partialling out")
  }
  obj$fit(); obj
}
influence <- function(obj, k = 1L, r = 1L) {
  ar <- obj$psi_a[, r, k]; br <- obj$psi_b[, r, k]
  (ar * (-mean(br) / mean(ar)) + br) / mean(ar)
}
clustV <- function(IF, cl) {
  n <- nrow(IF); G <- length(unique(cl))
  crossprod(rowsum(IF, cl)) / n^2 * (G / (G - 1))
}
tcrit <- function(G) qt(0.975, G - 1)

# ============================================================ TABLE 1 + FIG 1
message("Stage 2: descriptives")
tab1 <- cached("tab1_descriptives.csv", {
  vs <- c(infl_yoy_w = "Inflation, 12-month (%)",
          d_shock_12m = "ToT shock, 12-month (%)",
          shock_neg = "  Negative part", shock_pos = "  Positive part",
          d_log_fx_12m = "FX depreciation, 12-month (%)",
          d_log_imports_usd_12m = "Import growth, 12-month (%)",
          d_log_reserves_12m = "Reserve growth, 12-month (%)",
          d_log_exports_usd_12m = "Export growth, 12-month (%)")
  rbindlist(lapply(names(vs), function(v) {
    s <- d[[v]]; x <- s[!is.na(s)]
    data.table(Variable = vs[[v]], N = length(x), Mean = mean(x), SD = sd(x),
               P10 = quantile(x, .1), Median = median(x), P90 = quantile(x, .9),
               Coverage = 100 * mean(!is.na(s)))
  }))
})

bust_cut <- quantile(d$d_shock_12m, BUST_Q, na.rm = TRUE)
d[, bust := as.integer(d_shock_12m <= bust_cut)]

# ------------------------------------------------------------ figure theme
if (.Platform$OS.type == "windows")
  try(windowsFonts(Times = windowsFont("Times New Roman")), silent = TRUE)
FAM <- if (.Platform$OS.type == "windows") "Times" else "serif"

th <- theme_classic(base_size = 9, base_family = FAM) +
  theme(panel.grid.major.y = element_line(colour = "grey88", linewidth = .25),
        axis.line = element_line(linewidth = .35),
        legend.position = "bottom", legend.title = element_blank(),
        legend.key.width = unit(18, "pt"), legend.margin = margin(t = -4),
        plot.margin = margin(4, 8, 4, 4))

sv <- function(pl, name, w = 5.4, h = 3.2) {
  ggsave(file.path("output/figures", paste0(name, ".tiff")), pl, width = w, height = h,
         dpi = 600, compression = "lzw", bg = "white")
  ggsave(file.path("output/figures", paste0(name, ".png")), pl, width = w, height = h,
         dpi = 600, bg = "white")
}

sv(
  ggplot(d[!is.na(d_shock_12m) & abs(d_shock_12m) <= 30], aes(d_shock_12m)) +
    geom_histogram(bins = 70, fill = "grey78", colour = "grey35", linewidth = .15) +
    geom_vline(xintercept = bust_cut, linetype = "dashed", linewidth = .6) +
    geom_vline(xintercept = 0, linewidth = .35) +
    annotate("text", x = bust_cut - 11, y = Inf, vjust = 2, size = 2.6, family = FAM,
             label = sprintf("bust threshold\n%.2f%%", bust_cut)) +
    labs(x = "12-month change in the commodity terms of trade (%)",
         y = "Country-months") + th + theme(legend.position = "none"),
  "fig1", h = 2.8)

# ============================================================ TABLE 2 + FIG 2
message("Stage 3: baseline response of inflation")
base_raw <- cached("raw_baseline.csv", {
  rbindlist(lapply(H_MIN:H_MAX, function(h) {
    rbindlist(lapply(list(rf = lrn_rf, xgb = lrn_xgb), function(L) {
      cat(sprintf("  baseline h=%3d %s\n", h, L$id)); flush.console()
      w <- make_y(copy(d), "infl_yoy_w", h)
      xc <- c(CONTROLS_LEAN, "shock_pos",
              own_lags("infl_yoy_w"), own_lags("shock_neg"), own_lags("shock_pos"))
      w <- na.omit(w[, c("iso3","date","cl","y_h","shock_neg", xc), with = FALSE])
      if (nrow(w) < MIN_OBS) return(NULL)
      n <- nrow(w); G <- uniqueN(w$cl); nc <- uniqueN(w$iso3)
      w <- demean(w, c("y_h", "shock_neg", xc)); cl <- w$cl
      w[, c("iso3","date","cl") := NULL]
      o  <- fit_plr(as.data.frame(w), "y_h", "shock_neg", xc, cl, L, N_REP_BASE)
      IF <- sapply(seq_len(N_REP_BASE), function(r) influence(o, 1L, r))
      se <- median(sapply(seq_len(N_REP_BASE), function(r)
             sqrt(clustV(cbind(IF[, r]), cl)[1, 1])))
      th_ <- o$coef[[1]]
      data.table(h = h, learner = L$id, coef = th_, se_cl = se, t = th_/se,
                 sig = abs(th_/se) > tcrit(G), n_eff = n, n_country = nc, G = G)
    }), fill = TRUE)
  }), fill = TRUE)
})

tab2 <- base_raw[h %in% c(0,3,4,5,6,7,12,24)]
tab2 <- dcast(tab2, h ~ learner, value.var = c("coef","se_cl","t","sig"))
fwrite(tab2, "output/tables/tab2_baseline.csv")

tabA2 <- base_raw[h < 0, .(h, learner, coef, se_cl, t, sig)]
fwrite(tabA2, "output/tables/tabA2_placebo.csv")

b <- base_raw[learner == "regr.ranger" & h >= 0]
bx <- base_raw[learner == "regr.xgboost" & h >= 0]
G0 <- b$G[1]; tc0 <- tcrit(G0)
sv(
  ggplot(b, aes(h, coef)) +
    annotate("rect", xmin = 3, xmax = 7, ymin = -Inf, ymax = Inf, fill = "grey93") +
    geom_ribbon(aes(ymin = coef - tc0*se_cl, ymax = coef + tc0*se_cl),
                fill = "grey80", alpha = .85) +
    geom_hline(yintercept = 0, linewidth = .35) +
    geom_line(aes(linetype = "Random forest"), linewidth = .7) +
    geom_line(data = bx, aes(h, coef, linetype = "Gradient boosting"), linewidth = .5) +
    geom_point(data = b[sig == TRUE], size = 1.5) +
    scale_linetype_manual(values = c("Random forest" = "solid",
                                     "Gradient boosting" = "dashed")) +
    scale_x_continuous(breaks = seq(0, 24, 3)) +
    labs(x = "Months after the shock", y = "Response of inflation (pp)") + th,
  "fig2")

# ============================================================ TABLE 3 + FIG 3
message("Stage 4: symmetry test")
asym_raw <- cached("raw_symmetry.csv", {
  rbindlist(lapply(H_MIN:12L, function(h) {
    cat(sprintf("  symmetry h=%3d\n", h)); flush.console()
    w <- make_y(copy(d), "infl_yoy_w", h)
    xall <- c(CONTROLS_LEAN, "shock_neg", "shock_pos",
              own_lags("infl_yoy_w"), own_lags("shock_neg"), own_lags("shock_pos"))
    w <- na.omit(w[, c("iso3","date","cl","y_h", xall), with = FALSE])
    if (nrow(w) < MIN_OBS) return(NULL)
    G <- uniqueN(w$cl)
    w <- demean(w, c("y_h", xall)); cl <- w$cl
    w[, c("iso3","date","cl") := NULL]; df <- as.data.frame(w)
    one <- function(tr, ot) {
      xc <- c(CONTROLS_LEAN, ot, own_lags("infl_yoy_w"),
              own_lags("shock_neg"), own_lags("shock_pos"))
      o <- fit_plr(df, "y_h", tr, xc, cl, lrn_rf, N_REP)
      list(th = o$coef[[1]], inf = influence(o))
    }
    a <- one("shock_neg", "shock_pos"); bb <- one("shock_pos", "shock_neg")
    V <- clustV(cbind(a$inf, bb$inf), cl)
    dif <- a$th - bb$th; sed <- sqrt(V[1,1] + V[2,2] - 2*V[1,2])
    data.table(h = h, beta_neg = a$th, se_neg = sqrt(V[1,1]), t_neg = a$th/sqrt(V[1,1]),
               beta_pos = bb$th, se_pos = sqrt(V[2,2]), t_pos = bb$th/sqrt(V[2,2]),
               diff = dif, se_diff = sed, t_diff = dif/sed,
               p_diff = 2*pt(-abs(dif/sed), df = G-1),
               sig_diff = abs(dif/sed) > tcrit(G), mde = 2.80*sed, G = G)
  }), fill = TRUE)
})
tab3 <- asym_raw[h %in% c(3,4,5,6,7,12),
                 .(h, beta_neg, t_neg, beta_pos, t_pos, diff, p_diff)]
fwrite(tab3, "output/tables/tab3_symmetry.csv")

a3 <- melt(asym_raw[, .(h, `Negative shock` = beta_neg, `Positive shock` = beta_pos)],
           id.vars = "h", variable.name = "part", value.name = "coef")
sv(
  ggplot(a3, aes(h, coef, linetype = part)) +
    geom_ribbon(data = asym_raw, inherit.aes = FALSE,
                aes(h, ymin = beta_neg - tcrit(G)*se_neg,
                        ymax = beta_neg + tcrit(G)*se_neg),
                fill = "grey82", alpha = .8) +
    geom_hline(yintercept = 0, linewidth = .35) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = .4) +
    geom_line(linewidth = .65) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    scale_x_continuous(breaks = seq(-6, 12, 3)) +
    labs(x = "Months after the shock", y = "Response of inflation (pp)") + th,
  "fig3")

# ===================================================== SECTION 5.3 magnitude
message("Stage 5: size nonlinearity")
mag_raw <- cached("raw_magnitude.csv", {
  rbindlist(lapply(H_MIN:12L, function(h) {
    cat(sprintf("  magnitude h=%3d\n", h)); flush.console()
    w <- make_y(copy(d), "infl_yoy_w", h)
    xall <- c(CONTROLS_LEAN, "shock_big", "shock_small",
              own_lags("infl_yoy_w"), own_lags("shock_big"), own_lags("shock_small"))
    w <- na.omit(w[, c("iso3","date","cl","y_h", xall), with = FALSE])
    if (nrow(w) < MIN_OBS) return(NULL)
    G <- uniqueN(w$cl)
    w <- demean(w, c("y_h", xall)); cl <- w$cl
    w[, c("iso3","date","cl") := NULL]; df <- as.data.frame(w)
    one <- function(tr, ot) {
      xc <- c(CONTROLS_LEAN, ot, own_lags("infl_yoy_w"),
              own_lags("shock_big"), own_lags("shock_small"))
      o <- fit_plr(df, "y_h", tr, xc, cl, lrn_rf, N_REP)
      list(th = o$coef[[1]], inf = influence(o))
    }
    a <- one("shock_big", "shock_small"); bb <- one("shock_small", "shock_big")
    V <- clustV(cbind(a$inf, bb$inf), cl)
    dif <- a$th - bb$th; sed <- sqrt(V[1,1] + V[2,2] - 2*V[1,2])
    data.table(h = h, beta_big = a$th, beta_small = bb$th, diff = dif,
               se_diff = sed, t_diff = dif/sed,
               p_diff = 2*pt(-abs(dif/sed), df = G-1), mde = 2.80*sed)
  }), fill = TRUE)
})

# ============================================================ TABLE 4 + FIG 4
message("Stage 6: adjustment margins")
marg_raw <- cached("raw_margins.csv", {
  rbindlist(lapply(names(MARGINS), function(m) {
    yc <- MARGINS[[m]]
    rbindlist(lapply(H_MIN:H_MAX, function(h) {
      cat(sprintf("  margin %-9s h=%3d\n", m, h)); flush.console()
      w <- make_y(copy(d), yc, h)
      ctrl <- setdiff(CONTROLS_LEAN, yc)
      xc <- c(ctrl, "shock_pos", own_lags(yc), own_lags("shock_neg"))
      w <- na.omit(w[, c("iso3","date","cl","y_h","shock_neg", xc), with = FALSE])
      if (nrow(w) < MIN_OBS) return(NULL)
      G <- uniqueN(w$cl); nc <- uniqueN(w$iso3); n <- nrow(w)
      w <- demean(w, c("y_h","shock_neg", xc)); cl <- w$cl
      w[, c("iso3","date","cl") := NULL]
      o <- fit_plr(as.data.frame(w), "y_h", "shock_neg", xc, cl, lrn_rf, N_REP)
      th_ <- o$coef[[1]]; se <- sqrt(clustV(cbind(influence(o)), cl)[1,1])
      data.table(margin = m, h = h, coef = th_, se_cl = se, t = th_/se,
                 sig = abs(th_/se) > tcrit(G), n_eff = n, n_country = nc)
    }), fill = TRUE)
  }), fill = TRUE)
})

tab4 <- dcast(marg_raw[h %in% c(6,12,18)], margin ~ h, value.var = "coef")
setnames(tab4, c("6","12","18"), c("h=6","h=12","h=18"))
tab4[, `Significant horizons` := marg_raw[h >= 0 & sig == TRUE,
      .(s = paste(h, collapse = ", ")), by = margin][match(tab4$margin, margin), s]]
fwrite(tab4, "output/tables/tab4_margins.csv")

mm <- copy(marg_raw[h >= 0])
mm[, Margin := factor(margin, levels = c("reserves","imports","fx"),
                      labels = c("Reserves","Imports","Exchange rate"))]
sv(
  ggplot(mm, aes(h, coef, linetype = Margin)) +
    geom_hline(yintercept = 0, linewidth = .35) +
    geom_line(linewidth = .6) +
    geom_point(data = mm[sig == TRUE], aes(shape = Margin), size = 1.6, show.legend = FALSE) +
    scale_linetype_manual(values = c("solid","longdash","dotted")) +
    scale_x_continuous(breaks = seq(0, 24, 3)) +
    labs(x = "Months after the shock", y = "Response of margin (%)") + th,
  "fig4")

# ================================================= classification + FIG 5/6
message("Stage 7: margin classification")
for (m in names(MARGINS)) {
  yc <- MARGINS[[m]]
  d[, (paste0("fwd_", m)) := shift(get(yc), n = -H_COND, type = "lag"), by = iso3]
  d[, mo_f := shift(mo, n = -H_COND, type = "lag"), by = iso3]
  d[is.na(mo_f) | (mo_f - mo) != H_COND, (paste0("fwd_", m)) := NA_real_]
}
d[, mo_f := NULL]
ABS_SIGN <- c(imports = -1, reserves = -1, fx = +1)
ep <- d[bust == 1 & !is.na(fwd_imports) & !is.na(fwd_reserves) & !is.na(fwd_fx)]
for (m in names(MARGINS)) {
  v <- ep[[paste0("fwd_", m)]] * ABS_SIGN[[m]]
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
fwrite(ep[, .(iso3, date, d_shock_12m, z_imports, z_reserves, z_fx, margin_class_used)],
       "output/tables/raw_classification.csv")

d <- merge(d, ep[, .(iso3, date, margin_class_used)], by = c("iso3","date"), all.x = TRUE)
for (m in names(MARGINS))
  d[, (paste0("shk_", m)) := shock_neg *
      as.integer(!is.na(margin_class_used) & margin_class_used == m)]
INTER <- paste0("shk_", names(MARGINS))
for (v in INTER) for (L in seq_len(N_LAGS))
  d[, paste0(v, "_L", L) := shift(get(v), L), by = iso3]

message("Stage 8: margin-conditional inflation, horizon by horizon")
cond_raw <- cached("raw_conditional.csv", {
  rbindlist(lapply(0:H_MAX, function(h) {
    cat(sprintf("  conditional h=%3d\n", h)); flush.console()
    w <- make_y(copy(d), "infl_z", h)
    xall <- c(CONTROLS_LEAN, "shock_pos", INTER,
              own_lags("infl_z"), unlist(lapply(INTER, own_lags)))
    w <- na.omit(w[, c("iso3","date","cl","y_h", xall), with = FALSE])
    ntr <- sapply(INTER, function(v) sum(w[[v]] < 0))
    if (nrow(w) < MIN_OBS || any(ntr < MIN_TREAT)) return(NULL)
    G <- uniqueN(w$cl)
    w <- demean(w, c("y_h", xall)); cl <- w$cl
    w[, c("iso3","date","cl") := NULL]; df <- as.data.frame(w)
    fits <- lapply(INTER, function(v) {
      xc <- c(CONTROLS_LEAN, "shock_pos", setdiff(INTER, v),
              own_lags("infl_z"), unlist(lapply(INTER, own_lags)))
      o <- fit_plr(df, "y_h", v, xc, cl, lrn_rf, N_REP)
      list(th = o$coef[[1]], inf = influence(o))
    })
    V <- clustV(do.call(cbind, lapply(fits, `[[`, "inf")), cl)
    th_ <- sapply(fits, `[[`, "th"); se <- sqrt(diag(V))
    data.table(h = h, margin = names(MARGINS), coef = th_, se_cl = se,
               t = th_/se, sig = abs(th_/se) > tcrit(G), n_eff = nrow(df))
  }), fill = TRUE)
})

cc <- copy(cond_raw)
cc[, Margin := factor(margin, levels = c("fx","imports","reserves"),
     labels = c("Exchange rate absorbs","Imports absorb","Reserves absorb"))]
sv(
  ggplot(cc, aes(h, coef, linetype = Margin)) +
    annotate("rect", xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf, fill = "grey94") +
    geom_hline(yintercept = 0, linewidth = .35) +
    geom_line(linewidth = .65) +
    scale_linetype_manual(values = c("solid","longdash","dotted")) +
    scale_x_continuous(breaks = seq(0, 24, 3)) +
    labs(x = "Months after the shock",
         y = "Response of inflation (within-country SD)") + th,
  "fig5")

# ====================================================== TABLES 5, 6 + FIG 6
message("Stage 9: band averages")
band_out <- cached("raw_bands.csv", {
  rbindlist(lapply(names(BANDS), function(bn) {
    hs <- BANDS[[bn]]
    w <- copy(d); ycols <- character(0)
    for (h in hs) { nm <- sprintf("y_%02d", h); w <- make_y(w, "infl_z", h, nm)
                    ycols <- c(ycols, nm) }
    xall <- c(CONTROLS_LEAN, "shock_pos", INTER,
              own_lags("infl_z"), unlist(lapply(INTER, own_lags)))
    w <- na.omit(w[, c("iso3","date","cl", ycols, xall), with = FALSE])
    ntr <- sapply(INTER, function(v) sum(w[[v]] < 0))
    cat(sprintf("  band %-7s n=%d treated %s\n", bn, nrow(w), paste(ntr, collapse="/")))
    if (nrow(w) < MIN_OBS || any(ntr < MIN_TREAT)) return(NULL)
    n <- nrow(w); G <- uniqueN(w$cl)
    w <- demean(w, c(ycols, xall)); cl <- w$cl
    w[, c("iso3","date","cl") := NULL]; df <- as.data.frame(w)
    H <- length(hs); Mn <- length(INTER)
    TH <- matrix(NA_real_, H, Mn); IFs <- matrix(0, n, Mn)
    for (i in seq_along(hs)) { cat(sprintf("    h=%2d\n", hs[i])); flush.console()
      for (j in seq_along(INTER)) {
        v <- INTER[j]
        xc <- c(CONTROLS_LEAN, "shock_pos", setdiff(INTER, v),
                own_lags("infl_z"), unlist(lapply(INTER, own_lags)))
        o <- fit_plr(df, ycols[i], v, xc, cl, lrn_rf, N_REP)
        TH[i, j] <- o$coef[[1]]; IFs[, j] <- IFs[, j] + influence(o)
      } }
    V <- clustV(IFs / H, cl); tb <- colMeans(TH); se <- sqrt(diag(V)); tc <- tcrit(G)
    est <- data.table(band = bn, h_lo = min(hs), h_hi = max(hs),
                      margin = names(MARGINS), coef = tb, se_cl = se, t = tb/se,
                      p = 2*pt(-abs(tb/se), df = G-1),
                      lo = tb - tc*se, hi = tb + tc*se,
                      sig = abs(tb/se) > tc, n_eff = n, G = G,
                      n_treated = as.integer(ntr), type = "estimate")
    ct <- rbindlist(lapply(combn(Mn, 2, simplify = FALSE), function(pp) {
      i <- pp[1]; j <- pp[2]; dif <- tb[i] - tb[j]
      s2 <- sqrt(V[i,i] + V[j,j] - 2*V[i,j])
      data.table(band = bn, h_lo = min(hs), h_hi = max(hs),
                 margin = paste(names(MARGINS)[i], "-", names(MARGINS)[j]),
                 coef = dif, se_cl = s2, t = dif/s2,
                 p = 2*pt(-abs(dif/s2), df = G-1),
                 lo = dif - tc*s2, hi = dif + tc*s2, sig = abs(dif/s2) > tc,
                 n_eff = n, G = G, n_treated = NA_integer_, type = "contrast")
    }))
    R <- rbind(c(1,-1,0), c(0,1,-1)); Rb <- R %*% tb
    W <- as.numeric(t(Rb) %*% solve(R %*% V %*% t(R)) %*% Rb)
    jt <- data.table(band = bn, h_lo = min(hs), h_hi = max(hs),
                     margin = "JOINT equality", coef = W, se_cl = NA_real_,
                     t = NA_real_, p = 1 - pchisq(W, 2), lo = NA_real_, hi = NA_real_,
                     sig = (1 - pchisq(W, 2)) < .05, n_eff = n, G = G,
                     n_treated = NA_integer_, type = "joint")
    rbind(est, ct, jt)
  }), fill = TRUE)
})

tab5 <- band_out[type == "estimate", .(band, margin, coef, se_cl, t, lo, hi, p, n_eff)]
fwrite(tab5, "output/tables/tab5_bands.csv")
tab6 <- band_out[band == "long" & type %in% c("contrast","joint"),
                 .(contrast = margin, diff = coef, se_cl, t, p)]
fwrite(tab6, "output/tables/tab6_contrasts.csv")

bb6 <- band_out[type == "estimate"]
bb6[, Band := factor(band, levels = c("short","medium","long"),
      labels = c("Short\n(0-6)","Medium\n(7-11)","Long\n(12-24)"))]
bb6[, Margin := factor(margin, levels = c("imports","reserves","fx"),
      labels = c("Imports","Reserves","Exchange rate"))]
sv(
  ggplot(bb6, aes(Band, coef, shape = Margin)) +
    geom_hline(yintercept = 0, linewidth = .35) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = .12, linewidth = .45,
                  position = position_dodge(.55)) +
    geom_point(size = 2, position = position_dodge(.55), fill = "white") +
    scale_shape_manual(values = c(22, 21, 23)) +
    labs(x = NULL, y = "Response of inflation (within-country SD)") + th,
  "fig6", h = 3.0)

# ================================================================ APPENDIX A1
tabA1 <- d[, .(Months = .N,
               First = format(min(date), "%Y-%m"), Last = format(max(date), "%Y-%m"),
               `CPI obs` = sum(!is.na(infl_yoy_w))), by = .(ISO3 = iso3)][order(ISO3)]
fwrite(tabA1, "output/tables/tabA1_countries.csv")

# ============================================================== ALL TABLES
message("Stage 10: writing tables")
stars <- function(p) ifelse(is.na(p), "", ifelse(p < .01, "***",
                     ifelse(p < .05, "**", ifelse(p < .10, "*", ""))))
sink("output/tables/ALL_TABLES.txt")
cat("REPLICATION OUTPUT -- empirical section\n")
cat("generated", format(Sys.time()), "| R", R.version.string, "\n")
cat("fold fallback triggered:", FOLD_FALLBACK, "\n")
cat(strrep("=", 78), "\n\nTABLE 1. Descriptive statistics\n")
print(tab1[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
cat(sprintf("\nbust threshold (25th pct of the shock): %.3f\n", bust_cut))
cat("\n\nTABLE 2. Baseline response of inflation\n"); print(tab2)
cat("\n\nTABLE 3. Symmetry test\n")
print(tab3[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
cat("\n\nSECTION 5.3. Size nonlinearity\n")
print(mag_raw[h %in% 3:7, .(h, beta_big = round(beta_big, 4),
      beta_small = round(beta_small, 4), ratio = round(beta_big/beta_small, 2),
      p = round(p_diff, 3), mde = round(mde, 4))])
cat(sprintf("mean |shock|: large %.2f, small %.2f\n",
            d[big == 1, mean(abs(d_shock_12m), na.rm = TRUE)],
            d[big == 0, mean(abs(d_shock_12m), na.rm = TRUE)]))
cat("\n\nTABLE 4. Adjustment margins\n"); print(tab4)
cat("\n\nTABLE 5. Band-averaged conditional response\n")
print(tab5[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
cat("\n\nTABLE 6. Contrasts, long band\n")
print(tab6[, .(contrast, diff = round(diff, 4), se_cl = round(se_cl, 4),
               t = round(t, 2), p = round(p, 4), sig = stars(p))])
cat("\n\nEPISODES BY ABSORBING MARGIN\n")
print(ep[, .N, by = margin_class_used][order(-N)])
print(ep[!is.na(margin_class_used), .(countries = uniqueN(iso3)), by = margin_class_used])
cat("\n\nAPPENDIX A1. Countries\n"); print(tabA1)
cat("\n\nAPPENDIX A2. Pre-treatment placebo horizons\n")
print(tabA2[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
cat("\n\nTWO-WAY FIXED EFFECTS BENCHMARK (Section 5.1)\n")
print(summary(feols(infl_yoy_w ~ shock_neg + shock_pos + d_log_fx_12m | iso3 + date,
                    data = d, cluster = ~iso3)))
sink()

cat("\nDONE\n")
cat("  tables  -> out/ALL_TABLES.txt and out/tab*.csv\n")
cat("  figures -> fig/fig1..fig6 (.tiff 600 dpi and .png)\n")
cat("  raw estimates -> out/raw_*.csv\n")
cat("\nRe-running now rebuilds every table and figure from cache in seconds.\n")
cat("Set FORCE <- TRUE, or delete out/raw_*.csv, to re-estimate.\n")
