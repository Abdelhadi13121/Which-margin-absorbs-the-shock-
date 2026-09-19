# Which margin absorbs the shock?

Replication package for **"Which margin absorbs the shock? Commodity terms of trade,
external adjustment and inflation in resource-rich economies."**

Abdelhadi Benghalem, Maghnia University Center, Algeria
ORCID [0000-0003-4057-7448](https://orcid.org/0000-0003-4057-7448)

---

## What the paper does

A commodity exporter facing a deterioration in its terms of trade can absorb the shock
along three external margins: compress imports, run down reserves, or let the currency
depreciate. The paper asks whether that choice determines the inflationary consequences,
using a monthly panel of 49 resource-rich economies over 1990–2026 and double machine
learning applied to panel local projections.

Three results:

1. **Transmission is symmetric.** A direct test of the null — on the *difference*
   between the responses to the two partial sums, with the covariance between the
   estimators accounted for — fails to reject at every horizon.
2. **The margins move in sequence.** Reserves are drawn down between 11 and 15 months
   after a bust, import compression arrives from 15 months, the exchange rate moves last.
3. **The inflation response depends on which margin absorbs.** Over the second year,
   exchange-rate absorbers experience inflation higher by 0.042 within-country standard
   deviations; import- and reserve-absorbers show no response. The joint null of equality
   is rejected (W = 8.45, p = 0.015), and the differential *strengthens* when countries
   large enough to influence world commodity prices are excluded.

---

## Repository layout

```
.
├── code/
│   ├── 00_run_all.R                  master script — runs the pipeline end to end
│   ├── 01_build_master.py            builds the analysis panel from raw IMF extracts
│   ├── 02_replicate_empirical.R      every table and figure in the paper
│   ├── 03_robustness_pricesetters.R  dominant price-setter exclusion (Table 7)
│   ├── 04_robustness_learner_composition.R  learner and composition checks (Table 8)
│   └── 05_fix_twfe_band.R            corrected two-way FE band average
├── data/
│   ├── raw/                          raw IMF extracts (not redistributed — see below)
│   └── derived/
│       ├── view_dml_lp.csv           the analysis panel (21,560 × 33) — ships with the repo
│       ├── dim_variable.csv          data dictionary
│       └── dim_country.csv           country dimension table
├── output/
│   ├── tables/                       ALL_TABLES.txt and tab*.csv
│   ├── figures/                      Fig1–Fig6, 600 dpi
│   └── logs/                         coverage report and validation log
└── docs/
    ├── data_dictionary.md            variable definitions and construction
    └── results_summary.md            headline numbers, for checking a replication
```

---

## Quick start

The analysis panel ships with the repository, so you can reproduce every table and
figure without touching the raw data.

```r
# from the repository root
setwd("path/to/Which-margin-absorbs-the-shock-")
source("code/00_run_all.R")
```

`00_run_all.R` checks your package versions, runs the estimation pipeline, and writes
to `output/`. Expect **several hours on one core** for a cold run. Every stage caches its
results to `output/tables/raw_*.csv`; a second run rebuilds all tables and figures from
cache **in seconds**. Set `FORCE <- TRUE` in a script, or delete the relevant
`raw_*.csv`, to force re-estimation.

To run the robustness exercise separately:

```r
source("code/03_robustness_pricesetters.R")
```

---

## Requirements

R 4.5.0 or later. Python 3.10+ only if rebuilding the panel from raw extracts.

```r
install.packages(c("data.table", "DoubleML", "mlr3", "mlr3learners",
                   "ranger", "xgboost", "fixest", "ggplot2"))
```

Versions used in the published results:

| Package | Version |
|---|---|
| R | 4.5.0 (2025-04-11 ucrt) |
| data.table | 1.17.x |
| DoubleML | 1.0.x |
| mlr3 / mlr3learners | 1.x |
| ranger | 0.17.x |
| xgboost | 1.7.x |
| fixest | 0.12.x |
| ggplot2 | 3.5.x |

Results are seed-controlled (`SEED <- 20260914L`) and reproduce exactly on the same
package versions. Minor numerical differences across `ranger` versions are possible
because the random forest's internal RNG has changed between releases; the reported
coefficients are stable to three decimal places across the versions we have tested.

---

## Data

**Ships with the repository.** `data/derived/view_dml_lp.csv` — 21,560 country-month
observations, 49 countries, January 1990 to August 2026, 33 variables. This is the
only input `02_replicate_empirical.R` needs.

**Not redistributed.** The raw extracts underlying it come from the IMF
International Financial Statistics, the IMF Direction of Trade Statistics, and the IMF
Primary Commodity Price System, which are available to subscribers from the IMF and are
not redistributable here. `code/01_build_master.py` documents the construction in full
and will rebuild the panel from those extracts if you place them in `data/raw/`.

The commodity terms-of-trade index follows the fixed-weight construction of
Gruss and Kebhaj (2019, IMF WP 19/21).

Two measurement points are documented in `docs/data_dictionary.md` because they affect
comparability with earlier work:

- **Reserve adequacy** must be computed against the *monthly* import flow, not an
  annualised figure. The distinction changes the median buffer in this panel from
  70.3 to 6.06 months, and the low-buffer share from 2.4% to 23.4%.
- **Treatment and outcome must share a scale.** Mixing log differences with percentage
  log differences inflates the estimated elasticity by a factor of 100 without affecting
  signs, *t*-statistics, or the shape of the impulse response — which makes the error
  easy to miss.

---

## Method

Partially linear model (Robinson 1988) estimated by double machine learning
(Chernozhukov et al. 2018), applied to panel local projections (Jordà 2005).

- **Nuisance learners:** random forest (500 trees, max depth 8, min node size 10) and,
  for the baseline, gradient boosting (300 rounds, eta 0.05, max depth 5).
- **Cross-fitting:** 5 folds drawn **over countries, not over rows**. Row-level splitting
  would place observations from the same country in both training and evaluation samples,
  and serial dependence would leak across folds.
- **Inference:** computed directly from the orthogonal score, clustered at country level,
  with a G/(G−1) finite-cluster correction and critical values from t with G−1 degrees of
  freedom. At G ≈ 46 the two-sided 5% value is 2.013, not 1.96. Clustering roughly doubles
  the standard errors relative to the i.i.d. default.
- **Band averages:** the estimation sample is held fixed within each horizon band so the
  influence functions are conformable and the cross-horizon covariance enters the variance
  of the average directly.

---

## Reproducing specific results

| Paper object | Script | Output file |
|---|---|---|
| Table 1 (descriptives) | `02_replicate_empirical.R` | `output/tables/tab1_descriptives.csv` |
| Table 2 (baseline) | `02_replicate_empirical.R` | `output/tables/tab2_baseline.csv` |
| Table 3 (symmetry) | `02_replicate_empirical.R` | `output/tables/tab3_symmetry.csv` |
| Table 4 (margins) | `02_replicate_empirical.R` | `output/tables/tab4_margins.csv` |
| Tables 5–6 (bands, contrasts) | `02_replicate_empirical.R` | `output/tables/tab5_bands.csv`, `tab6_contrasts.csv` |
| Table 7 (price-setter exclusion) | `03_robustness_pricesetters.R` | `output/tables/ROBUSTNESS_pricesetters.txt` |
| Table 8 (learner, composition) | `04_robustness_learner_composition.R` | `output/tables/ROBUSTNESS_learner_composition.txt` |
| Table 8, FE column | `05_fix_twfe_band.R` | `output/tables/CORRECTION_twfe_band.txt` |
| Appendix A1–A2 | `02_replicate_empirical.R` | `output/tables/tabA1_countries.csv`, `tabA2_placebo.csv` |
| Figures 1–6 | `02_replicate_empirical.R` | `output/figures/Fig*.tiff` and `.png` |

Headline numbers for checking a replication are in `docs/results_summary.md`.

---

## Known limitations

Stated here as they are in the paper, because a replication package should not read
more confidently than the manuscript:

- The shock is a **constructed fixed-weight index**, not a structurally identified series.
  Excluding dominant producers weakens the path from domestic conditions to world prices
  without closing it.
- The **average response is not robust** to that exclusion, although the mechanism is.
  Both are reported.
- The margin classification **assigns countries, not episodes**. The exchange-rate group
  consists disproportionately of economies with weaker nominal anchors; standardising
  inflation within country removes the difference in scale but does not fully separate a
  channel effect from a composition effect.
- The margin differential is identified at **horizons beyond 12 months only**. At the
  horizons where the average response peaks, the joint test across margins is far from
  significance (p = 0.291).

---

## Citation

If you use this code or data, please cite the paper. See `CITATION.cff`, or:

> Benghalem, A. (2026). Which margin absorbs the shock? Commodity terms of trade,
> external adjustment and inflation in resource-rich economies. *Swiss Journal of
> Economics and Statistics* (under review).

---

## Licence

Code is released under the MIT Licence (see `LICENSE`). The derived panel in
`data/derived/` is provided for replication; the underlying IMF source data remain
subject to the IMF's own terms of use.

---

## Contact

Questions, or a replication that does not reproduce: open an issue, or write to
abdelhadi.benghalem@univ-tlemcen.dz.
