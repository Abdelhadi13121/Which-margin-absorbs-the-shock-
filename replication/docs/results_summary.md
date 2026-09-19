# Results summary

Headline numbers, for checking whether a replication has worked. If your run
reproduces these to three decimal places, the pipeline is behaving. All estimates
use the random forest nuisance learner unless stated otherwise; standard errors
are clustered at country level with critical values from t(G−1).

Seed: `20260914L`. Bust threshold (25th percentile of the shock): **−1.583**.

---

## Table 1 — descriptive statistics

| Variable | N | Mean | SD | P10 | Median | P90 | Coverage |
|---|---|---|---|---|---|---|---|
| Inflation, 12-month (%) | 11316 | 6.66 | 8.38 | 0.41 | 3.97 | 15.29 | 52.5% |
| ToT shock, 12-month (%) | 20513 | 0.41 | 6.93 | −5.63 | 0.18 | 7.08 | 95.1% |
| — negative part | 20513 | −1.85 | 4.32 | −5.63 | 0.00 | 0.00 | 95.1% |
| — positive part | 20513 | 2.25 | 4.58 | 0.00 | 0.18 | 7.08 | 95.1% |
| FX depreciation, 12-month (%) | 19156 | 9.15 | 48.00 | −6.84 | 0.46 | 22.54 | 88.9% |
| Import growth, 12-month (%) | 20508 | 6.90 | 37.44 | −30.89 | 6.97 | 43.37 | 95.1% |
| Reserve growth, 12-month (%) | 17575 | 9.84 | 42.23 | −22.02 | 6.70 | 45.67 | 81.5% |
| Export growth, 12-month (%) | 20508 | 6.25 | 54.62 | −45.16 | 6.80 | 55.60 | 95.1% |

---

## Table 2 — baseline inflation response

Estimation sample 10,231 observations, 47 countries. Critical value 2.013.

| h | RF coef | SE | t | XGB coef | t |
|---|---|---|---|---|---|
| 0 | −0.022 | 0.013 | −1.68 | −0.017 | −1.36 |
| 3 | −0.052 | 0.026 | **−2.05** | −0.024 | −1.02 |
| 4 | −0.056 | 0.028 | **−2.02** | −0.031 | −1.09 |
| 5 | −0.064 | 0.029 | **−2.18** | −0.037 | −1.20 |
| 6 | −0.062 | 0.030 | **−2.06** | −0.039 | −1.25 |
| 7 | −0.064 | 0.033 | −1.91 | −0.040 | −1.16 |
| 12 | −0.048 | 0.052 | −0.92 | −0.046 | −0.76 |
| 24 | −0.091 | 0.065 | −1.40 | −0.113 | −1.39 |

Significant block: **h = 3 to 6** (random forest). Gradient boosting reaches
significance at **no** horizon — this is reported in the paper as a limitation,
not as corroboration.

Two-way fixed-effects benchmark on the same sample: coefficient −0.007,
clustered SE 0.034, 10,489 observations, 47 country and 425 time fixed effects.

---

## Table 3 — symmetry test

| h | β⁻ | t | β⁺ | t | difference | p |
|---|---|---|---|---|---|---|
| 3 | −0.051 | −2.02 | −0.058 | −2.22 | 0.007 | 0.812 |
| 4 | −0.057 | −2.07 | −0.056 | −2.17 | −0.001 | 0.975 |
| 5 | −0.059 | −2.07 | −0.051 | −2.01 | −0.008 | 0.791 |
| 6 | −0.067 | −2.31 | −0.055 | −2.01 | −0.012 | 0.691 |
| 7 | −0.061 | −1.93 | −0.062 | −2.16 | 0.001 | 0.968 |
| 12 | −0.047 | −0.90 | −0.086 | −2.56 | 0.040 | 0.379 |

Symmetry is not rejected at any horizon from −6 to +12. Note that β⁻ is
individually significant at h = 3–6 and β⁺ at h = 3–7.

---

## Section 5.3 — size nonlinearity

Split at the 75th percentile of |shock| (= 5.083). Mean |shock|: 11.60 large,
1.60 small.

| h | β_large | β_small | ratio | p | MDE |
|---|---|---|---|---|---|
| 3 | −0.0488 | −0.0322 | 1.51 | 0.699 | 0.119 |
| 4 | −0.0486 | −0.0382 | 1.27 | 0.834 | 0.138 |
| 5 | −0.0476 | −0.0401 | 1.19 | 0.888 | 0.147 |
| 6 | −0.0475 | −0.0629 | 0.76 | 0.788 | 0.160 |
| 7 | −0.0478 | −0.0824 | 0.58 | 0.588 | 0.177 |

Inconclusive rather than null: minimum detectable differences are 0.119–0.177
against coefficients of roughly 0.05.

---

## Table 4 — adjustment margins

| Margin | h=6 | h=12 | h=18 | Countries | Significant horizons (5%) |
|---|---|---|---|---|---|
| Imports | −0.459 | 0.183 | 0.630 | 48 | 5, 6, 15–18, 22–24 |
| Reserves | 0.730 | 0.968 | 0.261 | 44 | 11–15 |
| Exchange rate | 0.116 | 0.092 | −0.169 | 47–48 | 24 |

---

## Tables 5 and 6 — band-averaged conditional response

Outcome: within-country standardised inflation. Ex ante classification,
leave-one-episode-out. Episode counts: reserves 1769, imports 1390, fx 1091,
unclassified 11. Countries: 16 reserves, 16 imports, 13 fx.

| Band | Margin | Coef | SE | t | p | n | G |
|---|---|---|---|---|---|---|---|
| Short (0–6) | imports | 0.0021 | 0.0063 | 0.33 | 0.744 | 9999 | 47 |
| | reserves | −0.0098 | 0.0116 | −0.85 | 0.402 | 9999 | 47 |
| | fx | −0.0157 | 0.0106 | −1.48 | 0.146 | 9999 | 47 |
| Medium (7–11) | imports | 0.0091 | 0.0112 | 0.81 | 0.420 | 9783 | 46 |
| | reserves | 0.0010 | 0.0133 | 0.07 | 0.942 | 9783 | 46 |
| | fx | −0.0214 | 0.0160 | −1.33 | 0.189 | 9783 | 46 |
| **Long (12–24)** | imports | 0.0133 | 0.0088 | 1.52 | 0.136 | 9208 | 46 |
| | reserves | 0.0131 | 0.0147 | 0.89 | 0.379 | 9208 | 46 |
| | **fx** | **−0.0417** | 0.0193 | **−2.17** | **0.036** | 9208 | 46 |

Contrasts, long band:

| Contrast | difference | SE | t | p |
|---|---|---|---|---|
| imports − fx | 0.0550 | 0.0195 | 2.82 | 0.007 |
| reserves − fx | 0.0548 | 0.0216 | 2.54 | 0.015 |
| imports − reserves | 0.0003 | 0.0154 | 0.02 | 0.986 |
| **Joint equality** | W = 8.448 | | | **0.015** |

---

## Table 7 — dominant price-setter exclusion

Conservative: drop SAU, RUS, CHL. Aggressive: those three plus IRN, IRQ, ARE,
KWT, PER, QAT, DZA.

| | Full | Conservative | Aggressive |
|---|---|---|---|
| Clusters | 46–47 | 43–44 | 37–38 |
| **Baseline**, coef h=6 | −0.067 | −0.067 | −0.046 |
| t | −2.31 | −2.04 | −1.24 |
| Significant horizons | 3, 4, 5, 6 | 5, 6 | none |
| **Long band**, fx coef | −0.0417 | −0.0464 | −0.0524 |
| t | −2.17 | −2.65 | −2.93 |
| p | 0.036 | 0.011 | 0.006 |
| imports − fx (p) | 0.007 | 0.009 | 0.004 |
| reserves − fx (p) | 0.015 | 0.008 | 0.003 |
| Joint W | 8.45 | 9.45 | 10.86 |
| Joint p | 0.015 | 0.009 | 0.004 |

The two results move in opposite directions under the same test. The average
response weakens; the mechanism strengthens. Note also that the aggressive
exclusion removes a disproportionate share of reserve-absorbing episodes
(954 → 572, against 824 → 669 for imports and 577 → 467 for fx), so the reserve
estimate in that column rests on a thinner base.

---

## Appendix A2 — pre-treatment placebos

No pre-treatment horizon is significant for either learner. Largest |t| is 1.80
(gradient boosting, h = −1). Coefficient magnitudes are 0.000–0.009, an order of
magnitude below the post-shock response.

---

## If your numbers differ

- **Small differences (fourth decimal):** likely a `ranger` version difference.
  The package's internal RNG has changed across releases.
- **A horizon flips significance near t = 2.013:** check that critical values come
  from t(G−1), not the normal, and that the G/(G−1) correction is applied.
- **Coefficients roughly 100× too large:** treatment and outcome are on different
  scales. See `docs/data_dictionary.md`.
- **Standard errors roughly half the reported size:** cross-fitting folds are being
  drawn over rows rather than over countries, or inference is not clustered.
- **Reserve buffer around 70 months rather than 6:** `reserves_months` is being
  computed against an annualised import figure.
