# Data dictionary

`data/derived/view_dml_lp.csv` — 21,560 rows (49 countries × 440 months), 33 columns.
Grain: one row per country-month. Primary key: `iso3` × `date`.
Span: 1990-01 to 2026-08. Vintage: 2026-08.

Sources: IMF International Financial Statistics (IFS), IMF Direction of Trade
Statistics (DOTS), IMF Primary Commodity Price System (PCPS). Nothing is
interpolated, spliced or back-cast; missing observations remain missing.

## Variables

| Column | Unit | Source | Definition |
|---|---|---|---|
| `iso3` | code | derived | ISO 3166-1 alpha-3 country code. Key. |
| `date` | date | derived | First day of month. Key. |
| `cpi` | index | IMF IFS | Consumer price index. Base year varies by country; use log differences only. Coverage 52.5%. |
| `fx` | LCU/USD | IMF IFS | Nominal exchange rate, local currency per US dollar. An increase is a depreciation. |
| `policy_rate` | per cent | IMF IFS | Central bank policy rate. Coverage 58% — used in the robustness control set only. |
| `money_market_rate` | per cent | IMF IFS | Money market rate. Not used in the reported specifications. |
| `exports_usd` | USD | IMF DOTS | Merchandise exports, monthly, current US dollars. |
| `imports_usd` | USD | IMF DOTS | Merchandise imports, monthly, current US dollars. |
| `reserves_total_usd` | USD | IMF IFS | Total reserves including gold. The analysis uses the ex-gold series instead. |
| `reserves_exgold_usd` | USD | IMF IFS | International reserves excluding gold. Margin outcome (as 12-month log change). |
| `xprice_ix` | index | IMF PCPS | Export commodity price index, pre-sample trade weights (Gruss & Kebhaj 2019). |
| `mprice_ix` | index | IMF PCPS | Import commodity price index, pre-sample trade weights. |
| `shock_ix` | index | derived | Commodity terms of trade: 100 × xprice_ix / mprice_ix. Equation (1). |
| `shock_rw_ix` | index | derived | Alternative rolling-weight variant. Not used in the reported results. |
| `vintage` | YYYY-MM | derived | Data vintage. All reported results use 2026-08. |
| `log_cpi` | log index | derived | Natural log of cpi. |
| `infl_yoy` | per cent | derived | 100 × (log cpi_t − log cpi_{t−12}). Unwinsorised. |
| `infl_mom` | per cent | derived | 100 × (log cpi_t − log cpi_{t−1}). |
| `infl_yoy_w` | per cent | derived | infl_yoy winsorised at the 1st and 99th percentiles. OUTCOME in the baseline. |
| `d_log_fx_12m` | per cent | derived | 100 × 12-month log change in fx. Positive = depreciation. Control and margin outcome. |
| `fx_vol_12m` | per cent | derived | Rolling 12-month standard deviation of monthly log fx changes. Control. |
| `log_shock` | log index | derived | Natural log of shock_ix. |
| `d_shock_12m` | per cent | derived | 100 × 12-month log change in shock_ix. Equation (2). The shock variable. |
| `shock_pos` | per cent | derived | max(d_shock_12m, 0). Positive partial sum, equation (3). |
| `shock_neg` | per cent | derived | min(d_shock_12m, 0). Negative partial sum — the TREATMENT. |
| `log_reserves` | log USD | derived | Natural log of reserves_exgold_usd. |
| `monthly_imports_approx` | USD | derived | Trailing 12-month mean of monthly imports_usd. |
| `reserves_months` | months | derived | reserves_exgold_usd / monthly_imports_approx. NOTE: computed against the MONTHLY flow, not an annualised figure — see README. |
| `state_peg` | 0/1 | derived | Indicator for a pegged arrangement. |
| `low_buffer` | 0/1 | derived | Indicator for reserves_months below the sample 25th percentile. Fires on 23.4% of observations. |
| `d_policy_rate_12m` | pp | derived | 12-month change in policy_rate. Robustness control set only. |
| `d_log_exports_usd_12m` | per cent | derived | 100 × 12-month log change in exports_usd. Control. |
| `d_log_imports_usd_12m` | per cent | derived | 100 × 12-month log change in imports_usd. Control and margin outcome. |

## Coverage

| Variable | Non-missing | Share |
|---|---|---|
| Trade flows | 20,508 | 95.1% |
| Terms-of-trade shock | 20,513 | 95.1% |
| Exchange rate | 19,156 | 88.9% |
| Reserves (ex gold) | 17,575 | 81.5% |
| Consumer prices | 11,316 | 52.5% |
| Policy rate | — | 58% |

Consumer prices are the binding constraint on the estimation sample, not the shock.
The effective sample for the inflation regressions is roughly 10,200 observations
across 46–47 countries. Papua New Guinea has no usable CPI series and Australia has
16 monthly observations, so the inflation results rest on 47 countries.

## Two measurement points

**Reserve adequacy.** `reserves_months` divides reserves by the *monthly* import flow.
Dividing by an annualised figure inflates the buffer twelvefold: the median in this
panel moves from 6.06 to 70.3 months and the low-buffer share from 23.4% to 2.4%.

**Scale consistency.** Treatment and outcome must both be in percentage log differences.
Mixing raw log differences with percentage log differences inflates the estimated
elasticity by a factor of 100 while leaving signs, *t*-statistics and the shape of the
impulse response unchanged — which is what makes the error easy to miss.

## Derived in the estimation scripts, not stored here

`d_log_reserves_12m`, `infl_z` (within-country standardised inflation), the size-split
terms `shock_big` / `shock_small`, the margin-interaction terms `shk_imports`,
`shk_reserves`, `shk_fx`, and six own lags of each. See `code/02_replicate_empirical.R`.
