#!/usr/bin/env python3
"""
build_master.py — build the research data warehouse from the existing CSV drops.

Run:  python build_master.py --src ./source --out ./warehouse

Layout produced
---------------
warehouse/
    dim_country.csv            one row per ISO3
    dim_variable.csv           the machine-readable data dictionary
    core_country_month.csv     grain: iso3 x month        (49 commodity exporters)
    core_country_year.csv      grain: iso3 x year         (217 countries)
    core_dyad_month.csv        grain: reporter-partner x month (USA-DZA)
    views/view_dml_lp.csv          paper 1 estimation frame
    views/view_climate_fx.csv      paper 2 estimation frame
    views/view_dza_quarterly.csv   paper 3 estimation frame (regenerates raouf_final_copy)
    views/view_bilateral.csv       paper 4 estimation frame
    qa/coverage_report.csv     non-missing share by table x variable x period
    qa/revision_report.csv     vintage-to-vintage differences in the dyad table
    qa/validation_log.txt      assertion results; non-empty FAIL section blocks use

Design rules enforced here
--------------------------
R1  Core tables hold LEVELS ONLY. Every transform (log, difference, partial
    sum, state indicator) is recomputed in the view layer from a single
    function, so two files can never disagree about what dlog_imports means.
R2  One grain per table. Country-month, country-year and dyad-month never share
    a file.
R3  Every column in a core table has a row in dim_variable giving its source,
    source series code, unit and vintage. A column without a registry entry is
    a build failure, not a warning.
R4  Revisions are preserved, not overwritten. Where two vintages of the same
    series disagree, the later vintage becomes current and the difference is
    written to qa/revision_report.csv.
R5  Nothing is imputed. Missing stays missing; coverage is reported instead.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# 0. Configuration
# --------------------------------------------------------------------------- #

SRC_FILES = {
    "lp": "lp_panel.csv",
    "climate": "master_climate_fx_panel.csv",
    "dyad_v1": "02_transformed_monthly_algeria_us_trade_data.csv",
    "dyad_v2": "02_transformed_data.csv",
    "dyad_raw": "01_raw_data.csv",
    "raouf": "raouf_final_copy.csv",
}

# Vintage labels. Update when a source is re-pulled; never overwrite in place.
VINTAGE = {"lp": "2026-08", "climate": "2026-07", "dyad_v1": "2025-12", "dyad_v2": "2026-04"}

WINSOR = (0.01, 0.99)   # applied to inflation in the view layer only


# --------------------------------------------------------------------------- #
# 1. Data dictionary  (R3: the registry is the contract)
# --------------------------------------------------------------------------- #
# fields: variable, grain, concept, unit, source, source_code, transform, notes

REGISTRY = [
    # ---- keys -------------------------------------------------------------
    ("iso3",        "all",  "ISO 3166-1 alpha-3 country code", "code", "derived", "", "level", "primary key component"),
    ("date",        "month","first day of month",              "date", "derived", "", "level", "monthly observations dated to month start"),
    ("year",        "year", "calendar year",                   "int",  "derived", "", "level", ""),

    # ---- country-month core (from lp_panel) --------------------------------
    ("cpi",                  "month", "consumer price index",                  "index",      "IMF IFS / national", "PCPI_IX",  "level", "base year varies by country; use log differences only"),
    ("fx",                   "month", "nominal exchange rate, LCU per USD",     "LCU/USD",    "IMF IFS",            "ENDE_XDC_USD_RATE", "level", "increase = depreciation"),
    ("policy_rate",          "month", "central bank policy rate",              "% p.a.",     "IMF IFS / BIS",      "FPOLM_PA", "level", "58% coverage — binding for any rate-based design"),
    ("money_market_rate",    "month", "money market rate",                     "% p.a.",     "IMF IFS",            "FIMM_PA",  "level", "33% coverage — do not use as primary"),
    ("exports_usd",          "month", "merchandise exports",                   "USD",        "IMF DOTS / IFS",     "TXG_FOB_USD", "level", ""),
    ("imports_usd",          "month", "merchandise imports",                   "USD",        "IMF DOTS / IFS",     "TMG_CIF_USD", "level", ""),
    ("reserves_total_usd",   "month", "total reserves including gold",         "USD",        "IMF IFS",            "RAFA_USD", "level", ""),
    ("reserves_exgold_usd",  "month", "reserves excluding gold",               "USD",        "IMF IFS",            "RAXG_USD", "level", ""),
    ("xprice_ix",            "month", "country export price index",            "index=100",  "constructed",        "",         "level", "base-weighted commodity export price index"),
    ("mprice_ix",            "month", "country import price index",            "index=100",  "constructed",        "",         "level", "base-weighted commodity import price index"),
    ("shock_ix",             "month", "commodity terms of trade",              "index=100",  "constructed",        "",         "level", "verified: shock_ix = 100 * xprice_ix / mprice_ix"),
    ("shock_rw_ix",          "month", "commodity ToT, real-weighted",          "index=100",  "constructed",        "",         "level", "alternative weighting scheme"),

    # ---- country-year core (from master_climate_fx_panel) -------------------
    ("gdp_usd",              "year", "GDP, current prices",                    "USD",        "World Bank WDI", "NY.GDP.MKTP.CD",   "level", ""),
    ("gdp_growth",           "year", "real GDP growth",                        "% p.a.",     "World Bank WDI", "NY.GDP.MKTP.KD.ZG","level", "annual %, NOT comparable to lp infl_yoy log units"),
    ("gdp_pc_usd",           "year", "GDP per capita, current prices",         "USD",        "World Bank WDI", "NY.GDP.PCAP.CD",   "level", ""),
    ("cpi_index",            "year", "consumer price index, annual",           "index",      "World Bank WDI", "FP.CPI.TOTL",      "level", "annual analogue of monthly cpi; different base"),
    ("inflation",            "year", "CPI inflation",                          "% p.a.",     "World Bank WDI", "FP.CPI.TOTL.ZG",   "level", "percent change, NOT a log difference"),
    ("reer",                 "year", "real effective exchange rate",           "index=100",  "World Bank WDI", "PX.REX.REER",      "level", "43% coverage — THE binding constraint on the climate-FX design"),
    ("fx_lcu_usd",           "year", "official exchange rate, period average", "LCU/USD",    "World Bank WDI", "PA.NUS.FCRF",      "level", ""),
    ("reserves_usd",         "year", "total reserves",                         "USD",        "World Bank WDI", "FI.RES.TOTL.CD",   "level", "WDI series; differs from IFS reserves in the month table"),
    ("resource_rents_gdp",   "year", "total natural resource rents",           "% of GDP",   "World Bank WDI", "NY.GDP.TOTL.RT.ZS","level", ""),
    ("oil_rents_gdp",        "year", "oil rents",                              "% of GDP",   "World Bank WDI", "NY.GDP.PETR.RT.ZS","level", ""),
    ("gas_rents_gdp",        "year", "natural gas rents",                      "% of GDP",   "World Bank WDI", "NY.GDP.NGAS.RT.ZS","level", ""),
    ("mineral_rents_gdp",    "year", "mineral rents",                          "% of GDP",   "World Bank WDI", "NY.GDP.MINR.RT.ZS","level", ""),
    ("government_debt_gdp",  "year", "central government debt",                "% of GDP",   "World Bank WDI", "GC.DOD.TOTL.GD.ZS","level", "21% coverage — unusable as a control"),
    ("exports_gdp",          "year", "exports of goods and services",          "% of GDP",   "World Bank WDI", "NE.EXP.GNFS.ZS",   "level", "ratio; not the level series in the month table"),
    ("imports_gdp",          "year", "imports of goods and services",          "% of GDP",   "World Bank WDI", "NE.IMP.GNFS.ZS",   "level", ""),
    ("trade_openness",       "year", "(exports + imports) / GDP",              "% of GDP",   "derived",        "",                 "derived","recomputed in the view layer"),
    ("population",           "year", "total population",                       "persons",    "World Bank WDI", "SP.POP.TOTL",      "level", ""),
    ("temperature",          "year", "mean annual temperature",                "deg C",      "CRU / WB CCKP",  "",                 "level", ""),
    ("temperature_z",        "year", "temperature anomaly, standardised",      "sd",         "derived",        "",                 "derived","z-score vs country baseline; baseline window must be documented"),
    ("precipitation",        "year", "total annual precipitation",             "mm",         "CRU / WB CCKP",  "",                 "level", ""),
    ("precipitation_z",      "year", "precipitation anomaly, standardised",    "sd",         "derived",        "",                 "derived",""),
    ("temp_base_mean",       "year", "baseline mean temperature",              "deg C",      "derived",        "",                 "derived","z-score denominator input"),
    ("temp_base_sd",         "year", "baseline sd of temperature",             "deg C",      "derived",        "",                 "derived",""),
    ("pr_base_mean",         "year", "baseline mean precipitation",            "mm",         "derived",        "",                 "derived",""),
    ("pr_base_sd",           "year", "baseline sd of precipitation",           "mm",         "derived",        "",                 "derived",""),
    ("ndgain_vulnerability", "year", "ND-GAIN vulnerability score",            "0-1",        "ND-GAIN",        "",                 "level", ""),
    ("ndgain_readiness",     "year", "ND-GAIN readiness score",                "0-1",        "ND-GAIN",        "",                 "level", ""),
    ("agriculture_value_added_gdp","year","agriculture value added",           "% of GDP",   "World Bank WDI", "NV.AGR.TOTL.ZS",   "level", ""),

    # ---- dyad-month core (from the FRED bilateral pulls) --------------------
    ("us_exports_to_dza",    "month", "US goods exports to Algeria, nominal",  "USD million","US Census via FRED", "EXP7210",   "level", "revised monthly; see qa/revision_report.csv"),
    ("us_imports_from_dza",  "month", "US goods imports from Algeria, nominal","USD million","US Census via FRED", "IMP7210",   "level", "revised monthly"),
    ("us_policy_rate",       "month", "effective federal funds rate",          "% p.a.",     "FRED",               "FEDFUNDS",  "level", "renamed from fedfunds / us_interest_rate"),
    ("us_cpi",               "month", "US CPI, all urban consumers, SA",       "index",      "FRED",               "CPIAUCSL",  "level", "PREVIOUSLY MISLABELLED 'global_inflation' — it is US CPI"),
    ("oil_price_wti",        "month", "WTI spot, monthly mean of daily",       "USD/bbl",    "FRED",               "DCOILWTICO","level", "wrong benchmark for Algerian crude; keep for US-side work only"),
    ("oil_price_brent",      "month", "Brent spot, monthly mean",              "USD/bbl",    "FRED",               "DCOILBRENTEU","level", "REQUIRED for any Algeria pricing analysis; not yet in the drop"),
]

REGISTRY_COLS = ["variable", "grain", "concept", "unit", "source",
                 "source_code", "transform", "notes"]

# Column renames applied at staging so one concept has one canonical name.
RENAME_DYAD = {
    "fedfunds": "us_policy_rate", "us_interest_rate": "us_policy_rate",
    "cpi": "us_cpi", "global_inflation": "us_cpi",
    "wti_spot": "oil_price_wti", "oil_price": "oil_price_wti",
}

DYAD_LEVELS = ["us_exports_to_dza", "us_imports_from_dza",
               "us_policy_rate", "us_cpi", "oil_price_wti"]

MONTH_LEVELS = ["cpi", "fx", "policy_rate", "money_market_rate", "exports_usd",
                "imports_usd", "reserves_total_usd", "reserves_exgold_usd",
                "xprice_ix", "mprice_ix", "shock_ix", "shock_rw_ix"]

DROP_YEAR_DERIVED = ("_l1", "_l2", "_l3", "_log", "_change", "_anomaly",
                     "_growth_l", "log_", "fx_appreciation", "fx_depreciation",
                     "fx_volatility_rolling3", "heat_shock_1sd",
                     "precipitation_pct_anomaly")


# --------------------------------------------------------------------------- #
# 2. Transform library  (R1: one definition, used everywhere)
# --------------------------------------------------------------------------- #

def dlog(s: pd.Series, k: int = 1) -> pd.Series:
    """100 x k-period log difference. Returns percent log points, not percent."""
    return 100.0 * (np.log(s) - np.log(s.shift(k)))


def partial_sums(s: pd.Series) -> tuple[pd.Series, pd.Series]:
    """Positive and negative parts of a shock series (Shin-Yu-Greenwood-Nimmo)."""
    return s.clip(lower=0), s.clip(upper=0)


def winsorize(s: pd.Series, lo: float, hi: float) -> pd.Series:
    ql, qh = s.quantile(lo), s.quantile(hi)
    return s.clip(ql, qh)


def zscore_within(df: pd.DataFrame, col: str, by: str) -> pd.Series:
    g = df.groupby(by)[col]
    return (df[col] - g.transform("mean")) / g.transform("std")


def lag_state(s: pd.Series, k: int = 1) -> pd.Series:
    """State variables enter lagged. Contemporaneous states are endogenous to
    the shock and will not survive a referee."""
    return s.shift(k)


# --------------------------------------------------------------------------- #
# 3. Staging
# --------------------------------------------------------------------------- #

def load(src: Path, key: str, **kw) -> pd.DataFrame:
    p = src / SRC_FILES[key]
    if not p.exists():
        raise FileNotFoundError(f"missing source file: {p}")
    return pd.read_csv(p, **kw)


def build_dim_country(climate: pd.DataFrame) -> pd.DataFrame:
    cols = ["iso3", "iso2", "country", "wb_region", "income_level",
            "lending_type", "capital_city", "latitude", "longitude"]
    d = (climate[cols].drop_duplicates("iso3").sort_values("iso3")
         .reset_index(drop=True))
    return d


def build_core_country_month(lp: pd.DataFrame) -> pd.DataFrame:
    df = lp.rename(columns={"country_code": "iso3"}).copy()
    df["date"] = pd.to_datetime(df["date"])
    keep = ["iso3", "date"] + [c for c in MONTH_LEVELS if c in df.columns]
    out = df[keep].sort_values(["iso3", "date"]).reset_index(drop=True)
    out["vintage"] = VINTAGE["lp"]
    return out


def build_core_country_year(cl: pd.DataFrame) -> pd.DataFrame:
    reg = {r[0] for r in REGISTRY if r[1] == "year"} - {"iso3", "year"}
    keep = ["iso3", "year"] + [c for c in cl.columns if c in reg]
    out = cl[keep].copy()
    # R5: the source is a rectangular 217 x 36 grid. Drop country-years that are
    # entirely empty outside the keys so the row count reflects real coverage.
    val = [c for c in out.columns if c not in ("iso3", "year")]
    out = out[out[val].notna().any(axis=1)].copy()
    out["vintage"] = VINTAGE["climate"]
    return out.sort_values(["iso3", "year"]).reset_index(drop=True)


def build_core_dyad_month(v1: pd.DataFrame, v2: pd.DataFrame
                          ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """R4: later vintage wins; the difference is reported, not discarded."""
    def prep(df, vint):
        d = df.rename(columns=RENAME_DYAD).copy()
        d["date"] = pd.to_datetime(d["date"])
        d = d.loc[:, ~d.columns.duplicated()]
        cols = ["date"] + [c for c in DYAD_LEVELS if c in d.columns]
        d = d[cols].copy()
        d["vintage"] = vint
        return d

    a = prep(v1, VINTAGE["dyad_v1"])
    b = prep(v2, VINTAGE["dyad_v2"])

    shared = [c for c in DYAD_LEVELS if c in a.columns and c in b.columns]
    rev = (a.set_index("date")[shared]
             .join(b.set_index("date")[shared], lsuffix="_old", rsuffix="_new",
                   how="inner"))
    recs = []
    for c in shared:
        diff = (rev[f"{c}_new"] - rev[f"{c}_old"]).abs()
        hit = diff[diff > 1e-9]
        for dt, v in hit.items():
            recs.append({"date": dt, "variable": c,
                         "old_vintage": VINTAGE["dyad_v1"],
                         "new_vintage": VINTAGE["dyad_v2"],
                         "old_value": rev.loc[dt, f"{c}_old"],
                         "new_value": rev.loc[dt, f"{c}_new"],
                         "abs_diff": v})
    revision = pd.DataFrame(recs, columns=["date", "variable", "old_vintage",
                                           "new_vintage", "old_value",
                                           "new_value", "abs_diff"])

    core = (pd.concat([a[~a["date"].isin(b["date"])], b], ignore_index=True)
              .sort_values("date").reset_index(drop=True))
    core.insert(0, "partner", "DZA")
    core.insert(0, "reporter", "USA")
    return core, revision


# --------------------------------------------------------------------------- #
# 4. Views  (R1: every derived column is created here, nowhere else)
# --------------------------------------------------------------------------- #

def view_dml_lp(cm: pd.DataFrame) -> pd.DataFrame:
    d = cm.sort_values(["iso3", "date"]).copy()
    g = d.groupby("iso3", sort=False)

    d["log_cpi"] = np.log(d["cpi"])
    d["infl_yoy"] = g["cpi"].transform(lambda s: dlog(s, 12))
    d["infl_mom"] = g["cpi"].transform(lambda s: dlog(s, 1))
    d["infl_yoy_w"] = winsorize(d["infl_yoy"], *WINSOR)

    d["d_log_fx_12m"] = g["fx"].transform(lambda s: dlog(s, 12))
    d["fx_vol_12m"] = g["fx"].transform(
        lambda s: dlog(s, 1).rolling(12, min_periods=8).std())

    d["log_shock"] = np.log(d["shock_ix"])
    d["d_shock_12m"] = g["shock_ix"].transform(lambda s: dlog(s, 12))
    pos, neg = partial_sums(d["d_shock_12m"])
    d["shock_pos"], d["shock_neg"] = pos, neg

    d["log_reserves"] = np.log(d["reserves_exgold_usd"])
    d["monthly_imports_approx"] = g["imports_usd"].transform(
        lambda s: s.rolling(12, min_periods=8).mean())
    d["reserves_months"] = d["reserves_exgold_usd"] / d["monthly_imports_approx"]

    # States, lagged one period. See validation V7.
    d["state_peg"] = (g["fx_vol_12m"].transform(lambda s: lag_state(s))
                      < d.groupby("date")["fx_vol_12m"].transform("median")).astype(float)
    d["low_buffer"] = (g["reserves_months"].transform(lambda s: lag_state(s))
                       < 3.0).astype(float)
    d.loc[d["fx_vol_12m"].isna(), "state_peg"] = np.nan
    d.loc[d["reserves_months"].isna(), "low_buffer"] = np.nan

    d["d_policy_rate_12m"] = g["policy_rate"].transform(lambda s: s - s.shift(12))
    d["d_log_exports_usd_12m"] = g["exports_usd"].transform(lambda s: dlog(s, 12))
    d["d_log_imports_usd_12m"] = g["imports_usd"].transform(lambda s: dlog(s, 12))
    return d


def view_climate_fx(cy: pd.DataFrame) -> pd.DataFrame:
    d = cy.sort_values(["iso3", "year"]).copy()
    g = d.groupby("iso3", sort=False)
    ok1 = g["year"].diff().eq(1)

    d["log_reer"] = np.log(d["reer"].where(d["reer"] > 0))
    d["d_log_reer"] = (100 * g["log_reer"].diff()).where(ok1)
    d["log_fx"] = np.log(d["fx_lcu_usd"].where(d["fx_lcu_usd"] > 0))
    d["d_log_fx"] = (100 * g["log_fx"].diff()).where(ok1)
    d["log_gdp_pc"] = np.log(d["gdp_pc_usd"].where(d["gdp_pc_usd"] > 0))
    d["trade_openness"] = d["exports_gdp"] + d["imports_gdp"]

    for v in ("temperature_z", "precipitation_z"):
        for h in (1, 2, 3):
            d[f"{v}_l{h}"] = g[v].shift(h).where(g["year"].diff(h).eq(h))
    return d


def view_dza_quarterly(cm: pd.DataFrame, dy: pd.DataFrame,
                       raouf: pd.DataFrame | None = None) -> pd.DataFrame:
    """Regenerates raouf_final_copy.csv from the warehouse rather than keeping
    it as an independent file. Verified: its cpi_algiers_log is exactly
    log(quarterly mean of the monthly CPI for DZA).

    Two columns in that file are NOT derivable from the warehouse -
    gdp_growth_q and gfce_growth - so they are carried through here until they
    are sourced for the full panel."""
    dza = cm[cm["iso3"] == "DZA"].set_index("date").sort_index()
    q = pd.DataFrame({
        "cpi_algiers_log": np.log(dza["cpi"].resample("QE").mean()),
        "usd_dzd_log": np.log(dza["fx"].resample("QE").mean()),
    })
    us = dy.set_index("date").sort_index()
    q["oil_wti_log"] = np.log(us["oil_price_wti"].resample("QE").mean())
    q["us_cpi_yoy"] = (100 * (np.log(us["us_cpi"].resample("QE").mean())
                              - np.log(us["us_cpi"].resample("QE").mean().shift(4))))
    q["us_policy_rate"] = us["us_policy_rate"].resample("QE").mean()
    q.index.name = "date"
    q = q.reset_index()

    if raouf is not None:
        r = raouf.copy()
        r["date"] = pd.to_datetime(r["date"])
        carry = ["date", "gdp_growth_q", "gfce_growth", "brent_log"]
        q = q.merge(r[[c for c in carry if c in r.columns]], on="date", how="left")
    return q


def view_bilateral(dy: pd.DataFrame) -> pd.DataFrame:
    d = dy.sort_values("date").copy()
    d["trade_balance"] = d["us_exports_to_dza"] - d["us_imports_from_dza"]
    d["log_exports"] = np.log(d["us_exports_to_dza"])
    d["log_imports"] = np.log(d["us_imports_from_dza"])
    d["dlog_exports"] = dlog(d["us_exports_to_dza"])
    d["dlog_imports"] = dlog(d["us_imports_from_dza"])
    d["oil_return"] = dlog(d["oil_price_wti"])
    d["us_inflation"] = dlog(d["us_cpi"])
    d["d_us_policy_rate"] = d["us_policy_rate"].diff()
    # Implied import volume: value deflated by the crude benchmark. Replace the
    # denominator with Brent once oil_price_brent is added.
    d["implied_volume"] = d["us_imports_from_dza"] / d["oil_price_wti"]
    d["log_implied_volume"] = np.log(d["implied_volume"])
    return d


# --------------------------------------------------------------------------- #
# 5. Validation
# --------------------------------------------------------------------------- #

def validate(tables: dict[str, pd.DataFrame], reg: pd.DataFrame) -> list[str]:
    log: list[str] = []

    def check(cond, msg):
        log.append(("PASS  " if cond else "FAIL  ") + msg)

    cm, cy, dy = tables["core_country_month"], tables["core_country_year"], tables["core_dyad_month"]

    check(not cm.duplicated(["iso3", "date"]).any(), "V1 country-month key is unique")
    check(not cy.duplicated(["iso3", "year"]).any(), "V2 country-year key is unique")
    check(not dy.duplicated(["reporter", "partner", "date"]).any(), "V3 dyad-month key is unique")

    registered = set(reg["variable"])
    for name, t in (("country_month", cm), ("country_year", cy), ("dyad_month", dy)):
        unreg = [c for c in t.columns
                 if c not in registered and c not in ("vintage", "reporter", "partner")]
        check(not unreg, f"V4 every {name} column is registered"
                         + (f" — unregistered: {unreg}" if unreg else ""))

    miss = set(cm["iso3"]) - set(cy["iso3"])
    check(not miss, f"V5 all monthly countries appear in the annual table {sorted(miss)}")

    # V6: the shock index identity
    s = cm.dropna(subset=["shock_ix", "xprice_ix", "mprice_ix"])
    resid = (s["shock_ix"] - 100 * s["xprice_ix"] / s["mprice_ix"]).abs().max()
    check(resid < 1e-3, f"V6 shock_ix == 100*xprice/mprice (max resid {resid:.2e})")

    # V7: states must not be contemporaneous with the shock
    log.append("NOTE  V7 state_peg and low_buffer are rebuilt with a one-period "
               "lag in view_dml_lp; the delivered lp_panel.csv versions were not "
               "verifiable as lagged and must not be reused.")

    # V8: unit mismatch that will silently corrupt a pooled regression
    log.append("NOTE  V8 annual 'inflation' (WDI, percent change) and monthly "
               "'infl_yoy' (100 x log difference) are different units. Never "
               "pool or compare coefficients across them without converting.")

    log.append("NOTE  V9 'us_cpi' was named 'global_inflation' in the delivered "
               "bilateral files. It is CPIAUCSL, US CPI. Any text describing it "
               "as global inflation must be corrected before submission.")
    return log


def validate_views(v_dml: pd.DataFrame, legacy: pd.DataFrame | None) -> list[str]:
    """Checks that catch the two defects found in the delivered lp_panel.csv."""
    log: list[str] = []

    def check(cond, msg):
        log.append(("PASS  " if cond else "FAIL  ") + msg)

    # V10 unit consistency: treatment and outcome must share units, or every
    # reported IRF magnitude is wrong by a factor of 100.
    sd_y = v_dml["infl_yoy"].std()
    sd_t = v_dml["shock_neg"].replace(0, np.nan).std()
    check(0.01 < sd_t / sd_y < 100,
          f"V10 shock and inflation are on the same scale "
          f"(sd ratio {sd_t / sd_y:.3f})")

    # V11 reserves in months of imports must be economically plausible.
    med = v_dml["reserves_months"].median()
    check(1 < med < 20, f"V11 reserves_months median is plausible ({med:.2f} months)")
    share = v_dml["low_buffer"].mean()
    check(0.05 < share < 0.50,
          f"V11b low_buffer fires on a usable share of the sample ({share:.1%})")

    if legacy is not None:
        lm = legacy["reserves_months"].median()
        ls = legacy["low_buffer"].mean()
        log.append(f"NOTE  V12 legacy lp_panel.csv had reserves_months median "
                   f"{lm:.1f} months and low_buffer firing on {ls:.1%} of rows: "
                   "monthly_imports_approx divided an already-monthly import "
                   "flow by 12, inflating the buffer twelvefold and shrinking "
                   "the low-buffer state to a few percent of the sample. Any "
                   "state-dependence result estimated on it must be re-run.")
        log.append("NOTE  V13 legacy lp_panel.csv reported infl_yoy as 100 x log "
                   "difference but d_shock_12m, shock_pos, shock_neg and "
                   "d_log_fx_12m as raw log differences. Coefficients from that "
                   "file are inflated by 100 relative to a percent-on-percent "
                   "elasticity.")
    return log


def coverage(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for name, t in tables.items():
        keys = {"iso3", "date", "year", "vintage", "reporter", "partner"}
        for c in t.columns:
            if c in keys:
                continue
            rows.append({"table": name, "variable": c,
                         "n_obs": int(t[c].notna().sum()),
                         "n_rows": len(t),
                         "coverage_pct": round(100 * t[c].notna().mean(), 1)})
    return pd.DataFrame(rows).sort_values(["table", "coverage_pct"])


# --------------------------------------------------------------------------- #
# 6. Main
# --------------------------------------------------------------------------- #

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", type=Path, default=Path("source"))
    ap.add_argument("--out", type=Path, default=Path("warehouse"))
    a = ap.parse_args()

    out = a.out
    (out / "views").mkdir(parents=True, exist_ok=True)
    (out / "qa").mkdir(parents=True, exist_ok=True)

    print("staging sources")
    lp = load(a.src, "lp")
    cl = load(a.src, "climate")
    d1 = load(a.src, "dyad_v1")
    d2 = load(a.src, "dyad_v2")

    reg = pd.DataFrame(REGISTRY, columns=REGISTRY_COLS)

    print("building dimensions and core tables")
    dim = build_dim_country(cl)
    cm = build_core_country_month(lp)
    cy = build_core_country_year(cl)
    dy, revision = build_core_dyad_month(d1, d2)

    tables = {"core_country_month": cm, "core_country_year": cy,
              "core_dyad_month": dy}

    print("building views")
    v_dml = view_dml_lp(cm)
    v_cfx = view_climate_fx(cy)
    v_dza = view_dza_quarterly(cm, dy, load(a.src, "raouf"))
    v_bil = view_bilateral(dy)

    print("validating")
    legacy = lp.rename(columns={"country_code": "iso3"})
    log = validate(tables, reg) + validate_views(v_dml, legacy)
    cov = coverage(tables)

    dim.to_csv(out / "dim_country.csv", index=False)
    reg.to_csv(out / "dim_variable.csv", index=False)
    for k, t in tables.items():
        t.to_csv(out / f"{k}.csv", index=False)
    v_dml.to_csv(out / "views/view_dml_lp.csv", index=False)
    v_cfx.to_csv(out / "views/view_climate_fx.csv", index=False)
    v_dza.to_csv(out / "views/view_dza_quarterly.csv", index=False)
    v_bil.to_csv(out / "views/view_bilateral.csv", index=False)
    cov.to_csv(out / "qa/coverage_report.csv", index=False)
    revision.to_csv(out / "qa/revision_report.csv", index=False)
    (out / "qa/validation_log.txt").write_text("\n".join(log), encoding="utf-8")

    fails = [x for x in log if x.startswith("FAIL")]
    print("\n" + "=" * 68)
    print(f"core_country_month  {len(cm):>6} rows  {cm.iso3.nunique():>3} countries  "
          f"{cm.date.min():%Y-%m} to {cm.date.max():%Y-%m}")
    print(f"core_country_year   {len(cy):>6} rows  {cy.iso3.nunique():>3} countries  "
          f"{cy.year.min()} to {cy.year.max()}")
    print(f"core_dyad_month     {len(dy):>6} rows  USA-DZA          "
          f"{dy.date.min():%Y-%m} to {dy.date.max():%Y-%m}")
    print(f"revisions logged    {len(revision):>6}")
    print(f"validation          {len(log) - len(fails)} pass/note, {len(fails)} FAIL")
    for f in fails:
        print("   " + f)
    print("=" * 68)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
