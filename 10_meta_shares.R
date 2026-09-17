#!/usr/bin/env Rscript
# ===========================================================================
# 10_meta_shares.R
# ---------------------------------------------------------------------------
# Meta MD + crosswalk -> one row per US county with folded shares o1/o2/o3,
# plus f0 and a coverage class.
#
# Output: meta_county_shares.csv, the only thing 11_fit_meta.R reads.
#
# Splitting this from the fitting is deliberate. The fold and the crosswalk
# are where the estimand is decided and where counties quietly disappear; a
# file you can open and count is worth more than the same logic buried at the
# top of a fitting script, which is where it used to live (in two copies,
# with different weights).
#
# Usage:  Rscript scripts/10_meta_shares.R
# ===========================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
})

source("config.R")
source("R/01_utils.R")
source("R/05_meta.R")

banner(sprintf("META SHARES  (scope: %s)", META_SCOPE))


# ---- inputs ----------------------------------------------------------------

md <- read_meta_md(META_MD_FILE)
say("MD rows: ", format(nrow(md), big.mark = ","),
    "  units: ", format(n_distinct(md$gadm_id), big.mark = ","),
    "  periods: ", n_distinct(md$ds))

xw <- load_crosswalk(XWALK_FILE)
say("crosswalk: ", nrow(xw), " counties, ",
    sum(xw$shared_meta_unit), " sharing a GADM unit with another county")

# The county frame. Everything is reported against this, so "3,143 counties,
# 2,980 with a Meta kernel" is a statement about the US rather than about
# whichever file happened to be joined.
frame <- read_keyed_csv(COUNTY_CENTROIDS_FILE, key_cols = "fips") |>
  select(fips, cen_pop) |>
  distinct(fips, .keep_all = TRUE)
say("county frame: ", nrow(frame), " counties")


# ---- shares ----------------------------------------------------------------

shares <- meta_county_shares(md, xw, scope = META_SCOPE)
say("counties with shares: ", nrow(shares))


# ---- coverage accounting ---------------------------------------------------
# Every county in the frame gets a reason. A county missing from the maps
# later should be explainable from this column alone.

ambiguous_fips <- read_keyed_csv(XWALK_FILE, key_cols = "fips") |>
  distinct(fips, gid_2) |>
  count(fips) |>
  filter(n > 1) |>
  pull(fips)

coverage <- frame |>
  left_join(select(shares, fips, o1, o2, o3, f0, n_periods, gid_2,
                   shared_meta_unit),
            by = "fips", relationship = "one-to-one") |>
  mutate(coverage = case_when(
    fips %in% ambiguous_fips           ~ "ambiguous_crosswalk",
    is.na(o1)  & is.na(gid_2)          ~ "no_crosswalk_row",
    is.na(o1)                          ~ "no_meta_unit",
    shared_meta_unit                   ~ "shared_meta_unit",
    TRUE                               ~ "ok"
  ))

say("coverage:")
coverage |> count(coverage, sort = TRUE) |> as.data.frame() |> print(row.names = FALSE)

pop_covered <- coverage |>
  filter(coverage %in% c("ok", "shared_meta_unit")) |>
  summarise(s = sum(as.numeric(cen_pop), na.rm = TRUE)) |>
  pull(s)
say(sprintf("population with a Meta kernel: %.1f%%",
            100 * pop_covered / sum(as.numeric(frame$cen_pop), na.rm = TRUE)))


# ===========================================================================
# CHECKS
# ---------------------------------------------------------------------------
# Hard stops for things that mean the code is wrong. Warnings for things that
# mean the data has moved.
# ===========================================================================

# 1. Shares are a distribution. If this fails the fold or the renormalisation
#    is wrong, and every fit downstream is meaningless.
bad <- shares |> filter(abs(o1 + o2 + o3 - 1) > 1e-8)
stopifnot("o1 + o2 + o3 does not sum to 1" = nrow(bad) == 0)
say("[ok] shares sum to 1 for all ", nrow(shares), " counties")

# 2. One row per county.
stopifnot("duplicate FIPS in shares" = !anyDuplicated(shares$fips))
say("[ok] one row per county")

# 3. Monotone decay. Meta's categories should fall away from home; a county
#    where o2 > o1 is either tiny or wrong, and it is worth seeing which.
non_mono <- shares |> filter(o2 > o1)
say(if (nrow(non_mono) == 0) "[ok] " else "[!]  ",
    nrow(non_mono), " counties with o2 > o1")

# 4. f0 near 0.346. Near-constant across US counties in every run so far
#    (sd 0.009). A drift here means the input file or the period set changed,
#    which is worth knowing before it shows up as a shifted map.
f0_med <- median(shares$f0, na.rm = TRUE)
f0_sd  <- sd(shares$f0, na.rm = TRUE)
say(sprintf("f0: median %.4f, sd %.4f  (expected ~0.346, sd ~0.009)",
            f0_med, f0_sd))
if (abs(f0_med - 0.346) > 0.02) {
  warning("f0 median has moved from its usual value -- check the period set ",
          "and that META_SCOPE is what you intended.")
}

# 5. Period count. Meta ships two-week periods; with one period there is no
#    averaging happening and nothing to average over.
say("periods per county: median ", median(shares$n_periods))


# ---- write -----------------------------------------------------------------

out <- shares |>
  select(fips, gid_2, f0, g1, g2, g3, o1, o2, o3, n_periods,
         shared_meta_unit, meta_scope) |>
  left_join(select(coverage, fips, coverage), by = "fips") |>
  arrange(fips) |>
  stamp_run()

write_atomic(out, META_SHARES_FILE)
say("wrote ", nrow(out), " rows to ", META_SHARES_FILE)

# Coverage for every county, including the ones with no shares, so the
# absences are documented rather than inferred.
write_atomic(select(coverage, fips, coverage, gid_2, cen_pop),
             file.path(TRACKB_DIR, "meta_coverage.csv"))
say("wrote coverage table")