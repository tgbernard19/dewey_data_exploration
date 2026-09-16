#!/usr/bin/env Rscript
# ===========================================================================
# build_veraset_observed.R  --  the empirical destination table
# ---------------------------------------------------------------------------
# Produces two files:
#
#   veraset_observed_pairs.csv
#     origin_fips, dest_fips, n_devices, n_device_days, n_visits, w_dwell,
#     d_km_median
#
#   veraset_observed_panel.csv
#     origin_fips, n_dest_counties, n_device_appearances, n_device_days,
#     n_visits, n_devices_panel, appearances_per_device, d_max_km, built_on
#
# This is NOT a fit. It is the observed Veraset destination distribution,
# which absolute_flows_compare.R scales to the population and treats as the
# benchmark the Meta model is trying to recapitulate.
#
# WHAT COUNTS
#   Non-home locations (is_home = FALSE) within D_MAX of the device's home,
#   from devices whose home is CBG-resolved.
#
# FOUR WEIGHTINGS, ONE PASS
#   n_devices      one count per device per destination county
#   n_device_days  device-days at that destination
#   n_visits       visit count
#   w_dwell        dwell minutes, = dwell_pos_sum + TAU * n_tau_visits
#
#   w_dwell is the one that matches the model. f0 is a share of PINGS, so
#   pop * (1 - f0) is the number of people away at a random instant, and
#   Meta's away categories are time-weighted for the same reason -- the
#   kernel says where people ARE, duration-weighted, not where they GO.
#   Benchmarking that against a device count compares prevalence against
#   incidence: a 15-minute stop and an 8-hour workday count equally on one
#   side and 32:1 on the other. This is the estimand mismatch that moved the
#   median Veraset-Meta gap from -0.111 to +0.007 in August, and the
#   resolution then was also the time-weighted quantity.
#
#   All four are emitted here because re-running this extract is the
#   expensive step; switching weighting downstream should be one line.
#
# THE CBG FILTER IS NOT OPTIONAL
#   ~23% of devices have geohash-5 homes carrying a median 3.467 km
#   displacement on rows whose true displacement is zero. They put 0.094 of
#   mass in [0,1) against 0.789 for CBG-resolved. Including them mixes two
#   populations and the near field becomes meaningless -- which is what made
#   six unrelated kernel families all fail at tv 0.10-0.18 before September.
#
# WHY THE TWO-STEP GROUP BY
#   Counting distinct devices per origin-destination pair directly is
#   COUNT(DISTINCT did) over millions of groups, which cannot stream and needs
#   a hash set per group. That is what stalled other_visits at 10% with 32 GB
#   spilled. Collapsing to (origin, dest, device) grain first turns the
#   distinct count into a plain row count. Get to the right grain, then count.
#
# THE CAVEAT THAT MATTERS MOST DOWNSTREAM
#   A device is observed only where there is a POI. So this table measures
#   visits-to-POIs, not trips, and it under-observes destinations with thin
#   POI coverage -- rural ones. The scaled-up flows inherit that. When the
#   Meta model over-predicts a rural destination, POI sparsity here is a live
#   explanation and has to be ruled out before it is called model error.
# ===========================================================================

suppressPackageStartupMessages({
  library(DBI); library(duckdb); library(dplyr); library(readr); library(stringr)
})

# ---- paths (edit for this machine) ----------------------------------------
KERNEL_DIR  <- "E:/dewey-june2025/kernel"
DEVDAY_GLOB <- file.path(KERNEL_DIR, "devday_parts_v2", "*.parquet")
HOMERES_PQ  <- file.path(KERNEL_DIR, "trackB", "device_home_resolution.parquet")
OUT_PAIRS   <- file.path(KERNEL_DIR, "trackB", "veraset_observed_pairs.csv")
OUT_PANEL   <- file.path(KERNEL_DIR, "trackB", "veraset_observed_panel.csv")

TMP_DIR   <- "E:/duckdb-tmp"
N_THREADS <- 32
MEM_LIMIT <- "400GB"

D_MAX <- 500   # km. must match D_MAX in absolute_flows_compare.R

# Minutes credited to a single-ping visit. minimum_dwell is last-ping minus
# first-ping, so a zero is a visit of unknown duration, not zero duration --
# ~33% of home rows, 42% of work, 73% of other. Deleting them would drop
# short visits preferentially. Profiles were invariant across tau = 1, 15, 60
# (largest mean band difference 0.0017), so 15 is a safe middle; the
# separately-carried columns make a different tau a re-query, not a re-scan.
TAU <- 15

# ---- home-resolution rule --------------------------------------------------
# device_home_resolution.parquet is device grain:
#   did, cbg_home (flag), n_home_cbg, n_home_geo
# n_home_cbg / n_home_geo count that device's is_home rows by how the location
# was resolved. A device with both is a mixture, and its away displacements are
# measured from an ambiguous home.
#
#  "pure_cbg"  n_home_cbg > 0 AND n_home_geo = 0   <- default, strictest
#  "majority"  n_home_cbg > n_home_geo
#  "any_cbg"   n_home_cbg > 0
#  "flag"      cbg_home = 1, whatever the builder meant by it
#
# Default is pure_cbg. The artifact is a 3.467 km median displacement on rows
# whose truth is zero, localised at 2.5-10 km -- precisely the range this
# comparison turns on -- so admitting mixed devices reintroduces a fraction of
# it for no gain in panel size worth having. The crosstab below tells you what
# that costs and what cbg_home actually encodes; check it before overriding.
HOME_RULE <- "pure_cbg"

HOME_WHERE <- switch(
  HOME_RULE,
  pure_cbg = "n_home_cbg > 0 AND n_home_geo = 0",
  majority = "n_home_cbg > n_home_geo",
  any_cbg  = "n_home_cbg > 0",
  flag     = "cbg_home = 1",
  stop("Unknown HOME_RULE: ", HOME_RULE))

# ---- connection ------------------------------------------------------------
# No on.exit() here. At the top level of a script, on.exit attaches to the
# current top-level expression rather than to a function frame, so it fires as
# soon as that line finishes and the connection is dead by the next statement
# ("Invalid connection ... rapi_prepare"). The disconnect is explicit at the
# bottom instead.
con <- dbConnect(duckdb::duckdb())
stopifnot("Connection did not open." = dbIsValid(con))

dbExecute(con, sprintf("SET threads=%d", N_THREADS))
dbExecute(con, sprintf("SET memory_limit='%s'", MEM_LIMIT))
dbExecute(con, sprintf("SET temp_directory='%s'", TMP_DIR))
dbExecute(con, "SET preserve_insertion_order=false")

# ---- 0. schema, footers only ----------------------------------------------
message("[schema] devday_parts_v2")
print(dbGetQuery(con, sprintf(
  "DESCRIBE SELECT * FROM read_parquet('%s') LIMIT 0", DEVDAY_GLOB)))
message("[schema] device_home_resolution")
hr_schema <- dbGetQuery(con, sprintf(
  "DESCRIBE SELECT * FROM read_parquet('%s') LIMIT 0", HOMERES_PQ))
print(hr_schema)

need_hr <- c("did", "n_home_cbg", "n_home_geo")
if (length(setdiff(need_hr, hr_schema$column_name)))
  stop("device_home_resolution is missing ",
       paste(setdiff(need_hr, hr_schema$column_name), collapse = ", "),
       ". Columns are: ", paste(hr_schema$column_name, collapse = ", "))

# ---- what does cbg_home actually mean? -------------------------------------
# Crosstab the flag against the three candidate definitions. If cbg_home = 1
# lines up exactly with pure_cbg, the flag is the strict rule and the choice
# is moot. If it lines up with any_cbg, the flag admits mixed devices and
# HOME_RULE is doing real work.
message("[home rule] what cbg_home encodes:")
print(dbGetQuery(con, sprintf("
  SELECT cbg_home,
         COUNT(*)                                                  AS n_devices,
         SUM(CASE WHEN n_home_cbg > 0 AND n_home_geo = 0 THEN 1 ELSE 0 END) AS pure_cbg,
         SUM(CASE WHEN n_home_cbg > 0 AND n_home_geo > 0 THEN 1 ELSE 0 END) AS mixed,
         SUM(CASE WHEN n_home_cbg = 0 THEN 1 ELSE 0 END)           AS no_cbg
  FROM read_parquet('%s')
  GROUP BY cbg_home ORDER BY cbg_home", HOMERES_PQ)))

# And how much panel each rule keeps, so the choice is made on numbers.
message("[home rule] devices kept by each rule:")
print(dbGetQuery(con, sprintf("
  SELECT COUNT(*) AS all_devices,
         SUM(CASE WHEN n_home_cbg > 0 AND n_home_geo = 0 THEN 1 ELSE 0 END) AS pure_cbg,
         SUM(CASE WHEN n_home_cbg > n_home_geo THEN 1 ELSE 0 END)           AS majority,
         SUM(CASE WHEN n_home_cbg > 0 THEN 1 ELSE 0 END)                    AS any_cbg,
         SUM(CASE WHEN cbg_home = 1 THEN 1 ELSE 0 END)                      AS flag
  FROM read_parquet('%s')", HOMERES_PQ)))

# ---- 1. device grain: CBG-resolved homes -----------------------------------
# NOTE ON did: the ids are hashed caid as UBIGINT and run past 1e19, well
# beyond the 2^53 an R double represents exactly. R prints them as <dbl> and
# will silently collide distinct devices if you dedupe or join on them in R.
# Every did operation here stays inside DuckDB for that reason. Do not pull
# the column into R and join it there.
message(sprintf("[1/4] device set  (rule: %s)", HOME_RULE))
dbExecute(con, sprintf("
  CREATE OR REPLACE TABLE dev AS
  SELECT DISTINCT did FROM read_parquet('%s')
  WHERE %s
", HOMERES_PQ, HOME_WHERE))

n_dev <- dbGetQuery(con, "SELECT COUNT(*) AS n FROM dev")$n
message(sprintf("  devices kept: %s", format(n_dev, big.mark = ",")))
if (n_dev == 0)
  stop("No devices matched HOME_RULE = '", HOME_RULE, "'. See the crosstab above.")

# ---- 2. origin x destination x device --------------------------------------
# One row per device per destination county it was seen in. This is the step
# that makes the device count cheap later.
message("[2/4] origin x destination x device")
# Indexed placeholders (%1$s, %2$f, %3$f) rather than bare %s / %f: the TAU
# term sits in the SELECT list, ahead of the file glob in the FROM, so plain
# positional formatting binds arguments to whichever placeholder comes first
# in the TEXT. Any edit that moves a term then silently rebinds every argument
# after it. Indexed placeholders make the query editable without that trap.
dbExecute(con, sprintf("
  CREATE OR REPLACE TABLE pair_dev AS
  SELECT
    v.home_county                AS origin_fips,
    v.county                     AS dest_fips,
    v.did,
    COUNT(*)                     AS n_loc_days,
    SUM(v.n_visits)              AS n_visits,
    SUM(v.dwell_pos_sum) + %2$f * SUM(v.n_tau_visits) AS w_dwell,
    MEDIAN(v.d_km)               AS d_km_median_dev
  FROM read_parquet('%1$s') v
  WHERE v.did IN (SELECT did FROM dev)
    AND NOT v.is_home
    AND v.d_km <= %3$f
    AND v.county IS NOT NULL
    AND v.home_county IS NOT NULL
  GROUP BY v.home_county, v.county, v.did
", DEVDAY_GLOB, TAU, D_MAX))
message(sprintf("  pair-device rows: %s",
                format(dbGetQuery(con, "SELECT COUNT(*) AS n FROM pair_dev")$n,
                       big.mark = ",")))

# ---- 3. pair grain ---------------------------------------------------------
message("[3/4] pair grain")
pairs <- dbGetQuery(con, "
  SELECT origin_fips, dest_fips,
         COUNT(*)                AS n_devices,
         SUM(n_loc_days)         AS n_device_days,
         SUM(n_visits)           AS n_visits,
         SUM(w_dwell)            AS w_dwell,
         MEDIAN(d_km_median_dev) AS d_km_median
  FROM pair_dev
  GROUP BY origin_fips, dest_fips") |>
  mutate(origin_fips = str_pad(as.character(origin_fips), 5, "left", "0"),
         dest_fips   = str_pad(as.character(dest_fips),   5, "left", "0"))

# ---- 4. panel totals per origin --------------------------------------------
# n_device_appearances is the denominator for the observed shares: a device
# seen in three destination counties contributes three times, so the shares
# partition APPEARANCES, not devices. Say so when reporting.
#
# n_devices_panel needs its own pass, because it must include devices that
# never left home and so never appear in pair_dev. It is a diagnostic only --
# nothing downstream depends on it.
message("[4/4] panel totals")
panel <- pairs |>
  group_by(origin_fips) |>
  summarise(n_dest_counties      = n(),
            n_device_appearances = sum(n_devices),
            n_device_days        = sum(n_device_days),
            n_visits             = sum(n_visits),
            w_dwell              = sum(w_dwell),
            .groups = "drop")

home_panel <- dbGetQuery(con, sprintf("
  SELECT home_county AS origin_fips, COUNT(*) AS n_devices_panel
  FROM (
    SELECT DISTINCT v.did, v.home_county
    FROM read_parquet('%s') v
    WHERE v.did IN (SELECT did FROM dev) AND v.home_county IS NOT NULL
  )
  GROUP BY home_county", DEVDAY_GLOB)) |>
  mutate(origin_fips = str_pad(as.character(origin_fips), 5, "left", "0"))

panel <- panel |>
  left_join(home_panel, by = "origin_fips") |>
  mutate(appearances_per_device = n_device_appearances / n_devices_panel,
         dwell_min_per_device = w_dwell / n_devices_panel,
         d_max_km = D_MAX, tau_minutes = TAU, home_rule = HOME_RULE,
         built_on = as.character(Sys.time())) |>
  arrange(desc(n_devices_panel))

message(sprintf("origins: %d   pairs: %s", nrow(panel),
                format(nrow(pairs), big.mark = ",")))
print(head(panel, 10))

# ---- checkpoint ------------------------------------------------------------
# The share going to the origin county itself. Expect it high -- most trips
# are short -- and expect it to fall with county size. If it is near 1 for a
# large county the destination county assignment is not working.
# The share going to the origin county itself, under both weightings. They
# will differ, and the size of the gap is how much the choice of estimand is
# worth. If it is large, the trips-versus-all decision is load-bearing and
# should be stated in the write-up rather than left as a default.
chk <- pairs |>
  group_by(origin_fips) |>
  summarise(self_dev   = sum(n_devices[dest_fips == origin_fips]) / sum(n_devices),
            self_dwell = sum(w_dwell[dest_fips == origin_fips]) / sum(w_dwell),
            .groups = "drop")
message("self-share across origins, device-weighted:")
print(summary(chk$self_dev))
message("self-share across origins, dwell-weighted:")
print(summary(chk$self_dwell))
message(sprintf("median |device - dwell| self-share gap: %.4f",
                median(abs(chk$self_dev - chk$self_dwell), na.rm = TRUE)))

# ---- atomic writes ---------------------------------------------------------
for (p in list(list(pairs, OUT_PAIRS), list(panel, OUT_PANEL))) {
  tmp <- paste0(p[[2]], ".tmp"); write_csv(p[[1]], tmp); file.rename(tmp, p[[2]])
  message("wrote ", p[[2]])
}

dbDisconnect(con, shutdown = TRUE)
message("done")

# ---------------------------------------------------------------------------
# RESIDUAL MISMATCH, WORTH ONE SENTENCE IN THE WRITE-UP
#   Even dwell-weighted, this is not Meta. Veraset's dwell is POI dwell only:
#   time in transit, and time at any non-POI location, is unallocated rather
#   than assigned. Meta counts every ping. Time at home is handled by f0 on
#   both sides, so the gap is specifically away-from-home time that Veraset
#   cannot see. Closer, not identical.
#
# EXERCISE
#   The four weightings disagree most where a few devices visit one
#   destination repeatedly. Write the query that measures it: per origin, the
#   total variation between the device-weighted and dwell-weighted destination
#   shares. Then check whether that gap correlates with rurality. If it does,
#   the choice of weighting is not neutral across the covariate the whole
#   comparison turns on.
# ---------------------------------------------------------------------------