# ---------------------------------------------------------------------------
# kernel_40counties.R
#
# Continues from kernel_walkthrough.Rmd. Two parts:
#
#   PART A  rebuild devday_parts as v2, carrying raw dwell components so that
#           tau, the home/non-home split, and the band edges all become
#           post-hoc choices. One 30-day loop over visit_parts (which already
#           exist on disk) — no rescan of the raw Veraset files.
#
#   PART B  county profiles for the full ~40 benchmark counties, under the
#           trips-outside-the-home estimand.
#
# Nothing here rescans home_visits/work_visits/other_visits. home_assignment,
# visit_parts and devday_parts were all built nationally with no origin
# filter, so the extra 30 counties need no new scanning at all.
# ---------------------------------------------------------------------------

library(DBI)
library(duckdb)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(glue)


# --- Config ----------------------------------------------------------------

ROOT    <- "E:/dewey-june2025"
TMP_DIR <- "E:/duckdb-tmp"
OUT     <- file.path(ROOT, "kernel")

MEM_LIMIT <- "400GB"   # measured: memory is not the constraint on this box
N_THREADS <- 32        # measured: peaks at 32, degrades past 64

dir.create(file.path(OUT, "devday_parts_v2"), showWarnings = FALSE, recursive = TRUE)

con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
dbExecute(con, glue("SET temp_directory = '{TMP_DIR}'"))
dbExecute(con, glue("SET memory_limit = '{MEM_LIMIT}'"))
dbExecute(con, glue("SET threads = {N_THREADS}"))
dbExecute(con, "SET preserve_insertion_order = false")

list_visit_days <- function() {
  list.files(file.path(OUT, "visit_parts"), pattern = "\\.parquet$") |>
    str_remove("\\.parquet$") |>
    sort()
}


# ===========================================================================
# PART A — rebuild device-day parts, with tau left free
#
# Three changes from the v1 devday_parts:
#
#  1. Carry `dwell_pos_sum` (sum of strictly positive dwell) and `n_tau_visits`
#     (count of visits with zero/null/negative dwell) instead of a pre-baked
#     w_dwell. Then for any tau:
#            w_dwell(tau) = dwell_pos_sum + tau * n_tau_visits
#     which makes the tau sensitivity a re-read rather than a re-scan.
#
#  2. Group by `is_home` as well as `cbg`. In v1 the grouping merged a home
#     row and a non-home row in the same block group into one location with
#     is_home = 1, which silently deleted trips made inside the home block
#     group. For a trips-outside-the-home estimand that is exactly the wrong
#     thing to lose.
#
#  3. Carry `home_county` from home_assignment, so a device-day does not need
#     to contain a home visit in order to have a home county. v1 derived the
#     home county from is_home rows per device-day, which dropped every
#     device-day where the person was away overnight — i.e. the long-distance
#     travel we most want to measure.
#
# No normalisation happens here. All of it moves to Part B, where it is cheap
# because we have already filtered to 40 counties.
# ===========================================================================

devday_v2_sql <- function(day) {
  glue("
    WITH v AS (
      SELECT
        p.did,
        p.day,
        p.cbg,
        p.county,
        p.loc_src,
        p.is_home,
        p.d_km,
        p.dwell,
        substr(h.home_cbg, 1, 5) AS home_county,
        h.home_cbg
      FROM read_parquet('{OUT}/visit_parts/{day}.parquet') p
      JOIN read_parquet('{OUT}/home_assignment.parquet')   h ON p.did = h.did
    )

    SELECT
      did,
      day,
      home_county,
      cbg,
      county,
      loc_src,
      is_home,
      MIN(d_km)                                             AS d_km,
      CASE WHEN cbg = MIN(home_cbg) THEN 1 ELSE 0 END       AS is_home_cbg,
      SUM(CASE WHEN dwell > 0 THEN dwell ELSE 0 END)        AS dwell_pos_sum,
      SUM(CASE WHEN dwell IS NULL OR dwell <= 0
               THEN 1 ELSE 0 END)                           AS n_tau_visits,
      COUNT(*)                                              AS n_visits
    FROM v
    GROUP BY did, day, home_county, cbg, county, loc_src, is_home
  ")
}

for (day in list_visit_days()) {
  out_path <- file.path(OUT, "devday_parts_v2", glue("{day}.parquet"))
  if (file.exists(out_path)) next
  
  message(glue("devday v2 {day} ..."))
  tmp_path <- paste0(out_path, ".tmp")
  
  dbExecute(con, glue("COPY ({devday_v2_sql(day)}) TO '{tmp_path}' (FORMAT PARQUET)"))
  file.rename(tmp_path, out_path)
}


# --- Checkpoint A ----------------------------------------------------------
# Compare v1 and v2 row counts for one day. v2 should have MORE rows, because
# splitting on is_home un-merges home/non-home locations that v1 collapsed.
# If the counts are identical, change 2 above did nothing and that is worth
# understanding.

probe_day <- list_visit_days()[3]

bind_rows(
  dbGetQuery(con, glue("
    SELECT 'v1' AS version, COUNT(*) AS n_rows
    FROM read_parquet('{OUT}/devday_parts/{probe_day}.parquet')
  ")),
  dbGetQuery(con, glue("
    SELECT 'v2' AS version, COUNT(*) AS n_rows
    FROM read_parquet('{OUT}/devday_parts_v2/{probe_day}.parquet')
  "))
) |>
  as_tibble()

# How many device-days have no home visit at all? These are the ones v1 was
# silently dropping. If the share is material, that alone explains some of the
# tail weirdness.
dbGetQuery(con, glue("
  WITH dd AS (
    SELECT did, day, MAX(is_home) AS has_home
    FROM read_parquet('{OUT}/devday_parts_v2/{probe_day}.parquet')
    GROUP BY did, day
  )
  SELECT
    COUNT(*)                                        AS n_devdays,
    SUM(CASE WHEN has_home = 0 THEN 1 ELSE 0 END)   AS n_no_home_visit
  FROM dd
")) |>
  as_tibble() |>
  mutate(share_no_home = n_no_home_visit / n_devdays)


# ===========================================================================
# PART B — the 40 benchmark counties
# ===========================================================================

# --- B1. Origin list -------------------------------------------------------
# Excel strips leading zeros from `fips`, so pad it back before anything else.

benchmark <- read_csv("E:/meta_movement_dist/benchmark_counties_gid2.csv",
                      show_col_types = FALSE) |>
  mutate(fips = str_pad(as.character(fips), 5, pad = "0"))

origins <- benchmark |>
  select(fips, county_name, state, rucc, urbanicity, pop, stratum, clean, gid_2)

dbWriteTable(con, "origins", select(origins, fips), overwrite = TRUE)

origins |>
  count(clean, name = "n_counties")


# --- B2. Settings, all post-hoc -------------------------------------------

TAU_MINUTES <- 60
BAND_EDGE   <- c(0, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000, Inf)

# Which visits count. "trips" is your estimand: locations outside the home,
# each device-day normalised across its own trips only.
#   "trips"     -> is_home = 0
#   "all"       -> everything, home included (the v1 behaviour, for comparison)
SCOPE <- "all"

# How a trip is weighted once it is in scope.
#   "dwell" -> minutes, with tau standing in for single-ping visits
#   "trip"  -> one per location, duration ignored
WEIGHTING <- "dwell"

band_case <- BAND_EDGE |>
  head(-1) |>
  seq_along() |>
  map_chr(function(i) {
    lo <- BAND_EDGE[i]
    hi <- BAND_EDGE[i + 1]
    if (is.infinite(hi)) {
      glue("WHEN d_km >= {lo} THEN {i}")
    } else {
      glue("WHEN d_km >= {lo} AND d_km < {hi} THEN {i}")
    }
  }) |>
  str_c(collapse = "\n        ")

scope_filter <- if (SCOPE == "trips") "is_home = 0" else "1 = 1"

weight_expr <- if (WEIGHTING == "dwell") {
  glue("dwell_pos_sum + {TAU_MINUTES} * n_tau_visits")
} else {
  "1"
}


# --- B3. The three-level normalisation, in SQL ----------------------------
#
# Done in the database because the national v2 table is ~1B rows and pulling
# it into R to filter afterwards is what made the v1 §6 chunk slow. The county
# filter is a semi-join pushed inside the scan; only the 40-county result
# comes back to R, which is a few hundred rows.
#
# The three levels, in order:
#   dd_share   weights normalised within device-day  (sums to 1 per device-day)
#   dev_band   device-days averaged within device    (sums to 1 per device)
#   county     devices averaged within county        (sums to 1 per county)
#
# COUNT(DISTINCT did) appears once, at county level, where there are 40 groups.
# That is safe — the rule is about distinct counts across hundreds of millions
# of groups, not about distinct counts as such.

profile_sql <- glue("
  WITH base AS (
    SELECT
      d.did,
      d.day,
      d.home_county,
      CASE
        {band_case}
        ELSE {length(BAND_EDGE)}
      END                              AS band_idx,
      {weight_expr}                    AS w
    FROM read_parquet('{OUT}/devday_parts_v2/*.parquet') d
    JOIN origins o ON d.home_county = o.fips
    WHERE {scope_filter}
      AND d.d_km IS NOT NULL
  ),

  dd_tot AS (
    SELECT did, day, SUM(w) AS tot
    FROM base
    GROUP BY did, day
  ),

  dd_band AS (
    SELECT
      b.did,
      b.day,
      b.home_county,
      b.band_idx,
      SUM(b.w) / MAX(t.tot) AS s
    FROM base b
    JOIN dd_tot t ON b.did = t.did AND b.day = t.day
    WHERE t.tot > 0
    GROUP BY b.did, b.day, b.home_county, b.band_idx
  ),

  dev_days AS (
    SELECT did, COUNT(*) AS n_days
    FROM (SELECT DISTINCT did, day FROM dd_band)
    GROUP BY did
  ),

  dev_band AS (
    SELECT
      b.did,
      b.home_county,
      b.band_idx,
      SUM(b.s) / MAX(d.n_days) AS s
    FROM dd_band b
    JOIN dev_days d ON b.did = d.did
    GROUP BY b.did, b.home_county, b.band_idx
  ),

  county_dev AS (
    SELECT home_county, COUNT(DISTINCT did) AS n_devices
    FROM dev_band
    GROUP BY home_county
  )

  SELECT
    b.home_county,
    b.band_idx,
    SUM(b.s) / MAX(c.n_devices) AS share,
    MAX(c.n_devices)            AS n_devices
  FROM dev_band b
  JOIN county_dev c ON b.home_county = c.home_county
  GROUP BY b.home_county, b.band_idx
")

county_profile <- dbGetQuery(con, profile_sql) |>
  as_tibble()


# --- B4. Complete the grid and label the bands ----------------------------
# complete() matters: a county with no mass in a band must appear as 0, not be
# absent. A missing band makes K and edges inconsistent and Stan fails with
# "failed to create the sampler", which surfaces later as log(NULL).

band_labels <- tibble(
  band_idx  = seq_len(length(BAND_EDGE) - 1),
  edge_low  = head(BAND_EDGE, -1),
  edge_high = tail(BAND_EDGE, -1)
) |>
  mutate(band = glue("[{edge_low},{edge_high})"))

county_profile_full <- county_profile |>
  complete(home_county, band_idx = band_labels$band_idx,
           fill = list(share = 0)) |>
  group_by(home_county) |>
  mutate(n_devices = max(n_devices, na.rm = TRUE)) |>
  ungroup() |>
  left_join(band_labels, by = "band_idx") |>
  arrange(home_county, band_idx)


# --- Checkpoint B ---------------------------------------------------------

county_profile_full |>
  group_by(home_county) |>
  summarise(total = sum(share), n_devices = first(n_devices), .groups = "drop") |>
  summarise(
    n_counties  = n(),
    min_total   = min(total),
    max_total   = max(total),
    min_devices = min(n_devices),
    med_devices = median(n_devices)
  )

# Each county's shares must sum to 1. min_devices is the county that will
# dominate your uncertainty — check it is not so small the fit is meaningless.

county_profile_full |>
  filter(home_county %in% c("06037", "30109", "19165")) |>
  select(home_county, band, share) |>
  pivot_wider(names_from = home_county, values_from = share) |>
  print(n = 20)


# --- B5. Write ------------------------------------------------------------

RUN_TAG <- glue("{SCOPE}_{WEIGHTING}_tau{TAU_MINUTES}")

profile_out <- county_profile_full |>
  left_join(select(origins, fips, county_name, state, rucc, pop, clean),
            by = c("home_county" = "fips")) |>
  mutate(
    run_tag     = RUN_TAG,
    scope       = SCOPE,
    weighting   = WEIGHTING,
    tau_minutes = TAU_MINUTES,
    built_on    = as.character(Sys.Date())
  ) |>
  relocate(run_tag, home_county, county_name, state, rucc, band,
           edge_low, edge_high, share, n_devices)

out_path <- file.path(OUT, glue("county_profile_40_{RUN_TAG}.csv"))
tmp_path <- paste0(out_path, ".tmp")

write_csv(profile_out, tmp_path)
file.rename(tmp_path, out_path)

out_path


# ===========================================================================
# PART C — the sensitivity grid
#
# Everything above is one query against v2, so the full grid is just a loop.
# Rebuild profile_sql inside it by re-sourcing B2/B3 with different settings,
# or wrap B2-B4 in a function and map over the grid.
#
# Suggested grid, in priority order:
#   tau        1, 15, 60    -- is the tail a tau artefact?
#   scope      trips, all   -- how much does excluding home change the shape?
#   weighting  dwell, trip  -- does duration matter once home is out?
#
# tau is the one that matters most for the Meta offset, because it moves rural
# profiles more than urban ones (a rural device-day is often home plus one
# stop, where tau sets the whole split). That means tau is confounded with the
# rural/urban gradient you are about to fit the offset against — so the range
# needs to be in the write-up, not just checked and forgotten.
# ===========================================================================

dbDisconnect(con, shutdown = TRUE)