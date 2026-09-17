#!/usr/bin/env Rscript
# ============================================================================
# NARROW PIPELINE -- 10 BENCHMARK ORIGINS, END TO END  (June 2025)
# ----------------------------------------------------------------------------
# WHY THIS IS FAST WHERE THE FULL VERSION WAS NOT
#   The full pipeline aggregated every device in the feed and discarded 99.9%
#   of the result afterwards. Here the panel -- devices whose HOME county is
#   one of ten -- is built first, and every scan is filtered to it before the
#   GROUP BY. The aggregate never grows large enough to spill, which is what
#   made other_visits crawl.
#
#   The filter is on the device's home county, never on the destination.
#   Destinations stay unrestricted: a device from Cook County that spends the
#   day in Lake County must still resolve to Lake County, or the off-diagonal
#   disappears by construction.
#
# WHAT COMES OUT
#   flow_sampled_all.csv             one row per (origin, destination, day),
#                                    every device-day of a panel device.
#   flow_sampled_excl_home_only.csv  the same, minus device-days whose only
#                                    rows came from home_visits.
#
#   Each carries BOTH a sampled integer count and an expected (fractional)
#   count. See section 6 -- the expected column is the same quantity without
#   Monte Carlo noise, and the gap between them measures how thin a cell is.
#
# THE ALLOCATION RULE
#   One county per device-day, drawn with probability proportional to dwell.
#   Where every destination has zero dwell, the draw is uniform over the
#   destinations actually visited. The draw is deterministic given
#   (caid_h, day, SEED), so reruns agree and a changed SEED gives a replicate.
#
# THE HOME-DWELL DECISION  -- still the most consequential parameter
#   Pure dwell weighting gives the home county weight 0 whenever any away
#   visit has positive dwell, because home rows carry little or no dwell, and
#   the diagonal collapses. HOME_DWELL_RULE = "residual" instead treats
#   unobserved time as time at home. Section 7 reports the diagonal under both.
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(duckdb)
  library(DBI)
  library(glue)
})

# ---- config ----------------------------------------------------------------
DATA_ROOT  <- "E:/dewey-data"
RAW_DIR    <- file.path(DATA_ROOT, "raw")
KERNEL_DIR <- file.path(DATA_ROOT, "kernel")
TMP_DIR    <- file.path(DATA_ROOT, "tmp")
OUT_DIR  <- file.path(KERNEL_DIR, "narrow10")
CACHE    <- file.path(OUT_DIR, "cache")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE,   showWarnings = FALSE, recursive = TRUE)

home_dir  <- file.path(RAW_DIR, "home_visits")
work_dir  <- file.path(RAW_DIR, "work_visits")
other_dir <- file.path(RAW_DIR, "other_visits")

home_glob <- file.path(home_dir, "*", "*.parquet")

dest_parts <- file.path(CACHE, "dest_parts")
flow_parts <- file.path(CACHE, "flow_parts")
dir.create(dest_parts, showWarnings = FALSE)
dir.create(flow_parts, showWarnings = FALSE)

home_lookup   <- file.path(KERNEL_DIR, "home_lookup.parquet")
BENCHMARK_CSV <- file.path(DATA_ROOT, "external", "meta_movement_dist",
                           "benchmark_counties_gid2.csv")
CENTROID_CSV  <- file.path(DATA_ROOT, "external", "county_centroids.csv")

# --- origin selection -------------------------------------------------------
N_RUCC1  <- 3     # counties drawn from rucc == 1
N_RANDOM <- 7     # drawn from everything else
ORIGIN_SEED <- 26082026L

# --- allocation -------------------------------------------------------------
HOME_DWELL_RULE <- "observed"   # "observed" or "residual"
DAY_BUDGET <- 1440              # 1440 if minimum_dwell is minutes, 86400 if seconds
MAX_KM <- Inf                   # set to 250 to apply the distance cutoff
SEED <- 20250601L

TS_DIVISOR <- 1
REBUILD <- FALSE

stopifnot(HOME_DWELL_RULE %in% c("observed", "residual"))

# Shared SQL fragments. The day expression must match every other stage.
day_sql  <- glue("CAST(CAST(to_timestamp(CAST(utc_timestamp AS DOUBLE) / {TS_DIVISOR}) AS TIMESTAMP) AS DATE)")
fips_sql <- "substr(lpad(CAST(census_block_group AS VARCHAR), 12, '0'), 1, 5)"

# ============================================================================
# 0. CONNECTION
# ============================================================================
con <- dbConnect(duckdb(), dbdir = file.path(CACHE, "narrow10.duckdb"))

dbExecute(con, "PRAGMA threads=8;")
dbExecute(con, "SET memory_limit='48GB';")
dbExecute(con, sprintf("SET temp_directory='%s';", gsub("\\\\", "/", TMP_DIR)))
dbExecute(con, "SET preserve_insertion_order=false;")
dbExecute(con, "PRAGMA enable_progress_bar;")

# ============================================================================
# 1. PICK THE TEN ORIGINS
# ----------------------------------------------------------------------------
# Three from rucc == 1 (the large metro core counties), then seven at random
# from everything else, so the set spans the urbanicity range rather than
# being all anchors. Seeded, and the chosen set is written out, so the same
# ten come back on every rerun and the selection is auditable.
#
# Excel strips leading zeros from FIPS, so the column arrives numeric and has
# to be padded back to five characters.
# ============================================================================
benchmark_raw <- read_csv(BENCHMARK_CSV, show_col_types = FALSE)

message("=== benchmark file columns ===")
print(names(benchmark_raw))

benchmark <- benchmark_raw %>%
  mutate(fips = str_pad(as.character(fips), 5, "left", "0"),
         pop = as.numeric(pop),
         rucc = as.integer(rucc)) %>%
  distinct(fips, .keep_all = TRUE)

message(glue("{nrow(benchmark)} benchmark counties, {sum(benchmark$rucc == 1)} at rucc 1"))

set.seed(ORIGIN_SEED)

rucc1_pool <- benchmark %>% filter(rucc == 1)
rucc1_pick <- rucc1_pool %>% slice_sample(n = min(N_RUCC1, nrow(rucc1_pool)))

rest_pool <- benchmark %>% filter(!fips %in% rucc1_pick$fips)
rest_pick <- rest_pool %>% slice_sample(n = min(N_RANDOM, nrow(rest_pool)))

origins <- bind_rows(rucc1_pick, rest_pick) %>%
  arrange(fips)

message("\n=== selected origins ===")
origins %>%
  select(fips, county_name, state, rucc, urbanicity, pop, stratum) %>%
  print(n = 20)

write_csv(origins, file.path(OUT_DIR, "selected_origins.csv"))

origin_vec  <- origins$fips
origin_list <- str_c("'", origin_vec, "'", collapse = ",")

# ============================================================================
# 2. THE PANEL
# ----------------------------------------------------------------------------
# Devices whose home county is one of the ten. This table is the filter that
# every subsequent scan uses, and it is small -- the whole point.
#
# If home_lookup.parquet is absent it is built here the way the probe script
# does: count each device's home rows by county and keep the modal county.
# That pass is over home_visits only, the smallest of the three sources.
# ============================================================================
if (!file.exists(home_lookup)) {
  message("\n[build] home_lookup.parquet (modal home county per device) ...")
  lookup_sql <- glue("
    COPY (
      WITH home_rows AS (
        SELECT caid, {fips_sql} AS home_fips
        FROM read_parquet('{home_glob}', union_by_name = true)
        WHERE census_block_group IS NOT NULL
          AND CAST(census_block_group AS VARCHAR) <> ''
      ),
      counted AS (
        SELECT caid, home_fips, COUNT(*) AS n,
               ROW_NUMBER() OVER (PARTITION BY caid ORDER BY COUNT(*) DESC, home_fips) AS rk
        FROM home_rows
        GROUP BY caid, home_fips
      )
      SELECT caid, home_fips FROM counted WHERE rk = 1
    ) TO '{home_lookup}' (FORMAT PARQUET, COMPRESSION zstd);")
  dbExecute(con, lookup_sql)
}

dbExecute(con, glue("
  CREATE OR REPLACE TABLE panel AS
  SELECT hash(caid) AS caid_h, caid, home_fips AS origin_fips
  FROM read_parquet('{home_lookup}')
  WHERE home_fips IN ({origin_list});"))

panel_size <- dbGetQuery(con, "
  SELECT origin_fips, COUNT(*) AS devices FROM panel GROUP BY 1 ORDER BY 1;")

message("\n=== panel devices by origin ===")
print(panel_size)
write_csv(panel_size, file.path(OUT_DIR, "panel_size.csv"))

n_panel <- sum(panel_size$devices)
message(glue("Panel total: {format(n_panel, big.mark = ',')} devices"))

if (n_panel == 0) {
  stop("panel is empty -- check that home_lookup's home_fips are zero-padded 5-char strings")
}

# ============================================================================
# 3. DESTINATION ROLLUPS, PANEL ONLY, ONE DAY AT A TIME
# ----------------------------------------------------------------------------
# Grain: (caid_h, day, dest_fips), with visits and dwell kept separate by
# source so downstream can still decide whether a work visit counts as being
# at that county.
#
# The semi-join against `panel` is what makes this cheap. It cannot be pushed
# into the parquet reader, so the rows are still read off disk, but the
# aggregate only ever holds panel devices and stays in memory.
#
# TRY_CAST because minimum_dwell arrives as a string. Zero dwell and missing
# dwell are different things and are counted separately.
# Rows with no usable CBG cannot resolve to a county and are dropped, counted
# in section 7 so the loss is visible.
#
# Atomic write: .tmp then rename, so a killed run leaves no half-written part.
# ============================================================================
if (REBUILD) {
  unlink(list.files(dest_parts, full.names = TRUE))
  unlink(list.files(flow_parts, full.names = TRUE))
}

home_day_dirs  <- list.dirs(home_dir,  recursive = FALSE)
work_day_dirs  <- list.dirs(work_dir,  recursive = FALSE)
other_day_dirs <- list.dirs(other_dir, recursive = FALSE)

day_tags <- sort(unique(c(basename(home_day_dirs),
                          basename(work_day_dirs),
                          basename(other_day_dirs))))

message(glue("\n{length(day_tags)} day partitions"))

for (this_tag in day_tags) {
  this_out <- file.path(dest_parts, str_c(this_tag, ".parquet"))
  this_tmp <- str_c(this_out, ".tmp")
  if (file.exists(this_out)) next
  
  message(glue("  [dest] {this_tag}"))
  
  home_glob_day  <- file.path(home_dir,  this_tag, "*.parquet")
  work_glob_day  <- file.path(work_dir,  this_tag, "*.parquet")
  other_glob_day <- file.path(other_dir, this_tag, "*.parquet")
  
  # Each source contributes its own columns and zeros elsewhere, then one
  # GROUP BY over the stack. Same join-free pattern as the coverage audit.
  this_sql <- glue("
    COPY (
      SELECT caid_h, day, dest_fips,
             SUM(n_visits_home)  AS n_visits_home,
             SUM(n_visits_work)  AS n_visits_work,
             SUM(n_visits_other) AS n_visits_other,
             SUM(dwell_home)     AS dwell_home,
             SUM(dwell_work)     AS dwell_work,
             SUM(dwell_other)    AS dwell_other,
             SUM(dwell_home + dwell_work + dwell_other) AS dwell_total,
             SUM(n_visits_home + n_visits_work + n_visits_other) AS n_visits_total
      FROM (
        SELECT hash(caid) AS caid_h, {day_sql} AS day, {fips_sql} AS dest_fips,
               COUNT(*) AS n_visits_home, 0 AS n_visits_work, 0 AS n_visits_other,
               CAST(0 AS DOUBLE) AS dwell_home,
               CAST(0 AS DOUBLE) AS dwell_work,
               CAST(0 AS DOUBLE) AS dwell_other
        FROM read_parquet('{home_glob_day}', union_by_name = true)
        WHERE census_block_group IS NOT NULL
          AND CAST(census_block_group AS VARCHAR) <> ''
          AND hash(caid) IN (SELECT caid_h FROM panel)
        GROUP BY 1, 2, 3
        UNION ALL
        SELECT hash(caid), {day_sql}, {fips_sql},
               0, COUNT(*), 0,
               CAST(0 AS DOUBLE),
               SUM(COALESCE(TRY_CAST(minimum_dwell AS DOUBLE), 0)),
               CAST(0 AS DOUBLE)
        FROM read_parquet('{work_glob_day}', union_by_name = true)
        WHERE census_block_group IS NOT NULL
          AND CAST(census_block_group AS VARCHAR) <> ''
          AND hash(caid) IN (SELECT caid_h FROM panel)
        GROUP BY 1, 2, 3
        UNION ALL
        SELECT hash(caid), {day_sql}, {fips_sql},
               0, 0, COUNT(*),
               CAST(0 AS DOUBLE),
               CAST(0 AS DOUBLE),
               SUM(COALESCE(TRY_CAST(minimum_dwell AS DOUBLE), 0))
        FROM read_parquet('{other_glob_day}', union_by_name = true)
        WHERE census_block_group IS NOT NULL
          AND CAST(census_block_group AS VARCHAR) <> ''
          AND hash(caid) IN (SELECT caid_h FROM panel)
        GROUP BY 1, 2, 3
      )
      GROUP BY 1, 2, 3
    ) TO '{this_tmp}' (FORMAT PARQUET, COMPRESSION zstd);")
  
  dbExecute(con, this_sql)
  file.rename(this_tmp, this_out)
}

dbExecute(con, glue("
  CREATE OR REPLACE VIEW devday_dest AS
  SELECT * FROM read_parquet('{dest_parts}/*.parquet');"))

dest_rows <- dbGetQuery(con, "
  SELECT COUNT(*) AS dest_rows,
         COUNT(DISTINCT caid_h) AS devices
  FROM devday_dest;")

message("\n=== destination table ===")
print(dest_rows)

# ============================================================================
# 4. DISTANCE LOOKUP
# ----------------------------------------------------------------------------
# Population-weighted county centroids (cenpop2020). Kept as a separate table
# so the 250km rule stays a filter applied at use time rather than baked in.
# Intra-county movement is distance zero, so the cutoff never touches the
# diagonal.
# ============================================================================
centroids <- read_csv(CENTROID_CSV, show_col_types = FALSE) %>%
  mutate(fips = str_pad(as.character(fips), 5, "left", "0")) %>%
  select(fips, lat, lon) %>%
  distinct(fips, .keep_all = TRUE)

dest_seen <- dbGetQuery(con, "SELECT DISTINCT dest_fips FROM devday_dest;") %>%
  as_tibble() %>%
  mutate(dest_fips = str_pad(as.character(dest_fips), 5, "left", "0"))

pairs <- expand_grid(origin_fips = origin_vec, dest_fips = dest_seen$dest_fips) %>%
  left_join(centroids, by = c("origin_fips" = "fips")) %>%
  rename(o_lat = lat, o_lon = lon) %>%
  left_join(centroids, by = c("dest_fips" = "fips")) %>%
  rename(d_lat = lat, d_lon = lon)

earth_km <- 6371
pair_distance <- pairs %>%
  mutate(o_lat_r = o_lat * pi / 180,
         o_lon_r = o_lon * pi / 180,
         d_lat_r = d_lat * pi / 180,
         d_lon_r = d_lon * pi / 180,
         dlat = d_lat_r - o_lat_r,
         dlon = d_lon_r - o_lon_r,
         a = sin(dlat / 2)^2 + cos(o_lat_r) * cos(d_lat_r) * sin(dlon / 2)^2,
         dist_km = 2 * earth_km * asin(sqrt(a))) %>%
  select(origin_fips, dest_fips, dist_km)

write_csv(pair_distance, file.path(OUT_DIR, "county_pair_distance.csv"))
dbWriteTable(con, "pair_distance", pair_distance, overwrite = TRUE)

n_no_centroid <- sum(is.na(pair_distance$dist_km))
message(glue("  {nrow(pair_distance)} origin-dest pairs, {n_no_centroid} with no centroid"))

use_distance <- is.finite(MAX_KM)
if (use_distance) {
  message(glue("  distance cutoff active at {MAX_KM} km"))
}

# ============================================================================
# 5. DWELL UNIT CHECK
# ----------------------------------------------------------------------------
# DAY_BUDGET only matters under the residual rule, but it has to be right.
# A median around 30-120 means minutes; in the thousands means seconds.
# ============================================================================
dwell_peek <- dbGetQuery(con, "
  SELECT MEDIAN(dwell_total) AS median_dwell,
         quantile_cont(dwell_total, 0.9) AS p90_dwell,
         MAX(dwell_total) AS max_dwell
  FROM devday_dest WHERE dwell_total > 0;")

message("\n=== dwell_total, positive values only ===")
print(dwell_peek)
message(glue("  DAY_BUDGET is {DAY_BUDGET}; confirm the unit matches"))

# ============================================================================
# 6. ALLOCATION -- SAMPLED AND EXPECTED, IN ONE PASS
# ----------------------------------------------------------------------------
# Because each device-day is an independent draw with probability
# proportional to dwell, a device-day's EXPECTED contribution to county j is
# exactly its weight share w_j / sum(w). Summing those fractions gives the
# expectation of the sampled table with no Monte Carlo noise -- a strictly
# better estimate of the same quantity for building a prior.
#
# Both are emitted. Use the expected columns for the prior; use the sampled
# columns if something downstream needs integer device-day counts. Where they
# disagree materially, the cell is thin and the apparent structure is noise.
#
#   totals    per-device-day sums via window functions, which attach a group
#             value to every row without collapsing it.
#   weighted  the all-zero-dwell fallback is tested FIRST, so a device-day
#             with no dwell anywhere goes uniform instead of the residual rule
#             handing everything to home.
#   cum       running weight in a fixed destination order, plus a uniform draw
#             hashed from (caid_h, day, SEED) -- deterministic and constant
#             within the device-day.
#   picked    first row whose running total reaches u x total. Inverse-CDF
#             sampling: each destination owns an interval proportional to its
#             weight.
# ============================================================================
if (HOME_DWELL_RULE == "residual") {
  home_weight_sql <- glue(
    "WHEN dest_fips = origin_fips
          THEN GREATEST(dwell_total, {DAY_BUDGET} - away_dwell, 0)")
} else {
  home_weight_sql <- ""
}

if (use_distance) {
  distance_clause <- glue("
    AND (d.dest_fips = p.origin_fips
         OR pd.dist_km IS NULL
         OR pd.dist_km <= {MAX_KM})")
} else {
  distance_clause <- ""
}

for (this_tag in day_tags) {
  this_dest <- file.path(dest_parts, str_c(this_tag, ".parquet"))
  this_out  <- file.path(flow_parts, str_c(this_tag, ".parquet"))
  this_tmp  <- str_c(this_out, ".tmp")
  if (!file.exists(this_dest)) next
  if (file.exists(this_out)) next
  
  message(glue("  [flow] {this_tag}"))
  
  this_sql <- glue("
    WITH joined AS (
      SELECT d.caid_h, p.origin_fips, d.dest_fips,
             d.dwell_total, d.n_visits_work, d.n_visits_other
      FROM read_parquet('{this_dest}') d
      JOIN panel p USING (caid_h)
      LEFT JOIN pair_distance pd
        ON pd.origin_fips = p.origin_fips AND pd.dest_fips = d.dest_fips
      WHERE d.dest_fips IS NOT NULL
        {distance_clause}
    ),
    totals AS (
      SELECT *,
             SUM(dwell_total) OVER (PARTITION BY caid_h) AS dd_dwell,
             SUM(CASE WHEN dest_fips <> origin_fips THEN dwell_total ELSE 0 END)
               OVER (PARTITION BY caid_h) AS away_dwell,
             MAX(CASE WHEN n_visits_work > 0 OR n_visits_other > 0 THEN 1 ELSE 0 END)
               OVER (PARTITION BY caid_h) AS has_away_source
      FROM joined
    ),
    weighted AS (
      SELECT *,
             CASE
               WHEN dd_dwell <= 0 THEN 1.0
               {home_weight_sql}
               ELSE dwell_total
             END AS w
      FROM totals
    ),
    cum AS (
      SELECT *,
             SUM(w) OVER (PARTITION BY caid_h ORDER BY dest_fips
                          ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW) AS cum_w,
             SUM(w) OVER (PARTITION BY caid_h) AS tot_w,
             (hash(CAST(caid_h AS VARCHAR) || '-{this_tag}-{SEED}') % 1000000) / 1000000.0 AS u
      FROM weighted
    ),
    expected AS (
      SELECT origin_fips, dest_fips,
             SUM(w / tot_w) AS exp_all,
             SUM(CASE WHEN has_away_source = 1 THEN w / tot_w ELSE 0 END) AS exp_excl
      FROM cum
      GROUP BY 1, 2
    ),
    ranked AS (
      SELECT caid_h, origin_fips, dest_fips, has_away_source, dd_dwell,
             ROW_NUMBER() OVER (PARTITION BY caid_h ORDER BY cum_w) AS rk
      FROM cum
      WHERE cum_w >= u * tot_w
    ),
    sampled AS (
      SELECT origin_fips, dest_fips,
             COUNT(*) AS flow_all,
             SUM(CASE WHEN has_away_source = 1 THEN 1 ELSE 0 END) AS flow_excl,
             SUM(CASE WHEN dd_dwell <= 0 THEN 1 ELSE 0 END) AS n_uniform_rule
      FROM ranked WHERE rk = 1
      GROUP BY 1, 2
    )
    SELECT origin_fips, dest_fips,
           DATE '{this_tag}' AS day,
           COALESCE(flow_all, 0)       AS flow_all,
           COALESCE(flow_excl, 0)      AS flow_excl,
           COALESCE(exp_all, 0)        AS exp_all,
           COALESCE(exp_excl, 0)       AS exp_excl,
           COALESCE(n_uniform_rule, 0) AS n_uniform_rule
    FROM sampled FULL OUTER JOIN expected USING (origin_fips, dest_fips)")
  
  dbExecute(con, glue("COPY ({this_sql}) TO '{this_tmp}' (FORMAT PARQUET, COMPRESSION zstd);"))
  file.rename(this_tmp, this_out)
}

flow_raw <- dbGetQuery(con, glue("
  SELECT * FROM read_parquet('{flow_parts}/*.parquet');"))

flow_tbl <- flow_raw %>%
  as_tibble() %>%
  mutate(origin_fips = str_pad(as.character(origin_fips), 5, "left", "0"),
         dest_fips   = str_pad(as.character(dest_fips), 5, "left", "0"),
         day = as.Date(day))

message(glue("\n{format(nrow(flow_tbl), big.mark = ',')} origin-dest-day cells"))

# ============================================================================
# 7. THE TWO OUTPUT FILES
# ----------------------------------------------------------------------------
# Shares are within origin-day, so each origin-day sums to 1 and the two files
# are comparable despite different denominators. `pop` rides along so the
# scaling step needs no second lookup.
# ============================================================================
origin_pop <- origins %>% select(origin_fips = fips, county_name, state, rucc, pop)

flow_all_out <- flow_tbl %>%
  filter(flow_all > 0 | exp_all > 0) %>%
  group_by(origin_fips, day) %>%
  mutate(dev_days_origin_day = sum(flow_all),
         share_sampled = flow_all / sum(flow_all),
         share_expected = exp_all / sum(exp_all)) %>%
  ungroup() %>%
  left_join(origin_pop, by = "origin_fips") %>%
  select(origin_fips, county_name, state, rucc, pop, dest_fips, day,
         flow = flow_all, flow_expected = exp_all,
         dev_days_origin_day, share_sampled, share_expected) %>%
  arrange(origin_fips, day, desc(flow_expected))

write_csv(flow_all_out, file.path(OUT_DIR, "flow_sampled_all.csv"))

flow_excl_out <- flow_tbl %>%
  filter(flow_excl > 0 | exp_excl > 0) %>%
  group_by(origin_fips, day) %>%
  mutate(dev_days_origin_day = sum(flow_excl),
         share_sampled = flow_excl / sum(flow_excl),
         share_expected = exp_excl / sum(exp_excl)) %>%
  ungroup() %>%
  left_join(origin_pop, by = "origin_fips") %>%
  select(origin_fips, county_name, state, rucc, pop, dest_fips, day,
         flow = flow_excl, flow_expected = exp_excl,
         dev_days_origin_day, share_sampled, share_expected) %>%
  arrange(origin_fips, day, desc(flow_expected))

write_csv(flow_excl_out, file.path(OUT_DIR, "flow_sampled_excl_home_only.csv"))

message(glue("Wrote flow_sampled_all.csv ({format(nrow(flow_all_out), big.mark = ',')} rows)"))
message(glue("Wrote flow_sampled_excl_home_only.csv ({format(nrow(flow_excl_out), big.mark = ',')} rows)"))

# ============================================================================
# 8. QA
# ============================================================================
totals <- flow_tbl %>%
  summarise(dev_days_all = sum(flow_all),
            dev_days_excl = sum(flow_excl),
            diagonal_all = sum(flow_all[origin_fips == dest_fips]),
            diagonal_excl = sum(flow_excl[origin_fips == dest_fips]),
            uniform_rule = sum(n_uniform_rule)) %>%
  mutate(pct_home_only_dropped = 100 * (dev_days_all - dev_days_excl) / dev_days_all,
         pct_diagonal_all = 100 * diagonal_all / dev_days_all,
         pct_diagonal_excl = 100 * diagonal_excl / dev_days_excl,
         pct_uniform_rule = 100 * uniform_rule / dev_days_all)

message("\n=== headline ===")
print(totals)
write_csv(totals, file.path(OUT_DIR, "qa_headline.csv"))

by_origin <- flow_tbl %>%
  group_by(origin_fips) %>%
  summarise(dev_days = sum(flow_all),
            n_dest = n_distinct(dest_fips),
            pct_diagonal = 100 * sum(flow_all[origin_fips == dest_fips]) / sum(flow_all),
            .groups = "drop") %>%
  left_join(origin_pop, by = "origin_fips") %>%
  mutate(dev_days_per_1k_pop = 1000 * dev_days / pop) %>%
  arrange(desc(dev_days))

message("\n=== by origin ===")
print(by_origin, n = 20)
write_csv(by_origin, file.path(OUT_DIR, "qa_by_origin.csv"))

diag_by_day <- flow_tbl %>%
  group_by(day) %>%
  summarise(dev_days = sum(flow_all),
            diagonal = sum(flow_all[origin_fips == dest_fips]),
            .groups = "drop") %>%
  mutate(dow = wday(day, label = TRUE),
         pct_diagonal = 100 * diagonal / dev_days)

message("\n=== diagonal share by day ===")
print(diag_by_day, n = 40)
write_csv(diag_by_day, file.path(OUT_DIR, "qa_diagonal_by_day.csv"))

# Sampling noise: how far the integer draw strays from its own expectation.
# Large gaps mark cells too thin to interpret.
noise <- flow_all_out %>%
  filter(flow_expected > 0) %>%
  mutate(abs_gap = abs(flow - flow_expected),
         rel_gap = abs_gap / flow_expected)

noise_summary <- noise %>%
  mutate(size_bucket = cut(flow_expected,
                           breaks = c(0, 1, 5, 20, 100, Inf),
                           labels = c("<1", "1-5", "5-20", "20-100", "100+"))) %>%
  group_by(size_bucket) %>%
  summarise(median_rel_gap = median(rel_gap),
            cells = n(),
            .groups = "drop")

message("\n=== sampled vs expected, by cell size ===")
print(noise_summary)
write_csv(noise_summary, file.path(OUT_DIR, "qa_sampling_noise.csv"))

message(glue("\nDone. Outputs in {OUT_DIR}"))
dbDisconnect(con, shutdown = TRUE)