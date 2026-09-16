### ---------------------------------------------------------------------
### Split-half profiles, all replicate splits in one pass
###
### Rebuilds the county band profile from disjoint random halves of the
### device panel, once per salt in REP_SALTS. Fitting each half
### independently and correlating across counties gives a reliability
### estimate that needs no distributional assumption.
###
### Why more than one split. A single split gives one correlation, but which
### devices landed in which half is itself a random draw. Repeating with
### different salts turns "reliability is 0.94" into "0.94, range 0.92-0.96
### across ten splits", which is a much harder number to argue with.
###
### The split is on the DEVICE, not the device-day. Devices are the
### independent unit; splitting on device-days would put the same person in
### both halves and inflate reliability.
###
### PERFORMANCE. An earlier version ran the whole scan once per salt, so ten
### replicates meant ten passes over devday_parts_v2. This version does one
### pass:
###
### 1. Scan, join and band ONCE, aggregating to home_county x did x
### band_idx. That collapses device-days to at most ten rows per
### device before anything else happens.
### 2. CROSS JOIN the small aggregate against the salt list, so all
### half-assignments are derived in one go.
### 3. Device counts come from a DISTINCT over the already-device-grain
### table rather than COUNT(DISTINCT did) over the raw scan. The
### latter is what stalled other_visits with 32 GB spilled.
###
### NOTE ON n_devices. Computed at salt x county x half grain in its own CTE
### and joined back on, NOT inside the band-level GROUP BY. A per-band count
### gives every band a different device count, breaks the balance check, and
### understates the effective sample size the fitting stage relies on.
###
### The aggregation here must mirror whatever built the full profile. If
### that used device-averaging rather than summed dwell weight, change THIS
### file to match - the halves have to be built the same way as the whole or
### they are not comparable. The pooled-vs-full check at the bottom tests
### exactly that.
### ---------------------------------------------------------------------

library(tidyverse)
library(DBI)
library(duckdb)

### CONSTANTS ###

DEVDAY_DIR <- "E:/dewey-june2025/kernel/devday_parts_v2"
HOMERES <- "E:/dewey-june2025/kernel/trackB/device_home_resolution.parquet"
PROFILE_FILE <- "E:/dewey-june2025/kernel/trackB/county_profile_cbghome_national_all_dwell_tau60.csv"
OUT_FILE <- "E:/dewey-june2025/kernel/trackB/county_profile_cbghome_splithalf_tau60.csv"

TAU <- 60 # minutes credited to a single-ping visit
N_THREADS <- 32
MEM_LIMIT <- "400GB"

# One salt per replicate split. hash(did || salt) is deterministic, so reruns
# reproduce the same splits without needing a seed. Start with two to check
# the balance and timing, then widen to the full set.
REP_SALTS <- c("s01", "s02", "s03", "s04", "s05",
               "s06", "s07", "s08", "s09", "s10")

BAND_EDGES <- tibble(
  band_idx = 1:10,
  edge_low = c(0, 1, 2.5, 5, 10, 25, 50, 100, 250, 500),
  edge_high = c(1, 2.5, 5, 10, 25, 50, 100, 250, 500, Inf)
)

### QUERY ###
#
# Reading it top to bottom:
#
# clean_devices devices whose home was resolved at CBG level. Geohash-5
# homes carry a median 3.5 km displacement on rows whose true
# displacement is zero, contaminating every band under 10 km.
# Not optional for sub-10 km work.
#
# rows device-day x location for clean devices, with the dwell
# weight. No salt here - the split has not happened yet.
#
# dev_band THE COLLAPSE. Sums weight to one row per device per band.
# Everything downstream operates on this much smaller table,
# which is what makes many salts affordable.
#
# salts the salt list as a one-column relation, so it can be
# cross-joined. Interpolated into the query text rather than
# bound: DuckDB's R driver reads a multi-element vector as
# multiple bind values, which fails against the scalar
# parameters. The salts are script constants, not input.
#
# assigned dev_band fanned out across salts, each row carrying its
# half assignment for that salt.
#
# dev_counts device count at salt x county x half grain. The inner
# DISTINCT runs over a table that is already device-grain per
# band, so it is cheap in a way COUNT(DISTINCT) over the raw
# scan is not.
#
# band_weights summed weight per band, divided by the salt x county x half
# total so each half's shares sum to 1 on their own.

sql_template <- "
WITH clean_devices AS (
SELECT did
FROM read_parquet(?)
WHERE cbg_home = 1
),

rows AS (
SELECT
d.home_county,
d.did,
d.d_km,
d.dwell_pos_sum + ? * d.n_tau_visits AS w
FROM read_parquet(?) AS d
INNER JOIN clean_devices AS c USING (did)
WHERE d.home_county IS NOT NULL
AND d.d_km IS NOT NULL
),

dev_band AS (
SELECT
home_county,
did,
CASE
WHEN d_km < 1.0 THEN 1
WHEN d_km < 2.5 THEN 2
WHEN d_km < 5.0 THEN 3
WHEN d_km < 10.0 THEN 4
WHEN d_km < 25.0 THEN 5
WHEN d_km < 50.0 THEN 6
WHEN d_km < 100.0 THEN 7
WHEN d_km < 250.0 THEN 8
WHEN d_km < 500.0 THEN 9
ELSE 10
END AS band_idx,
SUM(w) AS w
FROM rows
GROUP BY home_county, did, band_idx
),

salts AS (
SELECT unnest(__SALTS__) AS salt
),

assigned AS (
SELECT
b.home_county,
b.did,
b.band_idx,
b.w,
s.salt,
hash(CAST(b.did AS VARCHAR) || s.salt) % 2 AS half
FROM dev_band AS b
CROSS JOIN salts AS s
),

dev_counts AS (
SELECT salt, home_county, half, COUNT(*) AS n_devices
FROM (SELECT DISTINCT salt, home_county, half, did FROM assigned)
GROUP BY salt, home_county, half
),

band_weights AS (
SELECT
salt,
home_county,
half,
band_idx,
SUM(w) AS w_band,
SUM(w) / SUM(SUM(w)) OVER (PARTITION BY salt, home_county, half) AS share
FROM assigned
GROUP BY salt, home_county, half, band_idx
)

SELECT b.salt, b.home_county, b.half, b.band_idx, b.w_band, b.share, c.n_devices
FROM band_weights AS b
JOIN dev_counts AS c USING (salt, home_county, half)
"

salt_literal <- paste0("['", paste(REP_SALTS, collapse = "','"), "']")
sql <- sub("__SALTS__", salt_literal, sql_template, fixed = TRUE)

### RUN ###

con <- dbConnect(duckdb(), config = list(threads = as.character(N_THREADS),
                                         memory_limit = MEM_LIMIT))

dbExecute(con, "SET preserve_insertion_order=false")

message("scanning ", DEVDAY_DIR, " for ", length(REP_SALTS), " splits ...")
t0 <- Sys.time()

profile <- dbGetQuery(
  con, sql,
  params = list(HOMERES, TAU, file.path(DEVDAY_DIR, "*.parquet"))
)

message("query returned in ",
        signif(as.numeric(difftime(Sys.time(), t0, units = "mins")), 3), " min")

dbDisconnect(con, shutdown = TRUE)

# Sorting happens here rather than in SQL: the result is small, and an
# ORDER BY in the query would force a full sort of the aggregate for nothing.
profile <- profile |>
  as_tibble() |>
  mutate(home_county = str_pad(as.character(home_county), 5, pad = "0"),
         half = as.integer(half),
         n_devices = as.numeric(n_devices)) |>
  left_join(BAND_EDGES, by = "band_idx") |>
  arrange(salt, home_county, half, band_idx)

### CHECKS ###

# 0. Every salt came back.
message("salts returned: ", n_distinct(profile$salt), " of ", length(REP_SALTS))
stopifnot(n_distinct(profile$salt) == length(REP_SALTS))

# 1. Shares sum to 1 within every salt x county x half.
bad_shares <- profile |>
  group_by(salt, home_county, half) |>
  summarise(total = sum(share), .groups = "drop") |>
  filter(abs(total - 1) > 1e-8)
stopifnot(nrow(bad_shares) == 0)

# 2. n_devices constant within salt x county x half. If this fails, the
# device count leaked into the band-level GROUP BY.
bad_counts <- profile |>
  distinct(salt, home_county, half, n_devices) |>
  count(salt, home_county, half) |>
  filter(n > 1)
stopifnot(nrow(bad_counts) == 0)

# 3. The salts must actually differ. If the concatenation is not changing
# the hash - a did cast producing something unexpected, say - every
# replicate is the same split and the reported range is meaninglessly
# tight.
salt_variation <- profile |>
  filter(half == 0, band_idx == 1) |>
  group_by(home_county) |>
  summarise(n_distinct_counts = n_distinct(n_devices), .groups = "drop")
message("counties where all salts gave an identical half-0 device count: ",
        sum(salt_variation$n_distinct_counts == 1), " of ",
        nrow(salt_variation), " (should be few)")

# 4. Halves should be close to the same size. Ratio computed before any
# reshaping, so nothing depends on pivot_wider's column naming.
balance <- profile |>
  distinct(salt, home_county, half, n_devices) |>
  group_by(salt, home_county) |>
  summarise(n_total = sum(n_devices),
            n_half_0 = sum(n_devices[half == 0]),
            n_half_1 = sum(n_devices[half == 1]),
            ratio = n_half_0 / n_total,
            .groups = "drop")

message("balance across ", length(REP_SALTS), " splits:")
message(" ratio median ", signif(median(balance$ratio, na.rm = TRUE), 4),
        ", 1st pct ", signif(quantile(balance$ratio, 0.01, na.rm = TRUE), 3),
        ", 99th pct ", signif(quantile(balance$ratio, 0.99, na.rm = TRUE), 3))

# A large imbalance usually means very few devices rather than a broken
# hash, so report size alongside.
lopsided <- balance$ratio < 0.4 | balance$ratio > 0.6
message(" county-splits outside [0.4, 0.6]: ", sum(lopsided, na.rm = TRUE),
        " of ", nrow(balance),
        " (median n_total among them: ",
        signif(median(balance$n_total[lopsided], na.rm = TRUE), 4), ")")

### CHECK THE FULL PROFILE FOR THE SAME n_devices BUG ###
# If n_devices varies within a county in the full profile, it was built at
# band grain there too, and every n_devices >= 10000 filter downstream has
# been operating on the wrong quantity.

emp <- read_csv(PROFILE_FILE, show_col_types = FALSE) |>
  mutate(home_county = str_pad(as.character(home_county), 5, pad = "0"))

full_bad <- emp |>
  distinct(home_county, n_devices) |>
  count(home_county) |>
  filter(n > 1)

if (nrow(full_bad) > 0) {
  warning(nrow(full_bad), " counties have band-varying n_devices in ",
          basename(PROFILE_FILE),
          " - the full profile has the same bug and MIN_TRUST filters are wrong")
} else {
  message("full profile n_devices is constant within county - OK")
}

# Pooling the halves must reproduce the full profile. Any discrepancy means
# the aggregation here differs from the original build, in which case the
# reliability estimate would be measuring that difference rather than
# sampling noise.
pooled <- profile |>
  filter(salt == REP_SALTS[1]) |>
  group_by(home_county, band_idx) |>
  summarise(w_band = sum(w_band), .groups = "drop") |>
  group_by(home_county) |>
  mutate(share_pooled = w_band / sum(w_band)) |>
  ungroup() |>
  select(home_county, band_idx, share_pooled)

recon <- emp |>
  select(home_county, band_idx, share_full = share) |>
  inner_join(pooled, by = c("home_county", "band_idx")) |>
  mutate(d = abs(share_full - share_pooled))

message("pooled-vs-full share discrepancy: median ",
        signif(median(recon$d, na.rm = TRUE), 3),
        ", p99 ", signif(quantile(recon$d, 0.99, na.rm = TRUE), 3))

### WRITE ###

write_csv(profile, OUT_FILE)
message("wrote ", nrow(profile), " rows (",
        length(REP_SALTS), " splits) to ", OUT_FILE)
