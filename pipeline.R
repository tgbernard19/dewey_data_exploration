#!/usr/bin/env Rscript

# =============================================================================
# TRACK B — VERASET DISPLACEMENT KERNEL -> POOLED PRIOR ON p
# =============================================================================
#
# WHAT THIS SCRIPT DOES, END TO END
#
# Phase 0 Preflight. Checks packages, compiles Stan, validates the parquet
# schema. Fails here rather than at 4am.
# Phase 1 Builds a national device x band profile from devday_parts, one
# day partition at a time. Resumable.
# Phase 1b QC REPORT after the first day. This is the gate. Read it before
# letting the rest run.
# Phase 2 Combines the day parts into a county x band profile.
# Phase 3 Profile-level QC gates (hard stops).
# Phase 4 Per-county Stan fits on the 9 bands below 500 km.
# Phase 5 Fit diagnostics, including the posterior predictive check.
# Phase 6 Does p_v transfer? (the check that gates pooling)
# Phase 7 The pooled beta prior on p.
# Phase 8 Outputs and a run manifest.
#
# -----------------------------------------------------------------------------
# HOW TO RUN IT TONIGHT
#
# 1. Leave SMOKE_TEST <- TRUE. Run the script. It processes ONE day partition,
# prints the QC report, and stops.
# 2. Read the report. Every line has a stated expectation next to it. The two
# numbers that decide whether to proceed are the [0,1) share in check 5 and
# the distinct-day count in check 2.
# 3. Set SMOKE_TEST <- FALSE and run again. The first day is already cached
# on disk, so nothing is recomputed — it picks up at day two.
#
# Everything after the gate is defensive: the day loop is resumable, the fits
# are checkpointed every FIT_CHUNK counties, a county that fails to sample is
# recorded with its error message rather than killing the run, and a first
# chunk that fails wholesale aborts instead of burning the night.
#
# -----------------------------------------------------------------------------
# THE STATISTICAL OBJECT, IN ONE PARAGRAPH
#
# Displacement from a device's own home is modelled with a survival function
# that mixes two exponentials:
#
# S(d) = p * exp(-a1 * d) + (1 - p) * exp(-a2 * d), a2 < a1
#
# p is the weight on the short-range component, a1 its rate, a2 the tail rate.
# The probability of landing in a band [lo, hi) is S(lo) - S(hi).
#
# Meta's three away-categories pin exactly two numbers, S(10) and S(100). Three
# parameters, two constraints, so one direction is unidentified — a ridge. The
# global fitter currently picks a point on that ridge using an arbitrary
# beta(9,1) prior on p. Veraset's nine sub-500 km bands can see the shape that
# ridge coordinate controls, so this script estimates it and hands back a
# defensible prior in its place.
#
# The prior goes on p and the tail ONLY. Vague priors stay on the directions
# Meta's likelihood constrains, so Meta keeps determining S(10) and S(100).
# A county-level Veraset prior on the LEVEL would be roughly 5x tighter than
# Meta's entire cross-county spread and would overwrite Meta rather than inform
# it — and would smuggle in the finding-2 disagreement disguised as information.
#
# -----------------------------------------------------------------------------
# THE ESTIMAND, AND THE ONE CHANGE THIS REQUIRES IN THE GLOBAL FITTER
#
# SCOPE = "all": home dwell is included. Meta samples a device at a random
# moment and most random moments are at home, so a trip-weighted Veraset
# quantity is not the same object as Meta's ping fractions. Matching on this
# is what moved the median gap from -0.111 to +0.007.
#
# BUT the global fitter as written does NOT match that. Its conditional_period
# block computes o1 = b1/away with away = b1+b2+b3, so Meta's category "0" is
# stripped and f0 carried separately. Its (p, a1, a2) therefore describe
# displacement GIVEN you left the home tile. Veraset SCOPE = "all" describes
# displacement unconditionally. p is the short-component weight, so it is
# precisely the parameter that notices the difference.
#
# The fix is in the fitter, not here. Fold the home tile in:
#
# conditional_period <- per_period %>%
# transmute(gadm_id, ds, away = 1, o1 = b0 + b1, o2 = b2, o3 = b3)
#
# Meta's home tile is a few km across, unambiguously inside (0, 10). Folding
# costs no identification — three categories still pin S(10) and S(100), now of
# the unconditional survival — and it rescues the units currently QC'd out as
# "no_away_mass". It does change Map 1, and globally that is not a constant
# shift, because f0 is only near-constant within the US. The conditional map is
# recoverable afterwards by renormalising the kernel above a threshold.
#
# If you decide NOT to fold, do not use the prior this script produces. Build it
# from SCOPE = "trips" instead — accepting that this is the pairing that gave
# the -0.111 gap.
#
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DBI)
  library(duckdb)
  library(cmdstanr)
  library(posterior)
  library(furrr)
  library(arrow)
})


# =============================================================================
# PHASE 0 — CONFIGURATION AND PREFLIGHT
# =============================================================================

# ---- the gate ---------------------------------------------------------------

# TRUE = process one day, print the QC report, stop.
# FALSE = full run.
SMOKE_TEST <- TRUE

# ---- paths ------------------------------------------------------------------

DATA_ROOT <- "E:/dewey-data"
KERNEL_DIR <- file.path(DATA_ROOT, "kernel")
DEVDAY_DIR <- file.path(KERNEL_DIR, "devday_parts")
OUT_DIR <- file.path(KERNEL_DIR, "trackB")
TMP_DIR <- file.path(DATA_ROOT, "tmp")

# ---- estimand and weighting -------------------------------------------------

SCOPE <- "all" # "all" keeps home rows (matched to Meta); "trips" drops them
TAU <- 15 # minutes credited to a single-ping visit

# Why tau = 15: across a 60-fold range (1, 15, 60) the largest mean band
# difference was 0.0017. The profiles are effectively invariant, so this is a
# reported constant rather than a tuned parameter. It also falsified tau as the
# explanation for the far tail.
#
# Why the weight is dwell_pos_sum + tau * n_tau_visits: minimum_dwell is
# last-ping minus first-ping, so a zero means a single-ping visit of unknown
# duration, not zero time. Those visits get tau minutes rather than being
# deleted. Carrying the two components separately is what makes tau sensitivity
# a re-query rather than a re-scan.

# ---- band edges -------------------------------------------------------------

# 11 bands. The fit uses the 9 with an upper edge at or below TRUNC_KM.
BAND_EDGES <- c(0, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000, Inf)
TRUNC_KM <- 500

# Why truncate at 500: the [1000, Inf) bump is still unexplained, long distance
# will be handled separately by flights, and conditioning the likelihood on the
# observed range keeps the model from reading a cutoff as fast decay. Because
# the band probabilities are renormalised inside Stan, the fitted (p, a1, a2)
# still describe the same untruncated S(d).

# ---- compute ----------------------------------------------------------------

N_THREADS <- 32 # peaks here on this box; degrades at 64-128
MEM_LIMIT <- "400GB"
WORKERS <- 16 # counties in parallel during the fit phase
FIT_CHUNK <- 200 # checkpoint the fits every this many counties
SEED <- 20260831

# ---- Stan parameter space ---------------------------------------------------

# MUST MATCH the global fitter's A1_RANGE / A2_RANGE. A prior extracted in a
# different parameter space than the grid it is later applied on is not
# transportable. If one changes, change both.
A1_RANGE <- c(0.02, 2.0) # short-component median distance ~0.35 to 35 km
A2_RANGE <- c(1e-4, 0.5) # tail median distance ~1.4 to 6,900 km

# PPC replicate size is capped. n_devices runs into the millions for Los
# Angeles, and multinomial_rng at that size, times 1,500 draws, times two
# chains, is pointless expense. Capping makes the chi-square-style test less
# powerful for the largest counties, which is why a scale-free misfit measure
# (total variation between observed and fitted bands) is reported alongside it.
PPC_N_CAP <- 50000L

# ---- reporting --------------------------------------------------------------

# A flag, NOT an exclusion. Every county is fitted. Under the hierarchical
# extension this becomes the shrinkage story rather than a deletion rule.
MIN_REPORT <- 500

RUN_TAG <- sprintf("national_%s_dwell_tau%d", SCOPE, TAU)
PART_DIR <- file.path(OUT_DIR, sprintf("profile_parts_%s", RUN_TAG))
FIT_DIR <- file.path(OUT_DIR, sprintf("fit_parts_%s", RUN_TAG))

dir.create(PART_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

AREA_FILE <- file.path(OUT_DIR, "county_area.csv")

# ---- logging ----------------------------------------------------------------

LOG_FILE <- file.path(OUT_DIR, sprintf("trackB_%s.log",
                                       format(Sys.time(), "%Y%m%d_%H%M%S")))

say <- function(...) {
  line <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = LOG_FILE, append = TRUE)
}

banner <- function(txt) {
  say(strrep("=", 74))
  say(txt)
  say(strrep("=", 74))
}

banner("PHASE 0 — PREFLIGHT")
say("Log: ", LOG_FILE)
say("Mode: ", if (SMOKE_TEST) "SMOKE TEST (one day, then stop)" else "FULL RUN")

# ---- guard: packages used later, not at the top ------------------------------

# arrow::write_parquet runs inside the fit loop and sf/tigris in Phase 3. A
# missing one of those would otherwise surface hours in.
for (pkg in c("arrow", "sf", "tigris")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "arrow") {
      stop("arrow is required — the fit checkpoints are parquet.")
    }
    say("[!] ", pkg, " not installed. County areas will be skipped, and with ",
        "them the transfer check in Phase 6.")
  }
}

# ---- guard: pi --------------------------------------------------------------

# base R's pi is not locked, so a project-level `pi <- ...` shadows it silently.
# d_km is precomputed upstream, but a shadowed pi rescales every distance by a
# constant and the result still looks like a plausible distribution.
stopifnot(isTRUE(all.equal(pi, 3.141592653589793)))
say("pi is unshadowed.")

# ---- guard: inputs exist ----------------------------------------------------

stopifnot(dir.exists(DEVDAY_DIR))

day_files <- sort(list.files(DEVDAY_DIR, pattern = "\\.parquet$",
                             full.names = TRUE))
stopifnot(length(day_files) > 0)
say("Found ", length(day_files), " day partitions.")

# ---- guard: Stan compiles NOW, not after the long scan ----------------------

# Compiling here means a broken toolchain costs seconds instead of surfacing
# after six hours of DuckDB work.

stan_code <- "
data {
int<lower=2> K; // number of fitted bands
vector<lower=0>[K] edge_lo; // band lower edges, km
vector<lower=0>[K] edge_hi; // band upper edges, km
vector<lower=0>[K] share_raw; // observed band shares
real<lower=0> n_pseudo; // learning rate, NOT a row count
int<lower=1> n_rep; // replicate size for the PPC
real a1_lo; real a1_hi;
real a2_lo; real a2_hi;
}

transformed data {
vector[K] w = share_raw / sum(share_raw);
}

parameters {
real<lower=0, upper=1> p;
real<lower=a1_lo, upper=a1_hi> a1;
// Ordering by parameter-dependent upper bound, as in the global fitter.
// Enforced by the transform with the Jacobian handled automatically, so
// there is no -Inf wall mid-parameter-space the way target += log(a1 > a2)
// would create. Preferred to positive_ordered[2], which gives ordering but
// no upper bounds and would need explicit truncation.
real<lower=a2_lo, upper=fmin(a1, a2_hi)> a2;
}

transformed parameters {
vector[K] q; // model band probabilities
{
vector[K] num;
for (k in 1:K) {
num[k] = p * (exp(-a1 * edge_lo[k]) - exp(-a1 * edge_hi[k]))
+ (1 - p) * (exp(-a2 * edge_lo[k]) - exp(-a2 * edge_hi[k]));
}
// Dividing by the total is exactly 'condition the likelihood on the
// observed range'. Without it the model reads the 500 km cutoff as
// genuine fast decay.
q = num / sum(num);
}
}

model {
p ~ beta(1, 1); // vague; the nine bands should do the work
target += -log(a1) - log(a2); // bounded log-uniform, matching the fitter

// Gibbs / pseudo-posterior. n_pseudo is a learning rate in the sense of
// Bissiri, Holmes & Walker (2016), fed n_devices because device-days cluster
// ~13 per device, so the independent unit is the device.
target += n_pseudo * dot_product(w, log(q));
}

generated quantities {
real S10 = p * exp(-a1 * 10) + (1 - p) * exp(-a2 * 10);
real S100 = p * exp(-a1 * 100) + (1 - p) * exp(-a2 * 100);

// The comparison quantity from the band work: P(d > 10 | d < 100), logit.
real lo_v = logit((S10 - S100) / (1 - S100));

// Scale-free misfit: total variation between observed and fitted bands.
// Interpretable regardless of n_pseudo, unlike the PPC p-value.
real tv_dist = 0.5 * sum(abs(w - q));

// Boundary pressure, mirroring the fitter's edge_mass diagnostic. If this is
// widespread the grid is too narrow for what the sub-10 km bands can see and
// the extracted prior is truncated.
int at_edge = (a1 > 0.98 * a1_hi || a1 < 1.02 * a1_lo ||
a2 > 0.98 * a2_hi || a2 < 1.02 * a2_lo) ? 1 : 0;

// Posterior predictive check. Nine bands against three parameters leaves six
// degrees of freedom. The old A05 fit had three bands and three parameters,
// was saturated, and could never show misfit at all. This is the main
// scientific gain from the rebuild.
real T_obs = 0;
real T_rep = 0;
int ppc_exceed;
{
array[K] int y_rep = multinomial_rng(q, n_rep);
for (k in 1:K) {
if (w[k] > 0) T_obs += n_rep * w[k] * log(w[k] / q[k]);
if (y_rep[k] > 0) T_rep += y_rep[k] * log((y_rep[k] * 1.0 / n_rep) / q[k]);
}
ppc_exceed = (T_rep > T_obs) ? 1 : 0;
}
}
"

say("Compiling Stan model ...")
STAN_FILE <- write_stan_file(stan_code)
model <- cmdstan_model(STAN_FILE)

# The compiled executable's path, NOT the R6 object. cmdstan_model() returns an
# R6 wrapper around an external binary, and serialising that into multisession
# workers is unreliable — if it fails, every county lands in the try() and the
# night produces 3,100 failures and no prior. Each worker rebuilds a handle from
# this path instead, which is cheap because the binary already exists.
STAN_EXE <- model$exe_file()
say("Stan compiled OK -> ", STAN_EXE)

# ---- guard: schema ----------------------------------------------------------

con <- dbConnect(duckdb(), dbdir = file.path(OUT_DIR, "trackB_build.duckdb"))
dbExecute(con, sprintf("SET threads = %d", N_THREADS))
dbExecute(con, sprintf("SET memory_limit = '%s'", MEM_LIMIT))
dbExecute(con, sprintf("SET temp_directory = '%s'", TMP_DIR))
dbExecute(con, "SET preserve_insertion_order = false")

# DESCRIBE ... LIMIT 0 reads only the parquet footers. Instant, and it means a
# renamed column is caught before a multi-hour scan rather than during one.
schema <- dbGetQuery(con, sprintf(
  "DESCRIBE SELECT * FROM read_parquet('%s') LIMIT 0", day_files[1]))

needed <- c("did", "day", "home_county", "d_km", "is_home",
            "dwell_pos_sum", "n_tau_visits")
missing <- setdiff(needed, schema$column_name)

if (length(missing)) {
  stop("Missing columns in devday_parts: ", paste(missing, collapse = ", "),
       "\nPresent: ", paste(schema$column_name, collapse = ", "))
}
say("Schema OK. ", nrow(schema), " columns.")


# =============================================================================
# PHASE 1 — BUILD THE NATIONAL PROFILE, ONE DAY AT A TIME
# =============================================================================
#
# SQL LESSON — four constructs, stacked
#
# WITH name AS (...) A named intermediate result, a "common table
# expression". It lets the query read top to bottom
# instead of nesting inside-out. Purely for
# legibility; the planner inlines them.
#
# CASE WHEN x THEN y END An if-else ladder that returns a value rather
# than doing something. Used here to turn a
# continuous d_km into an integer band index.
# Conditions are evaluated in order, so the ladder
# must run smallest edge first.
#
# JOIN ... USING (a, b) Glues each row back to a row in another table
# matching on those columns. Here it attaches each
# visit's device-day total so the division that
# turns a weight into a share can happen row-wise.
#
# ANY_VALUE(x) Picks an arbitrary value from a group. Cheap.
# Correct only when the value is constant within
# the group — QC check 4 asserts that it is.
#
# NOTE what is absent: there is no COUNT(DISTINCT) anywhere in the build. Across
# hundreds of millions of groups it cannot stream, needs a hash set per group,
# and is what stalled other_visits at 10% with 32 GB spilled. It is not needed,
# because of the trick in the normalisation below. (It does appear in the QC
# block, but that runs against a single day partition, not the full month.)
#
# -----------------------------------------------------------------------------
# THE NORMALISATION, IN THREE STEPS
#
# 1. Within a device-day, band shares of w sum to 1.
# 2. Within a device, average those over ALL of its device-days.
# 3. Within a county, average over devices.
#
# Step 2 is the one that was previously wrong. The old code used
# n_distinct(day) grouped BY BAND, which counted only the days a device
# appeared in that band — so a device with one far-flung day divided it by 1.
# That inflated rare far bands and produced a nearly flat profile.
#
# The trick: because step 1 makes each device-day's shares sum to exactly 1,
# summing a device's shares across all bands and all days returns its
# device-day count for free. That is the denominator step 2 needs, obtained
# without a single distinct count.
# =============================================================================

banner("PHASE 1 — DAY LOOP")

band_case <- paste0(
  "CASE\n",
  paste(sprintf(" WHEN d_km < %s THEN %d", BAND_EDGES[2:11], 1:10),
        collapse = "\n"),
  "\n ELSE 11\n END"
)

scope_filter <- if (SCOPE == "trips") "AND is_home = FALSE" else ""

day_sql <- sprintf("
COPY (
WITH dd AS (
SELECT
did,
day,
home_county,
%s AS band_idx,
dwell_pos_sum + %f * n_tau_visits AS w
FROM read_parquet('%%s')
WHERE home_county IS NOT NULL
AND d_km IS NOT NULL
AND d_km >= 0
%s
),
pos AS (
-- A device-day with no positive weight has no defined share vector and is
-- dropped here rather than contributing a zero row.
SELECT * FROM dd WHERE w > 0
),
tot AS (
SELECT did, day, SUM(w) AS w_tot
FROM pos
GROUP BY did, day
)
SELECT
pos.did,
ANY_VALUE(pos.home_county) AS home_county,
pos.band_idx,
SUM(pos.w / tot.w_tot) AS share_sum
FROM pos
JOIN tot USING (did, day)
GROUP BY pos.did, pos.band_idx
) TO '%%s' (FORMAT PARQUET, COMPRESSION ZSTD)
", band_case, TAU, scope_filter)

first_new_day <- NA_character_
first_new_time <- NA_real_

for (f in day_files) {
  day <- tools::file_path_sans_ext(basename(f))
  dest <- file.path(PART_DIR, paste0(day, ".parquet"))
  
  if (file.exists(dest)) {
    say("[skip] ", day, " (cached)")
    next
  }
  
  tmp <- paste0(dest, ".tmp")
  t0 <- Sys.time()
  
  dbExecute(con, sprintf(day_sql, f, tmp))
  
  # Atomic write. A killed run must never leave a truncated file that the
  # file.exists() guard above later accepts as valid.
  file.rename(tmp, dest)
  
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  say(sprintf("[ok] %s %.1f min %.0f MB", day, mins,
              file.size(dest) / 1e6))
  
  if (is.na(first_new_day)) {
    first_new_day <- day
    first_new_time <- mins
    if (SMOKE_TEST) break
  }
}


# =============================================================================
# PHASE 1b — QC REPORT. THIS IS THE GATE.
# =============================================================================
# Every check states its expectation. Anything marked FAIL should be understood
# before the full run, because it will not get better at scale.
# =============================================================================

if (!is.na(first_new_day)) {
  
  banner("PHASE 1b — QC REPORT ON ONE DAY PARTITION")
  
  qc_part <- file.path(PART_DIR, paste0(first_new_day, ".parquet"))
  qc_src <- day_files[grepl(first_new_day, day_files, fixed = TRUE)][1]
  
  say("Day partition: ", first_new_day)
  
  # -- 1. output is non-empty -------------------------------------------------
  qc1 <- dbGetQuery(con, sprintf(
    "SELECT COUNT(*) AS n_rows, COUNT(DISTINCT did) AS n_dev
FROM read_parquet('%s')", qc_part))
  say(sprintf("1. Rows %s, devices %s. EXPECT: millions of devices, roughly",
              format(qc1$n_rows, big.mark = ","),
              format(qc1$n_dev, big.mark = ",")))
  say(" 3-8 rows per device (a device-day touches few bands).")
  stopifnot(qc1$n_rows > 0)
  
  # -- 2. shares sum to the device-day count ----------------------------------
  # The single most important check: it proves step 1 of the normalisation
  # closed, and therefore that the denominator trick in Phase 2 is valid.
  #
  # Counted against DEVICE-DAYS rather than devices on purpose. A partition
  # should hold one `day`, but the utc_timestamp day boundary resolving in the
  # session timezone is unaudited and is exactly the thing that would split a
  # stay across two dates. If that happens the identity still holds against
  # device-days — so this reports the finding instead of aborting on it.
  qc2 <- dbGetQuery(con, sprintf(
    "SELECT SUM(share_sum) AS total FROM read_parquet('%s')", qc_part))
  
  qc2b <- dbGetQuery(con, sprintf(
    "SELECT COUNT(DISTINCT day) AS n_days, COUNT(*) AS n_devday FROM (
SELECT did, day FROM read_parquet('%s')
WHERE home_county IS NOT NULL AND d_km IS NOT NULL AND d_km >= 0
AND dwell_pos_sum + %f * n_tau_visits > 0
GROUP BY did, day)", qc_src, TAU))
  
  gap2 <- abs(qc2$total - qc2b$n_devday) / qc2b$n_devday
  say(sprintf("2. sum(share_sum) = %.1f vs %s device-days, rel gap %.2e. %s",
              qc2$total, format(qc2b$n_devday, big.mark = ","), gap2,
              if (gap2 < 1e-7) "PASS" else "FAIL — normalisation is broken"))
  say(sprintf(" Distinct `day` values in this partition: %d. EXPECT 1.",
              qc2b$n_days))
  say(" More than 1 means the UTC day boundary is splitting stays — worth")
  say(" knowing, but not a reason to stop.")
  stopifnot(gap2 < 1e-7)
  
  # -- 3. no nulls, no bad bands ----------------------------------------------
  qc3 <- dbGetQuery(con, sprintf(
    "SELECT
SUM(CASE WHEN share_sum IS NULL OR NOT isfinite(share_sum)
THEN 1 ELSE 0 END) AS bad_share,
SUM(CASE WHEN home_county IS NULL THEN 1 ELSE 0 END) AS null_home,
SUM(CASE WHEN band_idx < 1 OR band_idx > 11 THEN 1 ELSE 0 END) AS bad_band
FROM read_parquet('%s')", qc_part))
  say(sprintf("3. bad_share %d, null_home %d, bad_band %d. EXPECT: all zero.",
              qc3$bad_share, qc3$null_home, qc3$bad_band))
  stopifnot(qc3$bad_share == 0, qc3$null_home == 0, qc3$bad_band == 0)
  
  # -- 4. one home county per device ------------------------------------------
  # ANY_VALUE(home_county) is only correct if home_county is constant within a
  # device. It comes from home_assignment at the device level, so it should be.
  qc4 <- dbGetQuery(con, sprintf(
    "SELECT COUNT(*) AS n_multi FROM (
SELECT did FROM read_parquet('%s')
GROUP BY did HAVING COUNT(DISTINCT home_county) > 1)", qc_part))
  say(sprintf("4. Devices with >1 home county: %d. EXPECT: 0. %s",
              qc4$n_multi,
              if (qc4$n_multi == 0) "" else "ANY_VALUE is unsafe — investigate."))
  stopifnot(qc4$n_multi == 0)
  
  # -- 5. the band profile itself ---------------------------------------------
  # The shape is the real check. Under SCOPE = "all" the [0,1) band carries home
  # dwell and should be HEAVY — well above the 0.202 seen under trips scope. If
  # it looks like the trips profile, the scope filter is misapplied.
  qc5 <- dbGetQuery(con, sprintf(
    "SELECT band_idx, SUM(share_sum) AS s FROM read_parquet('%s')
GROUP BY band_idx ORDER BY band_idx", qc_part)) |>
    mutate(share = s / sum(s),
           band = sprintf("[%g,%s)", BAND_EDGES[band_idx],
                          ifelse(is.infinite(BAND_EDGES[band_idx + 1]), "Inf",
                                 format(BAND_EDGES[band_idx + 1]))))
  say("5. Day band profile:")
  for (i in seq_len(nrow(qc5))) {
    say(sprintf(" %-12s %.4f", qc5$band[i], qc5$share[i]))
  }
  say(sprintf(" [0,1) share = %.3f. EXPECT under SCOPE='all': well above",
              qc5$share[1]))
  say(" 0.202 (the trips-scope value), because home dwell now counts.")
  if (SCOPE == "all" && qc5$share[1] < 0.25) {
    say(" *** WARNING: [0,1) looks like trips scope. Check is_home / the ***")
    say(" *** scope filter before running the full month. ***")
  }
  
  # -- 6. counties covered ----------------------------------------------------
  qc6 <- dbGetQuery(con, sprintf(
    "SELECT COUNT(DISTINCT home_county) AS n_cty FROM read_parquet('%s')",
    qc_part))
  say(sprintf("6. Counties present on this day: %s. EXPECT: near 3,100 —",
              format(qc6$n_cty, big.mark = ",")))
  say(" this is the national run, no origin filter anywhere upstream.")
  
  # -- 7. sanity against the source -------------------------------------------
  qc7 <- dbGetQuery(con, sprintf(
    "SELECT COUNT(*) AS n_src,
SUM(CASE WHEN home_county IS NULL THEN 1 ELSE 0 END) AS n_nohome,
SUM(CASE WHEN is_home THEN 1 ELSE 0 END) AS n_home_rows,
MEDIAN(d_km) AS med_d
FROM read_parquet('%s')", qc_src))
  say(sprintf("7. Source rows %s; no home county %s (%.1f%%); is_home rows %s.",
              format(qc7$n_src, big.mark = ","),
              format(qc7$n_nohome, big.mark = ","),
              100 * qc7$n_nohome / qc7$n_src,
              format(qc7$n_home_rows, big.mark = ",")))
  say(sprintf(" Median d_km in source: %.2f.", qc7$med_d))
  say(" EXPECT: is_home rows present and non-trivial under SCOPE='all'.")
  if (SCOPE == "all" && qc7$n_home_rows == 0) {
    say(" *** WARNING: no is_home rows at all. SCOPE='all' is a no-op. ***")
  }
  
  # -- 8. cost projection -----------------------------------------------------
  n_remaining <- sum(!file.exists(file.path(
    PART_DIR, paste0(tools::file_path_sans_ext(basename(day_files)),
                     ".parquet"))))
  est_hours <- first_new_time * n_remaining / 60
  est_gb <- file.size(qc_part) * length(day_files) / 1e9
  say(sprintf("8. %d days remaining at %.1f min each -> ~%.1f hours.",
              n_remaining, first_new_time, est_hours))
  say(sprintf(" Estimated parts footprint: ~%.1f GB on E:. Check free space.",
              est_gb))
  
  say("")
  say("QC report complete. If everything above reads sensibly, set")
  say("SMOKE_TEST <- FALSE and rerun. The cached day is not recomputed.")
}

if (SMOKE_TEST) {
  banner("STOPPING — SMOKE TEST MODE")
  dbDisconnect(con, shutdown = TRUE)
  # quit() would kill an RStudio session, so stop() when interactive.
  if (interactive()) stop("Smoke test complete — set SMOKE_TEST <- FALSE.")
  quit(save = "no", status = 0)
}


# =============================================================================
# PHASE 2 — COMBINE DAY PARTS INTO A COUNTY PROFILE
# =============================================================================

banner("PHASE 2 — COMBINE")

# Sum a device's band shares across every day it appears.
dbExecute(con, sprintf("
CREATE OR REPLACE TABLE dev_band AS
SELECT did, band_idx, SUM(share_sum) AS share_sum
FROM read_parquet('%s/*.parquet')
GROUP BY did, band_idx
", PART_DIR))

# One home county per device.
dbExecute(con, sprintf("
CREATE OR REPLACE TABLE dev_home AS
SELECT did, ANY_VALUE(home_county) AS home_county
FROM read_parquet('%s/*.parquet')
GROUP BY did
", PART_DIR))

# The denominator, for free: shares sum to 1 per device-day, so their total
# across bands and days is the device's device-day count. No distinct count.
dbExecute(con, "
CREATE OR REPLACE TABLE dev_denom AS
SELECT did, SUM(share_sum) AS n_devdays
FROM dev_band
GROUP BY did
")

say("Device tables built.")

county_long <- dbGetQuery(con, "
WITH dev_share AS (
SELECT h.home_county, b.band_idx,
b.share_sum / d.n_devdays AS dev_share
FROM dev_band b
JOIN dev_denom d USING (did)
JOIN dev_home h USING (did)
),
cty_n AS (
SELECT home_county, COUNT(*) AS n_devices
FROM dev_home
GROUP BY home_county
)
SELECT s.home_county, s.band_idx,
SUM(s.dev_share) / ANY_VALUE(n.n_devices) AS share,
ANY_VALUE(n.n_devices) AS n_devices
FROM dev_share s
JOIN cty_n n USING (home_county)
GROUP BY s.home_county, s.band_idx
")

say("County profile: ", n_distinct(county_long$home_county), " counties.")

band_key <- tibble(
  band_idx = 1:11,
  edge_low = BAND_EDGES[1:11],
  edge_high = BAND_EDGES[2:12]
) |>
  mutate(band = sprintf("[%g,%s)", edge_low,
                        ifelse(is.infinite(edge_high), "Inf",
                               format(edge_high))))

# complete() before anything downstream. A dropped empty band desynchronises K
# and edges, and Stan then fails with "failed to create the sampler", which
# surfaces much later as log(NULL).
profile <- county_long |>
  mutate(home_county = str_pad(as.character(home_county), 5, pad = "0")) |>
  complete(nesting(home_county, n_devices), band_idx = 1:11,
           fill = list(share = 0)) |>
  left_join(band_key, by = "band_idx") |>
  mutate(run_tag = RUN_TAG, scope = SCOPE, weighting = "dwell",
         tau_minutes = TAU, built_on = as.character(Sys.Date())) |>
  arrange(home_county, band_idx)


# =============================================================================
# PHASE 3 — PROFILE QC GATES (HARD STOPS)
# =============================================================================

banner("PHASE 3 — PROFILE QC")

chk <- profile |>
  group_by(home_county) |>
  summarise(tot = sum(share), k = n(), .groups = "drop")

say(sprintf("Counties with != 11 bands: %d (expect 0)", sum(chk$k != 11)))
say(sprintf("Max |sum(share) - 1|: %.2e (expect < 1e-6)",
            max(abs(chk$tot - 1))))
stopifnot(all(chk$k == 11), max(abs(chk$tot - 1)) < 1e-6)

nat <- profile |>
  group_by(band_idx, band) |>
  summarise(share = weighted.mean(share, n_devices), .groups = "drop")
say("National (device-weighted) profile:")
for (i in seq_len(nrow(nat))) {
  say(sprintf(" %-12s %.4f", nat$band[i], nat$share[i]))
}

write_csv(profile, file.path(OUT_DIR, sprintf("county_profile_%s.csv", RUN_TAG)))
say("Wrote county profile.")

# County areas, for the transfer check in Phase 6.
if (!file.exists(AREA_FILE) &&
    requireNamespace("tigris", quietly = TRUE) &&
    requireNamespace("sf", quietly = TRUE)) {
  area <- tigris::counties(cb = TRUE, year = 2020, progress_bar = FALSE) |>
    sf::st_drop_geometry() |>
    transmute(fips = GEOID, area_km2 = ALAND / 1e6)
  write_csv(area, AREA_FILE)
  say("Wrote county_area.csv (", nrow(area), " rows).")
}


# =============================================================================
# PHASE 4 — PER-COUNTY FITS
# =============================================================================

banner("PHASE 4 — FITTING")

fit_bands <- profile |>
  filter(edge_high <= TRUNC_KM) |>
  group_by(home_county) |>
  mutate(share_trunc = share / sum(share)) |>
  ungroup()

stopifnot(n_distinct(fit_bands$band_idx) == 9)

panel <- distinct(fit_bands, home_county, n_devices)

# The full threshold sweep, not just the chosen cut. Worth reading: if the SD
# never plateaus and the mean drifts as the threshold rises, n_devices is
# correlated with the estimand and trimming removes signal along with noise.
sweep <- tibble(threshold = c(100, 250, 500, 1000, 2000, 5000)) |>
  mutate(n_counties = map_int(threshold, \(t) sum(panel$n_devices >= t)))
say("Device-count sweep:")
for (i in seq_len(nrow(sweep))) {
  say(sprintf(" n >= %5d : %d counties", sweep$threshold[i],
              sweep$n_counties[i]))
}
say(sprintf("%d counties below %d devices — FLAGGED, NOT DROPPED.",
            sum(panel$n_devices < MIN_REPORT), MIN_REPORT))

# One county. The error message is kept rather than discarded, so a systematic
# problem is distinguishable from a genuinely hard county.
fit_one_county <- function(county_bands) {
  county_bands <- arrange(county_bands, band_idx)
  fips <- county_bands$home_county[1]
  n_dev <- county_bands$n_devices[1]
  
  stan_data <- list(
    K = nrow(county_bands),
    edge_lo = county_bands$edge_low,
    edge_hi = county_bands$edge_high,
    share_raw = county_bands$share_trunc,
    n_pseudo = n_dev,
    n_rep = min(max(as.integer(round(n_dev)), 50L), PPC_N_CAP),
    a1_lo = A1_RANGE[1], a1_hi = A1_RANGE[2],
    a2_lo = A2_RANGE[1], a2_hi = A2_RANGE[2]
  )
  
  res <- try({
    # Rebuild the model handle from the already-compiled binary. No recompile.
    mod <- cmdstan_model(exe_file = STAN_EXE)
    
    fit <- mod$sample(
      data = stan_data, seed = SEED, chains = 2, parallel_chains = 1,
      iter_warmup = 750, iter_sampling = 750, adapt_delta = 0.95,
      refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
    
    d <- fit$draws(c("p", "a1", "a2", "S10", "S100", "lo_v",
                     "tv_dist", "at_edge", "ppc_exceed"), format = "df")
    s <- fit$summary(c("p", "a1", "a2"))
    
    tibble(
      home_county = fips, n_devices = n_dev,
      p_v = mean(d$p), p_v_sd = sd(d$p),
      a1v = mean(d$a1), a1v_sd = sd(d$a1),
      a2v = mean(d$a2), a2v_sd = sd(d$a2),
      S10 = mean(d$S10), S100 = mean(d$S100),
      lo_v = mean(d$lo_v), lo_v_sd = sd(d$lo_v),
      tv_dist = mean(d$tv_dist),
      edge_mass = mean(d$at_edge),
      ppc_p = mean(d$ppc_exceed),
      max_rhat = max(s$rhat, na.rm = TRUE),
      min_ess = min(s$ess_bulk, na.rm = TRUE),
      n_divergent = sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent),
      status = "ok", err = NA_character_
    )
  }, silent = TRUE)
  
  if (inherits(res, "try-error")) {
    tibble(home_county = fips, n_devices = n_dev,
           status = "failed", err = as.character(res))
  } else {
    res
  }
}

# Checkpointed in chunks. A crash costs at most one chunk, and a rerun skips
# whatever is already on disk — the same discipline as the day loop.
county_list <- group_split(fit_bands, home_county)
chunks <- split(county_list, ceiling(seq_along(county_list) / FIT_CHUNK))

plan(multisession, workers = WORKERS)
t0 <- Sys.time()

for (i in seq_along(chunks)) {
  dest <- file.path(FIT_DIR, sprintf("chunk_%04d.parquet", i))
  if (file.exists(dest)) { say("[skip] chunk ", i); next }
  
  out <- future_map_dfr(chunks[[i]], fit_one_county,
                        .options = furrr_options(seed = TRUE,
                                                 packages = "cmdstanr"))
  tmp <- paste0(dest, ".tmp")
  arrow::write_parquet(out, tmp)
  file.rename(tmp, dest)
  
  n_bad <- sum(out$status != "ok")
  say(sprintf("[ok] chunk %d/%d (%d counties, %d failed, %.1f min elapsed)",
              i, length(chunks), nrow(out), n_bad,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))