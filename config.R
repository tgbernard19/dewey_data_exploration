# ===========================================================================
# config.R
# ---------------------------------------------------------------------------
# Sourced by every script. Two jobs:
#
#   1. Turn the machine-specific roots in config_local.R into the full set of
#      paths the pipeline uses, so no script contains an absolute path.
#   2. Hold every setting that defines constants and design decisions: e.g.
#      how long a dwell time do we infer for visits in Veraset without one?
#
# ===========================================================================


# ---- machine-specific settings ---------------------------------------------

if (!file.exists("config_local.R")) {
  stop("config_local.R not found. Copy config_local_template.R to ",
       "config_local.R and edit the paths for this machine.\n",
       "Working directory is: ", getwd(),
       "\n(Run scripts from the repo root, or open the .Rproj first.)")
}

source("config_local.R")

for (nm in c("VERASET_ROOT", "META_ROOT", "DERIVED_ROOT", "DUCKDB_TMP",
             "N_THREADS", "MEM_LIMIT")) {
  if (!exists(nm)) stop("config_local.R does not define ", nm,
                        ". Compare it against config_local_template.R.")
}


# ===========================================================================
# PATHS
# ---------------------------------------------------------------------------
# Naming: *_DIR is a directory, *_FILE is a single file, *_GLOB is a pattern
# passed to read_parquet(). Nothing below is machine-specific -- it is all
# built from the roots above.
# ===========================================================================

# ---- reference data, committed to the repo ---------------------------------
# Small, public, and needed to rebuild from raw. Kept in the repo rather than
# under DERIVED_ROOT so that a clone plus the licensed data is sufficient.
#
# NOTE: the folder is called data/processed for historical reasons. The
# contents are inputs, not outputs. Renaming it to data/reference means
# changing this line and three lines in .gitignore.
REFERENCE_DIR <- file.path("data", "processed")

CBG_CENTROIDS_FILE    <- file.path(REFERENCE_DIR, "cbg_centroids.csv")
COUNTY_CENTROIDS_FILE <- file.path(REFERENCE_DIR, "county_centroids.csv")
XWALK_FILE            <- file.path(REFERENCE_DIR, "county_gid2_crosswalk.csv")

# The crosswalk also lives under META_ROOT on the Windows box. The repo copy
# is the one of record; this is checked, not assumed, in 01_check_inputs.R.
XWALK_FILE_ALT <- file.path(META_ROOT, "county_gid2_crosswalk.csv")

# ---- licensed inputs -------------------------------------------------------

VERASET_SOURCES <- c("home_visits", "work_visits", "other_visits")

META_MD_FILE <- file.path(
  META_ROOT, "movement-distribution-1-june-2026_15-june-2026.csv")

# ---- Veraset build intermediates -------------------------------------------

HOME_PARTS_DIR  <- file.path(DERIVED_ROOT, "home_parts")
HOME_ASSIGN_PQ  <- file.path(DERIVED_ROOT, "home_assignment.parquet")
VISIT_PARTS_DIR <- file.path(DERIVED_ROOT, "visit_parts")
DEVDAY_DIR      <- file.path(DERIVED_ROOT, "devday_parts_v2")
DEVDAY_GLOB     <- file.path(DEVDAY_DIR, "*.parquet")

# ---- Track B outputs -------------------------------------------------------

TRACKB_DIR <- file.path(DERIVED_ROOT, "trackB")

HOMERES_FILE <- file.path(TRACKB_DIR, "device_home_resolution.parquet")

# Kernel parameters. Both are keyed on FIPS-like columns; see the schema
# expectations in 01_check_inputs.R.
VERASET_PARAMS_FILE <- file.path(TRACKB_DIR, "veraset_lognormal_params.csv")
META_PARAMS_FILE    <- file.path(TRACKB_DIR, "meta_lognormal_kernel-fit.csv")

# The Meta fit was previously written to the root of the data drive by
# kernel_gen-PPC-check.R. Checked as a fallback so an existing file is found
# rather than silently rebuilt.
META_PARAMS_FILE_ALT <- file.path(dirname(VERASET_ROOT),
                                  "meta_lognormal_kernel-fit.csv")

# Observed Veraset destinations -- the empirical benchmark.
OBSERVED_PAIRS_FILE <- file.path(TRACKB_DIR, "veraset_observed_pairs.csv")
OBSERVED_PANEL_FILE <- file.path(TRACKB_DIR, "veraset_observed_panel.csv")

# ---- flow stage ------------------------------------------------------------

FLOW_CACHE_DIR <- file.path(DERIVED_ROOT, "flow_cache")
GEODATA_DIR    <- file.path(FLOW_CACHE_DIR, "geodata")
FLOW_OUT_DIR   <- file.path(FLOW_CACHE_DIR, "county_flows")

# ---- figures and tables the repo produces ----------------------------------

OUT_DIR <- "outputs"


# ===========================================================================
# ESTIMAND SETTINGS
# ---------------------------------------------------------------------------
# These define what the numbers mean. Changing one invalidates every derived
# file that was built under the old value, which is why every output carries
# them as columns (see stamp_run() in R/01_utils.R).
# ===========================================================================

# ---- tau: minutes credited to a single-ping visit --------------------------
# minimum_dwell is last-ping minus first-ping, so a zero means a visit of
# unknown duration, not zero duration (33% of home rows, 42% of work, 73% of
# other). Deleting those rows would drop short visits preferentially.
#
# Profiles are effectively invariant across tau = 1, 15, 60 (largest mean
# band difference 0.0017), so this is a reported constant, not a tuned
# parameter. 60 is used because the completed national run used 60, and the
# clean profile on disk carries tau60 in its name.
TAU <- 60

# ---- home resolution rule --------------------------------------------------
# device_home_resolution.parquet is device grain, with n_home_cbg and
# n_home_geo counting a device's home rows by how the location resolved.
# Geohash-5 homes carry a median 3.467 km displacement on rows whose true
# displacement is zero, which contaminates every band under 10 km.
#
#   "flag"      cbg_home = 1: ANY home row resolved to a CBG
#   "pure_cbg"  n_home_cbg > 0 AND n_home_geo = 0: no geohash home rows
#   "majority"  n_home_cbg > n_home_geo
#   "any_cbg"   n_home_cbg > 0  (identical to "flag" by construction)
#
# "flag" is the lenient rule and is what the clean profile on disk was built
# with. "pure_cbg" is stricter and is what build_veraset_observed.R defaulted
# to. They must MATCH between the profile and the observed flows, or a
# difference in device populations shows up as model error.
HOME_RULE <- "flag"

# ---- scope -----------------------------------------------------------------
# "all"   home rows kept. Time-weighted and unconditional -- where people
#         ARE at a random moment. Matches Meta with the home tile folded into
#         (0, 10), and is the estimand the N = pop volume rule requires.
# "trips" home rows dropped. Where people GO. Needs a different volume rule.
#
# Matching this to the Meta side is what moved the median gap from -0.111 to
# +0.007 in the August work.
SCOPE <- "all"

# ---- band edges ------------------------------------------------------------
# 11 bands. The fit uses the 9 whose upper edge is at or below TRUNC_KM.
BAND_EDGES <- c(0, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000, Inf)

# Truncation for fitting. The [1000, Inf) bump is still unexplained, long
# distance will be handled separately by flights, and conditioning the
# likelihood on the observed range stops the model reading a cutoff as fast
# decay.
TRUNC_KM <- 500

# ---- day boundary ----------------------------------------------------------
# The device-day is the normalisation unit, so where the day is cut matters.
# "utc" takes the date from utc_timestamp, which is what every existing
# intermediate used: a Pacific evening visit lands on the following day.
# "local" would use local_timestamp instead.
#
# Kept at "utc" for consistency with the files already built. Recorded here
# so the choice is visible rather than implicit.
DAY_SOURCE <- "utc"

# ---- location resolution ---------------------------------------------------
# Block group centroid when the block group is small enough to be informative,
# otherwise the geohash-5 cell centre. Above this area the centroid can sit
# tens of km from the actual visit.
CBG_AREA_MAX <- 25   # km^2

# ---- flow stage ------------------------------------------------------------
# D_MAX is both the domain radius and the kernel support, kept identical so
# the spatial cut and the probabilistic cut are the same cut. It must match
# between the observed extractor and the flow comparison.
D_MIN <- 0
D_MAX <- 500

# Grid resolution for the population raster. The slice width used in the
# allocation is the realised cell size at the origin's latitude.
GRID_KM <- 5

# ---- Meta handling ---------------------------------------------------------
# "folded" folds the home tile into (0, 10), giving an unconditional
# distribution that matches SCOPE = "all". "away" drops it and renormalises,
# giving displacement conditional on having left home. Not interchangeable.
META_SCOPE <- "folded"

# Pseudo sample size for the Meta fit. Meta ships fractions with no sample
# size, so this scales the objective enough for BFGS to have a gradient.
#
# IMPORTANT: point estimates do not depend on it, but the covariance from the
# inverted Hessian scales as 1/META_OBS_WEIGHT. The uncertainty in the Meta
# fit file is therefore on an arbitrary scale and means nothing in absolute
# terms. It is fine for the maps, which use only mu and sigma; treat any
# interval derived from it as indicative.
META_OBS_WEIGHT <- 1000

# ---- analysis thresholds ---------------------------------------------------
# Reliable independent estimation of a Veraset kernel needs ~10,000 devices:
# runaway rates are 48% under 1,000, 7.7% at 2.5-5k, 0.6% at 5-10k, 0% above.
# A reliability flag, not a data floor -- hierarchical pooling would let
# smaller counties contribute with shrinkage.
MIN_DEVICES_KERNEL <- 10000

# Used for the observed side instead. Nothing is fitted there, so the
# threshold is about sampling noise in a destination share, not about whether
# an optimiser converges.
MIN_DEVICES_MAP <- 500

# ---- observed-flow weighting -----------------------------------------------
# "time"     dwell-weighted with one-user-one-vote normalisation, home rows
#            included. Matches the kernel estimand and the N = pop rule.
# "devices"  one count per device per destination county. Closest to "one
#            person away from home", which was the right match under the OLD
#            volume rule of pop * (1 - f0).
#
# The extractor emits several weightings in one pass, so switching is a
# re-read rather than a rescan.
OBS_WEIGHT <- "time"


# ===========================================================================
# RUN CONTROL
# ===========================================================================

# Tags every output built under the current settings, and is part of the
# filenames for profile parts.
RUN_TAG <- sprintf("national_%s_dwell_tau%d", SCOPE, TAU)

PROFILE_PARTS_DIR  <- file.path(DERIVED_ROOT,
                                sprintf("profile_parts_%s", RUN_TAG))
PROFILE_FILE       <- file.path(TRACKB_DIR,
                                sprintf("county_profile_%s.csv", RUN_TAG))
PROFILE_CLEAN_FILE <- file.path(
  TRACKB_DIR, "county_profile_cbghome_national_all_dwell_tau60.csv")

# Which Veraset build stage to rebuild from. Every stage at or after this
# number reruns; earlier ones are reused from disk. NA reuses everything that
# exists.
#
# This replaces the scattered file.exists() guards, which could not tell the
# difference between "cached" and "stale". Device ids are DuckDB hash(caid)
# values, and that function is not documented as stable across DuckDB
# versions -- so if DuckDB is upgraded, rebuild the whole chain (1) rather
# than one stage, or joins on did will quietly drop devices.
REBUILD_FROM <- NA_integer_

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)