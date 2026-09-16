# ============================================================================
# config.R -- shared machine-specific paths and constants
# ----------------------------------------------------------------------------
# Every script in this repo that reads/writes the Dewey/Veraset/Meta pipeline
# data sources this file instead of hardcoding its own copies of these paths.
# Edit the values below for your machine; the rest of each script should not
# need to change.
#
# Usage (each script assumes it is run with the repo root as the working
# directory, matching the .Rproj setup):
#
#   source("config.R")
#
# Environment-variable overrides are supported for the two base roots so the
# same scripts can run unmodified on a different machine/drive letter:
#   DEWEY_DATA_DIR  overrides DATA_DIR   (default "E:/dewey-june2025")
#   DEWEY_META_DIR  overrides META_DIR   (default "E:/meta_movement_dist")
# ============================================================================

# ---- base roots -------------------------------------------------------------
DATA_DIR <- Sys.getenv("DEWEY_DATA_DIR", unset = "E:/dewey-june2025")
META_DIR <- Sys.getenv("DEWEY_META_DIR", unset = "E:/meta_movement_dist")

KERNEL_DIR  <- file.path(DATA_DIR, "kernel")
TRACKB_DIR  <- file.path(KERNEL_DIR, "trackB")
DEVDAY_DIR  <- file.path(KERNEL_DIR, "devday_parts_v2")
DEVDAY_GLOB <- file.path(DEVDAY_DIR, "*.parquet")

# ---- raw visit feeds (narrow_national_pipeline.R) --------------------------
HOME_VISITS_DIR  <- file.path(DATA_DIR, "home_visits")
WORK_VISITS_DIR  <- file.path(DATA_DIR, "work_visits")
OTHER_VISITS_DIR <- file.path(DATA_DIR, "other_visits")

# ---- shared derived files ---------------------------------------------------
HOMERES_PQ      <- file.path(TRACKB_DIR, "device_home_resolution.parquet")
PROFILE_FILE    <- file.path(TRACKB_DIR, "county_profile_cbghome_national_all_dwell_tau60.csv")
CENTROID_CSV    <- file.path(DATA_DIR, "county_centroids.csv")
HOME_LOOKUP_PQ  <- file.path(DATA_DIR, "home_lookup.parquet")

MD_PATH          <- file.path(META_DIR, "movement-distribution-1-june-2026_15-june-2026.csv")
XWALK_FILE       <- file.path(META_DIR, "county_gid2_crosswalk.csv")
BENCHMARK_CSV    <- file.path(META_DIR, "benchmark_counties_gid2.csv")
META_LOGNORMAL_FIT <- "E:/meta_lognormal_kernel-fit.csv"
RUCC_CSV         <- "E:/rural_continuum/Ruralurbancontinuumcodes2023.csv"

# ---- scratch space -----------------------------------------------------------
DUCKDB_TMP_DIR <- "E:/duckdb-tmp"

# ---- shared constants --------------------------------------------------------
# Domain radius shared by the destination extract and the flow comparison --
# these two MUST agree, or the observed benchmark and the fitted kernel
# describe different truncations of the same distribution.
D_MAX <- 500   # km

# duckdb resource knobs used by the heavier scans (build_veraset_observed.R,
# split-half*.R).
N_THREADS <- 32
MEM_LIMIT <- "400GB"
