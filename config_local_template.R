# ===========================================================================
# config_local_template.R
# ---------------------------------------------------------------------------
# Copy this file to config_local.R and edit the paths for THIS machine.
#
#   file.copy("config_local_template.R", "config_local.R")
#
# config_local.R is git-ignored. It is the only file in the repo that should
# ever contain an absolute path or a machine-specific setting, which is what
# makes everything else portable.
#
# Nothing here is a scientific choice. Estimand settings -- tau, the home
# rule, band edges, the distance cutoff -- live in config.R, are the same on
# every machine, and are part of what a result means.
# ===========================================================================


# ---- where the licensed inputs live ---------------------------------------

# Veraset, as downloaded by scripts/00_download_veraset.py. Expected to hold
# home_visits/, work_visits/ and other_visits/, each partitioned by day.
VERASET_ROOT <- "E:/dewey-june2025"

# Meta Movement Distribution, plus the county <-> GADM crosswalk.
META_ROOT <- "E:/meta_movement_dist"


# ---- where the pipeline writes --------------------------------------------

# Every intermediate and output. Can be anywhere with room; it is never the
# repo, and never a synced folder (OneDrive/Dropbox will corrupt parquet
# parts mid-write).
DERIVED_ROOT <- file.path(VERASET_ROOT, "kernel")

# DuckDB spills here when a query exceeds MEM_LIMIT. Wants a fast local disk
# with tens of GB free.
DUCKDB_TMP <- "E:/duckdb-tmp"


# ---- machine capacity ------------------------------------------------------

# Measured on the Windows box: wall clock peaks at 32 threads and degrades
# past 64, because the hash aggregations are memory-bound and hyperthread
# siblings share cache. Re-benchmark rather than copy this onto new hardware.
N_THREADS <- 32

# Generous is fine -- DuckDB spills rather than failing -- but it must not
# exceed physical RAM, or the OS starts swapping and everything crawls.
MEM_LIMIT <- "400GB"


# ---- Veraset download (only needed to run scripts/00_download_veraset.py) --

# The Dewey API key is read from the DEWEY_API_KEY environment variable, not
# from this file, so that it cannot reach a commit even by accident. Set it
# in ~/.Renviron or the system environment.
#
# The month to pull. Everything downstream was built on June 2025.
DOWNLOAD_START <- "2025-06-01"
DOWNLOAD_END   <- "2025-06-30"