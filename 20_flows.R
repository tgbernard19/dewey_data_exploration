#!/usr/bin/env Rscript
# ===========================================================================
# 20_flows.R
# ---------------------------------------------------------------------------
# Kernels in, maps out. This is the end of the pipeline and the thing the
# whole repo exists to produce.
#
# Reads:  meta_lognormal_kernel-fit.csv        (required)
#         veraset_lognormal_params.csv         (optional)
#         veraset_observed_pairs.csv           (optional)
#         county_centroids.csv, county_gid2_crosswalk.csv
#         GADM + WorldPop, downloaded on first run
#
# Writes: outputs/maps/flows_<fips>.png, one per origin county
#         outputs/flow_summary.csv, one row per origin
#
# Usage:  Rscript scripts/20_flows.R
#         Rscript scripts/20_flows.R 06037 36061      # specific counties
# ===========================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(sf)
})

source("config.R")
source("R/01_utils.R")
source("R/04_lognormal.R")
source("R/07_flows.R")
source("R/08_maps.R")

banner("FLOWS")


# ---- the tables the flow functions use ------------------------------------
# These are module-level by inheritance from the original script. Every
# function takes them as arguments with these as defaults, so an explicit
# call can override them -- which is what a stale copy in the workspace
# would otherwise silently break.

COUNTY   <- load_flow_inputs()
XW       <- read_keyed_csv(XWALK_FILE, key_cols = "fips")
OBSERVED <- load_observed()

say("counties: ", nrow(COUNTY),
    " | with Meta kernel: ", sum(COUNTY$has_meta),
    " | with Veraset kernel: ", sum(COUNTY$has_veraset),
    " | observed flows: ", if (is.null(OBSERVED)) "no" else "yes")


# ---- which counties --------------------------------------------------------
# Command line if given, otherwise the largest counties that have everything
# needed. Starting with the big ones is deliberate: they have the thickest
# panels, so a disagreement there is about the model rather than about
# sampling noise.

args <- commandArgs(trailingOnly = TRUE)

fips_vec <- if (length(args) > 0) {
  pad_fips(args)
} else {
  COUNTY |>
    filter(has_meta, !is.na(gid_2)) |>
    arrange(desc(pop)) |>
    slice_head(n = 40) |>
    pull(fips)
}

say("running ", length(fips_vec), " counties")


# ---- maps ------------------------------------------------------------------

MAP_DIR <- file.path(OUT_DIR, "maps")
dir.create(MAP_DIR, showWarnings = FALSE, recursive = TRUE)

for (f in fips_vec) {
  
  say("  ", f, " ...")
  
  res <- tryCatch(
    compare_county(f, min_flow = MAP_MIN_FLOW),
    error = function(e) { message("    skipped: ", conditionMessage(e)); NULL })
  
  if (is.null(res)) next
  
  ggsave(file.path(MAP_DIR, sprintf("flows_%s.png", f)),
         res$plot, width = 13, height = 6, dpi = 150)
}

say("maps written to ", MAP_DIR)


# ---- cross-county summary --------------------------------------------------
# One row per origin. The three total-variation columns are the decomposition:
#
#   tv_obs_vs_veraset   error from the population allocation rule alone
#                       (same data, same kernel family, different spreading)
#   tv_veraset_vs_meta  difference between the kernels alone
#                       (same allocation rule, different kernel)
#   tv_obs_vs_meta      the two combined -- what the deliverable actually costs
#
# Read as shares of residents' time placed in a different county, so 0.05 is
# five percent of people in the wrong place.

summ <- run_all(fips_vec)

if (nrow(summ) > 0) {
  write_atomic(stamp_run(summ), file.path(OUT_DIR, "flow_summary.csv"))
  say("summary written for ", nrow(summ), " counties")
  
  say("")
  say("total variation, median across counties:")
  summ |>
    summarise(across(starts_with("tv_"), \(x) median(x, na.rm = TRUE))) |>
    as.data.frame() |>
    print(row.names = FALSE)
  
  say("")
  say("largest Meta vs observed disagreement:")
  summ |>
    filter(!is.na(tv_obs_vs_meta)) |>
    slice_max(tv_obs_vs_meta, n = 10) |>
    select(fips, pop, n_devices, tv_obs_vs_meta, tv_veraset_vs_meta,
           tv_obs_vs_veraset) |>
    as.data.frame() |>
    print(row.names = FALSE)
}