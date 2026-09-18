#!/usr/bin/env Rscript
# ===========================================================================
# 01_check_inputs.R
# ---------------------------------------------------------------------------
# Confirms that everything the pipeline reads exists, has the columns it is
# expected to have, and reads with keys as text rather than numbers.
#
# Run it first on any new machine, and after any rebuild. It writes nothing
# and takes seconds: parquet checks read only file footers, CSV checks read
# five rows.
#
# WHAT IT IS FOR
#   A missing file surfaces here rather than six hours into a scan. More
#   usefully, a file that exists under a slightly different name surfaces as
#   a list of candidates rather than as "file not found" -- the profile files
#   in particular exist in several permutations of scope, weighting and tau.
#
# WHAT IT DELIBERATELY DOES NOT DO
#   No scientific checks. Shares summing to 1, home rows sitting at zero
#   distance, f0 near 0.346 -- those belong to the scripts that build or
#   consume the files, where a failure means something about the data rather
#   than about the setup.
#
# Usage:  Rscript R/01_check_inputs.R
# ===========================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(purrr)
})

source("config.R")


# ===========================================================================
# 1. HELPERS
# ===========================================================================

results <- list()

record <- function(label, status, detail) {
  results[[length(results) + 1]] <<-
    tibble(label = label, status = status, detail = detail)
  invisible(NULL)
}

# Files whose names vary between runs are easier to find than to guess. When
# a path is missing, list what is actually in that directory so the config
# can be corrected in one edit.
candidates <- function(path, pattern = NULL) {
  dir <- dirname(path)
  if (!dir.exists(dir)) return(character(0))
  
  if (is.null(pattern)) {
    # Match on the first token of the filename: county_profile_x.csv finds
    # every county_profile_*.csv in the directory.
    stem <- basename(path) |> str_remove("\\.[^.]+$")
    pattern <- str_split(stem, "_")[[1]][1]
  }
  
  list.files(dir, pattern = pattern)
}

check_dir <- function(path, label, required = TRUE) {
  if (dir.exists(path)) {
    n <- length(list.files(path))
    record(label, "OK", sprintf("%d entries", n))
  } else {
    record(label, if (required) "MISSING" else "absent", path)
  }
  invisible(dir.exists(path))
}

# CSVs are read with every column as character, so nothing is coerced before
# it can be inspected. A FIPS column read as a number has already lost the
# leading zero of every Alaska and Alabama county by the time you see it.
check_csv <- function(path, label, need_cols = character(0),
                      key_cols = character(0), required = TRUE) {
  
  if (!file.exists(path)) {
    alt <- candidates(path)
    detail <- if (length(alt) > 0) {
      paste0("not found. In that directory: ",
             paste(head(alt, 8), collapse = ", "),
             if (length(alt) > 8) sprintf(" (+%d more)", length(alt) - 8))
    } else {
      paste0("not found: ", path)
    }
    record(label, if (required) "MISSING" else "absent", detail)
    return(invisible(NULL))
  }
  
  head_rows <- read_csv(path, n_max = 5, col_types = cols(.default = col_character()),
                        progress = FALSE)
  
  missing_cols <- setdiff(need_cols, names(head_rows))
  if (length(missing_cols) > 0) {
    record(label, "SCHEMA",
           paste0("missing column(s): ", paste(missing_cols, collapse = ", "),
                  " | present: ", paste(names(head_rows), collapse = ", ")))
    return(invisible(NULL))
  }
  
  # A key stored without its leading zeros is the single most common way a
  # join silently loses rows here. Excel strips them on every save.
  unpadded <- key_cols |>
    keep(~ .x %in% names(head_rows)) |>
    keep(~ any(!is.na(head_rows[[.x]]) &
                 nchar(head_rows[[.x]]) < 5 &
                 str_detect(head_rows[[.x]], "^[0-9]+$")))
  
  if (length(unpadded) > 0) {
    record(label, "WARN",
           paste0("key column(s) look unpadded: ",
                  paste(unpadded, collapse = ", "),
                  " -- read with col_character() and pad_fips() before joining"))
    return(invisible(NULL))
  }
  
  record(label, "OK", sprintf("%d columns", ncol(head_rows)))
  invisible(NULL)
}

# Parquet schema comes from the file footer, so this is effectively free even
# on a month of national data.
check_parquet <- function(path_or_glob, label, need_cols = character(0),
                          required = TRUE) {
  
  files <- if (str_detect(path_or_glob, "\\*")) {
    Sys.glob(path_or_glob)
  } else if (file.exists(path_or_glob)) {
    path_or_glob
  } else {
    character(0)
  }
  
  if (length(files) == 0) {
    record(label, if (required) "MISSING" else "absent", path_or_glob)
    return(invisible(NULL))
  }
  
  schema <- DBI::dbGetQuery(con, sprintf(
    "DESCRIBE SELECT * FROM read_parquet('%s') LIMIT 0", files[1]))
  
  missing_cols <- setdiff(need_cols, schema$column_name)
  if (length(missing_cols) > 0) {
    record(label, "SCHEMA",
           paste0("missing column(s): ", paste(missing_cols, collapse = ", "),
                  " | present: ", paste(schema$column_name, collapse = ", ")))
    return(invisible(NULL))
  }
  
  record(label, "OK", sprintf("%d file(s), %d columns",
                              length(files), nrow(schema)))
  invisible(NULL)
}


# ===========================================================================
# 2. CONNECTION
# ===========================================================================

con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")


# ===========================================================================
# 3. ROOTS AND SETTINGS
# ===========================================================================

cat("\n")
cat("PanDORA input check\n")
cat("  run tag      ", RUN_TAG, "\n")
cat("  tau          ", TAU, " minutes\n")
cat("  home rule    ", HOME_RULE, "\n")
cat("  scope        ", SCOPE, "\n")
cat("  D_MAX        ", D_MAX, " km\n")
cat("  obs weight   ", OBS_WEIGHT, "\n")
cat("  DuckDB       ", as.character(utils::packageVersion("duckdb")), "\n")
cat("\n")

check_dir(VERASET_ROOT, "Veraset root")
check_dir(META_ROOT,    "Meta root")
check_dir(DERIVED_ROOT, "Derived root")
check_dir(DUCKDB_TMP,   "DuckDB temp dir")
check_dir(REFERENCE_DIR, "Reference data dir")


# ===========================================================================
# 4. REFERENCE DATA (committed, needed to rebuild from raw)
# ===========================================================================

check_csv(CBG_CENTROIDS_FILE, "cbg_centroids.csv",
          need_cols = c("cbg", "cbg_lat", "cbg_lon", "area_km2"),
          key_cols  = "cbg")

check_csv(COUNTY_CENTROIDS_FILE, "county_centroids.csv",
          need_cols = c("fips", "lat", "lon", "cen_pop"),
          key_cols  = "fips")

check_csv(XWALK_FILE, "county_gid2_crosswalk.csv",
          need_cols = c("fips", "gid_2"),
          key_cols  = "fips")

# If the repo copy is missing but the Windows copy exists, say so explicitly
# rather than leaving the user to guess which is canonical.
if (!file.exists(XWALK_FILE) && file.exists(XWALK_FILE_ALT)) {
  record("crosswalk location", "WARN",
         paste0("found at ", XWALK_FILE_ALT,
                " but not in the repo. Copy it to ", REFERENCE_DIR))
}


# ===========================================================================
# 5. LICENSED INPUTS
# ===========================================================================

check_csv(META_MD_FILE, "Meta movement distribution",
          need_cols = c("gadm_id", "ds", "home_to_ping_distance_category",
                        "distance_category_ping_fraction"))

# Raw Veraset. Only needed to rebuild the chain from scratch.
for (src in VERASET_SOURCES) {
  check_dir(file.path(VERASET_ROOT, src), paste0("raw: ", src),
            required = FALSE)
}


# ===========================================================================
# 6. VERASET BUILD INTERMEDIATES
# ---------------------------------------------------------------------------
# Not required if the derived files below already exist -- the flow stage
# runs off small CSVs. Reported as absent rather than missing so that a
# clean analysis-machine setup does not look broken.
# ===========================================================================

check_parquet(HOME_ASSIGN_PQ, "home_assignment.parquet",
              need_cols = c("did", "home_cbg"), required = FALSE)

check_parquet(file.path(VISIT_PARTS_DIR, "*.parquet"), "visit_parts",
              need_cols = c("did", "day", "cbg", "county", "loc_src",
                            "is_home", "dwell", "d_km"),
              required = FALSE)

check_parquet(DEVDAY_GLOB, "devday_parts_v2",
              need_cols = c("did", "day", "home_county", "cbg", "county",
                            "loc_src", "is_home", "d_km", "dwell_pos_sum",
                            "n_tau_visits", "n_visits"),
              required = FALSE)

check_parquet(HOMERES_FILE, "device_home_resolution.parquet",
              need_cols = c("did", "cbg_home", "n_home_cbg", "n_home_geo"),
              required = FALSE)

# The strict home rules need the two count columns; the lenient one only
# needs the flag. Say so here rather than failing inside a query later.
if (HOME_RULE != "flag" && !file.exists(HOMERES_FILE)) {
  record("home rule inputs", "WARN",
         paste0("HOME_RULE = '", HOME_RULE,
                "' needs n_home_cbg / n_home_geo from ",
                basename(HOMERES_FILE)))
}


# ===========================================================================
# 7. WHAT THE FLOW STAGE ACTUALLY READS
# ===========================================================================

check_csv(PROFILE_CLEAN_FILE, "clean county profile",
          need_cols = c("home_county", "band_idx", "share", "n_devices"),
          key_cols  = "home_county",
          required  = FALSE)

check_csv(META_PARAMS_FILE, "Meta kernel parameters",
          need_cols = c("fips", "mu"),
          key_cols  = "fips")

if (!file.exists(META_PARAMS_FILE) && file.exists(META_PARAMS_FILE_ALT)) {
  record("Meta params location", "WARN",
         paste0("found at ", META_PARAMS_FILE_ALT, " -- move it to ",
                TRACKB_DIR, " so the path is derived, not special-cased"))
}

check_csv(VERASET_PARAMS_FILE, "Veraset kernel parameters",
          need_cols = c("home_county", "scope", "mu", "sigma"),
          key_cols  = "home_county",
          required  = FALSE)

check_csv(OBSERVED_PAIRS_FILE, "Veraset observed pairs",
          need_cols = c("origin_fips", "dest_fips", "n_devices", "w_dwell"),
          key_cols  = c("origin_fips", "dest_fips"),
          required  = FALSE)


# ===========================================================================
# 8. REPORT
# ===========================================================================

DBI::dbDisconnect(con, shutdown = TRUE)

report <- bind_rows(results)

cat("\n")
for (i in seq_len(nrow(report))) {
  cat(sprintf("  [%-7s] %-32s %s\n",
              report$status[i], report$label[i], report$detail[i]))
}
cat("\n")

n_missing <- sum(report$status %in% c("MISSING", "SCHEMA"))
n_warn    <- sum(report$status == "WARN")
n_absent  <- sum(report$status == "absent")

cat(sprintf("  %d ok, %d warning(s), %d absent (optional), %d problem(s)\n\n",
            sum(report$status == "OK"), n_warn, n_absent, n_missing))

if (n_missing > 0) {
  stop(n_missing, " required input(s) missing or wrong shape. ",
       "Fix the paths in config_local.R, or rebuild the stage that writes them.")
}

cat("Inputs OK.\n")