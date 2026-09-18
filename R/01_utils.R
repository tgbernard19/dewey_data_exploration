# ===========================================================================
# R/01_utils.R
# ---------------------------------------------------------------------------
# Small things used everywhere. Sourced after config.R.
#
# The rule for what belongs here: it is used by more than one script, it has
# no scientific content, and getting it subtly wrong would be silent. FIPS
# padding and the haversine both qualify -- neither errors when wrong, they
# just produce fewer rows or shorter distances.
# ===========================================================================


# ---- guards ----------------------------------------------------------------

# base R's pi is not locked, so a stray project-level `pi <- 3.14` shadows it
# and silently rescales every distance and the self-disc radius. The result
# still looks like a plausible distribution, which is what makes it worth an
# assertion rather than a comment.
stopifnot(isTRUE(all.equal(pi, 3.141592653589793)))


# ---- keys ------------------------------------------------------------------

# FIPS is a five-character string, not a number. Excel strips leading zeros
# on every save, so 02164 (Alaska) arrives as 2164 and quietly matches
# nothing. Pad defensively on both sides of every join rather than trusting
# whatever wrote the file.
pad_fips <- function(x) {
  stringr::str_pad(as.character(x), width = 5, side = "left", pad = "0")
}

# Block groups are twelve characters, and the county is the first five.
pad_cbg <- function(x) {
  stringr::str_pad(as.character(x), width = 12, side = "left", pad = "0")
}

county_of_cbg <- function(x) {
  substr(pad_cbg(x), 1, 5)
}


# ---- distance --------------------------------------------------------------

RAD <- 0.017453292519943295
EARTH_KM <- 6371.0088

# Great-circle distance in km. Vectorised over the second point, which is how
# it is used in the flow allocation: one origin cell against every
# destination cell.
hav_km <- function(lon0, lat0, lon, lat) {
  p0 <- lat0 * RAD
  p1 <- lat * RAD
  a <- sin((p1 - p0) / 2)^2 +
    cos(p0) * cos(p1) * sin((lon - lon0) * RAD / 2)^2
  2 * EARTH_KM * asin(pmin(1, sqrt(a)))
}

# One degree of longitude at the equator is about 111.19 km. A version with
# radians and degrees confused, or the arguments in the wrong order, fails
# this and nothing downstream would have.
stopifnot(abs(hav_km(0, 0, 1, 0) - 111.19) < 0.5)
stopifnot(hav_km(-122.4, 37.8, -122.4, 37.8) == 0)


# ---- logging ---------------------------------------------------------------
# Timestamped, because the useful question about a long run is usually "how
# long did that stage take", and printing the time is cheaper than
# instrumenting anything.

say <- function(...) {
  cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")
}

banner <- function(txt) {
  say(strrep("=", 74))
  say(txt)
  say(strrep("=", 74))
}


# ---- provenance ------------------------------------------------------------

# Short git hash, or NA outside a repo. Stamped onto every output so a file
# found on disk in six months can be traced back to the code that made it.
git_sha <- function() {
  out <- suppressWarnings(
    try(system2("git", c("rev-parse", "--short", "HEAD"),
                stdout = TRUE, stderr = FALSE),
        silent = TRUE))
  if (inherits(out, "try-error") || length(out) == 0) NA_character_ else out[1]
}

# Adds the settings that define what a file MEANS, as columns. A profile
# without its tau, scope and home rule is uninterpretable, and the filename
# carries only some of that.
#
# Takes the values from config.R rather than arguments, so a script cannot
# stamp one setting and run under another.
stamp_run <- function(df) {
  dplyr::mutate(
    df,
    run_tag    = RUN_TAG,
    scope      = SCOPE,
    tau        = TAU,
    home_rule  = HOME_RULE,
    d_max_km   = D_MAX,
    git_sha    = git_sha(),
    built_on   = as.character(Sys.time())
  )
}


# ---- atomic writes ---------------------------------------------------------

# Write to .tmp, then rename. A killed run must not leave a truncated file
# that a later file.exists() check accepts as complete -- that failure is
# invisible until the numbers come out wrong.
write_atomic <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(path, ".tmp")
  readr::write_csv(df, tmp)
  file.rename(tmp, path)
  invisible(path)
}


# ---- reading ---------------------------------------------------------------

# Every CSV read in this project wants its key columns as text. Wrapping it
# once means no script can forget.
read_keyed_csv <- function(path, key_cols = character(0)) {
  
  # key_cols may name columns the file does not have -- callers pass both the
  # current and the older name for a key. Parsers are intersected with the
  # real header first, because readr warns (noisily, and deferred to the end
  # of the run) about named parsers that match nothing.
  header <- names(readr::read_csv(path, n_max = 0, show_col_types = FALSE,
                                  progress = FALSE))
  present <- intersect(key_cols, header)
  
  spec <- rep(list(readr::col_character()), length(present))
  names(spec) <- present
  
  df <- readr::read_csv(path, col_types = do.call(readr::cols, spec),
                        progress = FALSE)
  
  for (k in present) {
    if (stringr::str_detect(k, "fips|county$")) df[[k]] <- pad_fips(df[[k]])
  }
  df
}


# ---- misc ------------------------------------------------------------------

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a
}