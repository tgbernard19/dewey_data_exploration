# ===========================================================================
# R/05_meta.R
# ---------------------------------------------------------------------------
# Meta Movement Distribution: raw file in, per-county shares out.
#
# WHAT META GIVES YOU
#   For each admin-2 unit and each two-week period, the fraction of pings
#   falling in four categories of distance from the user's home tile:
#
#     "0"           the home tile itself
#     "(0, 10)"     under 10 km
#     "[10, 100)"
#     "100+"
#
#   It is a share of PINGS, sampled at a random moment, so it is
#   time-weighted: it says where people ARE, not where they GO. That is why
#   the Veraset side is matched to it with SCOPE = "all".
#
# THE FOLD
#   Category "0" is an atom at a point. No continuous kernel can put mass on
#   a point, so it cannot be fitted as a fourth bin. Folding it into (0, 10)
#   is the only coherent way to include it -- and since the home tile is a
#   few km across, unambiguously inside 10 km, folding is arithmetically
#   exact rather than an approximation.
#
#   Folding leaves the family exactly identified: three categories still pin
#   two numbers, S(10) and S(100), now of the unconditional survival.
#
#   The alternative (META_SCOPE = "away") drops the home tile and
#   renormalises, giving displacement conditional on having left home. That
#   pairs with a trip-scope Veraset profile, not an all-scope one. Mixing
#   them is the estimand mismatch that produced a -0.111 median gap in the
#   August work; matching them moved it to +0.007.
#
# f0 IS CARRIED THROUGH EITHER WAY
#   It is diagnostic under the current volume rule (N = pop) rather than
#   load-bearing, but its near-constancy across US counties (~0.346, sd
#   0.009) is worth watching: if it moves, something upstream changed.
# ===========================================================================


# The four category labels, exactly as they appear in the file. Checked
# rather than assumed: if Meta changes the spacing in "(0, 10)", a pivot on
# these names silently produces a table of NAs and drop_na() then deletes
# every county. That failure looks like "no data" three steps later.
META_CATEGORIES <- c("0", "(0, 10)", "[10, 100)", "100+")


# ---- crosswalk -------------------------------------------------------------

# Two irregularities, handled differently, because they mean different things.
#
#   one gid_2 -> many fips   EXPECTED. GADM folds Virginia's independent
#                            cities into their surrounding county, so both
#                            counties legitimately receive the same Meta
#                            kernel. Kept, and flagged as shared_meta_unit so
#                            they can be excluded from cross-county
#                            regressions, where they are not independent
#                            observations.
#
#   one fips -> many gid_2   AMBIGUOUS. There is no principled way to pick,
#                            and a silent join would duplicate the county and
#                            split its population across the copies. Dropped.
#
# The key is always the zero-padded FIPS. Never `clean`: it is a lowercased
# county name with no state, and about thirty states have a Clark County.
load_crosswalk <- function(path = XWALK_FILE) {
  
  xw <- read_keyed_csv(path, key_cols = "fips") |>
    dplyr::select(fips, gid_2, dplyr::any_of(c("gadm_name", "state"))) |>
    dplyr::filter(!is.na(fips), !is.na(gid_2)) |>
    dplyr::distinct()
  
  ambiguous <- xw |>
    dplyr::distinct(fips, gid_2) |>
    dplyr::count(fips) |>
    dplyr::filter(n > 1) |>
    dplyr::pull(fips)
  
  shared <- xw |>
    dplyr::filter(!fips %in% ambiguous) |>
    dplyr::count(gid_2) |>
    dplyr::filter(n > 1) |>
    dplyr::pull(gid_2)
  
  xw |>
    dplyr::filter(!fips %in% ambiguous) |>
    dplyr::mutate(shared_meta_unit = gid_2 %in% shared)
}


# ---- raw file --------------------------------------------------------------

read_meta_md <- function(path = META_MD_FILE) {
  
  md <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  
  need <- c("gadm_id", "ds", "home_to_ping_distance_category",
            "distance_category_ping_fraction")
  missing <- setdiff(need, names(md))
  if (length(missing) > 0) {
    stop("Meta MD is missing column(s): ", paste(missing, collapse = ", "),
         "\nPresent: ", paste(names(md), collapse = ", "))
  }
  
  found <- unique(md$home_to_ping_distance_category)
  unexpected <- setdiff(META_CATEGORIES, found)
  if (length(unexpected) > 0) {
    stop("Meta MD category labels do not match what this code expects.\n",
         "  expected: ", paste(sprintf("'%s'", META_CATEGORIES), collapse = ", "),
         "\n  found:    ", paste(sprintf("'%s'", found), collapse = ", "),
         "\nUpdate META_CATEGORIES in R/05_meta.R rather than working around it.")
  }
  
  md
}


# ---- shares ----------------------------------------------------------------

# Per-county shares, averaged across periods.
#
# Two normalisations, in this order:
#
#   1. WITHIN a period, divide by that period's total. The four fractions
#      should already sum to 1; doing it explicitly means a unit with a
#      rounding shortfall does not drag its own average down.
#   2. ACROSS periods, take the unweighted mean. Meta ships no sample size,
#      so there is nothing to weight by.
#
# Units without all four categories in a period are dropped for that period
# rather than renormalised over three, which would inflate whichever
# categories survived.
meta_county_shares <- function(md, xw, scope = META_SCOPE) {
  
  stopifnot(scope %in% c("folded", "away"))
  
  periods <- md |>
    dplyr::inner_join(xw, by = dplyr::join_by(gadm_id == gid_2),
                      relationship = "many-to-many") |>
    dplyr::group_by(fips, ds) |>
    dplyr::filter(dplyr::n() == 4) |>
    dplyr::mutate(frac = distance_category_ping_fraction /
                    sum(distance_category_ping_fraction, na.rm = TRUE)) |>
    dplyr::ungroup()
  
  wide <- periods |>
    dplyr::group_by(fips, home_to_ping_distance_category) |>
    dplyr::summarise(frac = mean(frac, na.rm = TRUE),
                     n_periods = dplyr::n(), .groups = "drop") |>
    tidyr::pivot_wider(names_from = home_to_ping_distance_category,
                       values_from = c(frac, n_periods)) |>
    dplyr::rename(f0 = `frac_0`,
                  g1 = `frac_(0, 10)`,
                  g2 = `frac_[10, 100)`,
                  g3 = `frac_100+`,
                  n_periods = `n_periods_0`) |>
    dplyr::select(fips, f0, g1, g2, g3, n_periods) |>
    tidyr::drop_na(f0, g1, g2, g3)
  
  out <- if (scope == "folded") {
    # Home tile into the first away band. Shares already sum to 1.
    dplyr::mutate(wide, o1 = f0 + g1, o2 = g2, o3 = g3)
  } else {
    # Conditional on having left home.
    dplyr::mutate(wide,
                  o1 = g1 / (1 - f0),
                  o2 = g2 / (1 - f0),
                  o3 = g3 / (1 - f0))
  }
  
  out |>
    dplyr::left_join(dplyr::select(xw, fips, gid_2, shared_meta_unit),
                     by = "fips", relationship = "one-to-one") |>
    dplyr::mutate(meta_scope = scope)
}


# ---- one county, in the shape the fitters want -----------------------------

# Edges are 0 / 10 / 100 / Inf under either scope. Under "away" the first
# edge is still 0, because the renormalised g1 describes displacement from
# zero upward among those who left -- the home tile is removed from the
# population, not from the distance axis.
meta_bands_for <- function(shares, f) {
  r <- dplyr::filter(shares, fips == pad_fips(f))
  if (nrow(r) != 1) {
    stop("Expected exactly one row for ", f, ", found ", nrow(r))
  }
  tibble::tibble(
    edge_low  = c(0, 10, 100),
    edge_high = c(10, 100, Inf),
    w         = c(r$o1, r$o2, r$o3)
  )
}