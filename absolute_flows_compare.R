#!/usr/bin/env Rscript
# ===========================================================================
# absolute_flows_compare.R
# Absolute origin -> destination flows under the Veraset kernel and under
# the Meta kernel, on a shared population grid, and the log discrepancy
# between them mapped by destination county.
#
# Entry point: compare_county("06037")
# ---------------------------------------------------------------------------
#
# THE VOLUME RULE
#
# N = pop
#
# Both kernels are fitted with home included: Meta's o1 has the home tile
# folded in, and Veraset is SCOPE = "all". Each kernel therefore already
# describes where EVERYONE's time goes, with the at-home share sitting near
# zero distance. So the kernel is applied to the whole population.
#
# The previous version applied these kernels to n_leave = pop * (1 - f0).
# That removed the at-home share twice -- once in the kernel, once in the
# volume -- and scaled every flow by (1 - f0) ~ 0.654, a uniform log offset
# of about -0.42. The Veraset-vs-Meta panel could not see it because both
# sides got the same wrong N; it only shows against an external benchmark.
#
# f0 no longer enters the volume. It enters only through the fold,
# o1 = share(home tile) + share(0,10).
#
# UNITS
#
# A flow is the expected number of origin residents present in the
# destination at a random moment (person-days per day, time-weighted), not a
# trip count. The origin county receives most of it, because most time is
# spent at or near home.
#
# TRUNCATION
#
# N_source = pop * P_source(d <= D_MAX)
#
# Each source keeps the time its own kernel places inside the domain, then
# allocates it with the kernel renormalised over the domain. The in-domain
# totals differ between sources by exactly the difference in tail mass.
#
# KERNEL FAMILY
#
# Lognormal. On clean profiles it misplaces 4.2% of mass against the
# incumbent's 8.1%, with one fewer parameter, and Meta's three categories
# identify it exactly in closed form.
#
# ALLOCATION: CELL-WIDTH SLICES
#
# For each origin cell, distance is cut into slices one grid cell wide. The
# first slice is a disc with the same area as a cell, radius
# r0 = g / sqrt(pi), and holds only the origin cell itself. Slice k covers
# [r0 + (k-1) g, r0 + k g). The kernel's probability in each slice is split
# across the cells in that slice in proportion to their population.
#
# The previous 1 km bins failed because cells are ~4-5 km apart: from any
# cell, the bins between 1 km and the nearest neighbour hold no cell centres.
# Their mass was spread pro rata over every other bin, including ones 300 km
# away -- about 8% of the mass under the all-scope kernel in Los Angeles.
# A cell-wide slice always contains cells on land. The only slices left
# empty are genuinely empty (ocean, foreign), and their mass is returned to
# the populated slices: "conditional on being somewhere populated".
#
# For this to work the origin cells have to BE destination cells, so the
# origin's own cell sits at distance 0. Origins are therefore taken from the
# destination grid rather than from a separately masked raster.
# ===========================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr)
  library(stringr); library(tibble); library(terra); library(sf)
  library(ggplot2)
})

# `pi` is not locked in base R and a project-level assignment would silently
# rescale every distance and the self-disc radius.
stopifnot(isTRUE(all.equal(pi, 3.141592653589793)))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a

# ===========================================================================
# 1. CONFIGURATION
# ===========================================================================

# ---- paths (edit for this machine) ----------------------------------------
TRACKB_DIR <- "E:/dewey-june2025/kernel/trackB"
META_DIR <- "E:/meta_movement_dist"
CACHE_DIR <- "E:/dewey-june2025/kernel/flow_cache"
OUT_DIR <- file.path(CACHE_DIR, "county_flows")

VERASET_PARAMS <- file.path(TRACKB_DIR, "veraset_lognormal_params.csv")
META_PARAMS <- file.path("E:/meta_lognormal_kernel-fit.csv")
XWALK <- file.path(META_DIR, "county_gid2_crosswalk.csv")
CENTROIDS <- "E:/dewey-june2025/county_centroids.csv" # fips, lat, lon, cen_pop
GEODATA_CACHE <- file.path(CACHE_DIR, "geodata")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(GEODATA_CACHE, showWarnings = FALSE, recursive = TRUE)

# ---- modelling window ------------------------------------------------------
D_MIN <- 0 # km. Home is inside the kernel, so no distance floor
D_MAX <- 500 # km. domain radius and kernel support, kept identical so the
# spatial cut and the probabilistic cut are the same cut

# ---- grid ------------------------------------------------------------------
GRID_KM <- 5 # matches the geohash-5 resolution floor; the slice width is
# the realised cell size at the origin's latitude

# ---- which Veraset estimand -----------------------------------------------
# "all" is time-weighted and unconditional, home included (median
# displacement ~108 m -- the median moment of a day is at home). It is the
# estimand that matches Meta with the home tile folded in, and it is the
# one that goes with N = pop.
# "trip" is conditional on moving and trip-weighted. It would need a
# different volume rule, and it is not time-weighted like Meta, so it is
# not used here.
VERASET_SCOPE <- "all"

# ---- map scales -------------------------------------------------------------
# The log-ratio panel is fixed, in log2 units, so the same colour means the
# same thing in every county's map. 2 is a factor of 4 either way; anything
# past it is squished to the end of the scale rather than restretching the
# whole map around one outlier.
LOG2_LIM <- 4

# ---- cache ------------------------------------------------------------------
# Bumped whenever the stored object changes shape. Cached results from an
# older version are then simply not found, instead of being read back and
# silently missing fields. Bump it if you add anything to `res`.
CACHE_VERSION <- "v2"

# ===========================================================================
# 2. SMALL HELPERS
# ===========================================================================

pad_fips <- function(x) str_pad(as.character(x), 5, "left", "0")

RAD <- 0.017453292519943295
hav_km <- function(lon0, lat0, lon, lat) {
  p0 <- lat0 * RAD; p1 <- lat * RAD
  a <- sin((p1 - p0) / 2)^2 + cos(p0) * cos(p1) * sin((lon - lon0) * RAD / 2)^2
  2 * 6371.0088 * asin(pmin(1, sqrt(a)))
}
stopifnot(abs(hav_km(0, 0, 1, 0) - 111.19) < 0.5)

# Untruncated lognormal CDF. This is the object the volume rule needs, since
# N is scaled by the untruncated mass falling inside the domain.
lognormal_cdf <- function(d, mu, sigma) {
  pnorm((log(pmax(d, 1e-12)) - mu) / sigma)
}

# ...and the version renormalised over [D_MIN, D_MAX], which is what the
# allocation uses, so that both sources are spread over identical support and
# a heavier tail cannot make one look smaller everywhere inside the domain.
lognormal_cdf_trunc <- function(d, mu, sigma, lo = D_MIN, hi = D_MAX) {
  F_lo <- lognormal_cdf(lo, mu, sigma)
  F_hi <- lognormal_cdf(hi, mu, sigma)
  denom <- F_hi - F_lo
  if (!is.finite(denom) || denom <= 0)
    stop("Window carries no mass for mu=", mu, " sigma=", sigma)
  pmin(1, pmax(0, (lognormal_cdf(d, mu, sigma) - F_lo) / denom))
}

# Not used by the pipeline any more -- the Meta fit arrives as (mu,
# log_sigma). Kept because it is how to CHECK that fit: take a county's raw
# MD shares, fold the home tile into o1, and this returns the (mu, sigma)
# that reproduces them exactly. If the fitted kernel's F(10) sits near
# f0 + c1, the fold happened and N = pop is the right volume rule.
meta_lognormal_from_shares <- function(o1, o2) {
  z1 <- qnorm(o1)
  z2 <- qnorm(o1 + o2)
  sigma <- log(10) / (z2 - z1)
  mu <- log(10) - sigma * z1
  tibble(mu = mu, sigma = sigma)
}

# ---- cell-width slices ------------------------------------------------------
# Edges: 0, r0, r0 + g, r0 + 2g, ..., capped at D_MAX.
slice_edges <- function(g) {
  r0 <- g / sqrt(pi)
  K <- ceiling((D_MAX - r0) / g)
  pmin(c(0, r0 + g * (0:K)), D_MAX)
}

# 1-based slice index: 1 is the self disc, 2 the first ring, and so on.
slice_of <- function(d, g) {
  r0 <- g / sqrt(pi)
  K <- ceiling((D_MAX - r0) / g)
  ifelse(d < r0, 1L, as.integer(pmin(ceiling((d - r0) / g), K)) + 1L)
}

# Kernel probability per slice, window-normalised.
slice_mass <- function(k, edges) diff(lognormal_cdf_trunc(edges, k$mu, k$sigma))

# ===========================================================================
# 3. INPUTS
# ---------------------------------------------------------------------------
# Column contract. If a file on this machine names things differently, fix it
# here in one place rather than downstream.
#
# VERASET_PARAMS home_county (or fips), scope, mu, sigma;
# n_devices, at_bound, tv if present
# META_PARAMS fips, mu, log_sigma (or sigma), converged;
# f0 and the Hessian block if present, both diagnostic
# XWALK fips, gid_2 (join on padded fips, never on `clean`)
# CENTROIDS fips, lat, lon, cen_pop
#
# BOTH KERNEL FILES ARE KEYED ON FIPS, so neither is joined through the
# crosswalk. The crosswalk is still needed, but only for geometry: gid_2
# selects the origin polygon and labels destination cells.
#
# NEITHER FILE CARRIES o1 / o2 ANY MORE, and it does not matter. The shares
# were only ever a route to (mu, sigma), and (mu, sigma) is what arrives.
# What the shares would have let you verify is that the home tile was folded
# into (0,10) before fitting -- a folded fit satisfies F(10) = f0 + c1 from
# the raw MD file. That fold is assumed here, and it is the assumption
# N = pop rests on. Check it once against the MD file (see the note on
# meta_lognormal_from_shares below) and then leave it alone.
#
# f0 is NA in the current fit file. That costs nothing: f0 is diagnostic
# only under N = pop. It would be load-bearing under the old volume rule,
# which could not run on these files at all.
# ===========================================================================

load_inputs <- function() {
  
  # The Veraset file keys on `home_county`; older ones used `fips`.
  ver <- read_csv(VERASET_PARAMS, show_col_types = FALSE) |>
    rename(any_of(c(fips = "home_county"))) |>
    mutate(fips = pad_fips(fips)) |>
    filter(scope == VERASET_SCOPE) |>
    select(fips, mu_v = mu, sigma_v = sigma,
           any_of(c("n_devices", "at_bound", "tv")))
  stopifnot("Veraset params: no rows for this scope." = nrow(ver) > 0)
  stopifnot("Veraset params: duplicate FIPS within scope." = !anyDuplicated(ver$fips))
  if (!"n_devices" %in% names(ver)) ver$n_devices <- NA_real_
  if (!"at_bound" %in% names(ver)) ver$at_bound <- FALSE
  if (!"tv" %in% names(ver)) ver$tv <- NA_real_
  
  # The Meta fit is an MLE on (mu, log_sigma), keyed on FIPS, so sigma has to
  # be exponentiated. Getting that wrong is silent: a raw log_sigma of 0.70
  # read as sigma gives a kernel that is merely far too tight, with no error
  # anywhere. The check is that P(d < 10 km) lands near 0.8, not near 1.
  meta_raw <- read_csv(META_PARAMS, show_col_types = FALSE) |>
    mutate(fips = pad_fips(fips))
  stopifnot("Meta params: duplicate FIPS." = !anyDuplicated(meta_raw$fips))
  if (!"f0" %in% names(meta_raw)) meta_raw$f0 <- NA_real_
  if (!"converged" %in% names(meta_raw)) meta_raw$converged <- TRUE
  
  stopifnot("Meta params: need sigma or log_sigma." =
              any(c("sigma", "log_sigma") %in% names(meta_raw)))
  meta <- meta_raw |>
    mutate(sigma_m = if ("sigma" %in% names(meta_raw)) sigma else exp(log_sigma)) |>
    select(fips, f0, mu_m = mu, sigma_m, meta_converged = converged)
  
  bad <- meta |> filter(!meta_converged | !is.finite(mu_m) | !is.finite(sigma_m) |
                          sigma_m <= 0)
  if (nrow(bad) > 0)
    message(sprintf("[meta] %d counties with an unusable fit (not converged or non-finite)",
                    nrow(bad)))
  meta <- meta |> anti_join(bad, by = "fips")
  
  # Crosswalk: pad both sides, then dedupe. Seven counties map to more than
  # one GADM unit and are dropped rather than silently multiplied. One gid_2
  # to many fips is expected -- GADM folds Virginia independent cities into
  # their counties -- and is fine, since we join fips -> gid_2.
  xw_raw <- read_csv(XWALK, show_col_types = FALSE) |>
    mutate(fips = pad_fips(fips)) |>
    select(fips, gid_2, any_of(c("gadm_name", "state"))) |>
    filter(!is.na(fips), !is.na(gid_2)) |>
    distinct()
  if (!"gadm_name" %in% names(xw_raw)) xw_raw$gadm_name <- NA_character_
  
  multi <- xw_raw |> distinct(fips, gid_2) |> count(fips) |> filter(n > 1)
  if (nrow(multi) > 0)
    message(sprintf("[xwalk] dropping %d fips mapping to multiple GADM units",
                    nrow(multi)))
  xw <- xw_raw |> anti_join(multi, by = "fips")
  
  cen <- read_csv(CENTROIDS, show_col_types = FALSE) |>
    mutate(fips = pad_fips(fips)) |>
    select(fips, lat, lon, pop = cen_pop)
  
  # Both kernel files key on FIPS. The crosswalk contributes gid_2 for
  # geometry only, so a county missing from it has no polygon and cannot be
  # an origin -- though it can still be a destination, since destination
  # units come from the GADM layer itself.
  county <- cen |>
    inner_join(ver, by = "fips", relationship = "one-to-one") |>
    left_join(meta, by = "fips", relationship = "one-to-one") |>
    left_join(xw, by = "fips", relationship = "one-to-one")
  
  # Reliable independent estimation of a Veraset kernel needs roughly 10,000
  # devices (runaway rates 48% under 1,000, 27.5% at 1-2.5k, 7.7% at 2.5-5k,
  # 0.6% at 5-10k, 0% above). A reliability threshold, not a data floor --
  # hierarchical partial pooling would let smaller counties contribute with
  # shrinkage. n_devices now comes from the params file itself; previously it
  # came only from an optional volume file, so this flag was silently FALSE
  # everywhere when that file was missing.
  county <- county |>
    mutate(kernel_reliable = !is.na(n_devices) & n_devices >= 10000 &
             !coalesce(at_bound, FALSE))
  
  message(sprintf("[inputs] %d counties | %d with a Meta kernel | %d with a polygon | %d with a reliable Veraset kernel",
                  nrow(county), sum(!is.na(county$mu_m)), sum(!is.na(county$gid_2)),
                  sum(county$kernel_reliable, na.rm = TRUE)))
  message(sprintf("[inputs] P(d < 10 km): Meta median %.3f, Veraset median %.3f",
                  median(lognormal_cdf(10, county$mu_m, county$sigma_m), na.rm = TRUE),
                  median(lognormal_cdf(10, county$mu_v, county$sigma_v), na.rm = TRUE)))
  county
}

COUNTY <- load_inputs()

# ===========================================================================
# 4. THE SPATIAL DOMAIN FOR ONE ORIGIN COUNTY
# ---------------------------------------------------------------------------
# Destination cells: population within D_MAX of ANY point of the county,
# tagged with the destination county they fall in. The domain is padded by
# the county's own reach, so a cell on the far edge of a large county
# (San Bernardino) still sees a full D_MAX radius.
# Origin cells: the destination cells that fall inside the origin county,
# so each origin cell is its own destination at distance 0.
# Cached, because the domain depends on neither kernel nor volume.
# ===========================================================================

build_domain <- function(fips) {
  fips <- pad_fips(fips)
  cache <- file.path(CACHE_DIR, sprintf("domain_slices_%s_%dkm_%dkm.rds",
                                        fips, D_MAX, GRID_KM))
  if (file.exists(cache)) return(readRDS(cache))
  
  row <- COUNTY |> filter(fips == !!fips)
  stopifnot("Origin county not found in COUNTY." = nrow(row) == 1)
  
  adm2 <- geodata::gadm(country = "USA", level = 2, path = GEODATA_CACHE)
  wp <- geodata::population(year = 2020, res = 0.5, path = GEODATA_CACHE)
  
  poly <- adm2[adm2$GID_2 == row$gid_2, ]
  stopifnot("No GADM polygon for this county." = nrow(poly) > 0)
  
  c_lon <- row$lon; c_lat <- row$lat
  
  # how far the county itself reaches from its centroid
  e <- as.vector(terra::ext(poly))
  reach_km <- max(hav_km(c_lon, c_lat, e[c(1, 2, 1, 2)], e[c(3, 3, 4, 4)]))
  R_km <- D_MAX + reach_km
  
  pad_lat <- R_km / 111
  pad_lon <- R_km / (111 * cos(min(89, abs(c_lat) + pad_lat) * RAD))
  dom_ext <- terra::ext(c(c_lon - pad_lon, c_lon + pad_lon,
                          max(-90, c_lat - pad_lat), min(90, c_lat + pad_lat)))
  
  # people per cell, then aggregate to ~GRID_KM
  pop_d <- terra::crop(wp, dom_ext)
  pop_d <- pop_d * terra::cellSize(pop_d, unit = "km")
  fact <- max(1L, round(GRID_KM / (terra::res(pop_d)[1] * 111)))
  if (fact > 1L)
    pop_d <- terra::aggregate(pop_d, fact = fact, fun = "sum", na.rm = TRUE)
  
  # slice width: the geometric mean cell side at the origin's latitude
  res_deg <- terra::res(pop_d)
  slice_km <- sqrt(res_deg[2] * 111.32 * res_deg[1] * 111.32 * cos(c_lat * RAD))
  
  adm2_dom <- terra::crop(adm2, dom_ext)
  adm2_dom$county_idx <- seq_len(nrow(adm2_dom))
  cidx <- terra::rasterize(adm2_dom, pop_d, field = "county_idx")
  lut <- tibble(county_idx = adm2_dom$county_idx, gid_2 = adm2_dom$GID_2)
  
  dest <- terra::as.data.frame(pop_d, xy = TRUE, na.rm = TRUE) |>
    as_tibble() |>
    rename(lon = 1, lat = 2, pop = 3) |>
    filter(pop > 0) |>
    mutate(county_idx = terra::extract(cidx, cbind(lon, lat))[, 1]) |>
    filter(!is.na(county_idx)) |>
    left_join(lut, by = "county_idx") |>
    filter(!is.na(gid_2)) |>
    mutate(d_centre = hav_km(c_lon, c_lat, lon, lat)) |>
    filter(d_centre <= R_km) |>
    select(lon, lat, pop, gid_2)
  
  orig <- dest |> filter(gid_2 == row$gid_2) |> select(lon, lat, pop)
  
  stopifnot("Empty origin or destination grid." =
              nrow(orig) > 0 && nrow(dest) >= 10)
  
  dom <- list(fips = fips, gid_2 = row$gid_2, lon = c_lon, lat = c_lat,
              orig = orig, dest = dest, grid_km = GRID_KM, slice_km = slice_km,
              radius_km = D_MAX, reach_km = reach_km, built_on = Sys.time())
  saveRDS(dom, cache)
  dom
}

# ===========================================================================
# 5. EXPECTED ALLOCATION
# ---------------------------------------------------------------------------
# One pass over origin cells. Distances are shared across kernels, so both
# are evaluated inside the same loop.
#
# Per origin cell i and slice s:
# mass P_s from the window-normalised kernel
# split across cells in s in proportion to their population
#
# Cells further than D_MAX from the origin cell are excluded. Slices with no
# populated destination (ocean, foreign) return their mass to the pool: the
# per-origin weights are renormalised to sum to 1. empty_mass reports how
# much that was, population-weighted over origin cells; with cell-width
# slices it should be coastline and borders only.
# ===========================================================================

allocate_expected <- function(dom, kernels) {
  
  g <- dom$slice_km
  edges <- slice_edges(g)
  nb <- length(edges) - 1L
  
  bin_mass <- map(kernels, slice_mass, edges = edges)
  
  gid_levels <- sort(unique(dom$dest$gid_2))
  gidx <- match(dom$dest$gid_2, gid_levels)
  n_gid <- length(gid_levels)
  dpop <- dom$dest$pop
  
  acc <- set_names(rep(list(numeric(n_gid)), length(kernels)), names(kernels))
  prof <- set_names(rep(list(numeric(nb)), length(kernels)), names(kernels))
  empty <- set_names(numeric(length(kernels)), names(kernels))
  
  w_orig <- dom$orig$pop / sum(dom$orig$pop)
  
  for (i in seq_len(nrow(dom$orig))) {
    d <- hav_km(dom$orig$lon[i], dom$orig$lat[i], dom$dest$lon, dom$dest$lat)
    keep <- d <= D_MAX
    b <- slice_of(d, g)
    
    agg <- rowsum(dpop[keep], b[keep]) # population per occupied slice
    S <- numeric(nb)
    S[as.integer(rownames(agg))] <- as.numeric(agg)
    
    for (nm in names(kernels)) {
      pb <- bin_mass[[nm]]
      w <- ifelse(keep & S[b] > 0, pb[b] * dpop / S[b], 0)
      tot <- sum(w)
      if (!is.finite(tot) || tot <= 0) next
      
      empty[[nm]] <- empty[[nm]] + w_orig[i] * max(0, 1 - tot / sum(pb))
      w <- w * (w_orig[i] / tot)
      
      ac <- rowsum(w, gidx)
      ai <- as.integer(rownames(ac))
      acc[[nm]][ai] <- acc[[nm]][ai] + as.numeric(ac)
      
      pc <- rowsum(w, b)
      pci <- as.integer(rownames(pc))
      prof[[nm]][pci] <- prof[[nm]][pci] + as.numeric(pc)
    }
  }
  
  list(
    shares = map(acc, share_table, gid_levels = gid_levels),
    profile = map(prof, profile_table, edges = edges),
    empty_mass = empty
  )
}

share_table <- function(x, gid_levels) tibble(gid_2 = gid_levels, share = as.numeric(x))
profile_table <- function(x, edges) tibble(d_lo = edges[-length(edges)],
                                           d_hi = edges[-1], mass = as.numeric(x))

# ===========================================================================
# 6. ABSOLUTE FLOWS FOR ONE ORIGIN COUNTY
# ===========================================================================

REQUIRED_DIAGNOSTICS <- c("veraset_scope", "n_devices", "kernel_reliable", "tv_v",
                          "under_10km_v", "under_10km_m", "tail_beyond_v",
                          "tail_beyond_m", "away_v", "away_m", "slice_km",
                          "empty_mass_v", "empty_mass_m", "n_origin_cells",
                          "n_dest_cells")

county_flows <- function(fips, refresh = FALSE) {
  fips <- pad_fips(fips)
  out <- file.path(OUT_DIR, sprintf("flows_%s_%s_Npop_%s.rds",
                                    fips, VERASET_SCOPE, CACHE_VERSION))
  if (file.exists(out) && !refresh) {
    cached <- readRDS(out)
    # Belt and braces: a cache written by a version with a different set of
    # diagnostics is rebuilt rather than used.
    if (all(REQUIRED_DIAGNOSTICS %in% names(cached$diagnostics))) return(cached)
    message(fips, ": cached result predates the current diagnostics; rebuilding.")
  }
  
  row <- COUNTY |> filter(fips == !!fips)
  stopifnot("Origin county not found." = nrow(row) == 1)
  stopifnot("No Meta kernel for this county." = !is.na(row$mu_m))
  stopifnot("No gid_2 for this county, so no origin polygon." = !is.na(row$gid_2))
  if (!isTRUE(row$kernel_reliable))
    warning(sprintf("%s: n_devices = %s, below the ~10,000 needed to estimate a Veraset kernel independently.",
                    fips, format(row$n_devices, big.mark = ",")))
  
  dom <- build_domain(fips)
  
  kernels <- list(
    veraset = list(mu = row$mu_v, sigma = row$sigma_v),
    meta = list(mu = row$mu_m, sigma = row$sigma_m)
  )
  alloc <- allocate_expected(dom, kernels)
  
  # ---- LEVEL ---------------------------------------------------------------
  # N = pop. Each source keeps the share of everyone's time its own kernel
  # places inside the domain; the rest is the tail beyond D_MAX.
  in_window <- c(
    veraset = lognormal_cdf(D_MAX, row$mu_v, row$sigma_v) -
      lognormal_cdf(D_MIN, row$mu_v, row$sigma_v),
    meta = lognormal_cdf(D_MAX, row$mu_m, row$sigma_m) -
      lognormal_cdf(D_MIN, row$mu_m, row$sigma_m)
  )
  N <- row$pop * in_window
  
  flows <- alloc$shares$veraset |>
    rename(share_v = share) |>
    full_join(alloc$shares$meta |> rename(share_m = share), by = "gid_2") |>
    mutate(across(c(share_v, share_m), \(x) replace_na(x, 0)),
           flow_v = N[["veraset"]] * share_v,
           flow_m = N[["meta"]] * share_m,
           log_ratio_total = log(flow_v / flow_m),
           log_ratio_shape = log(share_v / share_m)) |>
    arrange(desc(flow_v))
  
  res <- list(
    fips = fips, gid_2 = row$gid_2, county_name = row$gadm_name, pop = row$pop,
    kernels = kernels,
    volume = list(
      rule = "N = pop * P(d <= D_MAX)",
      f0 = row$f0, # diagnostic only; not in the volume
      in_window = in_window, N = N
    ),
    diagnostics = list(
      veraset_scope = VERASET_SCOPE,
      n_devices = row$n_devices, kernel_reliable = row$kernel_reliable,
      tv_v = row$tv,
      # how much of each kernel sits under 10 km -- the range where Meta's
      # (0,10) category leaves the split unconstrained by construction, and
      # where 18.8 of the 21.2 points of misfit were found
      under_10km_v = lognormal_cdf(10, row$mu_v, row$sigma_v),
      under_10km_m = lognormal_cdf(10, row$mu_m, row$sigma_m),
      tail_beyond_v = 1 - in_window[["veraset"]],
      tail_beyond_m = 1 - in_window[["meta"]],
      # Share of the origin's residents each kernel places OUTSIDE the
      # origin county at a random moment. The in-domain shares sum to 1, so
      # this is 1 - in_window * self_share: everything not in the origin's
      # own polygon, the tail beyond D_MAX included.
      away_v = 1 - in_window[["veraset"]] *
        sum(alloc$shares$veraset$share[alloc$shares$veraset$gid_2 == row$gid_2]),
      away_m = 1 - in_window[["meta"]] *
        sum(alloc$shares$meta$share[alloc$shares$meta$gid_2 == row$gid_2]),
      slice_km = dom$slice_km,
      empty_mass_v = alloc$empty_mass[["veraset"]],
      empty_mass_m = alloc$empty_mass[["meta"]],
      n_origin_cells = nrow(dom$orig), n_dest_cells = nrow(dom$dest)
    ),
    flows = flows, profile = alloc$profile, built_on = Sys.time()
  )
  
  tmp <- paste0(out, ".tmp"); saveRDS(res, tmp); file.rename(tmp, out)
  res
}

# ===========================================================================
# 7. GEOMETRY FOR MAPPING
# ===========================================================================

.county_sf <- NULL
get_county_sf <- function() {
  if (!is.null(.county_sf)) return(.county_sf)
  adm2 <- geodata::gadm(country = "USA", level = 2, path = GEODATA_CACHE)
  sfd <- sf::st_as_sf(adm2) |>
    select(gid_2 = GID_2, gadm_name = NAME_2, state = NAME_1) |>
    sf::st_transform(5070) # Conus Albers, equal area
  .county_sf <<- sfd
  sfd
}

# `extra` is added before coord_sf; a geom_sf tacked on afterwards brings its
# own default coord and ggplot warns about replacing the one already there.
map_base <- function(g, bb, extra = NULL) {
  g + geom_sf(colour = NA) + extra +
    coord_sf(xlim = bb[c("xmin", "xmax")], ylim = bb[c("ymin", "ymax")],
             expand = FALSE) +
    theme_void(base_size = 10) +
    theme(legend.position = "bottom",
          legend.key.width = unit(22, "pt"),
          legend.key.height = unit(7, "pt"))
}

# ===========================================================================
# 8. THE ENTRY POINT
# ---------------------------------------------------------------------------
# compare_county("06037")
#
# Three panels:
# A absolute Veraset-kernel flows by destination county, log10
# B absolute Meta-kernel flows, same scale
# C log discrepancy. Both sides carry N = pop, so this is redistribution
# plus the tail-mass difference; the residual is the tail term.
# ===========================================================================

# min_flow is in people. Counties where neither source places at least this
# many are drawn grey rather than coloured: at 1e-3 the map is dominated by
# hundreds of distant counties receiving a fraction of a person, where the
# log ratio is extreme and meaningless, and a dozen counties holding all the
# people are outvoted by sheer polygon count.
compare_county <- function(fips, which = c("total", "shape"),
                           min_flow = 1, refresh = FALSE) {
  
  which <- match.arg(which)
  res <- county_flows(fips, refresh = refresh)
  sfd <- get_county_sf()
  d <- res$diagnostics
  v <- res$volume
  
  message(sprintf(
    "
origin %s (%s)
population %s
volume N = pop (home is inside the kernel, not the volume; f0 %s)
in-domain Veraset %s Meta %s
kernels Veraset mu %+.3f sigma %.3f | Meta mu %+.3f sigma %.3f
away from it Veraset %s (%.2f%%) Meta %s (%.2f%%) <- residents outside the origin county
%s places more away, by a factor of %.2f <- flow-weighted; individual far
counties can and do run the other way
under 10 km Veraset %.3f Meta %.3f <- Meta's (0,10) leaves this split free
beyond %d km Veraset %.4f Meta %.4f <- dropped from the domain
slices %.2f km wide empty-slice mass Veraset %.4f Meta %.4f
panel %s devices%s
",
    res$fips, res$county_name %||% "", format(round(res$pop), big.mark = ","),
    ifelse(is.na(v$f0 %||% NA_real_), "not in the fit file",
           sprintf("%.3f, diagnostic only", v$f0)),
    format(round(v$N[["veraset"]]), big.mark = ","),
    format(round(v$N[["meta"]]), big.mark = ","),
    res$kernels$veraset$mu, res$kernels$veraset$sigma,
    res$kernels$meta$mu, res$kernels$meta$sigma,
    format(round(res$pop * d$away_v), big.mark = ","), 100 * d$away_v,
    format(round(res$pop * d$away_m), big.mark = ","), 100 * d$away_m,
    ifelse(d$away_v > d$away_m, "Veraset", "Meta"),
    max(d$away_v, d$away_m) / min(d$away_v, d$away_m),
    d$under_10km_v, d$under_10km_m, D_MAX, d$tail_beyond_v, d$tail_beyond_m,
    d$slice_km, d$empty_mass_v, d$empty_mass_m,
    format(d$n_devices, big.mark = ","),
    ifelse(isTRUE(d$kernel_reliable), "",
           " [below ~10,000 -- kernel not reliably estimated alone]")))
  
  # Two views of the same disagreement, because they answer different
  # questions. The log ratio is scale-free: a county receiving 20 instead of
  # 10 reads the same as one receiving 20,000 instead of 10,000, so it shows
  # where the kernels disagree proportionally -- usually the sparse far
  # field. The absolute difference is in people, so it shows where the
  # disagreement is actually large enough to matter for a flow matrix, which
  # is nearly always near the origin. A county can be loud in one and
  # invisible in the other.
  fl <- res$flows |>
    filter(flow_v > min_flow | flow_m > min_flow) |>
    mutate(value_log = if (which == "total") log_ratio_total else log_ratio_shape,
           value_abs = if (which == "total") flow_v - flow_m else share_v - share_m) |>
    filter(is.finite(value_log), is.finite(value_abs)) |>
    # log2: one unit is a doubling, which is easier to read off a map than
    # natural log, and it is what LOG2_LIM is expressed in.
    mutate(value_log = value_log / log(2),
           is_origin = gid_2 == res$gid_2)
  
  # The origin county is drawn in its own flat shade rather than on either
  # scale. It holds most of the population, so on the absolute panel it
  # would compress every other county to white, and on the log panel it is
  # not really a comparison of the same kind -- it is where each kernel puts
  # the time it does not send anywhere.
  mp <- sfd |>
    inner_join(fl, by = "gid_2") |>
    mutate(negligible = pmax(flow_v, flow_m) < min_flow,
           across(c(value_log, value_abs),
                  \(x) if_else(is_origin | negligible, NA_real_, x)))
  bb <- sf::st_bbox(mp)
  
  # The map colours every county the same regardless of how many people it
  # receives, so it answers "where do the kernels disagree?" and NOT "which
  # kernel sends more?". Those can point opposite ways: Veraset's heavier
  # tail makes it larger in most far counties while Meta places more people
  # away in total. This is the flow-weighted direction, which is the one the
  # away counts in the header describe.
  wlog2 <- log2(sum(mp$flow_v[!mp$is_origin]) / sum(mp$flow_m[!mp$is_origin]))
  n_shown <- sum(!mp$is_origin & !mp$negligible)
  
  lim_abs <- mp |>
    sf::st_drop_geometry() |>
    pull(value_abs) |>
    quantile(c(0.02, 0.98), na.rm = TRUE) |>
    abs() |> max()
  
  origin_layer <- geom_sf(data = filter(mp, is_origin), fill = "grey55",
                          colour = "black", linewidth = 0.4)
  
  pA <- map_base(ggplot(mp, aes(fill = value_log)), bb, origin_layer) +
    scale_fill_gradient2(low = "#2166AC", mid = "grey94", high = "#B2182B",
                         midpoint = 0, limits = c(-LOG2_LIM, LOG2_LIM),
                         breaks = seq(-LOG2_LIM, LOG2_LIM, by = 1),
                         oob = scales::squish, name = NULL) +
    labs(title = sprintf("log2(Veraset / Meta), %s", which),
         subtitle = "per county, unweighted")
  
  pB <- map_base(ggplot(mp, aes(fill = value_abs)), bb, origin_layer) +
    scale_fill_gradient2(low = "#2166AC", mid = "grey94", high = "#B2182B",
                         midpoint = 0, limits = c(-lim_abs, lim_abs),
                         oob = scales::squish, name = NULL) +
    labs(title = sprintf("Veraset - Meta, %s", which),
         subtitle = if (which == "total") "people; where the weight is"
         else "share points; where the weight is")
  
  out <- if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(pA, pB, nrow = 1) +
      patchwork::plot_annotation(
        title = sprintf("Absolute flows from %s (%s)",
                        res$county_name %||% "", res$fips),
        subtitle = sprintf(
          "away from the origin county: Veraset %s (%.1f%%), Meta %s (%.1f%%), of %s residents\nflow-weighted overall: log2 %+.2f | %d counties coloured, origin and under %g person grey | log scale fixed +/-%g",
          format(round(res$pop * d$away_v), big.mark = ","), 100 * d$away_v,
          format(round(res$pop * d$away_m), big.mark = ","), 100 * d$away_m,
          format(round(res$pop), big.mark = ","), wlog2, n_shown, min_flow, LOG2_LIM),
        theme = ggplot2::theme(
          plot.title = element_text(face = "bold", size = 13)))
  } else {
    list(log_ratio = pA, difference = pB)
  }
  
  print(res$flows |> slice_head(n = 12) |>
          transmute(gid_2, flow_veraset = round(flow_v),
                    flow_meta = round(flow_m),
                    diff = round(flow_v - flow_m),
                    log_shape = round(log_ratio_shape, 3)))
  
  invisible(list(plot = out, result = res, mapped = mp))
}

# ===========================================================================
# 9. BATCH, AND THE CROSS-COUNTY SUMMARY
# ---------------------------------------------------------------------------
# Reduce each origin to one number: the total-variation distance between the
# two allocations, i.e. the share of residents' time the two kernels put in
# different counties. Same currency as the 4.2 / 8.1 / 21.2 percent figures.
# ===========================================================================

# A missing field would otherwise be length zero, and tibble() recycles
# every other column down to match it -- so one absent diagnostic silently
# returns no rows at all rather than raising anything. Force every entry to
# be exactly one value.
scalar <- function(x) if (length(x) != 1) NA_real_ else as.numeric(x)

summarise_county <- function(f, refresh = FALSE) {
  r <- tryCatch(county_flows(f, refresh = refresh),
                error = conditionMessage)
  if (is.character(r)) { message(f, ": ", r); return(NULL) }
  self_v <- r$flows$share_v[r$flows$gid_2 == r$gid_2]
  self_m <- r$flows$share_m[r$flows$gid_2 == r$gid_2]
  off <- r$flows |> filter(gid_2 != r$gid_2)
  out <- tibble(
    fips = f,
    pop = scalar(r$pop),
    f0 = scalar(r$volume$f0),
    tv_fit_v = scalar(r$diagnostics$tv_v),
    N_v = scalar(r$volume$N[["veraset"]]),
    N_m = scalar(r$volume$N[["meta"]]),
    n_devices = scalar(r$diagnostics$n_devices),
    kernel_reliable = isTRUE(r$diagnostics$kernel_reliable),
    mu_v = scalar(r$kernels$veraset$mu),
    sigma_v = scalar(r$kernels$veraset$sigma),
    mu_m = scalar(r$kernels$meta$mu),
    sigma_m = scalar(r$kernels$meta$sigma),
    under_10km_v = scalar(r$diagnostics$under_10km_v),
    under_10km_m = scalar(r$diagnostics$under_10km_m),
    tv_allocation = 0.5 * sum(abs(r$flows$share_v - r$flows$share_m)),
    self_share_v = scalar(self_v),
    self_share_m = scalar(self_m),
    away_v = scalar(r$diagnostics$away_v),
    away_m = scalar(r$diagnostics$away_m),
    away_n_v = scalar(r$pop) * scalar(r$diagnostics$away_v),
    away_n_m = scalar(r$pop) * scalar(r$diagnostics$away_m),
    away_log_ratio = log(scalar(r$diagnostics$away_v) / scalar(r$diagnostics$away_m)),
    max_abs_diff_off = if (nrow(off)) max(abs(off$flow_v - off$flow_m)) else NA_real_,
    empty_mass_m = scalar(r$diagnostics$empty_mass_m)
  )
  stopifnot("Summary row collapsed; a field was missing." = nrow(out) == 1)
  out
}

run_all <- function(fips_vec = NULL, refresh = FALSE) {
  if (is.null(fips_vec))
    fips_vec <- COUNTY |> filter(!is.na(mu_m), !is.na(mu_v), !is.na(gid_2)) |> pull(fips)
  out <- map(fips_vec, summarise_county, refresh = refresh) |> list_rbind()
  if (nrow(out) < length(fips_vec))
    message(sprintf("[run_all] %d of %d counties returned no row",
                    length(fips_vec) - nrow(out), length(fips_vec)))
  out
}

# ---------------------------------------------------------------------------
# usage
#
x <- compare_county("06037", min_flow = 50); x$plot # LA, level + shape
# compare_county("30109", which = "shape") # Wibaux, redistribution only
# summ <- run_all()
# summ |> arrange(desc(tv_allocation)) |> print(n = 40)
#
# summary(lm(tv_allocation ~ I(under_10km_v - under_10km_m), data = summ))
#
# If tv_allocation tracks the near-field difference and not the tail, the
# discrepancy is the sub-10 km split that Meta's (0,10) category leaves
# unconstrained -- not a disagreement about long-range travel.
# ---------------------------------------------------------------------------
