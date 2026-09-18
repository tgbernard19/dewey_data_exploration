# ===========================================================================
# R/07_flows.R
# ---------------------------------------------------------------------------
# From kernels to people on the map: how many of an origin county's residents
# are in each destination county at a random moment, under each source.
#
# THREE SOURCES, NOT TWO
#   meta      Meta's kernel, spread by the population rule. The deliverable:
#             it is the only one that exists outside the US.
#   veraset   Veraset's kernel, spread by the same rule. Optional.
#   obs       Veraset's actual destinations. No kernel, no allocation rule.
#
#   obs vs veraset   error from the allocation rule alone (same data)
#   veraset vs meta  difference between the kernels alone (same rule)
#   obs vs meta      the two combined -- what the deliverable costs
#
# THE VOLUME RULE: N = pop
#   Everyone is somewhere at every moment, including at home. The kernel is
#   unconditional and time-weighted -- SCOPE = "all" on the Veraset side, the
#   folded home tile on Meta's -- so the population at risk of being anywhere
#   is the whole population.
#
#   Each source keeps the share of that population its own kernel places
#   inside the domain, N = pop * P(d <= D_MAX); the rest is tail beyond the
#   window. Note this uses the UNTRUNCATED CDF while the allocation uses the
#   truncated one. Both are deliberate: the volume asks how much mass is in
#   the window, the allocation asks how in-window mass is distributed.
#
#   f0 plays no part. Under the previous rule, N = pop * (1 - f0), it was
#   load-bearing, and folding the home tile into the kernel while also
#   removing it from the volume subtracted the same people twice. That was
#   the -0.42 offset.
#
# ON THE TABLES
#   COUNTY, OBSERVED and XW are module-level, inherited from the original
#   script. Every function here takes them as arguments defaulting to those
#   globals, so an explicit call can override them -- which matters, because
#   a stale COUNTY in the workspace is otherwise invisible until a column
#   goes missing mid-run.
# ===========================================================================


# ===========================================================================
# 1. INPUTS
# ---------------------------------------------------------------------------
# Column contract. If a file names things differently on some machine, fix it
# here, in one place.
#
#   META_PARAMS_FILE      fips, mu, log_sigma (or sigma)        required
#   VERASET_PARAMS_FILE   home_county (or fips), scope, mu, sigma  optional
#   OBSERVED_PAIRS_FILE   origin_fips, dest_fips, weights          optional
#   XWALK_FILE            fips, gid_2
#   COUNTY_CENTROIDS_FILE fips, lat, lon, cen_pop
#
# Both kernel files are keyed on FIPS, so neither is joined through the
# crosswalk. The crosswalk is needed only for geometry: gid_2 picks the
# origin polygon and labels destination cells.
#
# The join is a LEFT join from the county frame, with Meta required and
# Veraset optional. An earlier version inner-joined Veraset first, which
# silently dropped every county without a Veraset kernel before Meta was
# even considered -- so a Meta-only run still produced a Veraset-shaped map.
#
# Optional columns are added to the frame BEFORE any transmute. Testing for
# them inside one with names(.data) does not work: .data is a pronoun, not
# the data frame, so every such test silently takes its else branch.
# ===========================================================================

load_flow_inputs <- function() {
  
  cen <- read_keyed_csv(COUNTY_CENTROIDS_FILE, key_cols = "fips") |>
    dplyr::select(fips, lat, lon, pop = cen_pop) |>
    dplyr::distinct(fips, .keep_all = TRUE)
  
  xw <- read_keyed_csv(XWALK_FILE, key_cols = "fips") |>
    dplyr::select(fips, gid_2, dplyr::any_of("gadm_name")) |>
    dplyr::distinct(fips, .keep_all = TRUE)
  
  # ---- Meta (required) -----------------------------------------------------
  meta_raw <- read_keyed_csv(META_PARAMS_FILE, key_cols = "fips")
  if (!"sigma" %in% names(meta_raw)) meta_raw$sigma <- exp(meta_raw$log_sigma)
  if (!"f0" %in% names(meta_raw))    meta_raw$f0 <- NA_real_
  
  meta <- meta_raw |>
    dplyr::transmute(fips = pad_fips(fips), mu_m = mu, sigma_m = sigma, f0) |>
    dplyr::distinct(fips, .keep_all = TRUE)
  
  county <- cen |>
    dplyr::left_join(xw, by = "fips") |>
    dplyr::left_join(meta, by = "fips")
  
  # ---- Veraset (optional) --------------------------------------------------
  if (file.exists(VERASET_PARAMS_FILE)) {
    
    ver <- read_keyed_csv(VERASET_PARAMS_FILE,
                          key_cols = c("home_county", "fips"))
    
    # Current files key on home_county; older ones on fips.
    key <- if ("home_county" %in% names(ver)) "home_county" else "fips"
    ver$fips <- pad_fips(ver[[key]])
    
    if ("scope" %in% names(ver)) ver <- dplyr::filter(ver, scope == SCOPE)
    for (nm in c("n_devices", "tv")) {
      if (!nm %in% names(ver)) ver[[nm]] <- NA_real_
    }
    if (!"at_bound" %in% names(ver)) ver$at_bound <- NA
    
    ver <- ver |>
      dplyr::transmute(fips, mu_v = mu, sigma_v = sigma,
                       n_devices, at_bound, tv_v = tv) |>
      dplyr::distinct(fips, .keep_all = TRUE)
    
    county <- dplyr::left_join(county, ver, by = "fips")
    
  } else {
    county <- dplyr::mutate(county, mu_v = NA_real_, sigma_v = NA_real_,
                            n_devices = NA_real_, at_bound = NA, tv_v = NA_real_)
  }
  
  county |>
    dplyr::mutate(
      kernel_reliable = !is.na(n_devices) & n_devices >= MIN_DEVICES_KERNEL &
        !isTRUE(at_bound),
      has_meta    = !is.na(mu_m) & !is.na(sigma_m),
      has_veraset = !is.na(mu_v) & !is.na(sigma_v)
    )
}

# The columns every downstream function assumes. Checked once, with a message
# that names what is missing, rather than warning per access and then failing
# somewhere unrelated.
COUNTY_REQUIRED <- c("fips", "gid_2", "pop", "mu_m", "sigma_m", "f0",
                     "mu_v", "sigma_v", "n_devices", "kernel_reliable",
                     "has_meta", "has_veraset")

check_county <- function(county) {
  missing <- setdiff(COUNTY_REQUIRED, names(county))
  if (length(missing) > 0) {
    stop("The county table is missing: ", paste(missing, collapse = ", "),
         "\nIt was probably built by an older load_flow_inputs(). ",
         "Run COUNTY <- load_flow_inputs() again.")
  }
  invisible(TRUE)
}


# ===========================================================================
# 2. SLICES
# ---------------------------------------------------------------------------
# The kernel is one-dimensional in distance; the map is two-dimensional. The
# bridge is a set of distance slices whose width is the grid cell size, so a
# slice is never thinner than the resolution at which population is known.
#
# The first slice is a disc of radius g/sqrt(pi): the radius of a circle with
# the same area as one grid cell. That is the "self" slice, and it carries
# the bulk of the mass under every kernel here.
# ===========================================================================

slice_edges <- function(g) {
  r0 <- g / sqrt(pi)
  K <- ceiling((D_MAX - r0) / g)
  pmin(c(0, r0 + g * (0:K)), D_MAX)
}

slice_of <- function(d, g) {
  r0 <- g / sqrt(pi)
  K <- ceiling((D_MAX - r0) / g)
  ifelse(d < r0, 1L, as.integer(pmin(ceiling((d - r0) / g), K)) + 1L)
}

slice_mass <- function(k, edges) {
  diff(ln_cdf_trunc(edges, k$mu, k$sigma, lo = D_MIN, hi = D_MAX))
}


# ===========================================================================
# 3. THE SPATIAL DOMAIN FOR ONE ORIGIN COUNTY
# ---------------------------------------------------------------------------
# Destination cells: population within D_MAX of ANY point of the county,
# tagged with the destination county they fall in. Padded by the county's own
# reach, so a cell on the far edge of a large county still sees a full D_MAX
# radius.
#
# Origin cells: the destination cells inside the origin county, so each
# origin cell is its own destination at distance zero.
#
# Cached: the domain depends on neither kernel nor volume, only on geometry.
# ===========================================================================

build_domain <- function(fips, county = COUNTY) {
  
  fips <- pad_fips(fips)
  cache <- file.path(FLOW_CACHE_DIR,
                     sprintf("domain_slices_%s_%dkm_%dkm.rds",
                             fips, D_MAX, GRID_KM))
  if (file.exists(cache)) return(readRDS(cache))
  
  dir.create(FLOW_CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
  
  row <- dplyr::filter(county, fips == !!fips)
  stopifnot("Origin county not found in the county table." = nrow(row) == 1)
  stopifnot("No gid_2 for this county, so no origin polygon." = !is.na(row$gid_2))
  
  adm2 <- geodata::gadm(country = "USA", level = 2, path = GEODATA_DIR)
  wp   <- geodata::population(year = 2020, res = 0.5, path = GEODATA_DIR)
  
  poly <- adm2[adm2$GID_2 == row$gid_2, ]
  stopifnot("No GADM polygon for this county." = nrow(poly) == 1)
  
  # Bounding box padded by the county's extent plus D_MAX, in degrees at this
  # latitude. Longitude degrees shrink with latitude, which matters north of
  # about 45.
  ext <- terra::ext(poly)
  lat_pad <- D_MAX / 111.32
  lon_pad <- D_MAX / (111.32 * cos(row$lat * pi / 180))
  box <- terra::ext(ext$xmin - lon_pad, ext$xmax + lon_pad,
                    ext$ymin - lat_pad, ext$ymax + lat_pad)
  
  pop <- terra::crop(wp, box)
  cells <- terra::as.data.frame(pop, xy = TRUE, na.rm = TRUE)
  names(cells) <- c("lon", "lat", "pop")
  cells <- dplyr::filter(cells, pop > 0)
  
  pts <- terra::vect(cells, geom = c("lon", "lat"), crs = terra::crs(adm2))
  hit <- terra::extract(adm2, pts)
  cells$gid_2 <- hit$GID_2
  
  dest <- cells |>
    dplyr::filter(!is.na(gid_2)) |>
    tibble::as_tibble()
  
  orig <- dplyr::filter(dest, gid_2 == row$gid_2)
  stopifnot("Origin county has no populated cells." = nrow(orig) > 0)
  
  # Realised cell width at this latitude, in km. Used as the slice width, so
  # slices and cells are at the same resolution.
  slice_km <- terra::res(pop)[1] * 111.32 * cos(row$lat * pi / 180)
  
  dom <- list(fips = fips, gid_2 = row$gid_2,
              dest = dest, orig = orig, slice_km = slice_km)
  
  tmp <- paste0(cache, ".tmp")
  saveRDS(dom, tmp)
  file.rename(tmp, cache)
  dom
}


# ===========================================================================
# 4. EXPECTED ALLOCATION
# ---------------------------------------------------------------------------
# One pass over origin cells. Distances are shared across kernels, so every
# kernel is evaluated inside the same loop.
#
# For origin cell i and slice s: take the mass P_s the kernel puts in that
# slice and split it across the cells in s in proportion to their population.
# That population-proportional split IS the allocation rule being tested
# against the observed data.
#
# Slices with no populated destination (ocean, across a border) return their
# mass to the pool: the per-origin weights are renormalised to sum to 1.
# empty_mass reports how much that was, population-weighted over origin
# cells, so a coastal county that quietly redistributes a third of its mass
# is visible rather than assumed away.
# ===========================================================================

allocate_expected <- function(dom, kernels) {
  
  g <- dom$slice_km
  edges <- slice_edges(g)
  nb <- length(edges) - 1L
  
  bin_mass <- purrr::map(kernels, slice_mass, edges = edges)
  
  gid_levels <- sort(unique(dom$dest$gid_2))
  gidx <- match(dom$dest$gid_2, gid_levels)
  n_gid <- length(gid_levels)
  dpop <- dom$dest$pop
  
  acc   <- purrr::set_names(rep(list(numeric(n_gid)), length(kernels)), names(kernels))
  prof  <- purrr::set_names(rep(list(numeric(nb)), length(kernels)), names(kernels))
  empty <- purrr::set_names(numeric(length(kernels)), names(kernels))
  
  w_orig <- dom$orig$pop / sum(dom$orig$pop)
  
  for (i in seq_len(nrow(dom$orig))) {
    
    d <- hav_km(dom$orig$lon[i], dom$orig$lat[i], dom$dest$lon, dom$dest$lat)
    keep <- d <= D_MAX
    b <- slice_of(d, g)
    
    agg <- rowsum(dpop[keep], b[keep])
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
    shares = purrr::map(acc, \(x) tibble::tibble(gid_2 = gid_levels,
                                                 share = as.numeric(x))),
    profile = purrr::map(prof, \(x) tibble::tibble(d_lo = edges[-length(edges)],
                                                   d_hi = edges[-1],
                                                   mass = as.numeric(x))),
    empty_mass = empty
  )
}


# ===========================================================================
# 5. THE OBSERVED LAYER
# ---------------------------------------------------------------------------
# No kernel and no allocation rule: where Veraset devices actually were.
#
# Destinations arrive as FIPS and are summed to gid_2, because GADM folds
# Virginia's independent cities into their surrounding counties and the
# predicted side is labelled by gid_2.
#
# ESTIMAND WARNING
#   The current file was built for the old volume rule: home rows excluded,
#   weights summed across devices rather than one user one vote. Under
#   N = pop those are not the same quantity as the predicted flows, and the
#   difference will look like model error. Until 05_observed_flows.R is
#   rewritten, the honest comparison is which = "shape" -- drop the origin
#   county and renormalise both sides, where the home rows never entered
#   either number.
# ===========================================================================

load_observed <- function(path = OBSERVED_PAIRS_FILE) {
  
  if (!file.exists(path)) return(NULL)
  
  obs <- read_keyed_csv(path, key_cols = c("origin_fips", "dest_fips"))
  
  if (!"includes_home" %in% names(obs) || !isTRUE(obs$includes_home[1])) {
    warning("Observed pairs were built WITHOUT home rows. Absolute levels are ",
            "not comparable to N = pop predictions; use which = \"shape\", ",
            "or rerun 05_observed_flows.R.", call. = FALSE)
  }
  
  weight_col <- dplyr::case_when(
    OBS_WEIGHT == "time"    & "w_time" %in% names(obs) ~ "w_time",
    OBS_WEIGHT == "time"                               ~ "w_dwell",
    OBS_WEIGHT == "devices"                            ~ "n_devices",
    TRUE                                               ~ NA_character_
  )
  stopifnot("No weight column for the requested OBS_WEIGHT" = !is.na(weight_col))
  
  obs |>
    dplyr::transmute(origin_fips, dest_fips,
                     w = .data[[weight_col]],
                     n_dev = n_devices)
}

# Observed destination shares for one origin, on the gid_2 labelling.
observed_shares <- function(obs, fips, xw) {
  
  if (is.null(obs)) return(NULL)
  
  o <- dplyr::filter(obs, origin_fips == pad_fips(fips))
  if (nrow(o) == 0) return(NULL)
  
  o |>
    dplyr::left_join(dplyr::select(xw, fips, gid_2),
                     by = dplyr::join_by(dest_fips == fips)) |>
    dplyr::filter(!is.na(gid_2)) |>
    dplyr::group_by(gid_2) |>
    dplyr::summarise(w = sum(w), n_dev = sum(n_dev), .groups = "drop") |>
    dplyr::mutate(share_obs = w / sum(w)) |>
    dplyr::select(gid_2, share_obs, n_dev)
}


# ===========================================================================
# 6. ONE COUNTY, ALL SOURCES
# ---------------------------------------------------------------------------
# CACHE_VERSION is part of the filename. Bump it whenever the stored object
# changes shape, so older results are simply not found rather than read back
# missing fields.
# ===========================================================================

CACHE_VERSION <- "v3"

county_flows <- function(fips, county = COUNTY, obs = OBSERVED, xw = XW,
                         refresh = FALSE) {
  
  fips <- pad_fips(fips)
  check_county(county)
  
  out <- file.path(FLOW_OUT_DIR,
                   sprintf("flows_%s_%s_Npop_%s.rds", fips, SCOPE, CACHE_VERSION))
  if (file.exists(out) && !refresh) return(readRDS(out))
  
  row <- dplyr::filter(county, fips == !!fips)
  stopifnot("Origin county not found." = nrow(row) == 1)
  stopifnot("No Meta kernel for this county." = isTRUE(row$has_meta))
  
  if (isTRUE(row$has_veraset) && !isTRUE(row$kernel_reliable)) {
    warning(sprintf("%s: n_devices = %s, below the %s needed to estimate a ",
                    fips, format(row$n_devices, big.mark = ","),
                    format(MIN_DEVICES_KERNEL, big.mark = ",")),
            "Veraset kernel independently.", call. = FALSE)
  }
  
  dom <- build_domain(fips, county)
  
  kernels <- list(meta = list(mu = row$mu_m, sigma = row$sigma_m))
  if (isTRUE(row$has_veraset)) {
    kernels$veraset <- list(mu = row$mu_v, sigma = row$sigma_v)
  }
  
  alloc <- allocate_expected(dom, kernels)
  
  # ---- level ---------------------------------------------------------------
  # Untruncated CDF: how much of everyone's time each kernel places inside
  # the window at all.
  in_window <- purrr::map_dbl(kernels, \(k)
                              ln_cdf(D_MAX, k$mu, k$sigma) - ln_cdf(D_MIN, k$mu, k$sigma))
  N <- row$pop * in_window
  
  flows <- purrr::imap(alloc$shares, \(s, nm)
                       dplyr::rename(s, !!paste0("share_", nm) := share)) |>
    purrr::reduce(dplyr::full_join, by = "gid_2")
  
  for (nm in names(kernels)) {
    flows[[paste0("flow_", nm)]] <-
      N[[nm]] * tidyr::replace_na(flows[[paste0("share_", nm)]], 0)
  }
  
  # ---- observed ------------------------------------------------------------
  obs_sh <- observed_shares(obs, fips, xw)
  
  if (!is.null(obs_sh)) {
    flows <- dplyr::full_join(flows, obs_sh, by = "gid_2") |>
      dplyr::mutate(share_obs = tidyr::replace_na(share_obs, 0),
                    flow_obs = row$pop * share_obs)
    
    # Sampling noise on an observed share, from the multinomial: the se of a
    # share p estimated from n devices is sqrt(p(1-p)/n). The z panel divides
    # the discrepancy by this, so a destination seen by six devices does not
    # dominate the colour scale purely by being thinly sampled.
    n_panel <- row$n_devices
    flows <- dplyr::mutate(
      flows,
      se_obs = sqrt(pmax(share_obs * (1 - share_obs), 0) / n_panel),
      z_meta = (share_meta - share_obs) / se_obs
    )
  }
  
  flows <- flows |>
    dplyr::mutate(dplyr::across(dplyr::starts_with("share_"),
                                \(x) tidyr::replace_na(x, 0)),
                  dplyr::across(dplyr::starts_with("flow_"),
                                \(x) tidyr::replace_na(x, 0))) |>
    dplyr::arrange(dplyr::desc(flow_meta))
  
  if ("flow_veraset" %in% names(flows)) {
    flows <- dplyr::mutate(flows,
                           log_ratio_total = log(flow_veraset / flow_meta),
                           log_ratio_shape = log(share_veraset / share_meta))
  }
  
  self <- purrr::set_names(
    purrr::map_dbl(names(kernels), \(nm)
                   sum(flows[[paste0("share_", nm)]][flows$gid_2 == row$gid_2])),
    names(kernels))
  
  res <- list(
    fips = fips, gid_2 = row$gid_2,
    county_name = row$gadm_name %||% NA_character_,
    pop = row$pop,
    kernels = kernels,
    flows = flows,
    profile = alloc$profile,
    volume = list(rule = "N = pop * P(d <= D_MAX)",
                  f0 = row$f0, in_window = in_window, N = N),
    diagnostics = list(
      scope = SCOPE,
      n_devices = row$n_devices,
      kernel_reliable = row$kernel_reliable,
      tv_v = row$tv_v,
      slice_km = dom$slice_km,
      n_origin_cells = nrow(dom$orig),
      n_dest_cells = nrow(dom$dest),
      empty_mass = alloc$empty_mass,
      # Share of residents each source places under 10 km -- the range Meta's
      # (0, 10) category leaves unconstrained by construction, and where most
      # of the misfit has been found.
      under_10km = purrr::map_dbl(kernels, \(k) ln_cdf(10, k$mu, k$sigma)),
      tail_beyond = 1 - in_window,
      # Residents outside the origin county at a random moment, tail included.
      away = 1 - in_window * self,
      has_observed = !is.null(obs_sh)
    ),
    built_on = Sys.time()
  )
  
  dir.create(FLOW_OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(out, ".tmp")
  saveRDS(res, tmp)
  file.rename(tmp, out)
  res
}


# ===========================================================================
# 7. CROSS-COUNTY SUMMARY
# ---------------------------------------------------------------------------
# One row per origin. The headline is total variation between two
# allocations: the share of residents' time they place in different counties.
# Same currency as the 4.2 / 8.1 percent fit figures.
# ===========================================================================

# A missing field would be length zero, and tibble() then recycles every
# other column down to nothing -- so one absent diagnostic silently returns
# no rows rather than raising. Force every entry to be one value.
scalar <- function(x) if (length(x) != 1) NA_real_ else as.numeric(x)

tv_between <- function(flows, a, b) {
  ca <- paste0("share_", a)
  cb <- paste0("share_", b)
  if (!all(c(ca, cb) %in% names(flows))) return(NA_real_)
  0.5 * sum(abs(flows[[ca]] - flows[[cb]]), na.rm = TRUE)
}

summarise_county <- function(f, county = COUNTY, obs = OBSERVED, xw = XW,
                             refresh = FALSE) {
  
  r <- tryCatch(county_flows(f, county = county, obs = obs, xw = xw,
                             refresh = refresh),
                error = conditionMessage)
  if (is.character(r)) { message(f, ": ", r); return(NULL) }
  
  d <- r$diagnostics
  off <- dplyr::filter(r$flows, gid_2 != r$gid_2)
  
  out <- tibble::tibble(
    fips = f,
    pop = scalar(r$pop),
    n_devices = scalar(d$n_devices),
    kernel_reliable = isTRUE(d$kernel_reliable),
    mu_m = scalar(r$kernels$meta$mu),
    sigma_m = scalar(r$kernels$meta$sigma),
    mu_v = scalar(r$kernels$veraset$mu),
    sigma_v = scalar(r$kernels$veraset$sigma),
    N_m = scalar(r$volume$N[["meta"]]),
    N_v = scalar(r$volume$N[["veraset"]]),
    under_10km_m = scalar(d$under_10km[["meta"]]),
    under_10km_v = scalar(d$under_10km[["veraset"]]),
    self_share_m = scalar(sum(r$flows$share_meta[r$flows$gid_2 == r$gid_2])),
    away_m = scalar(d$away[["meta"]]),
    away_v = scalar(d$away[["veraset"]]),
    # allocation rule error, kernel difference, and the two combined
    tv_obs_vs_veraset  = tv_between(r$flows, "obs", "veraset"),
    tv_veraset_vs_meta = tv_between(r$flows, "veraset", "meta"),
    tv_obs_vs_meta     = tv_between(r$flows, "obs", "meta"),
    max_abs_diff_off = if (nrow(off) && "flow_veraset" %in% names(off))
      max(abs(off$flow_veraset - off$flow_meta)) else NA_real_,
    empty_mass_m = scalar(d$empty_mass[["meta"]]),
    has_observed = isTRUE(d$has_observed)
  )
  
  stopifnot("Summary row collapsed; a field was missing." = nrow(out) == 1)
  out
}

run_all <- function(fips_vec = NULL, county = COUNTY, obs = OBSERVED, xw = XW,
                    refresh = FALSE) {
  
  check_county(county)
  
  if (is.null(fips_vec)) {
    fips_vec <- county |>
      dplyr::filter(has_meta, !is.na(gid_2)) |>
      dplyr::pull(fips)
  }
  
  out <- purrr::map(fips_vec,
                    \(f) summarise_county(f, county = county, obs = obs,
                                          xw = xw, refresh = refresh)) |>
    purrr::list_rbind()
  
  if (nrow(out) < length(fips_vec)) {
    message(sprintf("[run_all] %d of %d counties returned no row",
                    length(fips_vec) - nrow(out), length(fips_vec)))
  }
  out
}