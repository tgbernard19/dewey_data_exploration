### ---------------------------------------------------------------------
### Kernel family comparison, all counties
### VARIANT: configurable handling of Meta's home tile
###
### META_SCOPE controls what the Meta side measures:
###
###   "folded"  the home tile is folded into (0, 10). Meta then describes
###             where a device is at a random moment, INCLUDING at home,
###             which is the same quantity as Veraset SCOPE = "all".
###             Pair this with county_profile_cbghome_national_all_*.csv.
###
###   "away"    the home tile is dropped and the three away-categories are
###             renormalised, so Meta describes displacement CONDITIONAL on
###             having left home. Pair this with the trip-scope profile.
###
### These are not interchangeable. Running "away" against an all-scope
### Veraset profile is the estimand mismatch that produced a -0.111 median
### gap in the August work; matching them moved it to +0.007. The previous
### run of this script used "away" against an all-scope profile, which is
### why its Meta tv_near was 0.317 out of a 0.370 total.
###
### Meta's category "0" is an atom at a point. No continuous kernel can put
### mass on a point, so there is no way to fit it as a fourth bin - folding
### is the only coherent way to include it. The tile is a few km across, so
### folding it into (0, 10) is arithmetically exact rather than an
### approximation.
###
### f0 is carried through to the output either way. It is needed for the
### volume rule downstream, and its near-constancy across counties
### (~0.346) is itself worth watching.
### ---------------------------------------------------------------------

### LIBRARIES ###
library(tidyverse)

### CONSTANTS ###

META_SCOPE <- "folded"    # "folded" or "away" - see header

OUT_DIR      <- "E:/dewey-june2025/kernel/trackB"
PROFILE_FILE <- file.path(OUT_DIR, "county_profile_cbghome_national_all_dwell_tau60.csv")
META_FILE    <- "E:/meta_movement_dist/movement-distribution-1-june-2026_15-june-2026.csv"
XWALK_FILE   <- file.path("data", "processed", "county_gid2_crosswalk.csv")
FIT_FILE     <- file.path(OUT_DIR, paste0("kernel_family_fits_long_", META_SCOPE, ".csv"))

TRUNC_KM        <- 500
MIN_FIT         <- 500
MIN_TRUST       <- 10000
N_START         <- 10
START_SD        <- 0.75
META_OBS_WEIGHT <- 1000   # scales the Meta objective so BFGS has a gradient
SEED            <- 20260901

stopifnot(META_SCOPE %in% c("folded", "away"))
set.seed(SEED)
stopifnot(file.exists(PROFILE_FILE), file.exists(META_FILE), file.exists(XWALK_FILE))

### FAMILIES ###

families <- list(
  lognormal = list(
    npar   = 2,
    par0   = c(-3, log(4)),
    unpack = function(par) list(mu = par[1], sigma = exp(par[2])),
    surv   = function(d, th) 1 - pnorm((log(d) - th$mu) / th$sigma),
    dens   = function(d, th) dlnorm(d, th$mu, th$sigma)
  ),
  exp2 = list(
    npar   = 3,
    par0   = c(1, -1, -2),
    unpack = function(par) {
      a1 <- exp(par[2])
      list(p = plogis(par[1]), a1 = a1, a2 = a1 * plogis(par[3]))
    },
    surv   = function(d, th) th$p * exp(-th$a1 * d) + (1 - th$p) * exp(-th$a2 * d),
    dens   = function(d, th) th$p * th$a1 * exp(-th$a1 * d) +
      (1 - th$p) * th$a2 * exp(-th$a2 * d)
  ),
  exp3 = list(
    npar   = 5,
    par0   = c(0, 0, 0, -1, -1),
    unpack = function(par) {
      e  <- exp(c(0, par[1], par[2]))
      wt <- e / sum(e)
      a3 <- exp(par[3])
      a2 <- a3 * plogis(par[4])
      list(w1 = wt[1], w2 = wt[2], w3 = wt[3],
           a1 = a2 * plogis(par[5]), a2 = a2, a3 = a3)
    },
    surv   = function(d, th) th$w1 * exp(-th$a1 * d) + th$w2 * exp(-th$a2 * d) +
      th$w3 * exp(-th$a3 * d),
    dens   = function(d, th) th$w1 * th$a1 * exp(-th$a1 * d) +
      th$w2 * th$a2 * exp(-th$a2 * d) +
      th$w3 * th$a3 * exp(-th$a3 * d)
  )
)

FIT_PLAN <- list(
  veraset = c("lognormal", "exp2", "exp3"),
  meta    = c("lognormal", "exp2")
)

### FUNCTIONS ###

predict_bands <- function(fam, th, bands) {
  q <- fam$surv(bands$edge_low, th) - fam$surv(bands$edge_high, th)
  q / sum(q)
}

nll <- function(par, bands, fam) {
  q <- predict_bands(fam, fam$unpack(par), bands)
  if (any(!is.finite(q)) || any(q <= 0)) return(1e10)
  -bands$n_devices[1] * sum(bands$w * log(q))
}

# Closed-form lognormal from three shares with edges at 10 and 100 km.
# Works under either META_SCOPE: the identification only needs the first two
# cumulative shares and the fact that the edges are one decade apart.
lognormal_closed_form <- function(o1, o2) {
  z1    <- qnorm(o1)
  z2    <- qnorm(o1 + o2)
  sigma <- log(10) / (z2 - z1)
  list(mu = log(10) - sigma * z1, sigma = sigma)
}

surv_trunc <- function(fam, th, d, dmax = TRUNC_KM) {
  s_d <- fam$surv(d, th)
  s_m <- fam$surv(dmax, th)
  (s_d - s_m) / (1 - s_m)
}

functionals <- function(fam, th, dmax = TRUNC_KM) {
  med <- tryCatch(
    uniroot(function(d) surv_trunc(fam, th, d, dmax) - 0.5,
            interval = c(1e-6, dmax))$root,
    error = function(e) NA_real_
  )
  tibble(S1   = surv_trunc(fam, th, 1,   dmax),
         S10  = surv_trunc(fam, th, 10,  dmax),
         S100 = surv_trunc(fam, th, 100, dmax),
         median_km = med)
}

# Total variation: 0.5 * sum |observed - predicted| over the bands. Reads as
# the share of probability mass sitting in the wrong band.
tv_summary <- function(bands, q, cut_km = 10) {
  near <- bands$edge_high <= cut_km
  far  <- !near
  tibble(
    tv_all      = 0.5 * sum(abs(bands$w - q)),
    tv_near     = 0.5 * sum(abs(bands$w[near] - q[near])),
    tv_far      = 0.5 * sum(abs(bands$w[far]  - q[far])),
    tv_far_cond = 0.5 * sum(abs(bands$w[far] / sum(bands$w[far]) -
                                  q[far]        / sum(q[far])))
  )
}

score_theta <- function(fam, th, score_bands, convergence = 0L, nll_value = NA_real_) {
  q <- predict_bands(fam, th, score_bands)
  bind_cols(
    tibble(convergence = as.integer(convergence), nll_value = nll_value),
    as_tibble(th),
    functionals(fam, th),
    tv_summary(score_bands, q)
  )
}

fit_once <- function(fam, fit_bands, score_bands, par0) {
  o <- tryCatch(optim(par0, nll, bands = fit_bands, fam = fam, method = "BFGS"),
                error = function(e) NULL)
  if (is.null(o)) return(tibble(convergence = 99L, nll_value = NA_real_))
  score_theta(fam, fam$unpack(o$par), score_bands, o$convergence, o$value)
}

fit_starts <- function(fam, fit_bands, score_bands, n_start = N_START) {
  map_dfr(seq_len(n_start), function(i) {
    par0 <- if (i == 1) fam$par0 else fam$par0 + rnorm(fam$npar, sd = START_SD)
    fit_once(fam, fit_bands, score_bands, par0) |> mutate(start_id = i, .before = 1)
  })
}

### LOAD ###

emp <- read_csv(PROFILE_FILE, show_col_types = FALSE) |>
  mutate(home_county = str_pad(as.character(home_county), 5, pad = "0"))

meta_raw <- read_csv(META_FILE, show_col_types = FALSE)

xwalk <- read_csv(XWALK_FILE, show_col_types = FALSE) |>
  mutate(fips = str_pad(as.character(fips), 5, pad = "0")) |>
  distinct(gid_2, fips)

ambiguous <- xwalk |> count(fips) |> filter(n > 1) |> pull(fips)
message(length(ambiguous), " counties map to multiple GADM units - dropping")
xwalk_ok <- xwalk |> filter(!fips %in% ambiguous)

### VERASET BANDS ###

veraset_bands <- emp |>
  filter(edge_high <= TRUNC_KM) |>
  group_by(home_county) |>
  arrange(edge_low, .by_group = TRUE) |>
  mutate(w = share / sum(share)) |>
  ungroup() |>
  select(fips = home_county, band, edge_low, edge_high, w, n_devices)

### META BANDS ###
#
# Step 1: normalise the four categories within (gadm_id, ds), so each period
#         is a proper distribution regardless of how the raw file scales.
# Step 2: average across periods with equal weight. Equal weight rather than
#         volume weight because the periods are meant to be interchangeable
#         samples of the same fortnight, not a time series.
# Step 3: extract f0, then apply META_SCOPE.

META_ORDER <- c("(0, 10)", "[10, 100)", "100+")

meta_periods <- meta_raw |>
  inner_join(xwalk_ok, by = join_by(gadm_id == gid_2)) |>
  group_by(fips, ds) |>
  filter(n() == 4) |>
  mutate(frac = distance_category_ping_fraction /
           sum(distance_category_ping_fraction, na.rm = TRUE)) |>
  ungroup()

meta_shares <- meta_periods |>
  group_by(fips, home_to_ping_distance_category) |>
  summarise(frac = mean(frac, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = home_to_ping_distance_category, values_from = frac) |>
  rename(f0 = `0`, g1 = `(0, 10)`, g2 = `[10, 100)`, g3 = `100+`) |>
  drop_na(f0, g1, g2, g3)

meta_shares <- if (META_SCOPE == "folded") {
  # Home tile folded into (0, 10). Shares already sum to 1, no renormalising.
  meta_shares |> mutate(o1 = f0 + g1, o2 = g2, o3 = g3)
} else {
  # Conditional on having left home.
  meta_shares |> mutate(o1 = g1 / (1 - f0),
                        o2 = g2 / (1 - f0),
                        o3 = g3 / (1 - f0))
}

stopifnot(all(abs(meta_shares$o1 + meta_shares$o2 + meta_shares$o3 - 1) < 1e-8))

message("META_SCOPE = ", META_SCOPE,
        " | f0 median ",  signif(median(meta_shares$f0), 4),
        ", sd ",          signif(sd(meta_shares$f0), 3))
message("  o1 median ",   signif(median(meta_shares$o1), 4),
        ", sd ",          signif(sd(meta_shares$o1), 3))

meta_bands_all <- meta_shares |>
  select(fips, o1, o2, o3) |>
  pivot_longer(c(o1, o2, o3), names_to = "cat", values_to = "w") |>
  mutate(edge_low  = rep(c(0,  10,  100), length.out = n()),
         edge_high = rep(c(10, 100, Inf), length.out = n()),
         n_devices = META_OBS_WEIGHT) |>
  select(fips, edge_low, edge_high, w, n_devices)

### THE COUNTY SET ###

county_n <- veraset_bands |> distinct(fips, n_devices)

counties <- veraset_bands |>
  distinct(fips) |>
  semi_join(meta_bands_all, by = "fips") |>
  inner_join(county_n, by = "fips") |>
  filter(n_devices >= MIN_FIT) |>
  arrange(fips) |>
  pull(fips)

trust_fips <- county_n |> filter(n_devices >= MIN_TRUST) |> pull(fips)

message(length(counties), " counties with both sources and n_devices >= ", MIN_FIT)
message(sum(counties %in% trust_fips), " of them at or above ", MIN_TRUST)

### FIT ###

v_split <- veraset_bands |> filter(fips %in% counties) |> group_split(fips)
m_split <- meta_bands_all |> filter(fips %in% counties) |> group_split(fips)
names(v_split) <- map_chr(v_split, ~ .x$fips[1])
names(m_split) <- map_chr(m_split, ~ .x$fips[1])

f0_lookup <- meta_shares |> select(fips, f0) |> deframe()

fits <- map_dfr(seq_along(counties), function(i) {
  f  <- counties[i]
  vb <- v_split[[f]]
  mb <- m_split[[f]]
  if (i %% 100 == 0) message("  fitted ", i, " / ", length(counties))
  
  th_cf <- lognormal_closed_form(mb$w[1], mb$w[2])
  meta_ln <- score_theta(families$lognormal, th_cf, vb) |>
    mutate(start_id = 1L, source = "meta", family = "lognormal",
           method = "closed_form", .before = 1)
  
  bind_rows(
    map_dfr(FIT_PLAN$veraset, function(nm) {
      fit_starts(families[[nm]], vb, vb) |>
        mutate(source = "veraset", family = nm, method = "optim")
    }),
    meta_ln,
    map_dfr(setdiff(FIT_PLAN$meta, "lognormal"), function(nm) {
      fit_starts(families[[nm]], mb, vb) |>
        mutate(source = "meta", family = nm, method = "optim")
    })
  ) |>
    mutate(fips = f, n_devices = vb$n_devices[1],
           f0 = f0_lookup[[f]], meta_scope = META_SCOPE, .before = 1)
})

fits <- fits |> relocate(source, family, method, start_id, .after = n_devices)

### CLOSED-FORM VERIFICATION ###

verify <- meta_bands_all |>
  filter(fips %in% trust_fips) |>
  group_split(fips) |>
  map_dfr(function(mb) {
    cf <- lognormal_closed_form(mb$w[1], mb$w[2])
    o  <- optim(c(cf$mu, log(cf$sigma)), nll, bands = mb,
                fam = families$lognormal, method = "BFGS")
    th <- families$lognormal$unpack(o$par)
    tibble(fips = mb$fips[1],
           d_mu    = abs(th$mu    - cf$mu),
           d_sigma = abs(th$sigma - cf$sigma))
  })

message("closed-form verification on ", nrow(verify), " counties")
message("  |d_mu|    median ", signif(median(verify$d_mu), 3),
        ", p99 ",             signif(quantile(verify$d_mu, 0.99), 3))
message("  |d_sigma| median ", signif(median(verify$d_sigma), 3),
        ", p99 ",             signif(quantile(verify$d_sigma, 0.99), 3))

### WRITE ###

write_csv(fits, FIT_FILE)
message("wrote ", nrow(fits), " rows to ", FIT_FILE)

### QUICK LOOK ###

best <- fits |>
  filter(convergence == 0) |>
  group_by(fips, source, family) |>
  slice_min(nll_value, n = 1, with_ties = FALSE, na_rm = FALSE) |>
  ungroup()

# Veraset rows do not depend on META_SCOPE and should be identical across
# runs: lognormal ~0.042, exp2 ~0.081, exp3 ~0.017. If they move, something
# other than the Meta handling has changed.
best |>
  filter(n_devices >= MIN_TRUST) |>
  group_by(source, family) |>
  summarise(n = n(),
            tv_all      = median(tv_all,      na.rm = TRUE),
            tv_near     = median(tv_near,     na.rm = TRUE),
            tv_far      = median(tv_far,      na.rm = TRUE),
            tv_far_cond = median(tv_far_cond, na.rm = TRUE),
            .groups = "drop") |>
  print()

spread <- fits |>
  filter(convergence == 0, method == "optim") |>
  group_by(fips, source, family) |>
  summarise(sd_S1         = sd(S1, na.rm = TRUE),
            sd_log_median = sd(log(median_km), na.rm = TRUE),
            n_ok = n(),
            .groups = "drop")

spread |>
  filter(fips %in% trust_fips) |>
  group_by(source, family) |>
  summarise(med_sd_S1         = median(sd_S1, na.rm = TRUE),
            p90_sd_S1         = quantile(sd_S1, 0.9, na.rm = TRUE),
            med_sd_log_median = median(sd_log_median, na.rm = TRUE),
            .groups = "drop") |>
  print()