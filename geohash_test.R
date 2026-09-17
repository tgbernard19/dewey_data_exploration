#!/usr/bin/env Rscript

# =============================================================================
# TRACK B — THE GEOHASH FLOOR: SCOPE OF DAMAGE AND CLEAN REBUILD
# =============================================================================
# CONFIRMED by the audit:
#
#   is_home rows, median displacement      CBG 0.000 km    geohash 3.467 km
#   fraction within 100 m of zero          CBG 0.99999998  geohash 0.00047
#   share of home rows that are geohash    22.9%
#
# Home displacement is zero by construction, so that entire distribution is
# location error at the geohash-5 cell scale.
#
# Effect on the profile, conditioning both groups on being away from band 1 so
# the comparison is not a renormalisation artifact:
#
#   band        CBG-home   geohash-home
#   [1,2.5)       0.205       0.202
#   [2.5,5)       0.172       0.269      <- excess at the cell scale
#   [5,10)        0.164       0.194      <-
#   [10,25)       0.177       0.126      <- where it came from
#   [1000,Inf)    0.032       0.030         far tail untouched
#
# So the national profile is a MIXTURE of two populations: ~77% of devices with
# mass spiked at zero, ~23% with that spike smeared across 1-10 km. That is why
# six unrelated families all failed in the same range, and why every one of them
# under-predicted [2.5,5) by about +0.047 while over-predicting [10,25).
#
# -----------------------------------------------------------------------------
# WHY THIS SCRIPT EXISTS
#
# The family question is now secondary. The urgent question is whether the
# Track B headline survives: if geohash resolution rates vary with rurality —
# and sparser CBG matching in rural areas is exactly what one would expect —
# then Veraset's RUCC gradient may be a resolution artifact rather than
# behaviour. That is finding 2, and it is what goes to Charlie.
#
# PHASES
#   1  Device-level home resolution table. One pass over is_home rows only.
#   2  Stratified county profiles, by JOINING to the cached profile parts.
#      No rebuild — did and band_idx are already there.
#   3  Does geohash share explain the cross-county gradients? The test that
#      decides whether finding 2 stands.
#   4  Refit M1 / M4 / M7 on CBG-home devices only. If the mixture was the
#      whole story, tv should collapse.
#
# SCOPE NOTE. This stratifies on HOME resolution only. A geohash-resolved
# destination also carries ~3.5 km of error, but it affects one row, whereas
# home error displaces every row belonging to that device, since d_km is
# measured from home. Home resolution is therefore the dominant term. Doing
# destination resolution as well would require carrying loc_src through the
# day aggregation, which is a rescan — worth it only if Phase 4 leaves
# meaningful misfit behind.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DBI)
  library(duckdb)
  library(cmdstanr)
  library(posterior)
  library(furrr)
  library(arrow)
})

# ---- configuration ----------------------------------------------------------

DATA_ROOT  <- "E:/dewey-data"
KERNEL_DIR <- file.path(DATA_ROOT, "kernel")
DEVDAY_DIR <- file.path(KERNEL_DIR, "devday_parts")
OUT_DIR    <- file.path(KERNEL_DIR, "trackB")
TMP_DIR    <- file.path(DATA_ROOT, "tmp")

SCOPE <- "all"
TAU   <- 60
RUN_TAG  <- sprintf("national_%s_dwell_tau%d", SCOPE, TAU)
PART_DIR <- file.path(OUT_DIR, sprintf("profile_parts_%s", RUN_TAG))

CBG_LABEL <- "cbg"          # confirmed from the loc_src vocabulary
GEO_LABEL <- "geohash"

TRUNC_KM    <- 500
N_SAMPLE    <- 300
MIN_DEVICES <- 500
WORKERS     <- 16
SEED        <- 20260901

N_THREADS <- 32
MEM_LIMIT <- "400GB"

A1_RANGE <- c(0.02, 2.0)
A2_RANGE <- c(1e-4, 0.5)
BAND_EDGES <- c(0, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000, Inf)

HOMERES_FILE <- file.path(OUT_DIR, "device_home_resolution.parquet")
AREA_FILE    <- file.path(OUT_DIR, "county_area.csv")
FITS_FILE    <- file.path(OUT_DIR, sprintf("veraset_band_fits_%s.csv", SCOPE))

# Optional: fips, rucc. The 40-county benchmark file has it; a national RUCC
# file is better if one is to hand. Phase 3 degrades gracefully without it.
# Same file, same path, as RUCC_FILE in kernel_gen-PPC-check.R -- keep them in
# sync if this ever moves.
RUCC_FILE <- file.path(DATA_ROOT, "external", "rural_continuum",
                       "Ruralurbancontinuumcodes2023.csv")

say <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n",
                         sep = "")

stopifnot(dir.exists(PART_DIR))

con <- dbConnect(duckdb(), dbdir = ":memory:")
dbExecute(con, sprintf("SET threads = %d", N_THREADS))
dbExecute(con, sprintf("SET memory_limit = '%s'", MEM_LIMIT))
dbExecute(con, sprintf("SET temp_directory = '%s'", TMP_DIR))
dbExecute(con, "SET preserve_insertion_order = false")


# =============================================================================
# PHASE 1 — DEVICE-LEVEL HOME RESOLUTION
# =============================================================================
# Home resolution is a property of the DEVICE, not of a device-day: it comes
# from home_assignment upstream. Computing it per day would mark a device
# geohash-resolved simply because it had no is_home row that day.
#
# One pass, is_home rows only, so this is far cheaper than the main build.
# A device is treated as CBG-resolved if ANY of its home rows resolved to a
# CBG — the conservative direction, since a CBG match is the reliable one.
# =============================================================================

say("PHASE 1 — device home resolution")

if (!file.exists(HOMERES_FILE)) {
  day_files <- sort(list.files(DEVDAY_DIR, pattern = "\\.parquet$",
                               full.names = TRUE))
  glob_all <- paste(sprintf("'%s'", day_files), collapse = ", ")
  
  tmp <- paste0(HOMERES_FILE, ".tmp")
  
  dbExecute(con, sprintf("
  COPY (
    SELECT
      did,
      MAX(CASE WHEN loc_src = '%s' THEN 1 ELSE 0 END) AS cbg_home,
      SUM(CASE WHEN loc_src = '%s' THEN 1 ELSE 0 END) AS n_home_cbg,
      SUM(CASE WHEN loc_src = '%s' THEN 1 ELSE 0 END) AS n_home_geo
    FROM read_parquet([%s])
    WHERE is_home
    GROUP BY did
  ) TO '%s' (FORMAT PARQUET, COMPRESSION ZSTD)
  ", CBG_LABEL, CBG_LABEL, GEO_LABEL, glob_all, tmp))
  
  file.rename(tmp, HOMERES_FILE)
  say("  wrote ", HOMERES_FILE)
} else {
  say("  cached")
}

res_summary <- dbGetQuery(con, sprintf(
  "SELECT cbg_home, COUNT(*) AS n_devices FROM read_parquet('%s')
   GROUP BY cbg_home ORDER BY cbg_home", HOMERES_FILE))
print(res_summary)


# =============================================================================
# PHASE 2 — STRATIFIED COUNTY PROFILES
# =============================================================================
# The cached profile parts already hold did, home_county, band_idx and
# share_sum for the whole month. Stratifying is a join, not a rebuild.
# =============================================================================

say("PHASE 2 — stratified profiles")

dbExecute(con, sprintf("
CREATE OR REPLACE TABLE dev_band AS
SELECT p.did, p.band_idx, SUM(p.share_sum) AS share_sum
FROM read_parquet('%s/*.parquet') p
GROUP BY p.did, p.band_idx
", PART_DIR))

dbExecute(con, sprintf("
CREATE OR REPLACE TABLE dev_home AS
SELECT p.did,
       ANY_VALUE(p.home_county) AS home_county,
       ANY_VALUE(COALESCE(r.cbg_home, 0)) AS cbg_home
FROM read_parquet('%s/*.parquet') p
LEFT JOIN read_parquet('%s') r USING (did)
GROUP BY p.did
", PART_DIR, HOMERES_FILE))

# Shares sum to 1 per device-day, so the total across bands and days is the
# device's device-day count. Same trick as the main build, no COUNT(DISTINCT).
dbExecute(con, "
CREATE OR REPLACE TABLE dev_denom AS
SELECT did, SUM(share_sum) AS n_devdays FROM dev_band GROUP BY did
")

county_strat <- dbGetQuery(con, "
WITH dev_share AS (
  SELECT h.home_county, h.cbg_home, b.band_idx,
         b.share_sum / d.n_devdays AS dev_share
  FROM dev_band b
  JOIN dev_denom d USING (did)
  JOIN dev_home  h USING (did)
),
cty_n AS (
  SELECT home_county, cbg_home, COUNT(*) AS n_devices
  FROM dev_home GROUP BY home_county, cbg_home
)
SELECT s.home_county, s.cbg_home, s.band_idx,
       SUM(s.dev_share) / ANY_VALUE(n.n_devices) AS share,
       ANY_VALUE(n.n_devices) AS n_devices
FROM dev_share s
JOIN cty_n n ON n.home_county = s.home_county AND n.cbg_home = s.cbg_home
GROUP BY s.home_county, s.cbg_home, s.band_idx
")

band_key <- tibble(
  band_idx = 1:11, edge_low = BAND_EDGES[1:11], edge_high = BAND_EDGES[2:12]
) |>
  mutate(band = sprintf("[%g,%s)", edge_low,
                        ifelse(is.infinite(edge_high), "Inf",
                               format(edge_high))))

strat <- county_strat |>
  mutate(home_county = str_pad(as.character(home_county), 5, pad = "0")) |>
  complete(nesting(home_county, cbg_home, n_devices), band_idx = 1:11,
           fill = list(share = 0)) |>
  left_join(band_key, by = "band_idx") |>
  arrange(home_county, cbg_home, band_idx)

write_csv(strat, file.path(OUT_DIR, sprintf("county_profile_bystrat_%s.csv",
                                            RUN_TAG)))

# The clean profile: CBG-resolved homes only.
profile_clean <- strat |> filter(cbg_home == 1)
write_csv(profile_clean,
          file.path(OUT_DIR, sprintf("county_profile_cbghome_%s.csv", RUN_TAG)))

say(sprintf("  %d counties in the clean profile; median n_devices %d",
            n_distinct(profile_clean$home_county),
            as.integer(median(unique(profile_clean$n_devices)))))


# =============================================================================
# PHASE 3 — DOES GEOHASH SHARE EXPLAIN THE CROSS-COUNTY GRADIENTS?
# =============================================================================
# The test that decides whether finding 2 survives.
#
# If geohash_share correlates strongly with rurality AND absorbs the RUCC
# coefficient on lo_v, then Veraset's rural gradient is a data-resolution
# artifact, not behaviour, and the Charlie framing needs rewriting.
#
# Note the direction to expect. Geohash homes move mass from [0,1) into
# 1-10 km, which is BELOW 10, so it depresses P(d > 10 | d < 100) and hence
# lo_v. The reported Veraset RUCC coefficient was POSITIVE (+0.063), so a
# naive resolution story predicts the wrong sign — which is a point in
# finding 2's favour and worth checking explicitly rather than assuming.
# =============================================================================

say("PHASE 3 — does resolution explain the gradient?")

geo_share <- strat |>
  distinct(home_county, cbg_home, n_devices) |>
  pivot_wider(names_from = cbg_home, values_from = n_devices,
              names_prefix = "n_", values_fill = 0) |>
  mutate(n_total = n_0 + n_1,
         geohash_share = n_0 / pmax(n_total, 1))

say(sprintf("  geohash_share across counties: median %.3f, IQR %.3f-%.3f",
            median(geo_share$geohash_share),
            quantile(geo_share$geohash_share, 0.25),
            quantile(geo_share$geohash_share, 0.75)))

cf <- read_csv(FITS_FILE, show_col_types = FALSE) |>
  mutate(home_county = str_pad(as.character(home_county), 5, pad = "0")) |>
  filter(status == "ok", max_rhat <= 1.01) |>
  select(home_county, lo_v, lo_v_sd, p_v, n_devices) |>
  inner_join(select(geo_share, home_county, geohash_share, n_total),
             by = "home_county", relationship = "one-to-one")

if (file.exists(AREA_FILE)) {
  cf <- cf |>
    left_join(read_csv(AREA_FILE, show_col_types = FALSE) |>
                mutate(fips = str_pad(as.character(fips), 5, pad = "0")) |>
                select(fips, area_km2),
              by = join_by(home_county == fips), relationship = "one-to-one")
}

if (file.exists(RUCC_FILE)) {
  cf <- cf |>
    left_join(read_csv(RUCC_FILE, show_col_types = FALSE) |>
                filter(Attribute == "RUCC_2023") |>
                mutate(fips = str_pad(as.character(FIPS), 5, pad = "0"),
                       rucc = Value) |>
                select(fips, rucc),
              by = join_by(home_county == fips), relationship = "one-to-one")
}

cf <- cf %>%
  select(-rucc.x) %>%
  mutate(rucc = rucc.y) %>%
  select(-rucc.y)

cf <- filter(cf, n_total >= MIN_DEVICES)

if ("rucc" %in% names(cf) && sum(!is.na(cf$rucc)) > 100) {
  say("  geohash_share against rurality:")
  print(summary(lm(geohash_share ~ rucc, data = cf)))
  
  say("  lo_v, with and without geohash_share:")
  m_base <- lm(lo_v ~ rucc + log(area_km2), data = cf)
  m_full <- lm(lo_v ~ rucc + log(area_km2) + geohash_share, data = cf)
  print(summary(m_base))
  print(summary(m_full))
  
  say(sprintf("  RUCC coefficient: %.4f -> %.4f when resolution enters.",
              coef(m_base)["rucc"], coef(m_full)["rucc"]))
  say("  A large attenuation means the rural gradient is a resolution")
  say("  artifact and finding 2 needs rewriting. Little change means it")
  say("  stands, and this becomes a robustness paragraph instead.")
} else {
  say("  No national RUCC file at ", RUCC_FILE, " — running area only.")
  if ("area_km2" %in% names(cf)) {
    print(summary(lm(lo_v ~ log(area_km2) + geohash_share, data = cf)))
  }
  say("  Supply fips,rucc to complete the decisive test.")
}

write_csv(cf, file.path(OUT_DIR, "diag_resolution_covariates.csv"))


# =============================================================================
# PHASE 4 — REFIT THE THREE FAMILIES ON CBG-HOME DEVICES ONLY
# =============================================================================
# If the two-population mixture was the whole story, tv should fall sharply and
# the [2.5,5) residual should vanish. Whatever remains after this is either
# destination-side quantisation or genuine kernel structure.
# =============================================================================

say("PHASE 4 — refit on clean profile")

bands_clean <- profile_clean |>
  filter(edge_high <= TRUNC_KM) |>
  group_by(home_county) |>
  mutate(w = share / sum(share)) |>
  ungroup() |>
  select(home_county, band_idx, edge_low, edge_high, w, n_devices) |>
  filter(n_devices >= MIN_DEVICES)

set.seed(SEED)
samp <- bands_clean |>
  distinct(home_county, n_devices) |>
  mutate(decile = ntile(n_devices, 10)) |>
  group_by(decile) |>
  slice_sample(n = ceiling(N_SAMPLE / 10)) |>
  ungroup() |>
  pull(home_county)

fit_input <- filter(bands_clean, home_county %in% samp)

data_block <- "
data {
  int<lower=3> K;
  vector<lower=0>[K] edge_lo;
  vector<lower=0>[K] edge_hi;
  vector<lower=0>[K] share_raw;
  real<lower=0> n_pseudo;
  real a1_lo; real a1_hi; real a2_lo; real a2_hi;
}
transformed data { vector[K] w = share_raw / sum(share_raw); }
"

metrics_block <- "
generated quantities {
  vector[K] q_out = q;
  real tv_all = 0.5 * sum(abs(w - q));
}
"

code_M1 <- paste0(data_block, "
parameters {
  real<lower=0, upper=1> p;
  real<lower=a1_lo, upper=a1_hi> a1;
  real<lower=a2_lo, upper=fmin(a1, a2_hi)> a2;
}
transformed parameters {
  vector[K] q;
  { vector[K] num;
    for (k in 1:K)
      num[k] = p * (exp(-a1*edge_lo[k]) - exp(-a1*edge_hi[k]))
             + (1-p) * (exp(-a2*edge_lo[k]) - exp(-a2*edge_hi[k]));
    q = num / sum(num); }
}
model {
  p ~ beta(1, 1);
  target += -log(a1) - log(a2);
  target += n_pseudo * dot_product(w, log(q));
}
", metrics_block)

code_M4 <- paste0(data_block, "
parameters { real mu; real<lower=0.05, upper=10> sigma; }
transformed parameters {
  vector[K] q;
  { vector[K] num;
    for (k in 1:K) {
      real lo = edge_lo[k] > 0 ? lognormal_cdf(edge_lo[k] | mu, sigma) : 0;
      num[k] = lognormal_cdf(edge_hi[k] | mu, sigma) - lo;
    }
    q = num / sum(num); }
}
model {
  mu ~ normal(0, 5); sigma ~ normal(0, 3);
  target += n_pseudo * dot_product(w, log(q));
}
", metrics_block)

code_M7 <- paste0(data_block, "
parameters {
  real<lower=0, upper=1> f;
  real<lower=0.01, upper=50> d0;
  real<lower=0.05, upper=5> beta_t;
}
transformed parameters {
  vector[K] q;
  { vector[K] num;
    for (k in 1:K)
      num[k] = (1 - f) * (pow(1 + edge_lo[k]/d0, -beta_t)
                        - pow(1 + edge_hi[k]/d0, -beta_t));
    num[1] += f;
    q = num / sum(num); }
}
model {
  f ~ beta(1, 1); d0 ~ lognormal(0, 2); beta_t ~ lognormal(0, 1);
  target += n_pseudo * dot_product(w, log(q));
}
", metrics_block)

EXE <- c(
  M1_exp2       = cmdstan_model(write_stan_file(code_M1))$exe_file(),
  M4_lognormal  = cmdstan_model(write_stan_file(code_M4))$exe_file(),
  M7_atom_lomax = cmdstan_model(write_stan_file(code_M7))$exe_file()
)

fit_clean <- function(cb) {
  cb <- arrange(cb, band_idx)
  fips <- cb$home_county[1]; n_dev <- cb$n_devices[1]
  
  dat <- list(K = nrow(cb), edge_lo = cb$edge_low, edge_hi = cb$edge_high,
              share_raw = cb$w, n_pseudo = n_dev,
              a1_lo = A1_RANGE[1], a1_hi = A1_RANGE[2],
              a2_lo = A2_RANGE[1], a2_hi = A2_RANGE[2])
  
  run_one <- function(label, exe) {
    r <- try({
      mod <- cmdstan_model(exe_file = exe)
      f <- mod$sample(data = dat, seed = SEED, chains = 4, parallel_chains = 1,
                      iter_warmup = 1500, iter_sampling = 1000,
                      adapt_delta = 0.99, max_treedepth = 12, refresh = 0,
                      show_messages = FALSE, show_exceptions = FALSE)
      tibble(home_county = fips, n_devices = n_dev, model = label,
             band_idx = cb$band_idx, w = cb$w,
             q = f$summary("q_out", mean)$mean,
             tv_all = mean(f$draws("tv_all", format = "df")$tv_all),
             max_rhat = max(f$summary()$rhat, na.rm = TRUE))
    }, silent = TRUE)
    if (inherits(r, "try-error")) NULL else r
  }
  
  imap_dfr(EXE, \(exe, label) run_one(label, exe))
}

plan(multisession, workers = WORKERS)
t0 <- Sys.time()

clean_fits <- fit_input |>
  group_split(home_county) |>
  future_map_dfr(fit_clean,
                 .options = furrr_options(seed = TRUE, packages = "cmdstanr"),
                 .progress = TRUE)

plan(sequential)
say(sprintf("  %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

verdict <- clean_fits |>
  distinct(home_county, model, tv_all, max_rhat) |>
  group_by(model) |>
  summarise(n = n(), median_tv = median(tv_all),
            p90_tv = quantile(tv_all, 0.9),
            pct_bad = 100 * mean(max_rhat > 1.01), .groups = "drop") |>
  arrange(median_tv)

say("  CLEAN profile, CBG-resolved homes only:")
print(verdict, width = Inf)

say("  Contaminated baseline was M7 0.100, M4 0.117, M1 0.175.")

resid_clean <- clean_fits |>
  mutate(r = w - q) |>
  group_by(model, band_idx) |>
  summarise(mean_resid = mean(r), .groups = "drop") |>
  left_join(select(band_key, band_idx, band), by = "band_idx")

say("  Residual by band — [2.5,5) was +0.047 on the contaminated profile:")
print(resid_clean, n = Inf, width = Inf)

write_csv(clean_fits, file.path(OUT_DIR, "diag_clean_fits_raw.csv"))
write_csv(verdict,    file.path(OUT_DIR, "diag_clean_verdict.csv"))

say("")
say("READING PHASE 4")
say("  tv collapses and the [2.5,5) residual disappears -> the mixture was")
say("  the whole story. Rebuild the prior on the clean profile and the")
say("  incumbent two-exponential family may well be adequate after all.")
say("  tv falls but the residual persists -> destination-side quantisation")
say("  is next, and that one needs loc_src carried through a rescan.")
say("  tv barely moves -> genuine kernel structure, and the misfit becomes a")
say("  reported limitation rather than something to fix.")

dbDisconnect(con, shutdown = TRUE)