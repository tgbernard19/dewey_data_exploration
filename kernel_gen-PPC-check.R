###runs a posterior predictive check between meta fit and empirical Veraset Data 

library(tidyverse)
library(dplyr)

DATA_ROOT <- "E:/dewey-data"
TRACKB_DIR <- file.path(DATA_ROOT, "kernel", "trackB")

veraset_profile <- read_csv(file.path(TRACKB_DIR, "county_profile_cbghome_national_all_dwell_tau60.csv")) %>%
  mutate(w = share)

n_draws <- 1000

run_ppc_for_county <- function(par, cov, veraset_bands, fam, n_draws = 1000) {
  # --- Guard Checks ---
  if (is.null(par) || any(is.na(par)) || is.null(cov) || any(is.na(cov))) return(NULL)
  # CHANGED: Allow any valid non-empty band structure (7 or 9)
  if (is.null(veraset_bands) || nrow(veraset_bands) < 1) return(NULL)
  
  # DYNAMIC: Detect actual number of bands passed (7 for pooled, 9 for original)
  n_bands <- nrow(veraset_bands)
  
  # Coerce cov vector into 2x2 matrix if passed flat
  if (is.vector(cov) && length(cov) == 4) {
    cov <- matrix(cov, nrow = 2, ncol = 2)
  }
  if (!is.matrix(cov) || any(dim(cov) != c(2, 2))) return(NULL)
  
  # Invert Hessian to get true Covariance Matrix (Sigma = H^-1)
  cov_mat <- tryCatch(solve(cov), error = function(e) NULL)
  if (is.null(cov_mat) || any(is.na(cov_mat))) return(NULL)
  
  # Extract county device count
  N_veraset <- round(as.numeric(veraset_bands$n_devices[1]))
  if (is.na(N_veraset) || N_veraset <= 0) return(NULL)
  
  # --- Step 1: Posterior Parameter Sampling ---
  param_draws <- tryCatch(
    MASS::mvrnorm(n = n_draws, mu = par, Sigma = cov_mat),
    error = function(e) NULL
  )
  if (is.null(param_draws)) return(NULL)
  
  # --- Step 2: Calculate Continuous Predictions Across Dynamic Bands ---
  # CHANGED: Pre-allocate columns dynamically using n_bands
  predicted_shares <- matrix(NA_real_, nrow = n_draws, ncol = n_bands)
  
  for (s in 1:n_draws) {
    th_s <- fam$unpack(param_draws[s, ])
    
    q_s <- tryCatch({
      p_hi <- plnorm(veraset_bands$edge_high, meanlog = th_s$mu, sdlog = th_s$sigma)
      p_lo <- plnorm(veraset_bands$edge_low, meanlog = th_s$mu, sdlog = th_s$sigma)
      p_hi - p_lo
    }, error = function(e) NULL)
    
    if (!is.null(q_s) && !any(is.na(q_s)) && sum(q_s) > 0) {
      predicted_shares[s, ] <- q_s / sum(q_s)
    }
  }
  
  valid_rows <- complete.cases(predicted_shares)
  if (sum(valid_rows) == 0) return(NULL)
  valid_shares <- predicted_shares[valid_rows, , drop = FALSE]
  
  # --- Step 3: Multinomial Sampling ---
  simulated_counts <- t(apply(valid_shares, 1, function(q) {
    rmultinom(n = 1, size = N_veraset, prob = q)
  }))
  simulated_shares <- simulated_counts / N_veraset
  
  # --- Step 4: Summary Statistics ---
  tibble(
    band_idx  = veraset_bands$band_idx,
    edge_low  = veraset_bands$edge_low,
    edge_high = veraset_bands$edge_high,
    w_veraset = veraset_bands$w,
    q_median  = apply(simulated_shares, 2, median),
    q_low95   = apply(simulated_shares, 2, quantile, probs = 0.025),
    q_high95  = apply(simulated_shares, 2, quantile, probs = 0.975)
  )
}

library(furrr)

# Set up parallel execution across available cores (e.g., 8 workers or availableCores() - 1)
plan(multisession, workers = 124)

ppc_results <- meta_lognorm %>%
  filter(converged) %>%
  inner_join(veraset_nested, by = c("fips" = "home_county")) %>%
  mutate(
    ppc = future_pmap(
      list(par = par, cov = cov, veraset_bands = veraset_bands),
      ~ run_ppc_for_county(..1, ..2, ..3, fam = families$lognormal),
      .options = furrr_options(seed = TRUE)
    )
  ) %>%
  dplyr::select(fips, ppc) %>%
  filter(!map_lgl(ppc, is.null)) %>% # Drop any counties that failed guards
  unnest(ppc)

# Reset workers back to standard single-threaded plan when done
plan(sequential)

# 1. Add coverage and residual flags
ppc_eval <- ppc_results %>%
  mutate(
    covered = w_veraset >= q_low95 & w_veraset <= q_high95,
    residual = q_median - w_veraset, # positive = model over-predicts share
    abs_err = abs(residual)
  )

# 2. National aggregate coverage (Target: ~95%)
overall_coverage <- mean(ppc_eval$covered)
message(sprintf("National 95%% PPC Coverage Rate: %.2f%%", overall_coverage * 100))

# 3. Band-specific breakdown to locate failure modes
band_summary <- ppc_eval %>%
  group_by(band_idx, edge_low, edge_high) %>%
  summarise(
    coverage = mean(covered),
    mean_bias = mean(residual),  # Systematic shift
    mae = mean(abs_err),         # Mean absolute error
    .groups = "drop"
  )

print(band_summary)

#0-5 km behavior doesn't matter much at this resolution, and Veraset struggles to recapitulate it because intra-CBG dists are often less than 5km
#producing bands that pool over 0-5 instead


# Collapse sub-5 km bands into a single 0-5 km band
veraset_nested_pooled <- veraset_profile %>%
  ungroup() %>%
  filter(edge_high <= 500) %>%
  mutate(
    # Re-map bands 1, 2, 3 into a single 0-5 km band
    band_idx_pooled = case_when(
      edge_high <= 5.0 ~ 1,
      TRUE ~ band_idx - 2 # Shift indices 4..9 down to 2..7
    ))#,
    edge_low_pooled = case_when(
      edge_high <= 5.0 ~ 0.0,
      TRUE ~ edge_low
    ),
    edge_high_pooled = case_when(
      edge_high <= 5.0 ~ 5.0,
      TRUE ~ edge_high
    )
  ) %>%
  group_by(home_county, band_idx_pooled, edge_low_pooled, edge_high_pooled) %>%
  summarise(
    w = sum(w), # Sum Veraset shares for 0-1, 1-2.5, 2.5-5 km
    n_devices = first(n_devices),
    .groups = "drop"
  ) %>%
  rename(
    band_idx = band_idx_pooled,
    edge_low = edge_low_pooled,
    edge_high = edge_high_pooled
  ) %>%
  dplyr::select(home_county, band_idx, edge_low, edge_high, w, n_devices) %>%
  nest(veraset_bands = -home_county)

# Sanity Check: Ensure every county now has EXACTLY 7 bands
stopifnot(nrow(veraset_nested_pooled$veraset_bands[[1]]) == 7)

plan(multisession, workers = 123)

ppc_results_pooled <- meta_lognorm %>%
  filter(converged) %>%
  inner_join(veraset_nested_pooled, by = c("fips" = "home_county")) %>%
  mutate(
    ppc = future_pmap(
      list(par = par, cov = cov, veraset_bands = veraset_bands),
      ~ run_ppc_for_county(..1, ..2, ..3, fam = families$lognormal),
      .options = furrr_options(seed = TRUE)
    )
  ) %>%
  dplyr::select(fips, ppc) %>%
  filter(!map_lgl(ppc, is.null)) %>%
  unnest(ppc)

plan(sequential)

ppc_eval_2 <- ppc_results_pooled %>%
  mutate(
    covered = w_veraset >= q_low95 & w_veraset <= q_high95,
    residual = q_median - w_veraset, # positive = model over-predicts share
    abs_err = abs(residual)
  )

# 2. National aggregate coverage (Target: ~95%)
overall_coverage <- mean(ppc_eval_2$covered)
message(sprintf("National 95%% PPC Coverage Rate: %.2f%%", overall_coverage * 100))

rucc <- read_csv(file.path(DATA_ROOT, "external", "rural_continuum",
                           "Ruralurbancontinuumcodes2023.csv")) %>%
  mutate(fips = FIPS)

rucc

ppc_eval_2 %>%
  filter(covered == FALSE) %>%
  count(band_idx)

failed_counties <- ppc_eval_2 %>%
  filter(covered == FALSE) %>%
  distinct(fips) %>%
  left_join(rucc, by = "fips") %>%
  count(Value)# %>%
#  filter(Value < 10)

View(failed_counties)

# 3. Band-specific breakdown to locate failure modes
band_summary <- ppc_eval_2 %>%
  group_by(band_idx, edge_low, edge_high) %>%
  summarise(
    coverage = mean(covered),
    mean_bias = mean(residual),  # Systematic shift
    mae = mean(abs_err),         # Mean absolute error
    .groups = "drop"
  )

print(band_summary)



### plotting
library(tidyverse)

# --- Step 1: JS Divergence Helper & Fit Summary ---
calc_jsd <- function(p, q, eps = 1e-12) {
  p <- p + eps; p <- p / sum(p)
  q <- q + eps; q <- q / sum(q)
  m <- 0.5 * (p + q)
  kl_pm <- sum(p * log(p / m))
  kl_qm <- sum(q * log(q / m))
  0.5 * (kl_pm + kl_qm)
}

# Measure county-level fit quality across 7 bands
county_jsd <- ppc_results_pooled %>%
  group_by(fips) %>%
  summarise(
    jsd = calc_jsd(w_veraset, q_median),
    n_devices = first(q_median), # surrogate check
    .groups = "drop"
  ) %>%
  arrange(jsd)

# --- Step 2: Sample 9 Representative Counties Across the Spectrum ---
n_total <- nrow(county_jsd)

# Select 3 Top Fit (10th-30th percentile), 3 Median Fit (45th-55th), 3 Lower-Bound (90th-95th)
sampled_indices <- c(
  round(n_total * 0.05), round(n_total * 0.15), round(n_total * 0.25), # Excellent
  round(n_total * 0.45), round(n_total * 0.50), round(n_total * 0.55), # Typical
  round(n_total * 0.88), round(n_total * 0.92), round(n_total * 0.96)  # Outer Bound
)

selected_counties <- county_jsd[sampled_indices, ] %>%
  mutate(
    fit_category = case_when(
      row_number() <= 3 ~ "Top 25% Fit (Low JSD)",
      row_number() <= 6 ~ "Median Fit (Typical)",
      TRUE ~ "Lower 10% Fit (Outer Bound)"
    ),
    facet_label = sprintf("%s (FIPS: %s)\nJSD: %.4f | %s", fips, fips, jsd, fit_category)
  )

# --- Step 3: Prepare Data for Ribbon Plotting ---
band_labels <- c("0–5", "5–10", "10–25", "25–50", "50–100", "100–250", "250–500")

plot_data <- ppc_results_pooled %>%
  inner_join(selected_counties, by = "fips") %>%
  mutate(
    band_name = factor(band_labels[band_idx], levels = band_labels),
    # Order facets from best fit to worst fit
    facet_label = factor(facet_label, levels = unique(selected_counties$facet_label))
  )

# --- Step 4: Generate Faceted Ribbon Figure ---
ggplot(plot_data, aes(x = band_name, group = 1)) +
  # 1. Meta 95% Posterior Predictive Interval Ribbon
  geom_ribbon(
    aes(ymin = q_low95, ymax = q_high95, fill = "Meta 95% PPC Interval"), 
    alpha = 0.35
  ) +
  # 2. Meta Predicted Median Line
  geom_line(
    aes(y = q_median, color = "Meta Predicted Median"), 
    linewidth = 0.9
  ) +
  # 3. Veraset Ground Truth Points & Line
  geom_line(
    aes(y = w_veraset, color = "Veraset Ground Truth"), 
    linewidth = 0.6, linetype = "dashed"
  ) +
  geom_point(
    aes(y = w_veraset, color = "Veraset Ground Truth"), 
    size = 2.2
  ) +
  # 4. Faceting & Styling
  facet_wrap(~ facet_label, ncol = 3, scales = "free_y") +
  scale_fill_manual(name = "", values = c("Meta 95% PPC Interval" = "#1877F2")) +
  scale_color_manual(
    name = "", 
    values = c("Meta Predicted Median" = "#0B4EA2", "Veraset Ground Truth" = "#D9381E")
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "National Posterior Predictive Check: Observed vs. Predicted Trip Shares",
    subtitle = "7-Band Pooled Model (0–5 km to 250–500 km) Across Representative County Fit Profiles",
    x = "Distance Band (km)",
    y = "Trip Share (% of Total County Trips)",
    caption = "PPC draws: 1,000 parameter samples per county from N(θ_hat, H^-1) x Multinomial observation noise."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#555555", size = 10, margin = margin(b = 12)),
    strip.background = element_rect(fill = "#F0F2F5", color = NA),
    strip.text = element_text(face = "bold", size = 8.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "top",
    legend.margin = margin(b = -5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "#E0E0E0", fill = NA, linewidth = 0.5)
  )


### saving output

# Unpack parameter vector & Hessian elements into scalar columns
meta_export <- meta_lognorm %>%
  drop_na() %>%
  mutate(
    fips = str_pad(as.character(fips), 5, pad = "0"),
    # Extract lognormal parameters: c(mu, log_sigma)
    mu = map_dbl(par, 1),
    log_sigma = map_dbl(par, 2),
    # Flatten Hessian matrix elements (H_11, H_12, H_21, H_22)
    h_11 = map_dbl(cov, 1),
    h_12 = map_dbl(cov, 2),
    h_21 = map_dbl(cov, 3),
    h_22 = map_dbl(cov, 4)
  ) %>%
  dplyr::select(fips, mu, log_sigma, h_11, h_12, h_21, h_22, converged)


# Ensure output directory exists and write CSV
# This fitted-kernel export (derived from the raw Meta movement-distribution CSV
# under external/meta_movement_dist) lives under trackB/ alongside the other
# derived kernel artifacts; absolute_flows_compare.R reads it from here.
dir.create(TRACKB_DIR, recursive = TRUE, showWarnings = FALSE)
write_csv(meta_export, file.path(TRACKB_DIR, "meta_lognormal_kernel-fit.csv"))

message("Successfully saved ", nrow(meta_export), " county fits to meta_lognormal_kernel-fit.csv")
