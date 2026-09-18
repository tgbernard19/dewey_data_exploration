#!/usr/bin/env Rscript
# ===========================================================================
# 11_fit_meta.R
# ---------------------------------------------------------------------------
# meta_county_shares.csv -> one lognormal kernel per county.
#
# THE CLOSED FORM IS THE ESTIMATOR
#   Meta's three categories give two independent constraints, F(10) = o1 and
#   F(100) = o1 + o2. A lognormal has two parameters. So it is exactly
#   identified and there is nothing to optimise:
#
#     z1 = qnorm(o1),  z2 = qnorm(o1 + o2)
#     sigma = log(10) / (z2 - z1)
#     mu    = log(10) - sigma * z1
#
#   optim runs anyway, on every county, as a check that the two agree. A
#   disagreement means one of them has the band edges wrong -- which is a
#   coding error, not a data finding, and it should be impossible to ship.
#
#   This is the whole reason the lognormal replaced the two-exponential
#   mixture. Three parameters against two constraints left one direction
#   unidentified, and that ridge was being pinned by an arbitrary beta(9,1)
#   prior on p. The ridge, the prior, and the Gibbs-posterior machinery all
#   existed only because the family was too big for the data.
#
# WHAT THE COVARIANCE MEANS
#   Very little in absolute terms. Meta ships fractions with no sample size,
#   so the Hessian is taken at a made-up n (META_OBS_WEIGHT) and the
#   covariance scales as 1/n. Point estimates do not depend on it. Keep the
#   columns for relative comparisons and for posterior predictive work; do
#   not read them as confidence intervals.
#
# Output: meta_lognormal_kernel-fit.csv -- the file the flow stage reads.
#
# Usage:  Rscript R/11_fit_meta.R
# ===========================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

source("config.R")
source("R/01_utils.R")
source("R/04_lognormal.R")
source("R/05_meta.R")

# The existing file is the parity baseline. FALSE writes alongside it and
# reports the differences; set TRUE only once those differences are
# understood.
OVERWRITE <- FALSE

banner("META KERNEL FITS")


# ---- input -----------------------------------------------------------------

shares <- read_keyed_csv(META_SHARES_FILE, key_cols = "fips")
say("counties with shares: ", nrow(shares))

stopifnot("shares file is not folded -- rerun 10_meta_shares.R" =
            all(shares$meta_scope == META_SCOPE))


# ---- 1. closed form, vectorised over every county --------------------------

cf <- shares |>
  mutate(ln_closed_form(o1, o2)) |>
  mutate(log_sigma = log(sigma))

say("closed form: ", sum(is.finite(cf$mu) & is.finite(cf$sigma)),
    " of ", nrow(cf), " counties give finite parameters")

# Non-finite happens when o1 is 0 or 1, i.e. qnorm hits an infinity. Those
# are degenerate inputs rather than failed fits, and they are dropped with
# their FIPS named rather than silently.
degenerate <- cf |> filter(!is.finite(mu) | !is.finite(sigma) | sigma <= 0)
if (nrow(degenerate) > 0) {
  say("[!] dropping ", nrow(degenerate), " degenerate counties: ",
      paste(head(degenerate$fips, 10), collapse = ", "))
}
cf <- cf |> filter(is.finite(mu), is.finite(sigma), sigma > 0)


# ---- 2. optim with the Hessian, per county ---------------------------------
# The parameters come from step 1. This adds the covariance and, more
# importantly, checks that an independent route reaches the same answer.

say("fitting ", nrow(cf), " counties with optim (for covariance and check) ...")

hess <- map(cf$fips, function(f) {
  b <- meta_bands_for(shares, f)
  fit_lognormal_hessian(b) |> mutate(fips = f)
}, .progress = TRUE) |>
  list_rbind()

say("optim converged: ", sum(hess$converged), " of ", nrow(hess))


# ---- 3. verification -------------------------------------------------------

check <- cf |>
  select(fips, mu_cf = mu, sigma_cf = sigma) |>
  inner_join(hess |> filter(converged) |>
               transmute(fips, mu_op = mu, sigma_op = exp(log_sigma)),
             by = "fips") |>
  mutate(d_mu = abs(mu_cf - mu_op), d_sigma = abs(sigma_cf - sigma_op))

say(sprintf("closed form vs optim: |d mu| median %.2e, p99 %.2e",
            median(check$d_mu), quantile(check$d_mu, 0.99)))
say(sprintf("                      |d sigma| median %.2e, p99 %.2e",
            median(check$d_sigma), quantile(check$d_sigma, 0.99)))

# Loose, because optim only has to land in the same place, not to the last
# decimal. A disagreement past this is structural.
stopifnot("closed form and optim disagree -- check the band edges in both" =
            quantile(check$d_mu, 0.99) < 0.01)
say("[ok] closed form and optim agree")

# The identity that says the fold happened: a folded fit must put F(10) back
# on o1, because that is the constraint it was built from. If F(10) instead
# lands near g1 alone, the home tile was dropped somewhere.
ident <- cf |>
  mutate(F10 = ln_cdf(10, mu, sigma),
         F100 = ln_cdf(100, mu, sigma),
         d10 = abs(F10 - o1),
         d100 = abs(F100 - (o1 + o2)))
stopifnot("F(10) does not reproduce o1" = max(ident$d10) < 1e-8)
stopifnot("F(100) does not reproduce o1 + o2" = max(ident$d100) < 1e-8)
say(sprintf("[ok] F(10) = o1 and F(100) = o1 + o2 (max error %.1e)",
            max(c(ident$d10, ident$d100))))

# The blunt sanity check from the flow script: most time is spent near home,
# so P(d < 10 km) should sit near 0.8. A median near 1.0 means a raw
# log_sigma was read somewhere as a sigma.
say(sprintf("P(d < 10 km): median %.3f  (expect ~0.8)", median(ident$F10)))
say(sprintf("median displacement exp(mu): %.3f km  (expect small -- this is ",
            median(exp(cf$mu))))
say("             time-weighted and includes time at home)")


# ---- 4. assemble -----------------------------------------------------------

out <- cf |>
  select(fips, gid_2, f0, o1, o2, o3, mu, log_sigma, sigma,
         n_periods, shared_meta_unit, coverage) |>
  left_join(hess |> select(fips, cov_11, cov_12, cov_21, cov_22, converged),
            by = "fips", relationship = "one-to-one") |>
  mutate(meta_obs_weight = META_OBS_WEIGHT,
         meta_scope = META_SCOPE) |>
  arrange(fips) |>
  stamp_run()

stopifnot(!anyDuplicated(out$fips))


# ---- 5. parity against the existing file -----------------------------------

if (file.exists(META_PARAMS_FILE)) {
  
  old <- read_keyed_csv(META_PARAMS_FILE, key_cols = "fips")
  
  cmp <- out |>
    select(fips, mu_new = mu, ls_new = log_sigma) |>
    inner_join(old |> select(fips, mu_old = mu, ls_old = log_sigma),
               by = "fips")
  
  say("parity against existing file: ", nrow(cmp), " counties in common, ",
      nrow(out) - nrow(cmp), " new, ",
      nrow(old) - nrow(cmp), " in the old file only")
  say(sprintf("  |d mu|        median %.2e, max %.2e",
              median(abs(cmp$mu_new - cmp$mu_old)),
              max(abs(cmp$mu_new - cmp$mu_old))))
  say(sprintf("  |d log_sigma| median %.2e, max %.2e",
              median(abs(cmp$ls_new - cmp$ls_old)),
              max(abs(cmp$ls_new - cmp$ls_old))))
  
  worst <- cmp |>
    mutate(d = abs(mu_new - mu_old)) |>
    slice_max(d, n = 5)
  say("  largest differences:")
  print(as.data.frame(worst), row.names = FALSE)
  
} else {
  say("no existing file to compare against")
}


# ---- 6. write --------------------------------------------------------------

target <- if (OVERWRITE || !file.exists(META_PARAMS_FILE)) {
  META_PARAMS_FILE
} else {
  paste0(META_PARAMS_FILE, ".new")
}

write_atomic(out, target)
say("wrote ", nrow(out), " rows to ", target)

if (!identical(target, META_PARAMS_FILE)) {
  say("")
  say("The existing file was NOT replaced. Once the differences above are")
  say("understood, set OVERWRITE <- TRUE at the top and rerun.")
}