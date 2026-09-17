# ===========================================================================
# R/04_lognormal.R
# ---------------------------------------------------------------------------
# THE lognormal kernel. One definition of the survival function, the band
# probabilities, the likelihood, and the two fitters.
#
# Before consolidation there were at least four copies of this material, and
# they had drifted: two different nll() signatures, one taking the device
# count as an argument and one reading it from the data, plus two different
# values of the Meta pseudo-weight. Everything that fits a lognormal now
# comes through this file.
#
# WHY LOGNORMAL
#   On clean CBG-home profiles it misplaces 4.2% of mass against the
#   two-exponential mixture's 8.1%, with one fewer parameter. More usefully,
#   Meta's three categories identify it exactly in closed form, which removes
#   the ridge, the arbitrary prior on p, and the whole Gibbs-posterior
#   apparatus that existed only because a three-parameter family was being
#   fitted to two-parameter data.
#
# PARAMETERISATION
#   Optimiser scale is c(mu, log_sigma); natural scale is (mu, sigma).
#   sigma must stay positive and the optimiser must not have to know that,
#   hence the log. Every function here takes one or the other explicitly --
#   never guess which one a column holds. A raw log_sigma of 0.70 read as a
#   sigma gives a kernel that is merely far too tight, with no error anywhere.
# ===========================================================================


# ---- survival and distribution --------------------------------------------

# P(D > d) for a lognormal. d = 0 gives 1, since log(0) is -Inf and
# pnorm(-Inf) is 0 -- worth checking rather than assuming, because the first
# band edge is always zero.
ln_surv <- function(d, mu, sigma) {
  1 - pnorm((log(pmax(d, 1e-12)) - mu) / sigma)
}

ln_cdf <- function(d, mu, sigma) {
  pnorm((log(pmax(d, 1e-12)) - mu) / sigma)
}

# The CDF renormalised over [lo, hi]. This is what the flow allocation uses,
# so that both sources spread their mass over identical support and a heavier
# tail cannot make one look smaller everywhere inside the domain.
ln_cdf_trunc <- function(d, mu, sigma, lo = D_MIN, hi = D_MAX) {
  f_lo <- ln_cdf(lo, mu, sigma)
  f_hi <- ln_cdf(hi, mu, sigma)
  denom <- f_hi - f_lo
  if (!is.finite(denom) || denom <= 0) {
    stop("Window [", lo, ", ", hi, "] carries no mass for mu = ", mu,
         ", sigma = ", sigma)
  }
  pmin(1, pmax(0, (ln_cdf(d, mu, sigma) - f_lo) / denom))
}


# ---- band probabilities ----------------------------------------------------

# Probability of landing in each [edge_low, edge_high), renormalised to sum
# to 1 across the bands supplied.
#
# That renormalisation IS the truncation: pass nine bands ending at 500 km
# and the mass beyond 500 is redistributed rather than ignored, so the fitted
# (mu, sigma) still describe the same untruncated kernel. Pass all eleven and
# nothing is truncated. The caller decides by choosing which rows to hand in.
predict_bands <- function(mu, sigma, edge_low, edge_high) {
  q <- ln_surv(edge_low, mu, sigma) - ln_surv(edge_high, mu, sigma)
  q / sum(q)
}


# ---- fit statistics --------------------------------------------------------

# Total variation: the share of probability mass sitting in the wrong band.
# The 0.5 is what makes that reading exact -- mass removed from one band
# lands in another, so without it every error is counted twice.
tv_distance <- function(w, q) {
  0.5 * sum(abs(w - q))
}

# Negative multinomial log-likelihood. n is the effective sample size: the
# DEVICE count, not the device-day count and not the row count, because
# device-days within a device are not independent draws.
#
# Returns a large finite number rather than Inf on a bad step, so the
# optimiser backs away instead of crashing.
ln_nll <- function(par, bands, n) {
  q <- predict_bands(par[1], exp(par[2]), bands$edge_low, bands$edge_high)
  if (any(!is.finite(q)) || any(q <= 0)) return(1e10)
  -n * sum(bands$w * log(q))
}


# ---- closed form from Meta's three categories ------------------------------

# Meta gives three shares with the home tile folded into the first:
# o1 on (0, 10), o2 on [10, 100), o3 beyond 100. Two parameters, two
# constraints -- exactly identified, no optimiser needed.
#
#   z1 = qnorm(o1)            so F(10)  = o1
#   z2 = qnorm(o1 + o2)       so F(100) = o1 + o2
#   sigma = log(100/10) / (z2 - z1)
#   mu    = log(10) - sigma * z1
#
# o3 is not used: it is determined by the other two, and it is the check.
ln_closed_form <- function(o1, o2) {
  z1 <- qnorm(o1)
  z2 <- qnorm(o1 + o2)
  sigma <- log(10) / (z2 - z1)
  mu <- log(10) - sigma * z1
  tibble::tibble(mu = mu, sigma = sigma)
}


# ---- bounded multi-start fit (the Veraset side) ----------------------------

# Bounds, and why these:
#   mu    in [-8, 7]        median displacement from 0.3 m to 1,100 km
#   sigma in [0.2, 8]       a sigma of 8 is already four orders of magnitude
#                           per standard deviation
#
# A fit landing on a bound has failed. It is FLAGGED, not silently accepted,
# because a runaway mu of 210 looks like a number until you ask what it means.
# Runaway rates track panel size, not family: 48% under 1,000 devices, 0%
# above 10,000.
LN_LOWER <- c(-8, log(0.2))
LN_UPPER <- c(7, log(8))
LN_INIT  <- c(-3, log(4))

# `spread` is the range of objective values across starts. A large spread
# means the surface has more than one basin and the reported fit depends on
# where the optimiser began, which is worth knowing even when it converged.
fit_lognormal <- function(bands, n_starts = 10, seed = NULL) {
  
  stopifnot(all(c("edge_low", "edge_high", "w") %in% names(bands)))
  stopifnot(abs(sum(bands$w) - 1) < 1e-8)
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- if ("n_devices" %in% names(bands)) bands$n_devices[1] else 1
  
  starts <- vector("list", n_starts)
  starts[[1]] <- LN_INIT
  for (i in seq_len(n_starts - 1)) {
    s <- LN_INIT + rnorm(2, sd = 0.75)
    starts[[i + 1]] <- pmin(pmax(s, LN_LOWER + 1e-6), LN_UPPER - 1e-6)
  }
  
  runs <- lapply(starts, function(s) {
    out <- try(optim(s, ln_nll, bands = bands, n = n, method = "L-BFGS-B",
                     lower = LN_LOWER, upper = LN_UPPER,
                     control = list(maxit = 500)),
               silent = TRUE)
    if (inherits(out, "try-error")) NULL else out
  })
  runs <- Filter(Negate(is.null), runs)
  
  if (length(runs) == 0) return(NULL)
  
  vals <- vapply(runs, function(r) r$value, numeric(1))
  best <- runs[[which.min(vals)]]
  
  tol <- 1e-3 * (LN_UPPER - LN_LOWER)
  q <- predict_bands(best$par[1], exp(best$par[2]),
                     bands$edge_low, bands$edge_high)
  
  tibble::tibble(
    mu        = best$par[1],
    sigma     = exp(best$par[2]),
    nll       = best$value,
    spread    = diff(range(vals)),
    at_bound  = any(best$par <= LN_LOWER + tol | best$par >= LN_UPPER - tol),
    converged = best$convergence == 0,
    tv        = tv_distance(bands$w, q),
    n_devices = n
  )
}


# ---- Hessian fit (the Meta side) -------------------------------------------

# Same likelihood, unbounded BFGS, keeping the Hessian so the parameter
# covariance comes back as solve(hessian).
#
# READ THE SCALE CAVEAT. Meta ships fractions with no sample size, so n comes
# from META_OBS_WEIGHT, a made-up number chosen to give BFGS a gradient. The
# point estimates do not depend on it; the covariance scales as 1/n. So the
# covariance is on an arbitrary scale and means nothing in absolute terms.
# Fine for the maps, which use only mu and sigma. Any interval built from it
# is indicative, not a confidence statement.
#
# The columns are named cov_*, not h_*, because they are entries of the
# inverse Hessian. An earlier version named them h_11..h_22, which invites
# someone to invert them a second time.
fit_lognormal_hessian <- function(bands, n = META_OBS_WEIGHT) {
  
  stopifnot(all(c("edge_low", "edge_high", "w") %in% names(bands)))
  
  fit <- try(optim(LN_INIT, ln_nll, bands = bands, n = n,
                   method = "BFGS", hessian = TRUE),
             silent = TRUE)
  
  if (inherits(fit, "try-error") || fit$convergence != 0) {
    return(tibble::tibble(mu = NA_real_, log_sigma = NA_real_,
                          cov_11 = NA_real_, cov_12 = NA_real_,
                          cov_21 = NA_real_, cov_22 = NA_real_,
                          converged = FALSE))
  }
  
  cov_mat <- try(solve(fit$hessian), silent = TRUE)
  if (inherits(cov_mat, "try-error")) {
    return(tibble::tibble(mu = fit$par[1], log_sigma = fit$par[2],
                          cov_11 = NA_real_, cov_12 = NA_real_,
                          cov_21 = NA_real_, cov_22 = NA_real_,
                          converged = FALSE))
  }
  
  tibble::tibble(
    mu        = fit$par[1],
    log_sigma = fit$par[2],
    cov_11    = cov_mat[1, 1],
    cov_12    = cov_mat[1, 2],
    cov_21    = cov_mat[2, 1],
    cov_22    = cov_mat[2, 2],
    converged = TRUE
  )
}


# ===========================================================================
# SELF-TEST
# ---------------------------------------------------------------------------
# source("R/04_lognormal.R"); lognormal_self_test()
#
# Four checks, in order of how much they would cost if they failed silently.
# ===========================================================================

lognormal_self_test <- function(verbose = TRUE) {
  
  ok <- TRUE
  note <- function(...) if (verbose) cat("  ", ..., "\n", sep = "")
  
  # 1. Survival at the edges. S(0) must be exactly 1 or every band
  #    probability is wrong by the same unnoticeable amount.
  s0 <- ln_surv(0, -3, 4)
  sinf <- ln_surv(Inf, -3, 4)
  test1 <- abs(s0 - 1) < 1e-12 && abs(sinf) < 1e-12
  ok <- ok && test1
  note("[", if (test1) "ok" else "FAIL", "] S(0) = 1, S(Inf) = 0")
  
  # 2. Parameter recovery. Generate exact band probabilities from known
  #    parameters, then fit them back. This is the test that the likelihood,
  #    the band probabilities and the optimiser agree with each other.
  edges <- BAND_EDGES[BAND_EDGES <= TRUNC_KM]
  bands <- tibble::tibble(
    edge_low  = head(edges, -1),
    edge_high = tail(edges, -1)
  )
  true_mu <- -1.2
  true_sigma <- 3.1
  bands$w <- predict_bands(true_mu, true_sigma, bands$edge_low, bands$edge_high)
  bands$n_devices <- 50000
  
  f <- fit_lognormal(bands, seed = 1)
  test2 <- !is.null(f) && abs(f$mu - true_mu) < 0.01 &&
    abs(f$sigma - true_sigma) < 0.01 && f$tv < 1e-6
  ok <- ok && test2
  note("[", if (test2) "ok" else "FAIL", "] recovers mu = ", true_mu,
       ", sigma = ", true_sigma,
       " (got ", signif(f$mu, 4), ", ", signif(f$sigma, 4),
       "; tv ", signif(f$tv, 3), ")")
  
  # 3. Closed form against the optimiser. The closed form is the estimator
  #    for the Meta side and optim is only the verification, so a
  #    disagreement means one of them is wrong about the band edges.
  cf <- ln_closed_form(o1 = 0.82, o2 = 0.14)
  meta_bands <- tibble::tibble(
    edge_low  = c(0, 10, 100),
    edge_high = c(10, 100, Inf),
    w         = c(0.82, 0.14, 0.04)
  )
  h <- fit_lognormal_hessian(meta_bands, n = 1000)
  test3 <- h$converged &&
    abs(h$mu - cf$mu) < 1e-3 &&
    abs(exp(h$log_sigma) - cf$sigma) < 1e-3
  ok <- ok && test3
  note("[", if (test3) "ok" else "FAIL", "] closed form matches optim ",
       "(mu ", signif(cf$mu, 4), " vs ", signif(h$mu, 4),
       "; sigma ", signif(cf$sigma, 4), " vs ", signif(exp(h$log_sigma), 4), ")")
  
  # 4. The closed form reproduces its own inputs. F(10) must come back as o1
  #    and F(100) as o1 + o2, which is what "exactly identified" means.
  test4 <- abs(ln_cdf(10, cf$mu, cf$sigma) - 0.82) < 1e-8 &&
    abs(ln_cdf(100, cf$mu, cf$sigma) - 0.96) < 1e-8
  ok <- ok && test4
  note("[", if (test4) "ok" else "FAIL", "] closed form reproduces F(10), F(100)")
  
  if (verbose) cat(if (ok) "\nlognormal: all checks passed\n"
                   else "\nlognormal: SOMETHING FAILED\n")
  invisible(ok)
}