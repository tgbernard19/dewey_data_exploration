#!/usr/bin/env Rscript
# ===========================================================================
# cpc.R
# Common Part of Commuters between the Veraset-kernel and Meta-kernel
# flows, per origin county and averaged.
#
# source("absolute_flows_compare.R")
# source("cpc.R")
# out <- cpc_all(c("06037", "48201", "25015", "37009", "30109"))
# cpc_report(out)
# ---------------------------------------------------------------------------
#
# CPC(T, T') = 2 * sum_j min(T_j, T'_j) / (sum_j T_j + sum_j T'_j)
#
# 1 is identical, 0 is no overlap. It is the share of people both kernels
# agree about, so it reads directly: CPC 0.92 means the two allocations
# place 92% of residents the same way.
#
# TWO NUMBERS, AND THE SECOND IS THE INFORMATIVE ONE
#
# cpc_all every destination, the origin county included. Most people are
# at or near home, and both kernels know that, so this is high
# (0.9+) almost everywhere and mostly measures the agreement
# about the diagonal.
# cpc_off the origin county dropped and each side renormalised over what
# is left. This is agreement about where travellers GO, which is
# what the kernels are actually being compared on.
#
# Report both. cpc_off will be markedly lower, and that is not a problem
# with it -- it is the diagonal no longer doing the work.
#
# THE TAIL BEYOND D_MAX
#
# Each source leaves pop * (1 - P(d <= D_MAX)) unallocated, and the two
# leave different amounts (Veraset's tail is the heavier one). That is real
# disagreement, so it enters as one extra destination called "BEYOND"
# rather than being dropped. Set INCLUDE_BEYOND to FALSE to see what it was
# contributing.
# ===========================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble); library(ggplot2)
})

stopifnot("Run source('absolute_flows_compare.R') first." = exists("county_flows"))

INCLUDE_BEYOND <- TRUE

# The formula, on two vectors of the same length and in the same order.
cpc <- function(a, b) {
  stopifnot(length(a) == length(b))
  denom <- sum(a) + sum(b)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  2 * sum(pmin(a, b)) / denom
}

# ---- one county -------------------------------------------------------------
cpc_county <- function(f, refresh = FALSE) {
  
  r <- tryCatch(county_flows(f, refresh = refresh), error = conditionMessage)
  if (is.character(r)) { message(f, ": ", r); return(NULL) }
  
  fl <- r$flows |> select(gid_2, flow_v, flow_m)
  
  if (INCLUDE_BEYOND)
    fl <- fl |>
    add_row(gid_2 = "BEYOND",
            flow_v = r$pop * r$diagnostics$tail_beyond_v,
            flow_m = r$pop * r$diagnostics$tail_beyond_m)
  
  # Off-diagonal: drop the origin, then renormalise each side over its own
  # remaining total. Without that, a difference in how many people travel at
  # all would show up as disagreement about where they went.
  off <- fl |> filter(gid_2 != r$gid_2)
  
  tibble(
    fips = f,
    pop = r$pop,
    cpc_all = cpc(fl$flow_v, fl$flow_m),
    cpc_off = cpc(off$flow_v / sum(off$flow_v), off$flow_m / sum(off$flow_m)),
    away_v = r$diagnostics$away_v,
    away_m = r$diagnostics$away_m,
    n_dest = nrow(off)
  )
}

# ---- many counties ----------------------------------------------------------
cpc_all <- function(fips_vec, refresh = FALSE) {
  out <- map(pad_fips(fips_vec), cpc_county, refresh = refresh) |> list_rbind()
  message(sprintf("[cpc] %d of %d counties returned a row",
                  nrow(out), length(fips_vec)))
  out
}

# ---- display ----------------------------------------------------------------
# The mean is unweighted, so it is the typical county rather than the
# typical person; population-weighted, a handful of metros would be the
# whole answer. Both are printed.
cpc_report <- function(out) {
  
  summary_tbl <- out |>
    summarise(across(c(cpc_all, cpc_off),
                     list(mean = \(x) mean(x, na.rm = TRUE),
                          median = \(x) median(x, na.rm = TRUE),
                          min = \(x) min(x, na.rm = TRUE),
                          max = \(x) max(x, na.rm = TRUE)))) |>
    pivot_longer(everything(), names_to = c("measure", "stat"),
                 names_pattern = "(cpc_all|cpc_off)_(.*)") |>
    pivot_wider(names_from = stat, values_from = value) |>
    mutate(pop_weighted_mean = c(weighted.mean(out$cpc_all, out$pop),
                                 weighted.mean(out$cpc_off, out$pop)))
  print(summary_tbl)
  
  plot_df <- out |>
    select(fips, pop, cpc_all, cpc_off) |>
    pivot_longer(c(cpc_all, cpc_off), names_to = "measure", values_to = "cpc") |>
    mutate(measure = recode(measure,
                            cpc_all = "all destinations",
                            cpc_off = "excluding the origin county"))
  
  p <- ggplot(plot_df, aes(reorder(fips, cpc), cpc, colour = measure)) +
    geom_point(size = 2) +
    coord_flip(ylim = c(0, 1)) +
    scale_colour_manual(values = c("all destinations" = "grey60",
                                   "excluding the origin county" = "#B2182B")) +
    labs(x = NULL, y = "CPC", colour = NULL,
         title = "Agreement between the Veraset and Meta allocations",
         subtitle = sprintf("mean %.3f all destinations, %.3f excluding the origin; %d counties",
                            mean(out$cpc_all, na.rm = TRUE),
                            mean(out$cpc_off, na.rm = TRUE), nrow(out))) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  
  print(p)
  invisible(list(summary = summary_tbl, plot = p))
}

# ---------------------------------------------------------------------------
# usage
#
#source("absolute_flows_compare.R")
# source("cpc.R")

top_40_fips <- c(
  "06037", "17031", "48201", "04013", "06073", "06059", "12086", "48113",
  "36047", "06065", "32003", "53033", "36081", "48439", "06071", "48029",
  "12011", "06085", "26163", "25017", "36061", "06001", "06067", "12099",
  "12057", "42101", "36103", "12095", "36005", "48453", "36059", "39049",
  "48085", "26125", "27053", "37183", "37119", "39035", "49035", "42003"
)

out <- cpc_all(top_40_fips)
cpc_report(out)
#
# # CPC and the total-variation distance from run_all() are the same
# # quantity when both sides carry equal totals: CPC = 1 - TV. They differ
# # here only because N_v and N_m differ by the tail mass.
left_join(out, run_all(out$fips), by = "fips") |>
mutate(check = cpc_all + tv_allocation) |>
select(fips, cpc_all, tv_allocation, check)
# ---------------------------------------------------------------------------