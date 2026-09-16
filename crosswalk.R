#!/usr/bin/env Rscript
# ============================================================================
# Crosswalk: (county name, state) -> GADM gid_2, derived from the MD file
# ----------------------------------------------------------------------------
# The MD file gives gadm_id + gadm_name (county name only, NO state), so name
# matching alone is ambiguous -- there are 30+ Washington Counties nationally.
#
# The fix: gadm_id looks like "USA.<state_idx>.<county_idx>_1", and <state_idx>
# identifies the state. Rather than hardcode that mapping (a wrong guess shifts
# every subsequent state by one and fails SILENTLY), we derive it: for each
# state_idx, compare its set of county names against every real US state's
# county list and take the best overlap. Alabama's 67 counties don't
# accidentally match Alaska's, so the pairing is unambiguous and self-checking.
#
# Outputs a crosswalk CSV, and -- importantly -- a list of counties that failed
# to match, so gaps are visible rather than silently dropped.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(stringr); library(tidyr)
})

# ---- config -----------------------------------------------------------------
source("config.R")

BENCH_PATH <- "E:/pandora-out"   # the benchmark file itself (machine-specific)
OUT_XWALK  <- XWALK_FILE
OUT_BENCH  <- BENCHMARK_CSV

# ---- name normalisation ----------------------------------------------------
# Aggressive but reversible-in-spirit: the goal is that "St. Louis City",
# "St Louis city" and "Saint Louis" all collapse to the same token.
norm_county <- function(x) {
  x %>%
    as.character() %>%
    iconv(to = "ASCII//TRANSLIT") %>%          # Dona Ana -> Dona Ana
    tolower() %>%
    str_replace_all("[.'`]", "") %>%           # drop periods / apostrophes
    str_replace_all("\\bst\\b", "saint") %>%   # st -> saint (after periods gone)
    str_replace_all("\\bste\\b", "sainte") %>%
    # strip trailing type words (may repeat, e.g. "city and borough")
    str_remove_all("\\b(county|parish|borough|census area|municipality|city and borough|city|town|municipio)\\b") %>%
    str_replace_all("[^a-z0-9]+", " ") %>%     # hyphens/underscores -> space
    str_squish()
}

# ---- 1. distinct GADM counties from the MD file ----------------------------
md <- fread(MD_PATH, select = c("gadm_id", "gadm_name", "country"),
            showProgress = FALSE)

gadm <- md %>%
  filter(country == "USA") %>%
  distinct(gadm_id, gadm_name) %>%
  mutate(state_idx = as.integer(str_match(gadm_id, "^USA\\.(\\d+)\\.")[, 2]),
         gadm_clean = norm_county(gadm_name)) %>%
  filter(!is.na(state_idx))

message(sprintf("GADM counties in MD file: %d across %d state indices",
                nrow(gadm), n_distinct(gadm$state_idx)))

# ---- 2. reference county list per real state (no download needed) ----------
ref <- tigris::fips_codes %>%
  transmute(state, state_name,
            fips = paste0(state_code, county_code),
            ref_clean = norm_county(county)) %>%
  filter(!is.na(ref_clean), ref_clean != "")

# ---- 3. derive state_idx -> state, by county-set overlap --------------------
# For each (state_idx, candidate state) count how many county names coincide.
overlap <- gadm %>%
  select(state_idx, gadm_clean) %>%
  inner_join(ref %>% select(state, ref_clean), by = c("gadm_clean" = "ref_clean"),
             relationship = "many-to-many") %>%
  count(state_idx, state, name = "n_match")

state_map <- overlap %>%
  group_by(state_idx) %>%
  arrange(desc(n_match), .by_group = TRUE) %>%
  summarise(state       = first(state),
            best        = first(n_match),
            runner_up   = if (n() > 1) nth(n_match, 2) else 0L,
            .groups     = "drop") %>%
  mutate(margin = best - runner_up)

# Sanity: the best match should beat the runner-up decisively, and each real
# state should be claimed exactly once.
weak <- state_map %>% filter(margin < 3)
if (nrow(weak) > 0) {
  message("[warn] state indices with a weak/ambiguous match (inspect these):")
  print(weak)
}
dupes <- state_map %>% count(state) %>% filter(n > 1)
if (nrow(dupes) > 0) {
  message("[warn] a state was claimed by more than one index:")
  print(dupes)
}
message(sprintf("Mapped %d state indices; median match margin = %d",
                nrow(state_map), median(state_map$margin)))

# ---- 4. full crosswalk: gid_2 <-> (state, clean county name) <-> fips -------
xwalk <- gadm %>%
  left_join(state_map %>% select(state_idx, state), by = "state_idx") %>%
  left_join(ref %>% select(state, ref_clean, fips, state_name),
            by = c("state" = "state", "gadm_clean" = "ref_clean")) %>%
  select(gid_2 = gadm_id, gadm_name, state, state_name, fips, clean = gadm_clean)

# Within-state duplicate clean names (the Virginia independent-city problem:
# "Richmond County" and "Richmond city" both normalise to "richmond").
ambig <- xwalk %>% count(state, clean) %>% filter(n > 1)
if (nrow(ambig) > 0) {
  message(sprintf("[warn] %d within-state name collisions (likely VA cities). ",
                  nrow(ambig)),
          "These gid_2 <-> fips pairings may be swapped; check if any are in your draw:")
  print(head(ambig, 20))
}

fwrite(xwalk, OUT_XWALK)
message(sprintf("Wrote crosswalk -> %s (%d rows, %d with FIPS)",
                OUT_XWALK, nrow(xwalk), sum(!is.na(xwalk$fips))))

# ---- 5. attach gid_2 to the benchmark draw ---------------------------------
# fread reads any delimited text regardless of extension; if your draw is an
# RDS instead, swap in readRDS().
bench <- fread(BENCH_PATH)

# Benchmark has `state` as a 2-letter abbreviation and `county_name`.
bench_x <- bench %>%
  mutate(clean = norm_county(county_name)) %>%
  left_join(xwalk %>% select(gid_2, state, clean, gadm_name),
            by = c("state", "clean"))

missed <- bench_x %>% filter(is.na(gid_2))
message(sprintf("\nBenchmark counties matched: %d / %d",
                sum(!is.na(bench_x$gid_2)), nrow(bench_x)))

if (nrow(missed) > 0) {
  message("[unmatched] fix these by hand (nearest GADM names shown):")
  for (i in seq_len(nrow(missed))) {
    cand <- xwalk %>% filter(state == missed$state[i])
    near <- cand$gadm_name[agrep(missed$clean[i], cand$clean, max.distance = 0.3)]
    message(sprintf("  %s, %s  ->  %s",
                    missed$county_name[i], missed$state[i],
                    if (length(near)) paste(head(near, 4), collapse = " | ") else "(no near match)"))
  }
}

fwrite(bench_x, OUT_BENCH)
message(sprintf("Wrote benchmark+gid_2 -> %s", OUT_BENCH))