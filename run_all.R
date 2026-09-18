#!/usr/bin/env Rscript
# ===========================================================================
# run_all.R
# ---------------------------------------------------------------------------
# The whole pipeline, in order, from a fresh clone with the data in place.
#
#   Rscript run_all.R              everything not already built
#   Rscript run_all.R --force      rebuild everything
#   Rscript run_all.R --from 11    from stage 11 onward
#   Rscript run_all.R --only 20    one stage
#
# Stages are skipped when their output already exists, so a rerun after a
# crash picks up where it stopped. --force overrides that.
#
# NOT INCLUDED: the Veraset build (stages 00-07). Those need the licensed
# data and hours of scanning, and their outputs are treated as inputs here.
# See README for how to regenerate them.
#
# NOTE: run_all() -- with brackets -- is a different thing: the cross-county
# summary function in R/07_flows.R, called by stage 20.
# ===========================================================================

source("config.R")
source("R/01_utils.R")

args <- commandArgs(trailingOnly = TRUE)

FORCE <- "--force" %in% args
FROM  <- if ("--from" %in% args) as.numeric(args[which(args == "--from") + 1]) else 0
ONLY  <- if ("--only" %in% args) as.numeric(args[which(args == "--only") + 1]) else NA


# ---- the stages ------------------------------------------------------------
# `produces` is what the stage writes. If it exists and --force was not
# given, the stage is skipped.

STAGES <- list(
  list(n = 1,  script = "R/01_check_inputs.R",
       label = "check inputs",     produces = NULL),
  list(n = 10, script = "R/10_meta_shares.R",
       label = "Meta shares",      produces = META_SHARES_FILE),
  list(n = 11, script = "R/11_fit_meta.R",
       label = "Meta kernel fits", produces = META_PARAMS_FILE),
  list(n = 20, script = "R/20_flows.R",
       label = "flows and maps",   produces = file.path(OUT_DIR, "flow_summary.csv"))
)


# ---- run -------------------------------------------------------------------

banner(sprintf("PanDORA  (%s)", RUN_TAG))
say("git ", git_sha() %||% "not a repo",
    " | scope ", SCOPE, " | tau ", TAU, " | home rule ", HOME_RULE)

t_start <- Sys.time()
ran <- character(0)

for (st in STAGES) {
  
  if (!is.na(ONLY) && st$n != ONLY) next
  if (is.na(ONLY) && st$n < FROM) next
  
  # Stage 1 writes nothing and is always worth running: it is the check that
  # the inputs are where config says they are.
  skip <- !FORCE && !is.null(st$produces) && file.exists(st$produces)
  
  if (skip) {
    say(sprintf("[%2d] %-18s skipped (output exists)", st$n, st$label))
    next
  }
  
  say("")
  banner(sprintf("[%d] %s", st$n, st$label))
  t0 <- Sys.time()
  
  # Each stage runs in its own environment, so one leaving a variable behind
  # cannot change how the next one behaves. That is the difference between a
  # pipeline and a long console session, and it is what makes a rerun mean
  # the same thing as a first run.
  env <- new.env(parent = globalenv())
  result <- try(sys.source(st$script, envir = env), silent = TRUE)
  
  if (inherits(result, "try-error")) {
    say("")
    say("FAILED at stage ", st$n, ": ", st$script)
    say(conditionMessage(attr(result, "condition")))
    quit(status = 1)
  }
  
  say(sprintf("[%2d] %s done in %.1f min", st$n, st$label,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  ran <- c(ran, st$label)
}

say("")
banner(sprintf("finished in %.1f min",
               as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

if (length(ran) > 0) {
  say("ran: ", paste(ran, collapse = ", "))
} else {
  say("nothing to do -- every output already exists. Use --force to rebuild.")
}

say("maps:    ", file.path(OUT_DIR, "maps"))
say("summary: ", file.path(OUT_DIR, "flow_summary.csv"))