# ===========================================================================
# R/08_maps.R
# ---------------------------------------------------------------------------
# The panels. Ported from plot_flow_difference.R, with one change: the
# comparison takes a PAIR of sources rather than assuming Veraset and Meta.
#
#   compare_county("06037")                      meta vs observed (the goal)
#   compare_county("06037", b = "veraset")       meta vs Veraset's kernel
#   compare_county("06037", a = "veraset", b = "obs")   allocation-rule error
#
# Everything else -- the fixed log2 scale, the grey origin, the negligible
# cutoff, the two-panel layout -- is unchanged, so maps made now are
# comparable with the ones already saved.
# ===========================================================================

SOURCE_LABEL <- c(meta = "Meta", veraset = "Veraset kernel", obs = "Veraset observed")


# ---- geometry --------------------------------------------------------------
# Conus Albers (EPSG 5070): equal area, so a county's size on the page is
# proportional to its size on the ground. Cached in the closure because
# reading GADM is slow and it never changes.

.county_sf <- NULL

get_county_sf <- function() {
  if (!is.null(.county_sf)) return(.county_sf)
  adm2 <- geodata::gadm(country = "USA", level = 2, path = GEODATA_DIR)
  sfd <- sf::st_as_sf(adm2) |>
    dplyr::select(gid_2 = GID_2, gadm_name = NAME_2, state = NAME_1) |>
    sf::st_transform(5070)
  .county_sf <<- sfd
  sfd
}

# `extra` is added before coord_sf; a geom_sf tacked on afterwards brings its
# own default coord and ggplot warns about replacing the one already there.
map_base <- function(g, bb, extra = NULL) {
  g + ggplot2::geom_sf(colour = NA) + extra +
    ggplot2::coord_sf(xlim = bb[c("xmin", "xmax")], ylim = bb[c("ymin", "ymax")],
                      expand = FALSE) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(legend.position = "bottom",
                   legend.key.width = ggplot2::unit(22, "pt"),
                   legend.key.height = ggplot2::unit(7, "pt"))
}

DIVERGING <- function(limits, breaks = ggplot2::waiver()) {
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "grey94", high = "#B2182B",
                                midpoint = 0, limits = limits, breaks = breaks,
                                oob = scales::squish, name = NULL)
}


# ===========================================================================
# THE ENTRY POINT
# ---------------------------------------------------------------------------
# Two panels, because they answer different questions and a county can be
# loud in one and invisible in the other:
#
#   log ratio    scale-free. A county getting 20 instead of 10 reads the same
#                as one getting 20,000 instead of 10,000, so it shows where
#                the sources disagree proportionally -- usually the sparse
#                far field.
#   difference   in people. Shows where the disagreement is large enough to
#                matter for a flow matrix, which is nearly always near home.
#
# A third panel appears when the comparison is against observed data: the
# discrepancy divided by its sampling standard deviation. Without it, a
# destination seen by six devices produces the loudest colour on the map
# purely by being thinly sampled.
#
# min_flow is in people. Counties where neither source places at least this
# many are grey: otherwise hundreds of distant counties receiving a fraction
# of a person, where the log ratio is extreme and meaningless, outvote the
# dozen counties holding all the people.
# ===========================================================================

compare_county <- function(fips,
                           a = "meta", b = NULL,
                           which = c("total", "shape"),
                           min_flow = MAP_MIN_FLOW,
                           county = COUNTY, obs = OBSERVED, xw = XW,
                           refresh = FALSE) {
  
  which <- match.arg(which)
  res <- county_flows(fips, county = county, obs = obs, xw = xw, refresh = refresh)
  d <- res$diagnostics
  
  # Default: compare against observed if it exists, else the Veraset kernel.
  if (is.null(b)) b <- if (isTRUE(d$has_observed)) "obs" else "veraset"
  
  ca <- paste0(if (which == "total") "flow_" else "share_", a)
  cb <- paste0(if (which == "total") "flow_" else "share_", b)
  if (!all(c(ca, cb) %in% names(res$flows))) {
    stop("This county has no ", if (!ca %in% names(res$flows)) a else b,
         " layer. Available: ",
         paste(sub("^flow_", "", grep("^flow_", names(res$flows), value = TRUE)),
               collapse = ", "))
  }
  
  sfd <- get_county_sf()
  
  message(sprintf(
    "
origin        %s (%s)
population    %s
volume        N = pop  (home is inside the kernel, not the volume)
kernels       %s
away          %s
under 10 km   %s   <- Meta's (0, 10) leaves this split free
beyond %d km  %s   <- dropped from the domain
slices        %.2f km wide, empty-slice mass %s
panel         %s devices%s
",
    res$fips, res$county_name %||% "", format(round(res$pop), big.mark = ","),
    paste(sprintf("%s mu %+.3f sigma %.3f", names(res$kernels),
                  purrr::map_dbl(res$kernels, "mu"),
                  purrr::map_dbl(res$kernels, "sigma")), collapse = " | "),
    paste(sprintf("%s %.2f%%", names(d$away), 100 * d$away), collapse = "  "),
    paste(sprintf("%s %.3f", names(d$under_10km), d$under_10km), collapse = "  "),
    D_MAX,
    paste(sprintf("%s %.4f", names(d$tail_beyond), d$tail_beyond), collapse = "  "),
    d$slice_km,
    paste(sprintf("%s %.4f", names(d$empty_mass), d$empty_mass), collapse = "  "),
    format(d$n_devices, big.mark = ","),
    if (isTRUE(d$kernel_reliable)) "" else
      sprintf(" [below %s -- Veraset kernel not reliably estimated alone]",
              format(MIN_DEVICES_KERNEL, big.mark = ","))))
  
  fl <- res$flows |>
    dplyr::mutate(va = .data[[ca]], vb = .data[[cb]]) |>
    dplyr::filter(va > min_flow | vb > min_flow) |>
    dplyr::mutate(
      # log2 rather than natural log: one unit is a doubling, which is easier
      # to read off a map, and it is what LOG2_LIM is expressed in.
      value_log = log2(va / vb),
      value_abs = va - vb,
      is_origin = gid_2 == res$gid_2
    ) |>
    dplyr::filter(is.finite(value_log), is.finite(value_abs))
  
  # The origin county gets its own flat shade rather than a place on either
  # scale. It holds most of the population, so on the absolute panel it would
  # compress every other county to white; and it is not the same kind of
  # comparison -- it is where each source puts the time it sends nowhere.
  mp <- sfd |>
    dplyr::inner_join(fl, by = "gid_2") |>
    dplyr::mutate(negligible = pmax(va, vb) < min_flow,
                  dplyr::across(c(value_log, value_abs),
                                \(x) dplyr::if_else(is_origin | negligible,
                                                    NA_real_, x)))
  bb <- sf::st_bbox(mp)
  
  # Each county is coloured the same regardless of how many people it
  # receives, so the map answers "where do they disagree?" and NOT "which
  # sends more?". Those can point opposite ways. This is the flow-weighted
  # direction, which is the one the away figures describe.
  wlog2 <- log2(sum(mp$va[!mp$is_origin]) / sum(mp$vb[!mp$is_origin]))
  n_shown <- sum(!mp$is_origin & !mp$negligible)
  
  lim_abs <- mp |>
    sf::st_drop_geometry() |>
    dplyr::pull(value_abs) |>
    quantile(c(0.02, 0.98), na.rm = TRUE) |>
    abs() |> max()
  
  origin_layer <- ggplot2::geom_sf(data = dplyr::filter(mp, is_origin),
                                   fill = "grey55", colour = "black",
                                   linewidth = 0.4)
  
  pA <- map_base(ggplot2::ggplot(mp, ggplot2::aes(fill = value_log)), bb,
                 origin_layer) +
    DIVERGING(c(-LOG2_LIM, LOG2_LIM), seq(-LOG2_LIM, LOG2_LIM, by = 1)) +
    ggplot2::labs(title = sprintf("log2(%s / %s), %s",
                                  SOURCE_LABEL[[a]], SOURCE_LABEL[[b]], which),
                  subtitle = "per county, unweighted")
  
  pB <- map_base(ggplot2::ggplot(mp, ggplot2::aes(fill = value_abs)), bb,
                 origin_layer) +
    DIVERGING(c(-lim_abs, lim_abs)) +
    ggplot2::labs(title = sprintf("%s - %s, %s",
                                  SOURCE_LABEL[[a]], SOURCE_LABEL[[b]], which),
                  subtitle = if (which == "total") "people; where the weight is"
                  else "share points; where the weight is")
  
  panels <- list(pA, pB)
  
  # ---- z panel, only against observed ---------------------------------------
  # (predicted share - observed share) / se(observed share), where the se is
  # multinomial: sqrt(p(1-p)/n) on the origin's panel size. Reads as "how many
  # standard deviations of sampling noise is this discrepancy worth", so a
  # thin destination has to be badly wrong to show colour.
  if (b == "obs" && "z_meta" %in% names(res$flows) && a == "meta") {
    mz <- sfd |>
      dplyr::inner_join(dplyr::select(fl, gid_2, is_origin, negligible = va),
                        by = "gid_2") |>
      dplyr::left_join(dplyr::select(res$flows, gid_2, z_meta), by = "gid_2") |>
      dplyr::mutate(z = dplyr::if_else(is_origin, NA_real_, z_meta))
    
    pZ <- map_base(ggplot2::ggplot(mz, ggplot2::aes(fill = z)), bb,
                   origin_layer) +
      DIVERGING(c(-5, 5), seq(-5, 5, by = 2.5)) +
      ggplot2::labs(title = "z: (Meta - observed) / sampling sd",
                    subtitle = "discrepancy in units of multinomial noise")
    panels <- c(panels, list(pZ))
  }
  
  out <- if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(panels, nrow = 1) +
      patchwork::plot_annotation(
        title = sprintf("Flows from %s (%s): %s vs %s",
                        res$county_name %||% "", res$fips,
                        SOURCE_LABEL[[a]], SOURCE_LABEL[[b]]),
        subtitle = sprintf(
          "of %s residents | flow-weighted overall: log2 %+.2f | %d counties coloured, origin and under %g person grey | log scale fixed +/-%g",
          format(round(res$pop), big.mark = ","), wlog2, n_shown, min_flow, LOG2_LIM),
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 13)))
  } else {
    panels
  }
  
  print(res$flows |>
          dplyr::slice_head(n = 12) |>
          dplyr::transmute(gid_2,
                           a = round(.data[[ca]]),
                           b = round(.data[[cb]]),
                           diff = round(.data[[ca]] - .data[[cb]]),
                           log2 = round(log2(.data[[ca]] / .data[[cb]]), 3)) |>
          dplyr::rename(!!a := a, !!b := b))
  
  invisible(list(plot = out, result = res, mapped = mp))
}