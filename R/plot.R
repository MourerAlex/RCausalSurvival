# ============================================================================
# Plot methods for causal_survival_* S3 classes
# ----------------------------------------------------------------------------
# Adapted from `separable_effects/R/plot.R` (627 lines) for the
# CausalSurvival v0.1.0 binary point-treatment scope. The SE-specific
# four-arm display, decomposition annotations
# (`build_contrast_annotations()`), and the contrast / diagnostic plot
# bodies have been dropped or stubbed per spec §3.4 line 769
# ("placeholders for contrast and diagnostic, message-only in v0.1.0").
# The `risk_table` panel arg is preserved in the signature but errors
# with a deferral message until step 8b lands the underlying
# `causal_risk_table()` accessor.
# ============================================================================

#' Default Okabe-Ito palette for binary treatment arms
#'
#' Returns a named character vector of hex colors keyed by the
#' character representation of the two treatment levels in
#' `fit$treatment_levels`. Used as the per-arm color default by
#' [plot.causal_survival_risk()]. Users override individual entries
#' via the `arm_colors` argument.
#'
#' @param levels_vec Length-2 numeric vector (typically `c(0, 1)`).
#' @return Named character vector of hex strings, length 2.
#' @family internal
#' @keywords internal
default_arm_palette <- function(levels_vec) {
  if (length(levels_vec) != 2L) {
    stop("default_arm_palette requires exactly two levels (v0.1.0 ",
         "binary scope).", call. = FALSE)
  }
  # Lower level = blue (#0072B2 control); higher = vermillion
  # (#D55E00 treated). Both are colorblind-safe Okabe-Ito picks.
  out <- c("#0072B2", "#D55E00")
  names(out) <- as.character(c(min(levels_vec), max(levels_vec)))
  out
}


#' Default linetype mapping for treatment arms
#'
#' Returns a named character vector mapping each arm value to a linetype.
#' Used by [plot.causal_survival_risk()] so the curves remain
#' distinguishable in B&W or colorblind viewing. Users override individual
#' entries via the `arm_linetypes` argument.
#'
#' @param levels_vec Length-2 numeric vector (typically `c(0, 1)`).
#' @return Named character vector of linetypes, length 2.
#' @family internal
#' @keywords internal
default_arm_linetypes <- function(levels_vec) {
  if (length(levels_vec) != 2L) {
    stop("default_arm_linetypes requires exactly two levels (v0.1.0 ",
         "binary scope).", call. = FALSE)
  }
  # Lower level = solid; higher = longdash. Same lower/higher convention
  # as default_arm_palette() so color and linetype stay aligned per arm.
  out <- c("solid", "longdash")
  names(out) <- as.character(c(min(levels_vec), max(levels_vec)))
  out
}


#' Shared ggplot2 theme for causal_survival plot panels
#'
#' Single source of truth for the visual language shared between the
#' curves panel and the risk-table panel(s). Built on `theme_minimal()`
#' with: no minor gridlines, no vertical major gridlines (the x-axis
#' ticks already mark cut times), bold titles / axis titles, legend at
#' the bottom. Panel-specific tweaks are layered on top of this base
#' via further `theme(...)` calls.
#'
#' @param base_size Numeric. Base font size passed to `theme_minimal()`.
#' @return A `ggplot2` theme object.
#' @family internal
#' @keywords internal
.causal_survival_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position    = "bottom",
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title         = element_text(face = "bold"),
      plot.subtitle      = element_text(face = "plain"),
      axis.title.x       = element_text(face = "bold"),
      axis.title.y       = element_text(face = "bold"),
      legend.title       = element_text(face = "bold")
    )
}


#' Plot cumulative incidence / survival curves
#'
#' Plots the counterfactual curve under each arm of the binary
#' treatment as a step function on the reporting grid `fit$cut_times`,
#' with an optional bootstrap CI ribbon (step-transformed to match the
#' curve). Pairs with [causal_risk()] — pass a `"causal_survival_risk"`
#' object.
#'
#' @param x A `"causal_survival_risk"` object from [causal_risk()].
#' @param arms Numeric vector. Subset of `fit$treatment_levels` to
#'   draw. `NULL` (default) draws all levels.
#' @param arm_colors Named character vector overriding arm hex
#'   colors. Names must be the character representation of treatment
#'   levels. Defaults to Okabe-Ito via [default_arm_palette()]; only
#'   the arms you pass are overridden.
#' @param arm_labels Named character vector overriding the legend
#'   labels for each arm. Defaults to `c("0" = "Control", "1" =
#'   "Treated")` when the levels are `c(0, 1)`; otherwise the level
#'   value is used directly.
#' @param arm_linetypes Named character vector overriding the linetype
#'   for each arm. Defaults to `c("solid", "longdash")` via
#'   [default_arm_linetypes()] so curves remain distinguishable in B&W
#'   or colorblind viewing; only the arms you pass are overridden.
#' @param curves Logical. `TRUE` (default) draws the survival /
#'   incidence curves panel; `FALSE` suppresses the curves panel
#'   entirely so only the requested `risk_table` panel(s) render
#'   (table-only mode). Requires a non-NULL `risk_table`.
#' @param risk_table NULL (default — no risk-table panel) or a
#'   character vector of one or more entries in `c("at_risk", "events_y",
#'   "events_y_interval", "censored")`. When a single string, one
#'   labelled count-per-(arm, cut_time) panel renders below the curves
#'   panel. When a vector, one panel per requested count renders stacked
#'   in the supplied order. `"events_y"` is the ggsurvfit-style
#'   cumulative Y-events count; `"events_y_interval"` is the
#'   per-interval count.
#' @param risk_table_height Numeric. Per-panel height ratio relative to
#'   the curves panel (which is `1`). When `risk_table` is a vector of
#'   length `N`, the total bottom area is `N * risk_table_height`.
#'   Default `0.23`.
#' @param cut_times Controls the x-axis tick positions on both the curves
#'   panel and the risk-table panel (when shown). `0` is always included
#'   (auto-prepended if the user forgets it).
#'   \describe{
#'     \item{`NULL` (default)}{Display `c(0, fit$cut_times)`.}
#'     \item{Single positive integer `N`}{Display `c(0, ...)` where
#'       `...` is `N` indices equidistant along `fit$cut_times`.}
#'     \item{Numeric vector (length >= 2)}{Explicit subset, values must
#'       be in `c(0, fit$cut_times)`. `0` is auto-prepended if absent.
#'       Note: a length-1 numeric is interpreted as a count (above), not
#'       a subset — to display a single specific cut time, pair it with
#'       0, e.g. `cut_times = c(0, 5)`.}
#'   }
#' @param title,subtitle Character or NULL. Plot title / subtitle.
#' @param x_label,y_label Axis labels. Defaults reflect the selected
#'   `scale` (incidence vs survival).
#' @param base_size Numeric. Base font size for `theme_minimal()`.
#' @param linewidth Numeric. Width of the step lines. Default 0.8.
#' @param ribbon_alpha Numeric in `[0, 1]`. Transparency of the CI
#'   ribbons. Default 0.15.
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot2 object.
#' @family plot
#' @export
plot.causal_survival_risk <- function(x,
                                      arms = NULL,
                                      arm_colors = NULL,
                                      arm_labels = NULL,
                                      arm_linetypes = NULL,
                                      curves = TRUE,
                                      risk_table = NULL,
                                      risk_table_height = 0.23,
                                      cut_times = NULL,
                                      title = NULL,
                                      subtitle = NULL,
                                      x_label = "Time",
                                      y_label = NULL,
                                      base_size = 11,
                                      linewidth = 0.8,
                                      ribbon_alpha = 0.15,
                                      ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot.causal_survival_risk requires the 'ggplot2' package. ",
         "Install via install.packages('ggplot2').", call. = FALSE)
  }
  if (!is.null(risk_table)) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("risk_table panel requires the 'patchwork' package. ",
           "Install via install.packages('patchwork').",
           call. = FALSE)
    }
    if (is.null(x$pt_data) || is.null(x$id_col) ||
        is.null(x$treatment_col) || is.null(x$cut_times)) {
      stop("risk_table panel needs person-time references on the ",
           "causal_survival_risk object — refit with ",
           "`causal_survival(..., keep_data = TRUE)`.", call. = FALSE)
    }
  }

  # --- Slice ---
  risk_df    <- x$risk
  have_bands <- any(!is.na(risk_df$lower)) && any(!is.na(risk_df$upper))

  # --- Resolve arms ---
  avail_arms <- sort(unique(risk_df$treatment))
  if (is.null(arms)) arms <- avail_arms
  bad <- setdiff(arms, avail_arms)
  if (length(bad) > 0L) {
    stop("Unknown arm(s): ", paste(bad, collapse = ", "),
         ". Available: ", paste(avail_arms, collapse = ", "),
         call. = FALSE)
  }

  # --- Palette + labels (Okabe-Ito default, user-overridable per-arm) ---
  default_colors <- default_arm_palette(avail_arms)
  default_labels <- setNames(as.character(avail_arms),
                                    as.character(avail_arms))
  if (identical(sort(avail_arms), c(0, 1))) {
    default_labels[["0"]] <- "Control"
    default_labels[["1"]] <- "Treated"
  }
  if (!is.null(arm_colors)) {
    bad_c <- setdiff(names(arm_colors), names(default_colors))
    if (length(bad_c) > 0L) {
      stop("arm_colors has unknown entries: ",
           paste(bad_c, collapse = ", "), call. = FALSE)
    }
    default_colors[names(arm_colors)] <- arm_colors
  }
  if (!is.null(arm_labels)) {
    bad_l <- setdiff(names(arm_labels), names(default_labels))
    if (length(bad_l) > 0L) {
      stop("arm_labels has unknown entries: ",
           paste(bad_l, collapse = ", "), call. = FALSE)
    }
    default_labels[names(arm_labels)] <- arm_labels
  }
  arm_colors_v <- default_colors
  arm_labels_v <- default_labels

  # Linetypes per arm (so the curves remain distinguishable in B&W /
  # colorblind viewing). Same merge semantics as arm_colors / arm_labels.
  default_linetypes <- default_arm_linetypes(avail_arms)
  if (!is.null(arm_linetypes)) {
    bad_lt <- setdiff(names(arm_linetypes), names(default_linetypes))
    if (length(bad_lt) > 0L) {
      stop("arm_linetypes has unknown entries: ",
           paste(bad_lt, collapse = ", "), call. = FALSE)
    }
    default_linetypes[names(arm_linetypes)] <- arm_linetypes
  }
  arm_linetypes_v <- default_linetypes

  arms_chr <- as.character(arms)

  # --- Shared x-axis ticks (subset of c(0, cut_times) per `cut_times` arg) ---
  fit_cuts        <- sort(unique(risk_df$time))
  display_cuts    <- .resolve_display_cuts(cut_times, fit_cuts)
  shared_x_limits <- c(0, max(fit_cuts))
  shared_x_ticks  <- display_cuts

  # --- Plot data: origin row per arm so curves start at (t = 0, value = 0)
  # for incidence or (t = 0, value = 1) for survival. ---
  origin_value <- if (x$scale == "incidence") 0 else 1
  body_rows <- risk_df[risk_df$treatment %in% arms,
                       c("time", "treatment", "value"), drop = FALSE]
  origin_rows <- data.frame(
    time      = 0,
    treatment = arms,
    value     = origin_value,
    stringsAsFactors = FALSE
  )
  plot_data <- rbind(origin_rows, body_rows)
  plot_data$arm_chr   <- as.character(plot_data$treatment)
  plot_data$arm_label <- arm_labels_v[plot_data$arm_chr]

  # --- Default y-label per scale ---
  if (is.null(y_label)) {
    y_label <- if (x$scale == "incidence")
      "Cumulative incidence" else "Survival"
  }

  # --- Base plot ---
  p <- ggplot(plot_data, aes(
    x = .data$time, y = .data$value,
    color = .data$arm_label, linetype = .data$arm_label,
    group = .data$arm_label
  )) +
    geom_step(linewidth = linewidth) +
    scale_color_manual(
      values = setNames(arm_colors_v[arms_chr],
                               arm_labels_v[arms_chr])
    ) +
    scale_linetype_manual(
      values = setNames(arm_linetypes_v[arms_chr],
                               arm_labels_v[arms_chr])
    ) +
    scale_x_continuous(
      breaks = shared_x_ticks,
      limits = shared_x_limits,
      expand = expansion(mult = 0.02)
    ) +
    scale_y_continuous(
      labels = function(p) paste0(round(100 * p), "%"),
      limits = c(0, 1),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    labs(
      x        = x_label,
      y        = y_label,
      color    = "Arm",
      linetype = "Arm",
      title    = title,
      subtitle = subtitle
    ) +
    .causal_survival_theme(base_size = base_size)

  # --- CI ribbons (step-transformed to match geom_step) ---
  if (have_bands) {
    body_bands <- risk_df[risk_df$treatment %in% arms,
                          c("time", "treatment", "lower", "upper"),
                          drop = FALSE]
    origin_bands <- data.frame(
      time      = 0,
      treatment = arms,
      lower     = origin_value,
      upper     = origin_value,
      stringsAsFactors = FALSE
    )
    ribbon_data <- rbind(origin_bands, body_bands)

    # Step-transform per arm: duplicate each interior point so the ribbon
    # stays flat between cut times and jumps at each cut. Matches the
    # discrete-time hazard step structure.
    ribbon_data <- do.call(rbind, lapply(arms, function(a) {
      d <- ribbon_data[ribbon_data$treatment == a, , drop = FALSE]
      d <- d[order(d$time), , drop = FALSE]
      n <- nrow(d)
      if (n < 2L) return(d)
      time_step  <- c(d$time[1],
                      rep(d$time[2:n], each = 2L))
      lower_step <- c(rep(d$lower[1:(n - 1L)], each = 2L), d$lower[n])
      upper_step <- c(rep(d$upper[1:(n - 1L)], each = 2L), d$upper[n])
      data.frame(
        time      = time_step,
        lower     = lower_step,
        upper     = upper_step,
        treatment = a,
        stringsAsFactors = FALSE
      )
    }))
    ribbon_data$arm_chr   <- as.character(ribbon_data$treatment)
    ribbon_data$arm_label <- arm_labels_v[ribbon_data$arm_chr]

    p <- p + geom_ribbon(
      data = ribbon_data,
      aes(
        x = .data$time, ymin = .data$lower, ymax = .data$upper,
        fill = .data$arm_label
      ),
      alpha = ribbon_alpha, inherit.aes = FALSE
    ) +
      scale_fill_manual(
        values = setNames(arm_colors_v[arms_chr],
                                 arm_labels_v[arms_chr]),
        guide  = "none"
      )
  }

  # --- Risk table panel(s) (one per requested count, stacked via patchwork) ---
  if (!is.null(risk_table)) {
    valid_counts <- c("at_risk", "events_y", "events_y_interval", "censored")
    if (!is.character(risk_table)) {
      stop("`risk_table` must be NULL or a character vector with entries in: ",
           paste(shQuote(valid_counts), collapse = ", "), ".",
           call. = FALSE)
    }
    bad_rt <- setdiff(risk_table, valid_counts)
    if (length(bad_rt) > 0L) {
      stop("`risk_table` has unknown entries: ",
           paste(shQuote(bad_rt), collapse = ", "),
           ". Valid: ", paste(shQuote(valid_counts), collapse = ", "), ".",
           call. = FALSE)
    }

    n_tbl <- length(risk_table)
    tbl_plots <- lapply(seq_along(risk_table), function(i) {
      build_risk_table_plot(
        pt_data      = x$pt_data,
        id_col       = x$id_col,
        trt_col      = x$treatment_col,
        cut_times    = x$cut_times,
        count        = risk_table[[i]],
        base_size    = base_size,
        x_breaks     = shared_x_ticks,
        x_limits     = shared_x_limits,
        display_cuts = display_cuts,
        # Bottom-most panel gets the cut-time tick labels; the ones
        # above stay clean (curves panel / upper tables carry no
        # x-axis text).
        show_x_axis  = (i == n_tbl)
      )
    })

    if (isTRUE(curves)) {
      # curves on top + risk-table panel(s) below. `risk_table_height` is
      # the per-panel ratio (curves panel = 1); total bottom = N * h.
      p <- wrap_plots(
        c(list(p), tbl_plots), ncol = 1,
        heights = c(1, rep(risk_table_height, length(tbl_plots)))
      )
    } else {
      # Table-only: drop curves panel, stack the table panels equally.
      p <- wrap_plots(
        tbl_plots, ncol = 1,
        heights = rep(1, length(tbl_plots))
      )
    }
  } else if (!isTRUE(curves)) {
    stop("`curves = FALSE` requires a non-NULL `risk_table` ",
         "(nothing else would render).", call. = FALSE)
  }

  p
}


#' Build the risk-table panel for stacking below the curves
#'
#' Renders the output of [risk_table_internal()] as a minimalist
#' ggplot (text-on-tick) suitable for patchwork stacking below a
#' `plot.causal_survival_risk()` curves plot. Aligns its x-axis to
#' the curves panel via the `x_breaks` / `x_limits` arguments.
#'
#' @param pt_data Person-time data.frame.
#' @param id_col,trt_col Character column names.
#' @param cut_times Numeric vector of cut times.
#' @param count One of `"at_risk"`, `"events_y"`, `"censored"`.
#' @param base_size Base font size (inherited from the parent plot).
#' @param x_breaks,x_limits Tick positions / axis limits from the
#'   curves panel, so both panels align under `wrap_plots()`.
#' @return A ggplot2 object.
#' @family internal
#' @keywords internal
build_risk_table_plot <- function(pt_data, id_col, trt_col,
                                  cut_times, count, base_size = 11,
                                  x_breaks = NULL, x_limits = NULL,
                                  display_cuts = NULL,
                                  show_x_axis = FALSE) {
  # Explicit %in% (not match.arg) to avoid partial matching, e.g.
  # "events" -> "events_y" silently.
  valid_counts <- c("at_risk", "events_y", "events_y_interval", "censored")
  if (length(count) != 1L || !count %in% valid_counts) {
    stop("`risk_table` must be one of: ",
         paste(shQuote(valid_counts), collapse = ", "),
         ". Got: ", shQuote(count), ".", call. = FALSE)
  }

  tbl <- risk_table_internal(pt_data, id_col, trt_col, cut_times, count)
  arm_cols <- setdiff(names(tbl), "k")

  # display_cuts is the user-chosen subset of c(0, cut_times) to label.
  # NULL = show all (baseline + every cut_time).
  display_cuts <- display_cuts %||% c(0, cut_times)

  # X-axis breaks: reuse parent's for visual alignment with the curves
  # panel above; otherwise use display_cuts.
  axis_breaks  <- x_breaks  %||% display_cuts
  scale_limits <- x_limits %||% c(0, max(cut_times))

  # Baseline value at t = 0 per count type:
  #   at_risk  -> N per arm (everyone at risk at study start)
  #   events_y -> 0 (no events have occurred yet, by definition)
  #   censored -> 0 (no one censored yet, by definition)
  if (count == "at_risk") {
    n_at_zero <- vapply(arm_cols, function(ac) {
      a <- sub(paste0("^", trt_col, "_"), "", ac)
      length(unique(pt_data[[id_col]][as.character(pt_data[[trt_col]]) == a]))
    }, integer(1), USE.NAMES = FALSE)
  } else {
    n_at_zero <- rep(0L, length(arm_cols))
  }

  # Per-arm long row: for each position in display_cuts, look up the
  # count (baseline at t = 0, else tbl[[ac]][which(tbl$k == pos)]).
  # Out-of-grid values are pre-rejected by .resolve_display_cuts().
  long <- do.call(rbind, lapply(seq_along(arm_cols), function(i) {
    ac <- arm_cols[i]
    vals <- vapply(display_cuts, function(pos) {
      if (pos == 0) return(n_at_zero[i])
      as.integer(tbl[[ac]][which(tbl$k == pos)])
    }, integer(1))
    data.frame(
      k   = display_cuts,
      arm = ac,
      n   = vals,
      stringsAsFactors = FALSE
    )
  }))
  long$arm <- factor(long$arm, levels = rev(arm_cols))

  # hjust per position: left-align at x = 0 (so "100" doesn't get
  # clipped by the panel edge), right-align at the trailing display
  # position, center elsewhere.
  long$hjust <- 0.5
  long$hjust[long$k == 0]                 <- 0
  long$hjust[long$k == max(display_cuts)] <- 1

  title_map <- c(
    at_risk           = "Number at risk",
    events_y          = "Cumulative Y-events",
    events_y_interval = "Y-events per interval",
    censored          = "Number censored"
  )
  title_text <- title_map[count] %||% count

  ggplot(long, aes(x = .data$k, y = .data$arm)) +
    geom_text(
      aes(label = .data$n, hjust = .data$hjust),
      size  = base_size * 0.32
    ) +
    scale_x_continuous(
      breaks = axis_breaks,
      labels = axis_breaks,
      limits = scale_limits,
      expand = expansion(mult = 0.02)
    ) +
    scale_y_discrete(expand = expansion(add = 0.5)) +
    labs(
      x     = NULL,
      y     = trt_col,
      title = title_text
    ) +
    .causal_survival_theme(base_size = base_size) +
    theme(
      # Risk-table panel: drop horizontal gridlines (rows are arm labels,
      # not numeric levels), hide x-axis text by default (curves panel /
      # upper tables carry the ticks), smaller title, black bottom +
      # left borders to frame the count grid.
      panel.grid.major.y  = element_blank(),
      axis.line.x.bottom  = element_line(color     = "black",
                                                  linewidth = 0.4),
      axis.line.y.left    = element_line(color     = "black",
                                                  linewidth = 0.4),
      axis.ticks.x        = if (show_x_axis) {
        element_line(color = "black", linewidth = 0.3)
      } else {
        element_blank()
      },
      axis.text.x         = if (show_x_axis) {
        element_text(size = base_size * 0.85)
      } else {
        element_blank()
      },
      axis.title.y        = element_text(face = "bold",
                                                  angle = 90,
                                                  size  = base_size * 0.85),
      axis.text.y         = element_text(face = "bold",
                                                  size = base_size * 0.85),
      plot.title          = element_text(face = "bold",
                                                  size = base_size * 0.9),
      plot.title.position = "plot"
    )
}


#' Forest plot of counterfactual contrasts
#'
#' Renders the contrast object as a horizontal forest plot, one row per
#' `(method, contrast)` pair, with the point estimate, confidence
#' interval, and a vertical dashed reference line at the null value
#' (0 for difference, 1 for ratio).
#'
#' Two side-by-side panels — RD on a linear x-axis, RR on a log x-axis
#' — so that protective and harmful ratio effects appear equidistant
#' from the null (forest-plot convention since Cochrane meta-analyses).
#' Marker size is proportional to precision (`1 / (upper - lower)`) so
#' tighter CIs receive larger squares; rows with no bootstrap
#' (`lower = NA`) render as a point only.
#'
#' Per the lean-API design: the function expects the contrast object
#' to be method-scoped already (one method's rows; the typical case
#' since `causal_survival()` fits one method at a time). When multiple
#' methods are present the forest still renders one row per method,
#' but the design target is one-method-at-a-time.
#'
#' @param x A `"causal_survival_contrast"` object from
#'   [causal_contrast()].
#' @param base_size Numeric. Base font size for `theme_minimal()`.
#'   Default `11`.
#' @param ... Additional arguments (currently unused).
#' @return A `ggplot`/`patchwork` object.
#' @family plot
#' @export
plot.causal_survival_contrast <- function(x, base_size = 11, ...) {
  ctr <- x$contrasts
  if (is.null(ctr) || nrow(ctr) == 0L) {
    stop("No contrasts to plot — `x$contrasts` is empty.", call. = FALSE)
  }

  # Per-row label, null reference, precision weight.
  ctr$null_val <- ifelse(ctr$contrast == "difference", 0, 1)
  width_ci     <- ctr$upper - ctr$lower
  prec_raw     <- ifelse(is.finite(width_ci) & width_ci > 0,
                         1 / width_ci, NA_real_)
  ctr$prec     <- if (all(is.na(prec_raw))) {
    rep(1, nrow(ctr))
  } else {
    prec_raw / max(prec_raw, na.rm = TRUE)
  }

  # Build one panel per contrast type (RD / RR). Each panel uses its
  # own x-axis: linear for RD, log for RR.
  panels <- list()
  for (ct in c("difference", "ratio")) {
    rows <- ctr[ctr$contrast == ct, , drop = FALSE]
    if (nrow(rows) == 0L) next
    panels[[length(panels) + 1L]] <- .build_contrast_forest_panel(
      rows, ct, base_size
    )
  }
  if (length(panels) == 0L) {
    stop("No RD or RR rows to plot.", call. = FALSE)
  }

  # Caption (scale, alpha, time, identification reminder).
  scale_label <- if (!is.null(x$scale) && x$scale == "survival")
    "survival S^a" else "cumulative incidence F^a"
  alpha_label <- if (!is.null(x$alpha))
    sprintf("%.0f%% CIs", (1 - x$alpha) * 100) else "point estimates only"
  caption <- sprintf("%s  |  scale: %s  |  t = %g",
                     alpha_label, scale_label, x$time)

  wrap_plots(panels, ncol = length(panels)) +
    plot_annotation(
      title   = "Counterfactual contrasts",
      caption = caption,
      theme   = theme(
        plot.title   = element_text(
                         face = "bold", size = base_size * 1.1),
        plot.caption = element_text(
                         size = base_size * 0.78, hjust = 0,
                         color = "grey30")
      )
    )
}


#' Build a single RD or RR forest panel
#'
#' Internal helper for [plot.causal_survival_contrast()]. Renders one
#' contrast type (difference or ratio) as a horizontal forest: dashed
#' null reference line, error bars, point markers sized by precision.
#'
#' @keywords internal
.build_contrast_forest_panel <- function(rows, contrast_type, base_size) {
  null_val <- rows$null_val[1]

  # X-axis limits: include null + small margin on each side.
  vals <- c(rows$lower, rows$upper, rows$estimate, null_val)
  vals <- vals[is.finite(vals)]
  if (contrast_type == "ratio") {
    vals <- vals[vals > 0]
    if (length(vals) == 0L) {
      x_lim <- c(0.5, 2)
    } else {
      log_lo <- log10(min(vals))
      log_hi <- log10(max(vals))
      pad    <- max(0.1, 0.15 * (log_hi - log_lo))
      x_lim  <- 10 ^ c(log_lo - pad, log_hi + pad)
    }
  } else {
    if (length(vals) == 0L) {
      x_lim <- c(-1, 1)
    } else {
      span <- diff(range(vals))
      pad  <- max(0.05, 0.15 * span)
      x_lim <- c(min(vals) - pad, max(vals) + pad)
    }
  }

  axis_label  <- if (contrast_type == "difference") "RD" else "RR (log scale)"
  panel_title <- if (contrast_type == "difference")
    "Risk Difference (RD)" else "Risk Ratio (RR)"

  g <- ggplot(rows, aes(y = .data$method)) +
    geom_vline(xintercept = null_val, linetype = "dashed",
                        color = "grey50", linewidth = 0.4) +
    geom_errorbarh(
      aes(xmin = .data$lower, xmax = .data$upper),
      height = 0.2, linewidth = 0.6, color = "black", na.rm = TRUE
    ) +
    geom_point(
      aes(x = .data$estimate, size = .data$prec),
      shape = 15, color = "black"
    ) +
    scale_size_continuous(range = c(3, 6), guide = "none") +
    scale_y_discrete(limits = rev) +
    labs(x = axis_label, y = NULL, title = panel_title) +
    .causal_survival_theme(base_size = base_size) +
    theme(
      plot.title         = element_text(
                             face = "bold", size = base_size * 0.95),
      panel.grid.major.y = element_blank(),
      axis.title.y       = element_blank()
    )

  if (contrast_type == "ratio") {
    g + scale_x_log10(limits = x_lim)
  } else {
    g + scale_x_continuous(limits = x_lim)
  }
}


#' Plot a causal_survival_diagnostic object (placeholder)
#'
#' Per spec §3.4 line 769: ships as a message-only placeholder in
#' v0.1.0. Implementation deferred to v0.2.
#'
#' @param x A `"causal_survival_diagnostic"` object.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`. Emits an informational `message()`.
#' @export
plot.causal_survival_diagnostic <- function(x, ...) {
  message(
    "plot.causal_survival_diagnostic() is a placeholder in v0.1.0. ",
    "Inspect `x$model_checks` and `x$weight_summary` directly; a ",
    "graphical view ships in v0.2 (see dev/TODO.md)."
  )
  invisible(x)
}


#' Resolve the `cut_times` plot argument to a vector of x-axis positions
#'
#' Accepts:
#' - `NULL`: display all of `c(0, fit_cuts)`.
#' - single positive integer `N`: display `c(0, fit_cuts[idx])` with
#'   `idx <- round(seq(1, length(fit_cuts), length.out = N))`.
#' - numeric vector (length >= 2): explicit subset, must be a subset of
#'   `c(0, fit_cuts)`; `0` is auto-prepended if absent (so the baseline
#'   column is always shown).
#'
#' A length-1 numeric is always interpreted as a count, not a subset.
#' To display a single specific cut time, pair it with 0, e.g.
#' `cut_times = c(0, 5)`.
#'
#' @keywords internal
.resolve_display_cuts <- function(cut_times, fit_cuts) {
  if (is.null(cut_times)) return(c(0, fit_cuts))
  if (!is.numeric(cut_times)) {
    stop("`cut_times` must be NULL, a single positive integer (count), ",
         "or a numeric vector (length >= 2) subset of c(0, fit$cut_times).",
         call. = FALSE)
  }
  if (length(cut_times) == 1L) {
    n <- as.integer(cut_times)
    if (n < 1L || n > length(fit_cuts)) {
      stop("`cut_times` (count) must be between 1 and ",
           length(fit_cuts), " (number of intervals in the fit). ",
           "Got: ", cut_times, ".", call. = FALSE)
    }
    idx <- round(seq(1, length(fit_cuts), length.out = n))
    return(c(0, fit_cuts[idx]))
  }
  bad <- setdiff(cut_times, c(0, fit_cuts))
  if (length(bad) > 0L) {
    stop("`cut_times` (vector) must be a subset of c(0, fit$cut_times). ",
         "Unrecognized values: ", paste(bad, collapse = ", "), ".",
         call. = FALSE)
  }
  sort(unique(c(0, cut_times)))
}
