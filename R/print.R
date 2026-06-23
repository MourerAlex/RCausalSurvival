#' Print a causal_survival_fit object
#'
#' One-screen summary: estimand framing, method (+ IPW engine when
#' relevant), cohort size, time grid, per-arm cumulative incidence at
#' the final cut time, and the warning / model-check tallies.
#'
#' @param x A `"causal_survival_fit"` object.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.causal_survival_fit <- function(x, ...) {
  cat("Counterfactual survival fit (causal_survival_fit)\n")
  cat("-------------------------------------------------\n")
  if (x$method == "km") {
    cat("Estimand: E[Y_k^{a, c = 0}] under MARGINAL exchangeability\n",
        sep = "")
  } else {
    cat("Estimand: E[Y_k^{a, c = 0}] under exchangeability ",
        "conditional on baseline L_0\n", sep = "")
  }
  cat("See `causal_assumptions(fit)` for the full identification block.\n")

  estimator_str <- if (x$method == "ipw") {
    paste0(" (estimator: ", x$ipw_estimator, ")")
  } else if (x$method == "km" && isTRUE(x$ipcw)) {
    " (with IPCW)"
  } else {
    ""
  }
  cat("Method: ", x$method, estimator_str, "\n", sep = "")
  cat("N subjects: ", length(unique(x$pt_data[[x$id_col]])), "\n", sep = "")
  cat("Cut times: ", length(x$cut_times),
      " (T_max = ", max(x$cut_times), ")\n", sep = "")

  # Per-arm cumulative incidence at final cut time
  est <- x$cumulative_incidence[[x$method]]
  if (!is.null(est)) {
    K_max <- max(est$k)
    last  <- est[est$k == K_max, , drop = FALSE]
    cat("\nCumulative incidence at k = ", K_max,
        " (t = ", max(x$cut_times), "):\n", sep = "")
    for (i in seq_len(nrow(last))) {
      cat(sprintf("  a = %s : F^a = %.4f\n",
                  format(last$treatment[i]), last$inc[i]))
    }
  }

  # Model-checks tally
  if (!is.null(x$model_checks)) {
    issues <- 0L
    for (chk in x$model_checks) {
      if (is.null(chk)) next
      if (!isTRUE(chk$converged)) issues <- issues + 1L
      if (length(chk$glm_warnings) > 0L) issues <- issues + 1L
    }
    if (issues > 0L) {
      cat("\nModel checks: ", issues,
          " issue(s) - use `fit$model_checks` to inspect.\n", sep = "")
    }
  }

  # Warnings - spec §3.5 line 537
  if (length(x$warnings) > 0L) {
    cat("\nFit completed with ", length(x$warnings),
        " warning(s) (see fit$warnings).\n", sep = "")
  }

  cat("\nUse causal_risk(), causal_contrast(), bootstrap() to extract components.\n")
  invisible(x)
}


#' Summary of a causal_survival_fit object
#'
#' Per-arm cumulative incidence at the selected cut time, the
#' model-check tally, and (when a bootstrap is supplied) the RD/RR
#' contrast at the same time. The identification block is delegated
#' to [causal_assumptions()] - not inlined here.
#'
#' @param object A `"causal_survival_fit"` object.
#' @param ci Optional. A `"causal_survival_bootstrap"` object from
#'   [bootstrap()]. When supplied, the contrast row at the selected
#'   `time` is printed beneath the per-arm risk block.
#' @param time `NULL` (default, resolves to the final cut time per
#'   spec §3.4) or a numeric scalar resolved on the reporting grid
#'   via [snap_time()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns the cumulative incidence rows at the
#'   selected time.
#' @export
summary.causal_survival_fit <- function(object, ci = NULL,
                                        time = NULL,
                                        scale = c("incidence", "survival"),
                                        ...) {
  if (!is.null(ci)) {
    stopifnot(inherits(ci, "causal_survival_bootstrap"))
  }
  scale <- match.arg(scale)
  if (is.null(object$pt_data)) {
    stop("summary() requires fit$pt_data. Refit with `keep_data = TRUE`.",
         call. = FALSE)
  }
  t_at  <- snap_time(time, object$cut_times)
  k_at  <- which(object$cut_times == t_at)

  # Cohort context (baseline N per arm, cut-grid stats, estimator label)
  trt_col    <- object$treatment_col
  id_col     <- object$id_col
  pt         <- object$pt_data
  arm_levels <- sort(unique(pt[[trt_col]]))
  n_per_arm  <- vapply(arm_levels, function(a) length(unique(
    pt[[id_col]][as.character(pt[[trt_col]]) == as.character(a)])),
    integer(1))
  K_max   <- length(object$cut_times)
  T_max   <- max(object$cut_times)
  est_str <- if (object$method == "ipw") {
    paste0("ipw (", object$ipw_estimator, ")")
  } else if (object$method == "km" && isTRUE(object$ipcw)) {
    "km (with IPCW)"
  } else {
    object$method
  }

  # 1. Banner
  .render_summary_banner()

  # 2. Method / cohort info box
  .render_summary_info_box(est_str, sum(n_per_arm), K_max, T_max)

  # 3. Baseline section (N per arm + visual bars at t = 0)
  .render_summary_baseline(arm_levels, n_per_arm)

  # 4. Counterfactual risk at the snapped time (F^a and S^a side by side)
  est <- object$cumulative_incidence[[object$method]]
  row <- est[est$k == k_at, , drop = FALSE]
  .render_summary_risk(row, t_at)

  # 5. Contrasts (only when bootstrap attached) — table + proportional bar
  if (!is.null(ci)) {
    ctr <- causal_contrast(object, ci = ci, time = t_at,
                           scale = scale)$contrasts
    .render_summary_contrasts(ctr, t_at, ci, scale)
  }

  # 6. Footer pointers + model-checks tally
  .render_summary_footer(ci)
  .render_summary_model_checks(object$model_checks)

  invisible(row)
}


# ----------------------------------------------------------------------------
# Internal helpers for summary.causal_survival_fit(): tile-by-tile renderers
# for the box-drawn console layout (banner, info box, baseline bars, risk
# table, contrasts + proportional CI bar, footer, model-checks tally).
# ----------------------------------------------------------------------------

#' Render the top banner of `summary.causal_survival_fit()`
#' @keywords internal
.render_summary_banner <- function() {
  cat("\n")
  # ASCII-only banner: every char is exactly 1 column wide in any
  # monospace font, so the right edge always aligns with the border.
  title    <- "COUNTERFACTUAL SURVIVAL - SUMMARY"
  inner_w  <- 70L
  side_pad <- (inner_w - nchar(title)) %/% 2L
  middle   <- paste0(strrep(" ", side_pad), title,
                     strrep(" ", inner_w - side_pad - nchar(title)))
  bar      <- strrep("=", inner_w)
  cat("+", bar, "+\n", sep = "")
  cat("|", middle, "|\n", sep = "")
  cat("+", bar, "+\n\n", sep = "")
}

#' Render the method / cohort info box
#' @keywords internal
.render_summary_info_box <- function(est_str, N, K_max, T_max) {
  method_line <- sprintf("Method: %s", est_str)
  stats_line  <- sprintf("N = %d  |  K = %d  |  T_max = %g",
                         N, K_max, T_max)
  # ASCII-only box (every char exactly 1 col wide). Inner width =
  # longest content line + 2 spaces of breathing room each side so the
  # right edge closes flush just past "T_max = 10".
  content_w <- max(nchar(method_line), nchar(stats_line))
  inner_w   <- content_w + 4L
  pad_to    <- function(s) paste0(s, strrep(" ", content_w - nchar(s)))
  indent    <- strrep(" ", 11L)
  bar       <- strrep("-", inner_w)
  cat(indent, "+", bar, "+\n", sep = "")
  cat(indent, "|  ", pad_to(method_line), "  |\n", sep = "")
  cat(indent, "|  ", pad_to(stats_line),  "  |\n", sep = "")
  cat(indent, "+", bar, "+\n\n", sep = "")
}

#' Render the baseline (t = 0) section with N per arm and bar visuals
#' @keywords internal
.render_summary_baseline <- function(arm_levels, n_per_arm) {
  cat("  ", strrep("-", 66), "\n", sep = "")
  cat("  BASELINE   at t = 0\n")
  cat("  ", strrep("-", 66), "\n\n", sep = "")
  bar_max <- 20L
  for (i in seq_along(arm_levels)) {
    bar_len <- max(1L, round(n_per_arm[i] / max(n_per_arm) * bar_max))
    cat(sprintf("      arm %s  %s%s   N = %3d\n",
                format(arm_levels[i]),
                strrep("#", bar_len),
                strrep(" ", bar_max - bar_len),
                n_per_arm[i]))
  }
  cat("\n")
}

#' Render the counterfactual risk table (F^a and S^a side by side)
#' @keywords internal
.render_summary_risk <- function(row, t_at) {
  cat("  ", strrep("-", 66), "\n", sep = "")
  cat(sprintf("  COUNTERFACTUAL RISK   at t = %g\n", t_at))
  cat("  ", strrep("-", 66), "\n\n", sep = "")
  cat("      arm     F^a(t)    S^a(t)\n")
  for (i in seq_len(nrow(row))) {
    cat(sprintf("       %s      %.3f     %.3f\n",
                format(row$treatment[i]), row$inc[i], row$surv[i]))
  }
  cat("\n")
}

#' Render the contrasts section: per-row table + proportional CI bar
#'
#' Each contrast row prints: a header line with `RD/RR`, formula on the
#' chosen scale, and the null reference; an estimate-and-CI line; then a
#' proportional bar where `[` and `]` mark the CI endpoints, `|` marks
#' the null, `*` marks the point estimate (`<>` when null = estimate).
#' The bar's display range spans `[min(lower, null), max(upper, null)]`,
#' so the null is always visible even when it falls outside the CI.
#'
#' @keywords internal
.render_summary_contrasts <- function(ctr, t_at, ci, scale) {
  scale_label <- if (scale == "incidence")
    "cumulative incidence F^a" else "survival S^a"
  cat("  ", strrep("-", 66), "\n", sep = "")
  cat(sprintf("  CONTRASTS on %s (t = %g)\n", scale_label, t_at))
  cat(sprintf("  (%.0f%% CIs from %d bootstrap replicates)\n",
              (1 - ci$alpha) * 100, ci$n_boot_effective))
  cat("  ", strrep("-", 66), "\n\n", sep = "")

  formula_for <- function(contrast_type, scale) {
    base <- if (scale == "incidence") "F" else "S"
    sym  <- if (contrast_type == "difference") "-" else "/"
    sprintf("%s^1 %s %s^0", base, sym, base)
  }

  any_includes_null <- FALSE
  for (i in seq_len(nrow(ctr))) {
    r <- ctr[i, ]
    null_val      <- if (r$contrast == "difference") 0 else 1
    includes_null <- r$lower <= null_val && null_val <= r$upper
    if (includes_null) any_includes_null <- TRUE
    label         <- if (r$contrast == "difference") "RD" else "RR"

    cat(sprintf("      %s   formula: %-12s   null = %g\n",
                label, formula_for(r$contrast, scale), null_val))
    cat(sprintf("           estimate = %7.3f   CI [%7.3f, %7.3f]\n",
                r$estimate, r$lower, r$upper))

    bar_str <- .build_contrast_bar(r$lower, r$upper, r$estimate, null_val,
                                   contrast = r$contrast)
    cat("           ", bar_str, "\n", sep = "")

    annot <- if (includes_null) {
      sprintf("null = %g in CI  -> CI does not exclude null *", null_val)
    } else if (null_val < r$lower) {
      sprintf("null = %g < CI   -> CI excludes null at alpha", null_val)
    } else {
      sprintf("null = %g > CI   -> CI excludes null at alpha", null_val)
    }
    cat("           ", strrep(" ", 30L), annot, "\n\n", sep = "")
  }

  if (any_includes_null) {
    cat("      * inconclusive at the chosen alpha level\n",
        "        (conditional on identifying assumptions -",
        " see causal_assumptions(fit))\n\n", sep = "")
  }
}

#' Build a single proportional contrast bar
#'
#' Returns a string of the form `"<lo>   [....|.*....]   <hi>"` where:
#' `[` / `]` are the CI bounds, `|` is the null reference, `*` is the
#' point estimate, and `-` is filler. When `null = estimate` the marker
#' becomes `X`. Bar width fixed at 14 chars; display range spans
#' `[min(lower, null), max(upper, null)]` so the null is always shown.
#'
#' For `contrast = "ratio"` the bar geometry uses `log(value)` so that
#' RR = 0.5 and RR = 2.0 are equidistant from the null = 1 (forest-plot
#' convention since Cochrane). Printed endpoint labels stay on the
#' linear scale so the displayed numbers remain interpretable.
#'
#' Edge guards: non-finite bounds, non-positive values for ratios, and
#' broken ordering (`lower > estimate` or `estimate > upper` from
#' bootstrap pathologies) bail out to a numbers-only string.
#'
#' @keywords internal
.build_contrast_bar <- function(lower, upper, estimate, null_val,
                                contrast, width = 14L) {

  # Edge guards: non-finite values or broken ordering.
  if (!all(is.finite(c(lower, upper, estimate, null_val))) ||
      lower > estimate || estimate > upper) {
    return(sprintf("%7.3f   (non-finite or ill-ordered)   %7.3f",
                   lower, upper))
  }

  # Ratio contrasts: position markers on log scale so symmetric
  # multiplicative effects map symmetrically across the null.
  if (contrast == "ratio") {
    if (lower <= 0 || estimate <= 0 || upper <= 0 || null_val <= 0) {
      return(sprintf("%7.3f   (non-positive ratio)   %7.3f",
                     lower, upper))
    }
    geom_lower <- log(lower)
    geom_upper <- log(upper)
    geom_est   <- log(estimate)
    geom_null  <- log(null_val)
  } else {
    geom_lower <- lower
    geom_upper <- upper
    geom_est   <- estimate
    geom_null  <- null_val
  }

  display_lo <- min(geom_lower, geom_null)
  display_hi <- max(geom_upper, geom_null)
  span       <- display_hi - display_lo
  if (span <= 0) {
    return(sprintf("%7.3f   (point CI = null)   %7.3f", lower, upper))
  }

  pos <- function(v) {
    max(1L, min(width,
                round((v - display_lo) / span * (width - 1L)) + 1L))
  }
  i_lo   <- pos(geom_lower)
  i_hi   <- pos(geom_upper)
  i_null <- pos(geom_null)
  i_est  <- pos(geom_est)

  bar <- rep("-", width)
  bar[i_lo] <- "["
  bar[i_hi] <- "]"
  if (i_null != i_lo && i_null != i_hi) {
    bar[i_null] <- "|"
  }
  if (i_est != i_lo && i_est != i_hi) {
    bar[i_est] <- if (i_est == i_null) "X" else "*"
  }

  # Endpoint labels: keep linear-scale values so printed numbers stay
  # interpretable as the original quantities (especially RR).
  end_lo <- min(lower, null_val)
  end_hi <- max(upper, null_val)
  sprintf("%7.3f   %s   %7.3f",
          end_lo, paste(bar, collapse = ""), end_hi)
}

#' Render the footer block (identification + bootstrap pointers)
#' @keywords internal
.render_summary_footer <- function(ci) {
  cat("  ", strrep("-", 66), "\n", sep = "")
  cat("           identification -> causal_assumptions(fit)\n")
  if (!is.null(ci)) {
    cat(sprintf(
      "           bootstrap detail -> fit$ci, %d replicates @ %.0f%%\n",
      ci$n_boot_effective, (1 - ci$alpha) * 100))
  } else {
    cat("           for CIs: boot <- bootstrap(fit, n_boot = 500); ",
        "summary(fit, ci = boot)\n", sep = "")
  }
  cat("  ", strrep("-", 66), "\n")
}

#' Render the model-checks tally (warns on non-converged models)
#' @keywords internal
.render_summary_model_checks <- function(model_checks) {
  if (is.null(model_checks)) return(invisible())
  non_converged <- 0L
  for (chk in model_checks) {
    if (is.null(chk)) next
    if (!isTRUE(chk$converged)) non_converged <- non_converged + 1L
  }
  if (non_converged > 0L) {
    cat(sprintf(
      "\n  ! Model checks: %d non-converged model(s). See `fit$model_checks`.\n",
      non_converged))
  }
  invisible()
}


#' Print a causal_survival_risk object
#'
#' Per-method per-arm value at the final cut time on the requested
#' scale; indicates whether bootstrap bands are attached.
#'
#' @param x A `"causal_survival_risk"` object from [causal_risk()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.causal_survival_risk <- function(x, ...) {
  cat("Counterfactual ", x$scale,
      " curves (causal_survival_risk)\n", sep = "")
  cat("------------------------------------------------\n")

  methods_avail <- unique(x$risk$method)
  cat("Methods: ", paste(methods_avail, collapse = ", "), "\n", sep = "")
  cat("Bootstrap CIs: ",
      if (any(!is.na(x$risk$lower))) "yes" else "no", "\n\n", sep = "")

  value_label <- if (x$scale == "incidence") "F^a" else "S^a"

  for (m in methods_avail) {
    sub <- x$risk[x$risk$method == m, , drop = FALSE]
    K_max <- max(sub$k)
    last  <- sub[sub$k == K_max, , drop = FALSE]
    cat(sprintf("[%s] at final time (k = %d, t = %g):\n",
                m, K_max, last$time[1]))
    for (i in seq_len(nrow(last))) {
      band <- if (is.na(last$lower[i])) "" else
        sprintf("  [%.4f, %.4f]", last$lower[i], last$upper[i])
      cat(sprintf("  a = %s : %s = %.4f%s\n",
                  format(last$treatment[i]), value_label,
                  last$value[i], band))
    }
    cat("\n")
  }
  cat("Use plot(causal_risk(fit)) to visualize.\n")
  invisible(x)
}


#' Print a causal_survival_contrast object
#'
#' Renders the contrast table at the selected time, with the method
#' and significance level.
#'
#' @param x A `"causal_survival_contrast"` object from
#'   [causal_contrast()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.causal_survival_contrast <- function(x, ...) {
  cat("Counterfactual contrasts (causal_survival_contrast)\n")
  cat("---------------------------------------------------\n")
  scale_label <- if (x$scale == "incidence")
    "cumulative incidence F^a" else "survival S^a"
  cat("Method: ", x$method, "  |  Scale: ", scale_label, "\n", sep = "")
  if (!is.null(x$alpha)) {
    cat(sprintf("Significance level: %g (%.0f%% CIs)\n\n",
                x$alpha, (1 - x$alpha) * 100))
  } else {
    cat("Significance level: - (no bootstrap supplied)\n\n")
  }

  cat(sprintf("At t = %g:\n", x$time))
  out <- x$contrasts[, c("name", "contrast", "estimate", "lower", "upper"),
                     drop = FALSE]
  # Reference value under the null of no causal effect:
  #   difference -> 0
  #   ratio      -> 1
  out$null_value <- ifelse(out$contrast == "difference", 0, 1)
  out <- out[, c("name", "contrast", "null_value",
                 "estimate", "lower", "upper"),
             drop = FALSE]
  print(out, row.names = FALSE)
  cat("\nReading: if the [lower, upper] interval includes the ",
      "`null_value` (0 for difference, 1 for ratio), the CI does not ",
      "exclude the null at the chosen alpha level (inconclusive). ",
      "If the interval excludes the `null_value`, the CI excludes the ",
      "null at the chosen alpha level. Both conclusions are ",
      "conditional on identifying assumptions - ",
      "see causal_assumptions(fit).\n",
      sep = "")
  invisible(x)
}


#' Print a causal_survival_bootstrap object
#'
#' Replicate count (requested vs effective), significance level,
#' failed-replicate count, and a pointer to the accessor that pairs
#' the bands with a fit.
#'
#' @param x A `"causal_survival_bootstrap"` object from
#'   [bootstrap()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.causal_survival_bootstrap <- function(x, ...) {
  cat("Bootstrap confidence bands (causal_survival_bootstrap)\n")
  cat("------------------------------------------------------\n")
  cat("Replicates requested: ", x$n_boot_requested, "\n", sep = "")
  cat("Replicates effective: ", x$n_boot_effective, "\n", sep = "")
  if (length(x$failed_reps) > 0L) {
    cat("Failed replicates: ", length(x$failed_reps),
        " (see $failed_reps for indices)\n", sep = "")
  }
  if (!is.na(x$warnings_count) && x$warnings_count > 0L) {
    cat("Warnings inside replicates: ", x$warnings_count, "\n", sep = "")
  }
  cat(sprintf("Significance level: %g (%.0f%% CIs)\n",
              x$alpha, (1 - x$alpha) * 100))
  cat("\nUse `causal_contrast(fit, ci = <this>)` for contrast bands,\n")
  cat("or `plot(causal_risk(fit, ci = <this>))` for curve bands.\n")
  invisible(x)
}


#' Confidence intervals for causal_survival_fit (intentionally not provided)
#'
#' This method exists to redirect callers to the supported pattern.
#' Confidence intervals are not stored on the fit object; they live
#' on a separate `"causal_survival_bootstrap"` object that pairs with
#' the fit at accessor time.
#'
#' @param object A `"causal_survival_fit"` object.
#' @param parm,level,... Unused.
#' @return Always errors - use the [bootstrap()] + [causal_contrast()]
#'   pattern instead.
#' @export
confint.causal_survival_fit <- function(object, parm = NULL,
                                        level = 0.95, ...) {
  stop(
    "Confidence intervals are not stored on `fit` in this package. ",
    "Pair the fit with a bootstrap object:\n  ",
    "boot <- bootstrap(fit, n_boot = 500)\n  ",
    "causal_contrast(fit, ci = boot)",
    call. = FALSE
  )
}
