#' Null-Coalescing Operator
#'
#' Returns `x` if not NULL, otherwise `y`.
#'
#' @param x,y Values to coalesce.
#' @return `x` if not NULL, else `y`.
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


#' Clone Baseline Across Interval Indices
#'
#' Broadcast each subject's baseline row across all `K_max` interval
#' indices, setting the treatment column to a fixed value `a`. Used by
#' the g-formula and IPW-MSM workers to predict counterfactual hazards
#' at every (subject, k) regardless of the subject's observed
#' event/censoring time.
#'
#' The `k` column on the clone holds the integer interval index
#' (`1..K_max`), matching the fit-time encoding of the Y-hazard / Y-MSM
#' model. The time grid `cut_times` is preserved on the calling side
#' (e.g. via `attr(pt_data, "cut_times")`) for report-time alignment.
#'
#' @param baseline data.frame. One row per subject (typically
#'   `pt_data[pt_data$k == 1, ]`).
#' @param cut_times Numeric vector of interval-end times `t_1, ..., T_max`.
#'   Used here only for its length `K_max`.
#' @param treatment_col Character. Treatment column name.
#' @param a Numeric (0 or 1). Counterfactual treatment value.
#' @return data.frame with `nrow(baseline) * K_max` rows.
#' @export
#' @keywords internal
make_clone <- function(baseline, cut_times, treatment_col, a) {
  n <- nrow(baseline)
  K <- length(cut_times)
  clone <- baseline[rep(seq_len(n), each = K), , drop = FALSE]
  clone$k                <- rep(seq_len(K), times = n)
  clone[[treatment_col]] <- a
  rownames(clone) <- NULL
  clone
}


#' Build the long-format estimates data.frame
#'
#' Reshape a list of per-arm CIF vectors (`cif_by_arm[[1]]` for `a = 0`,
#' `cif_by_arm[[2]]` for `a = 1`) into the canonical long-format
#' data.frame used by all `causal_survival_fit$cumulative_incidence`
#' slots: `2 * K_max` rows, columns `treatment`, `k`, `time`, `surv`,
#' `inc`. Shared by `fit_gformula()`, `fit_ipw_msm()`, and `fit_ipw_km()`.
#'
#' @param cif_by_arm List of length 2 with per-arm CIF vectors of length
#'   `K_max`.
#' @param cut_times Numeric vector of interval-end times.
#' @return data.frame with `2 * length(cut_times)` rows.
#' @keywords internal
make_estimates_long <- function(cif_by_arm, cut_times) {
  K_max <- length(cut_times)
  data.frame(
    treatment = rep(c(0, 1), each = K_max),
    k         = rep(seq_len(K_max), times = 2),
    time      = rep(cut_times, times = 2),
    surv      = c(1 - cif_by_arm[[1]], 1 - cif_by_arm[[2]]),
    inc       = c(    cif_by_arm[[1]],     cif_by_arm[[2]])
  )
}


#' Resolve a user-supplied time to the reporting grid
#'
#' Used by [causal_contrast()] and [summary.causal_survival_fit()] to
#' interpret the optional `time` argument. `NULL` resolves to the
#' final cut time. A numeric value outside `[min(cut_times),
#' max(cut_times)]` is rejected with a hard error (spec §3.5). An
#' in-bounds value is snapped to the nearest entry of `cut_times`;
#' a `message()` is emitted if snapping changes it.
#'
#' @param time `NULL` or a numeric scalar.
#' @param cut_times Numeric vector of available cut times.
#' @return A single numeric value drawn from `cut_times`.
#' @family internal
#' @export
#' @keywords internal
snap_time <- function(time, cut_times) {
  if (is.null(time)) return(max(cut_times))
  if (!is.numeric(time) || length(time) != 1L || is.na(time)) {
    stop("`time` must be NULL or a single non-missing numeric.",
         call. = FALSE)
  }
  if (time > max(cut_times) || time < min(cut_times)) {
    stop(sprintf(
      "`time = %g` is outside the reporting grid [%g, %g].",
      time, min(cut_times), max(cut_times)
    ), call. = FALSE)
  }
  idx <- which.min(abs(cut_times - time))
  k_at <- cut_times[idx]
  if (!isTRUE(all.equal(k_at, time))) {
    message(sprintf("`time = %g` snapped to nearest cut time: %g.",
                    time, k_at))
  }
  k_at
}
