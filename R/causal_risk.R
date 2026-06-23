# ----------------------------------------------------------------------------
# Public accessor: causal_risk() — counterfactual cumulative incidence /
# survival under each arm.
#
# Exposes:
#   - causal_risk()       public S3 accessor on a causal_survival_fit
#   - build_risk_long()   internal long-format pivot used by causal_risk()
# ----------------------------------------------------------------------------


#' Cumulative incidence / survival under each arm
#'
#' Extract the counterfactual cumulative incidence
#' \eqn{F^a(k) = E[Y_k^{a, c = 0}]} (or the survival
#' \eqn{S^a(k) = 1 - F^a(k)}) under each level `a` of the binary
#' treatment, indexed by the reporting grid `fit$cut_times`. Optionally
#' folds bootstrap confidence bands into the same table.
#'
#' @param fit A `"causal_survival_fit"` object from
#'   [causal_survival()].
#' @param scale One of `"incidence"` (default, returns
#'   \eqn{F^a(k)}) or `"survival"` (returns \eqn{S^a(k)}).
#' @param ci Optional. A `"causal_survival_bootstrap"` object from
#'   [bootstrap()]. When provided, the `lower` / `upper` columns of
#'   `$risk` are populated; otherwise they are `NA_real_`.
#'
#' @return An S3 object of class `"causal_survival_risk"` with:
#'   \describe{
#'     \item{risk}{Long-format data.frame with columns `method`,
#'       `treatment`, `k`, `time`, `value`, `lower`, `upper`. One row
#'       per `(method, treatment, k)` triple.}
#'     \item{scale}{The selected scale.}
#'     \item{replicates, alpha}{Carried from `ci` for per-contrast
#'       bands in [plot.causal_survival_risk()] (or `NULL`).}
#'     \item{pt_data, id_col, treatment_col, cut_times}{References to
#'       the fit, used by the plot method's optional risk-table panel.}
#'   }
#'
#' @seealso [causal_contrast()], [plot.causal_survival_risk()]
#' @family accessors
#' @export
causal_risk <- function(fit, scale = c("incidence", "survival"),
                        ci = NULL) {

  # 1. Validate / canonicalize args.
  stopifnot(inherits(fit, "causal_survival_fit"))
  scale <- match.arg(scale)
  if (!is.null(ci)) {
    stopifnot(inherits(ci, "causal_survival_bootstrap"))
  }

  # 2. Pivot fit$cumulative_incidence into the long-format risk table.
  #    One row per (method, treatment, k); folds in bootstrap CI bands
  #    when `ci` is supplied. (See build_risk_long() below.)
  risk <- build_risk_long(fit, ci, scale)

  # 3. Assemble S3 output. Fit references (pt_data, id_col, treatment_col,
  #    cut_times) travel through so plot.causal_survival_risk() can render
  #    the optional risk-table panel; replicates + alpha carry the
  #    bootstrap context for per-contrast bands.
  structure(
    list(
      risk          = risk,
      scale         = scale,
      replicates    = if (!is.null(ci)) ci$replicates else NULL,
      alpha         = if (!is.null(ci)) ci$alpha      else NULL,
      pt_data       = fit$pt_data,
      id_col        = fit$id_col,
      treatment_col = fit$treatment_col,
      cut_times     = fit$cut_times
    ),
    class = "causal_survival_risk"
  )
}


#' Build the long-format `$risk` data.frame
#'
#' Pivot `fit$cumulative_incidence` (per-method long table with
#' columns `treatment`, `k`, `time`, `surv`, `inc`) and optionally
#' the bootstrap bands `ci$bands` (per-method long data.frame with
#' columns `treatment`, `k`, `lower`, `upper`) into the canonical
#' accessor shape used by [causal_risk()].
#'
#' @param fit A `"causal_survival_fit"` object.
#' @param ci A `"causal_survival_bootstrap"` object or `NULL`.
#' @param scale One of `"incidence"` or `"survival"`.
#' @return Long-format data.frame with columns `method`, `treatment`,
#'   `k`, `time`, `value`, `lower`, `upper`. `lower` / `upper` are
#'   `NA_real_` when `ci` is `NULL`.
#' @family internal
#' @keywords internal
build_risk_long <- function(fit, ci, scale) {
  cum_inc_list <- fit$cumulative_incidence
  value_col    <- if (scale == "incidence") "inc" else "surv"

  rows <- list()
  for (m in names(cum_inc_list)) {
    est <- cum_inc_list[[m]]
    if (is.null(est)) next  # method not fit this run

    out <- data.frame(
      method    = m,
      treatment = est$treatment,
      k         = est$k,
      time      = est$time,
      value     = est[[value_col]],
      lower     = NA_real_,
      upper     = NA_real_,
      stringsAsFactors = FALSE,
      row.names = NULL
    )

    # Join in bootstrap bands when present. The bootstrap object stores
    # a single `bands` data.frame with columns (treatment, k, lower,
    # upper) for the method that produced the fit; match by (treatment, k).
    if (!is.null(ci)) {
      key_out  <- paste(out$treatment, out$k, sep = "|")
      key_band <- paste(ci$bands$treatment, ci$bands$k, sep = "|")
      idx <- match(key_out, key_band)
      out$lower <- ci$bands$lower[idx]
      out$upper <- ci$bands$upper[idx]
      # When `scale = "survival"`, the cumulative incidence bands
      # map to survival bands via `S = 1 - F` — flipped order.
      if (scale == "survival") {
        flipped_lower <- 1 - out$upper
        flipped_upper <- 1 - out$lower
        out$lower <- flipped_lower
        out$upper <- flipped_upper
      }
    }
    rows[[length(rows) + 1L]] <- out
  }
  if (length(rows) == 0L) {
    return(data.frame(
      method = character(), treatment = integer(),
      k = integer(), time = numeric(),
      value = numeric(), lower = numeric(), upper = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}
