# ----------------------------------------------------------------------------
# Public accessor: causal_contrast() — counterfactual contrasts on the
# reporting grid (risk difference "-" and risk ratio "/" between treatment
# arms; binary point treatment in v0.1.0).
#
# Exposes:
#   - causal_contrast()         public S3 accessor on a causal_survival_fit
#   - compute_contrast_table()  internal long-format contrast builder;
#                                also reshapes bootstrap replicates into
#                                per-replicate matrices for CI bands
#   - contrast_op()             internal tiny difference / ratio dispatcher
#
# Uses snap_time() from R/utils.R (shared with summary.causal_survival_fit
# in R/print.R).
# ----------------------------------------------------------------------------


#' Counterfactual contrasts on the reporting grid
#'
#' Compute pairwise contrasts of the counterfactual cumulative
#' incidence between treatment arms on the reporting grid. The risk
#' difference `RD(k) = F^{a}(k) - F^{a'}(k)` (operator `"-"`) and the
#' risk ratio `RR(k) = F^{a}(k) / F^{a'}(k)` (operator `"/"`) are the
#' two default operators.
#'
#' Under the v0.1.0 binary point-treatment scope there is a single
#' arm pair `(max, min)` of `fit$treatment_levels`. The `reference` /
#' `contrasts` arguments are carried in the signature for forward
#' compatibility with multi-arm (v2) and separable-effects (downstream
#' package) scopes; in v0.1.0 supplying a non-`NULL` `contrasts`
#' raises an error.
#'
#' @param fit A `"causal_survival_fit"` object from
#'   [causal_survival()].
#' @param reference Optional reference arm (one of the levels in
#'   `fit$treatment_levels`). `NULL` defaults to the lowest level.
#' @param contrasts Optional named list of custom contrasts (deferred
#'   to v2 multi-arm; see `dev/unused_code/resolve_contrast_pairs_v2.R`).
#'   Must be `NULL` in v0.1.0.
#' @param ci Optional. A `"causal_survival_bootstrap"` object from
#'   [bootstrap()]. When `NULL`, point estimates only are returned
#'   with `lower` / `upper` set to `NA_real_`. A loud `warning()` is
#'   emitted in that case: contrasts without confidence intervals are
#'   not meaningful and their interpretation is strongly discouraged.
#' @param time `NULL` (default, resolves to the final cut time per
#'   spec §3.4) or a numeric scalar resolved on the reporting grid
#'   via [snap_time()]. Out-of-bounds values raise an error.
#' @param scale One of `"incidence"` (default — RD and RR computed on
#'   `F^a(k)`) or `"survival"` (RD and RR computed on
#'   `S^a(k) = 1 - F^a(k)`). The null reference values (`0` for RD,
#'   `1` for RR) are the same on both scales; only the magnitudes
#'   (and the sign of RD) change. RD on survival is the sign-flipped
#'   RD on incidence; RR on survival is a different ratio entirely.
#'
#' @return An S3 object of class `"causal_survival_contrast"` with:
#'   \describe{
#'     \item{contrasts}{Long-format data.frame with columns `method`,
#'       `name`, `treatment_a`, `treatment_b`, `contrast`, `k`, `time`,
#'       `estimate`, `lower`, `upper`.}
#'     \item{method}{The fit's estimation method.}
#'     \item{alpha}{Bootstrap significance level (or `NULL`).}
#'     \item{time}{The resolved cut time.}
#'     \item{scale}{Either `"incidence"` or `"survival"`.}
#'   }
#'
#' @seealso [causal_risk()], [plot.causal_survival_contrast()],
#'   [bootstrap()]
#' @family accessors
#' @export
causal_contrast <- function(fit, reference = NULL, contrasts = NULL,
                            ci = NULL, time = NULL,
                            scale = c("incidence", "survival")) {

  # 1. Validate / canonicalize args. (See .validate_causal_contrast_args()
  #    below — includes the loud `ci = NULL` warning, the v0.1.0
  #    contrasts-must-be-NULL guard, the binary-treatment check, the
  #    reference-arm default, and the scale resolution.)
  args <- .validate_causal_contrast_args(fit, reference, contrasts, ci, scale)

  # 2. Build the binary arm-pair × contrast list. Under v0.1.0 there is
  #    a single comparison: comparator vs reference, on both "difference"
  #    (risk difference) and "ratio" (risk ratio).
  comparator <- setdiff(args$levels_vec, args$reference)
  pairs <- list(
    name        = rep(paste0(comparator, "_vs_", args$reference), 2L),
    treatment_a = rep(comparator, 2L),
    treatment_b = rep(args$reference, 2L),
    contrast    = c("difference", "ratio")
  )

  # 3. Compute the long-format contrast table across methods on the
  #    requested scale. Bootstrap bands are produced per replicate
  #    inside compute_contrast_table() (transformed to the requested
  #    scale before quantile aggregation) and reduced to
  #    alpha/2 + 1-alpha/2 quantiles. (See below.)
  out_df <- compute_contrast_table(args$fit, pairs, args$ci, args$scale)

  # 4. Snap user-supplied `time` to the reporting grid; filter to that
  #    single k. `time = NULL` resolves to the final cut time per spec §3.4.
  t_at <- snap_time(time, args$fit$cut_times)
  out_df <- out_df[out_df$time == t_at, , drop = FALSE]

  # 5. Assemble S3 output.
  structure(
    list(
      contrasts = out_df,
      method    = args$fit$method,
      alpha     = if (!is.null(args$ci)) args$ci$alpha else NULL,
      time      = t_at,
      scale     = args$scale
    ),
    class = "causal_survival_contrast"
  )
}


# ----------------------------------------------------------------------------
# Internal helpers for causal_contrast(): validation. Plumbing lives here
# so the public API body stays a 5-tile pipeline.
# ----------------------------------------------------------------------------

#' Validate and canonicalize `causal_contrast()` arguments
#'
#' Runs the four validation steps the public function used to do
#' inline: fit class check; `ci` handling (emits the loud "contrasts
#' without uncertainty are not meaningful" warning when `ci = NULL`,
#' otherwise checks its class); v0.1.0 guard on the not-yet-supported
#' `contrasts` custom list; binary-treatment check + reference-arm
#' default (`NULL` resolves to `min(treatment_levels)`). Returns a
#' canonical args list with `fit`, `ci`, the resolved `reference`,
#' and `levels_vec` (used by the public function to derive the
#' comparator arm).
#'
#' @keywords internal
.validate_causal_contrast_args <- function(fit, reference, contrasts, ci,
                                           scale) {
  stopifnot(inherits(fit, "causal_survival_fit"))
  scale <- match.arg(scale, c("incidence", "survival"))

  if (is.null(ci)) {
    warning(
      "`ci = NULL`: causal_contrast() is returning point estimates ",
      "without confidence intervals. Contrasts without uncertainty ",
      "are not meaningful — their interpretation is strongly ",
      "discouraged. Compute a bootstrap first:\n  ",
      "boot <- bootstrap(fit, n_boot = 500)\n  ",
      "Then: causal_contrast(fit, ci = boot)",
      call. = FALSE
    )
  } else {
    stopifnot(inherits(ci, "causal_survival_bootstrap"))
  }

  if (!is.null(contrasts)) {
    stop(
      "`contrasts` (custom list) is not supported in v0.1.0 ",
      "(binary point treatment only). The drafted helper is parked ",
      "at dev/unused_code/resolve_contrast_pairs_v2.R; see TODO ",
      "'Multi-arm treatment'.",
      call. = FALSE
    )
  }

  levels_vec <- fit$treatment_levels
  if (length(levels_vec) != 2L) {
    stop("v0.1.0 supports binary treatment only; ",
         "fit$treatment_levels has length ", length(levels_vec), ".",
         call. = FALSE)
  }
  if (is.null(reference)) reference <- min(levels_vec)
  if (!reference %in% levels_vec) {
    stop("`reference` must be one of the fit's treatment levels: ",
         paste(levels_vec, collapse = ", "), ".", call. = FALSE)
  }

  list(fit = fit, ci = ci, reference = reference,
       levels_vec = levels_vec, scale = scale)
}


#' Build the long-format `$contrasts` data.frame
#'
#' Iterate over the methods present in `fit$cumulative_incidence` and
#' the arm pairs in `pairs`, computing the per-`k` contrast under the
#' requested operator. Bootstrap bands (when supplied) are obtained
#' by computing the contrast per replicate and taking the requested
#' lower/upper quantile.
#'
#' @param fit A `"causal_survival_fit"` object.
#' @param pairs A list with parallel vectors `name`, `treatment_a`,
#'   `treatment_b`, `contrast` — one element per (pair x contrast) cell.
#' @param ci A `"causal_survival_bootstrap"` object or `NULL`.
#' @return Long-format data.frame with columns `method`, `name`,
#'   `treatment_a`, `treatment_b`, `contrast`, `k`, `time`, `estimate`,
#'   `lower`, `upper`. `lower` / `upper` are `NA_real_` when `ci` is
#'   `NULL`.
#' @family internal
#' @keywords internal
compute_contrast_table <- function(fit, pairs, ci, scale) {
  alpha <- if (!is.null(ci)) ci$alpha else NA_real_
  reps  <- if (!is.null(ci)) ci$replicates else NULL

  out <- list()
  for (m in names(fit$cumulative_incidence)) {
    est <- fit$cumulative_incidence[[m]]
    if (is.null(est)) next  # method not fit this run

    for (j in seq_along(pairs$name)) {
      a             <- pairs$treatment_a[j]
      b             <- pairs$treatment_b[j]
      contrast_type <- pairs$contrast[j]
      nm            <- pairs$name[j]

      # Pull arm-level values on the requested scale. Bootstrap replicates
      # are always stored on the incidence scale (per spec §4.3); we
      # transform to survival here when requested.
      F_a <- est$inc[est$treatment == a]
      F_b <- est$inc[est$treatment == b]
      val_a <- if (scale == "survival") 1 - F_a else F_a
      val_b <- if (scale == "survival") 1 - F_b else F_b

      estimate <- contrast_op(val_a, val_b, contrast_type)

      if (!is.null(reps)) {
        # `reps` is long-format `data.frame(boot_id, treatment, k,
        # value)` per spec §4.3. Reshape to a `[k, boot_id]` matrix
        # per arm, transform to the requested scale if needed, compute
        # the contrast per (k, boot_id), then take quantiles along the
        # boot_id axis.
        wide_a <- tapply(
          reps$value[reps$treatment == a],
          list(reps$k[reps$treatment == a],
               reps$boot_id[reps$treatment == a]),
          identity
        )
        wide_b <- tapply(
          reps$value[reps$treatment == b],
          list(reps$k[reps$treatment == b],
               reps$boot_id[reps$treatment == b]),
          identity
        )
        if (scale == "survival") {
          wide_a <- 1 - wide_a
          wide_b <- 1 - wide_b
        }
        per_rep <- contrast_op(wide_a, wide_b, contrast_type)
        lower   <- apply(per_rep, 1L, quantile,
                         probs = alpha / 2,       na.rm = TRUE)
        upper   <- apply(per_rep, 1L, quantile,
                         probs = 1 - alpha / 2,   na.rm = TRUE)
      } else {
        lower <- rep(NA_real_, length(estimate))
        upper <- rep(NA_real_, length(estimate))
      }

      out[[length(out) + 1L]] <- data.frame(
        method      = m,
        name        = nm,
        treatment_a = a,
        treatment_b = b,
        contrast    = contrast_type,
        k           = est$k[est$treatment == a],
        time        = est$time[est$treatment == a],
        estimate    = estimate,
        lower       = lower,
        upper       = upper,
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  }
  if (length(out) == 0L) {
    return(data.frame(
      method = character(), name = character(),
      treatment_a = numeric(), treatment_b = numeric(),
      contrast = character(), k = integer(), time = numeric(),
      estimate = numeric(), lower = numeric(), upper = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, out)
}


#' Tiny operator dispatcher used by [compute_contrast_table()]
#'
#' Applies the contrast operator (`"-"` or `"/"`) to either a pair of
#' numeric vectors (point estimate) or a pair of matrices indexed by
#' `(k, boot_id)` (per-replicate bands). Kept as a single function so
#' the two shapes share one implementation.
#'
#' @param x,y Numeric vector or matrix.
#' @param contrast_type One of `"difference"` or `"ratio"`.
#' @return Numeric of the same shape as `x` / `y`.
#' @family internal
#' @keywords internal
contrast_op <- function(x, y, contrast_type) {
  switch(contrast_type,
    difference = x - y,
    ratio      = x / y,
    stop("Unknown contrast type: ", contrast_type, call. = FALSE)
  )
}
