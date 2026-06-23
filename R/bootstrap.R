#' Bootstrap confidence bands
#'
#' S3 generic. Subject-level percentile bootstrap pairing a fit with
#' confidence bands. CausalSurvival ships the
#' `"causal_survival_fit"` method; CausalCompetingRisks registers a
#' method for its `"causal_competing_risks_fit"` fits, so `bootstrap()`
#' dispatches correctly whichever package masks the other.
#'
#' Each replicate samples unique subject IDs with replacement, stitches
#' the corresponding person-time rows back together with synthetic IDs,
#' and re-evaluates `fit$call` with the resampled `pt_data` (spec
#' line 767). Failed replicates (errors during refit) are tracked in
#' `$failed_reps`; the effective B is `n_boot - length(failed_reps)`.
#'
#' @param fit A fitted object. For this package: a
#'   `"causal_survival_fit"` from [causal_survival()].
#' @param n_boot Positive integer. Number of bootstrap replicates.
#' @param alpha Two-sided significance level in `(0, 1)`.
#' @param seed Optional integer RNG seed for reproducibility.
#' @param ... Passed on to methods.
#'
#' @return An S3 object of class `"causal_survival_bootstrap"` with
#'   the shape documented in spec §4.3.
#'
#' @export
bootstrap <- function(fit, ...) UseMethod("bootstrap")


#' @rdname bootstrap
#' @export
bootstrap.default <- function(fit, ...) {
  stop("No bootstrap() method for class '",
       paste(class(fit), collapse = "/"),
       "'. Expected a causal_survival_fit (CausalSurvival) or a ",
       "causal_competing_risks_fit (CausalCompetingRisks).", call. = FALSE)
}


#' @rdname bootstrap
#' @export
bootstrap.causal_survival_fit <- function(fit, n_boot = 500, alpha = 0.05,
                                          seed = NULL, ...) {

  # 1. Validate / canonicalize args (incl. set.seed, the keep_data = TRUE
  #    refit hint when fit$pt_data is NULL, and cross-checks that the
  #    attached pt_data's cut grid + treatment encoding + treatment column
  #    match the fit's — silent-failure paths if they don't).
  args <- .validate_bootstrap_args(fit, n_boot, alpha, seed)

  # 2. Set up the resampling state: split pt_data by subject id, cache
  #    unique ids and the sample size n.
  state <- bootstrap_init_state(args$pt_data, args$id_col)

  # 3. Replicate loop: report progress, resample subjects, refit fit$call
  #    on the bootstrap sample, record per-replicate cumulative incidence
  #    (or note the failure). Replicate-level warnings (glm non-
  #    convergence, IPCW cliffs, propensity extremes, ...) are counted
  #    here via a calling handler — they are MUFFLED at the source so
  #    the loop output stays clean, and the running total surfaces on
  #    the final object as `warnings_count`.
  reps_long       <- list()
  failed_reps     <- integer()
  n_warnings      <- 0L
  report_progress <- bootstrap_progress_reporter(args$n_boot)
  warning_handler <- function(w) {
    n_warnings <<- n_warnings + 1L
    invokeRestart("muffleWarning")
  }
  for (b in seq_len(args$n_boot)) {
    report_progress(b)
    boot_data <- bootstrap_resample_pt(state, args$pt_data, args$id_col)
    call_b <- args$fit$call
    call_b$pt_data <- boot_data
    res <- tryCatch(
      withCallingHandlers(eval(call_b), warning = warning_handler),
      error = function(e) NULL
    )
    if (is.null(res)) {
      failed_reps <- c(failed_reps, b)
      next
    }
    est <- res$cumulative_incidence[[args$method]]
    reps_long[[length(reps_long) + 1L]] <- data.frame(
      boot_id   = b,
      treatment = est$treatment,
      k         = est$k,
      value     = est$inc,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }

  # 4. Aggregate per-replicate CIFs into the long-format replicates table.
  replicates <- if (length(reps_long) == 0L) {
    data.frame(boot_id = integer(), treatment = numeric(),
               k = integer(), value = numeric(),
               stringsAsFactors = FALSE)
  } else do.call(rbind, reps_long)

  # 5. Compute (alpha/2, 1 - alpha/2) percentile bands per (treatment, k).
  bands <- bootstrap_percentile_bands(replicates, args$alpha)

  # 6. Assemble S3 output.
  structure(
    list(
      fit_call         = args$fit$call,
      n_boot_requested = as.integer(args$n_boot),
      n_boot_effective = as.integer(args$n_boot - length(failed_reps)),
      alpha            = args$alpha,
      replicates       = replicates,
      bands            = bands,
      failed_reps      = failed_reps,
      warnings_count   = as.integer(n_warnings)
    ),
    class = "causal_survival_bootstrap"
  )
}


# ----------------------------------------------------------------------------
# Internal helpers for bootstrap(): validation, state init, progress
# reporter, resampling, percentile bands. Plumbing lives here so the public
# API body stays a 6-tile pipeline.
# ----------------------------------------------------------------------------

#' Validate and canonicalize `bootstrap()` arguments
#'
#' Checks `fit` class, the integer-positivity of `n_boot`, the
#' `(0, 1)` range of `alpha`, sets the RNG seed when supplied, and
#' verifies `fit$pt_data` is present (refit with `keep_data = TRUE`,
#' or attach manually via `to_person_time()` — the error message gives
#' a copy-pasteable reconstruction snippet). Also cross-checks that
#' the attached `pt_data`'s cut grid, `treatment_levels`, and
#' `treatment_col` match the corresponding fields on `fit` — three
#' silent-failure paths (wrong replicate alignment, swapped arms,
#' wrong column dispatched) that aren't caught by R's normal errors.
#' Returns a canonical args list with `fit`, `n_boot`, `alpha`,
#' `pt_data`, `id_col`, and `method`.
#'
#' @keywords internal
.validate_bootstrap_args <- function(fit, n_boot, alpha, seed) {
  stopifnot(inherits(fit, "causal_survival_fit"))
  if (!is.numeric(n_boot) || length(n_boot) != 1L ||
      n_boot < 1 || n_boot != round(n_boot)) {
    stop("`n_boot` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be in (0, 1).", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  if (is.null(fit$pt_data)) {
    stop(
      "`fit$pt_data` is NULL. Either refit with\n",
      "  causal_survival(..., keep_data = TRUE)\n",
      "or attach the person-time data manually before bootstrapping:\n",
      "  fit$pt_data <- to_person_time(\n",
      "      data,\n",
      "      id         = fit$id_col,\n",
      "      treatment  = fit$treatment_col,\n",
      "      covariates = fit$covariates,\n",
      "      cut_points = fit$cut_times[-length(fit$cut_times)],\n",
      "      T_max      = max(fit$cut_times),\n",
      "      ...)\n",
      "(time, status, and — for the per-subject `ipcw` classification ",
      "that drives cond_indep_cens / indep_cens, not the model-level ipcw flag — ",
      "must match the original subject-level data).",
      call. = FALSE
    )
  }

  # Cross-check: attached pt_data must use the same cut grid, treatment
  # encoding, and treatment column the fit was built on. Otherwise the
  # bootstrap silently produces wrong-grid replicates, swapped-arm
  # contrasts, or wrong-column propensity / outcome fits.
  pt_cuts <- attr(fit$pt_data, "cut_times")
  if (is.null(pt_cuts) || !isTRUE(all.equal(pt_cuts, fit$cut_times))) {
    stop(
      "`fit$pt_data` cut grid does not match `fit$cut_times`. ",
      "The attached pt_data was built on a different cut spec. Rebuild:\n",
      "  fit$pt_data <- to_person_time(data, ...,\n",
      "      cut_points = fit$cut_times[-length(fit$cut_times)],\n",
      "      T_max      = max(fit$cut_times))",
      call. = FALSE
    )
  }
  pt_levels <- attr(fit$pt_data, "treatment_levels")
  if (is.null(pt_levels) || !identical(pt_levels, fit$treatment_levels)) {
    stop(
      "`fit$pt_data` treatment_levels do not match the fit's ",
      "treatment_levels (fit: ",
      paste(fit$treatment_levels, collapse = " / "),
      "; attached: ",
      paste(pt_levels %||% "<missing>", collapse = " / "), "). ",
      "The attached pt_data uses a different encoding (e.g., factor ",
      "levels swapped). Verify `data[[fit$treatment_col]]` matches the ",
      "original values / factor level ordering, then rebuild:\n",
      "  fit$pt_data <- to_person_time(data, ...,\n",
      "      treatment  = fit$treatment_col,\n",
      "      cut_points = fit$cut_times[-length(fit$cut_times)],\n",
      "      T_max      = max(fit$cut_times))",
      call. = FALSE
    )
  }
  pt_treatment_col <- attr(fit$pt_data, "treatment_col")
  if (is.null(pt_treatment_col) ||
      !identical(pt_treatment_col, fit$treatment_col)) {
    stop(
      "`fit$pt_data` treatment_col (",
      pt_treatment_col %||% "<missing>",
      ") does not match `fit$treatment_col` (", fit$treatment_col, "). ",
      "Rebuild pt_data with `treatment = fit$treatment_col`.",
      call. = FALSE
    )
  }

  list(
    fit     = fit,
    n_boot  = n_boot,
    alpha   = alpha,
    pt_data = fit$pt_data,
    id_col  = fit$id_col,
    method  = fit$method
  )
}


#' Set up the bootstrap resampling state
#'
#' Splits `pt_data` into a per-subject list (`pt_by_id`) keyed by
#' subject id, returns the vector of unique ids and the sample size
#' `n` (used by `sample()` for each replicate). Pre-computing this
#' outside the loop avoids repeated `split()` calls.
#'
#' @export
#' @keywords internal
bootstrap_init_state <- function(pt_data, id_col) {
  unique_ids <- unique(pt_data[[id_col]])
  list(
    unique_ids = unique_ids,
    pt_by_id   = split(pt_data, pt_data[[id_col]]),
    n          = length(unique_ids)
  )
}


#' Build a stateful bootstrap progress reporter
#'
#' Returns a closure that, when called with the current replicate
#' index `b`, emits the cadenced progress messages (first replicate
#' announce, every 10 of the first 50, a one-shot ETA at `b = 50`
#' when `n_boot > 50`, then every 100 thereafter). The `tic` clock
#' and the "estimate already announced" latch live inside the closure
#' so the public [bootstrap()] body doesn't need to manage them.
#'
#' @export
#' @keywords internal
bootstrap_progress_reporter <- function(n_boot) {
  tic <- proc.time()[["elapsed"]]
  estimate_announced <- FALSE
  fmt_duration <- function(sec) {
    m <- as.integer(floor(sec / 60))
    s <- as.integer(round(sec - 60 * m))
    if (m > 0) sprintf("%d min %02d sec", m, s) else sprintf("%d sec", s)
  }
  function(b) {
    if (b == 1 || (b <= 50 && b %% 10 == 0)) {
      message("Bootstrap replicate ", b, "/", n_boot)
    }
    if (!estimate_announced && b == 50 && n_boot > 50) {
      elapsed   <- proc.time()[["elapsed"]] - tic
      per_rep   <- elapsed / 50
      remaining <- (n_boot - 50) * per_rep
      message(sprintf(
        "Bootstrap: 50 replicates done in %s. Estimated remaining: %s (%d more replicates).",
        fmt_duration(elapsed), fmt_duration(remaining), n_boot - 50
      ))
      estimate_announced <<- TRUE
    }
    if (b > 50 && b %% 100 == 0) {
      message("Bootstrap replicate ", b, "/", n_boot)
    }
  }
}


#' Draw a single bootstrap sample of person-time rows
#'
#' Samples `n` unique subject ids with replacement, rebuilds the
#' corresponding person-time rows with synthetic ids (so refit code
#' sees independent subjects even after duplicate sampling), and
#' restores the `person_time` class and attributes that `rbind` would
#' otherwise drop.
#'
#' @export
#' @keywords internal
bootstrap_resample_pt <- function(state, pt_data, id_col) {
  sampled <- sample(state$unique_ids, size = state$n, replace = TRUE)
  boot_data <- do.call(rbind, lapply(seq_along(sampled), function(i) {
    rows <- state$pt_by_id[[as.character(sampled[i])]]
    rows[[id_col]] <- paste0(sampled[i], "_", i)
    rows
  }))
  # Preserve class + attributes lost by rbind
  class(boot_data) <- class(pt_data)
  for (a in names(attributes(pt_data))) {
    if (a %in% c("names", "row.names", "class")) next
    attr(boot_data, a) <- attr(pt_data, a)
  }
  boot_data
}


#' Percentile confidence bands from bootstrap replicates
#'
#' Groups the long-format `replicates` table by `by` and, within each
#' group, takes the lower and upper percentiles of the `value` column.
#' Returns a single long data.frame: the `by` columns plus `lower` and
#' `upper`.
#'
#' @param replicates Long data.frame with a `value` column and the
#'   grouping columns named in `by`.
#' @param alpha Two-sided significance level in `(0, 1)`; bands are the
#'   `alpha/2` and `1 - alpha/2` percentiles.
#' @param by Grouping column names. Default `c("treatment", "k")` (CS);
#'   CCR passes `c("method", "arm", "k")`.
#' @return Data.frame: the `by` columns plus `lower`, `upper`. Zero rows
#'   (correct columns) when `replicates` is empty.
#' @export
#' @keywords internal
bootstrap_percentile_bands <- function(replicates, alpha,
                                       by = c("treatment", "k")) {

  lower_prob <- alpha / 2
  upper_prob <- 1 - alpha / 2

  # No replicates survived (e.g. every bootstrap fit failed): return an
  # empty table that still carries the by + lower/upper columns.
  if (nrow(replicates) == 0L) {
    empty <- replicates[, by, drop = FALSE]
    empty$lower <- numeric(0)
    empty$upper <- numeric(0)
    return(empty)
  }

  # Group replicate rows by the `by` columns, then within each group take
  # the lower / upper percentile of `value`.
  groups <- as.list(replicates[, by, drop = FALSE])
  lower  <- aggregate(replicates["value"], by = groups,
                      FUN = function(v) unname(quantile(v, lower_prob, na.rm = TRUE)))
  upper  <- aggregate(replicates["value"], by = groups,
                      FUN = function(v) unname(quantile(v, upper_prob, na.rm = TRUE)))

  # Both share the `by` columns in identical order; bind the value cols.
  bands <- lower
  names(bands)[names(bands) == "value"] <- "lower"
  bands$upper <- upper$value
  bands
}
