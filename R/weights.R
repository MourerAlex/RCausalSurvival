#' Inverse-Probability Weights (IPW)
#'
#' Helpers for constructing IPW weights from fitted models.
#'
#' * [ipw()]            — generic core: per-row inverse of cumulative
#'                        probability of the observed history (cumprod
#'                        baked in). Used for time-varying mechanisms.
#' * [ipw_cens()]       — censoring weights (IPCW), routes through
#'                        [ipw()].
#' * [ipw_static_trt()] — point-treatment weights (per-subject scalar
#'                        broadcast across the subject's rows).
#'
#' **Output convention (read before adding new helpers):** every
#' function in this file returns a *row-ready* numeric vector of
#' length `nrow(pt_data)`, suitable for direct use as `glm(weights = .)`.
#' Combine multiple weight vectors by row-wise multiplication, NOT
#' cumprod — cumulative logic is already baked into the producer.
#'
#' Stabilization is structural (two-cumprod ratio, Robins/Hernán form).
#' Truncation lives downstream in [apply_weight_truncation()] and is
#' applied to the final combined IPW weight (Cole & Hernán 2008 AJE),
#' not to each component.
#'
#' @keywords internal
#' @name weights
NULL


# Internal: validate a probability vector. Errors on NA, out-of-(0,1],
# or non-numeric input. Used by ipw() for the denominator and (if
# supplied) the numerator probability vector.
check_prob_vec <- function(x) {
  stopifnot(
    is.numeric(x),
    !any(is.na(x)),
    all(x > 0 & x <= 1)
  )
  invisible(x)
}


#' Inverse-Probability Weight (Generic Core)
#'
#' Per-row inverse of the cumulative probability of the observed history:
#'
#'   W_{i,k} = prod_{j <= k} prob_num_{ij} / prod_{j <= k} prob_denom_{ij}
#'
#' With `prob_num = NULL` the numerator collapses to 1 (unstabilized form).
#' Caller is responsible for constructing the per-row probabilities; this
#' function is agnostic to which causal mechanism is being weighted.
#'
#' Cole SR, Hernán MA (2008) AJE; Hernán & Robins, *Causal Inference: What
#' If* (ch. 12, 17).
#'
#' @param prob_denom Numeric vector in (0, 1]. P(observed history | full
#'   conditioning), one entry per row of the person-time data. Zeros
#'   indicate positivity violation — diagnose upstream.
#' @param id Long-format subject id, **same length as `prob_denom`** (one
#'   entry per person-time row, repeated within subject). Used to group
#'   the cumprod within subject.
#' @param prob_num Optional numeric vector in (0, 1], same length as
#'   `prob_denom`, or NULL. P(observed history | reduced conditioning).
#'   NULL = unstabilized.
#'
#' @return Numeric weight vector, length `length(prob_denom)`.
#' @export
#' @keywords internal
ipw <- function(prob_denom, id, prob_num = NULL) {
  check_prob_vec(prob_denom)
  stopifnot(length(prob_denom) == length(id))
  cum_denom <- ave(prob_denom, id, FUN = cumprod)
  cum_num   <- if (is.null(prob_num)) 1 else {
    check_prob_vec(prob_num)
    stopifnot(length(prob_num) == length(prob_denom))
    ave(prob_num, id, FUN = cumprod)
  }
  cum_num / cum_denom
}


#' Inverse-Probability Weight for Static (Point) Treatment
#'
#' IPT weights at baseline only.
#'
#' Predicts P(A = 1 | L) on baseline rows only (efficient: n_subjects
#' predictions, not n_rows), picks P(A_obs | L) per subject, inverts
#' (or takes the stabilized ratio with a numerator model), then
#' broadcasts to the long-format pt_data via subject id.
#'
#' Stabilization (Robins/Hernán form): provide `model_num` fit with
#' reduced conditioning, e.g. `A ~ 1` (marginal stabilization).
#'
#' @param model_full Fitted glm of A on full baseline conditioning set.
#' @param pt_data Person-time data frame.
#' @param treatment_col Character. Treatment column name. Values must
#'   be in `{0, 1}` (standardized upstream in [to_person_time()]).
#' @param id_col Character. Subject id column name.
#' @param time_col Character. Interval index column name (default "k").
#' @param model_num Optional glm of A with reduced conditioning, for
#'   stabilization. NULL = unstabilized.
#'
#' @return Numeric weight vector, length `nrow(pt_data)`. Constant
#'   within subject (broadcast from baseline).
#' @export
#' @keywords internal
ipw_static_trt <- function(model_full, pt_data, treatment_col, id_col,
                           time_col = "k", model_num = NULL) {
  baseline_idx <- pt_data[[time_col]] == 1
  baseline     <- pt_data[baseline_idx, ]

  p_full     <- predict(model_full, newdata = baseline, type = "response")
  p_obs_full <- ifelse(baseline[[treatment_col]] == 1, p_full, 1 - p_full)

  w_per_subj <- if (is.null(model_num)) {
    1 / p_obs_full
  } else {
    p_n     <- predict(model_num, newdata = baseline, type = "response")
    p_obs_n <- ifelse(baseline[[treatment_col]] == 1, p_n, 1 - p_n)
    p_obs_n / p_obs_full
  }

  # Broadcast per-subject scalar to long-format via id match
  w_per_subj[match(pt_data[[id_col]], baseline[[id_col]])]
}


#' Inverse Probability of Censoring Weights (IPCW)
#'
#' Thin wrapper around [ipw()] for censoring weights. Predicts the
#' censoring hazard on observed covariates, builds `prob = 1 - hazard`
#' per row, and delegates.
#'
#' Stabilization (Hernán & Robins §12.6, Technical Point 12.2):
#' `SW^C = Pr[C=0 | A] / Pr[C=0 | L, A]`. Supply `model_num` fit with
#' reduced conditioning (typically `c ~ A`) for the stabilized form.
#' `model_num = NULL` gives the unstabilized weight.
#'
#' @param model_full Fitted hazard glm with full conditioning.
#' @param pt_data Person-time data frame.
#' @param id_col Character. Subject id column name.
#' @param model_num Optional fitted hazard glm with reduced conditioning
#'   for stabilization. NULL = unstabilized.
#'
#' @return Numeric weight vector, length `nrow(pt_data)`.
#' @export
#' @keywords internal
ipw_cens <- function(model_full, pt_data, id_col, model_num = NULL) {
  haz_full   <- predict(model_full, newdata = pt_data, type = "response")
  prob_denom <- 1 - haz_full
  prob_num   <- if (is.null(model_num)) NULL else {
    1 - predict(model_num, newdata = pt_data, type = "response")
  }
  ipw(prob_denom, pt_data[[id_col]], prob_num)
}


#' Apply Weight Truncation
#'
#' Symmetric percentile truncation (Cole & Hernán 2008 AJE) applied to
#' the IPW weight columns on `pt_data`. Operates on `w_cens` and `w_a`
#' (whichever are present).
#'
#' Flagged rows are recorded BEFORE clipping so the log carries the raw
#' values. A single warning summarizes how many subjects and rows were
#' affected.
#'
#' @param pt_data Person-time data frame with raw weights attached.
#' @param id_col Character. Subject id column name.
#' @param truncate Either NULL (no truncation) or a length-2 numeric
#'   vector of percentile bounds, e.g. `c(0.001, 0.999)`. Setting the
#'   lower bound to 0 reduces to upper-tail-only truncation.
#'
#' @return A list with:
#'   \describe{
#'     \item{pt_data}{Adjusted person-time data with weights clipped.}
#'     \item{flagged_ids}{Vector of unique subject IDs whose weights
#'       exceeded the bounds across any weight column.}
#'     \item{flagged_log}{Long-format data.frame with one row per
#'       flagged `(subject, interval, weight column)` triple.}
#'   }
#' @export
#' @keywords internal
apply_weight_truncation <- function(pt_data, id_col, truncate = NULL,
                                    w_cols = c("w_cens", "w_a")) {

  # 1. Build the empty flagged_log scaffold used in both the no-op and
  #    no-flag-fired return paths.
  empty_log <- data.frame(
    id     = pt_data[[id_col]][0],
    weight = character(0),
    k      = numeric(0),
    value  = numeric(0),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # 2. Early return when truncation is not requested.
  if (is.null(truncate)) {
    return(list(
      pt_data     = pt_data,
      flagged_ids = integer(0),
      flagged_log = empty_log
    ))
  }

  # 3. Per-column truncation loop. For each present weight column:
  #    compute the percentile bounds, find rows whose weight is outside,
  #    log them BEFORE clipping (so the log carries raw extreme values),
  #    then symmetric-clip (Cole & Hernán 2008).
  w_cols <- intersect(w_cols, names(pt_data))
  log_rows <- list()
  for (wc in w_cols) {
    qs <- quantile(pt_data[[wc]], probs = truncate, na.rm = TRUE)
    lower_val <- qs[1]
    upper_val <- qs[2]

    extreme_rows <- which(
      pt_data[[wc]] < lower_val | pt_data[[wc]] > upper_val
    )

    if (length(extreme_rows) > 0) {
      log_rows[[wc]] <- data.frame(
        id     = pt_data[[id_col]][extreme_rows],
        weight = wc,
        k      = pt_data$k[extreme_rows],
        value  = pt_data[[wc]][extreme_rows],
        stringsAsFactors = FALSE,
        row.names = NULL
      )
      pt_data[[wc]] <- pmin(pmax(pt_data[[wc]], lower_val), upper_val)
    }
  }

  # 4. Aggregate the per-column logs; compute flagged subject id set
  #    and total flagged row count.
  flagged_log <- if (length(log_rows) == 0) empty_log
                 else do.call(rbind, log_rows)
  flagged_ids   <- unique(flagged_log$id)
  n_flagged_rows <- nrow(flagged_log)

  # 5. Emit a single grouped warning summarizing the impact.
  if (length(flagged_ids) > 0) {
    warning(
      length(flagged_ids), " subject(s) had weights truncated at ",
      "p=[", truncate[1], ", ", truncate[2], "] (",
      n_flagged_rows, " row(s) affected). ",
      "See $flagged_log for details.",
      call. = FALSE
    )
  }

  list(
    pt_data     = pt_data,
    flagged_ids = flagged_ids,
    flagged_log = flagged_log
  )
}


#' Summarize Weight Distributions
#'
#' Reports distributional statistics for every IPW weight column on
#' `pt_data`, for both raw (`*_raw`) and truncated versions when both
#' exist. Recognized columns: `w_cens(_raw)`, `w_a(_raw)`.
#'
#' @param pt_data Data frame with weight columns.
#' @return Data frame with one row per weight column (raw and truncated
#'   for each base name when present) and columns
#'   `weight, n_nonNA, mean, median, min, p001, p01, p99, p999, max`.
#'   NULL if no weight columns are present.
#' @export
#' @keywords internal
summarize_weights <- function(pt_data, w_base = c("w_cens", "w_a")) {
  candidates <- c(paste0(w_base, "_raw"), w_base)
  w_cols <- intersect(candidates, names(pt_data))
  if (length(w_cols) == 0) return(NULL)

  do.call(rbind, lapply(w_cols, function(wc) {
    w <- pt_data[[wc]]
    data.frame(
      weight  = wc,
      n_nonNA = sum(!is.na(w)),
      mean    = mean(w, na.rm = TRUE),
      median  = median(w, na.rm = TRUE),
      min     = min(w, na.rm = TRUE),
      p001    = quantile(w, 0.001, na.rm = TRUE),
      p01     = quantile(w, 0.01,  na.rm = TRUE),
      p99     = quantile(w, 0.99,  na.rm = TRUE),
      p999    = quantile(w, 0.999, na.rm = TRUE),
      max     = max(w, na.rm = TRUE),
      row.names = NULL
    )
  }))
}
