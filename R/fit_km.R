#' Unadjusted Kaplan-Meier per arm (with optional IPCW)
#'
#' Discrete-time Kaplan-Meier per arm. For each arm `a in {0, 1}`, filter
#' to the at-risk rows and compute the pooled-hazard
#' \deqn{\hat\lambda^a_k \;=\;
#'   \frac{\sum_i W_i\, 1\{Y_{ik}=1,\, A_i = a\}}
#'        {\sum_i W_i\, 1\{\text{at risk at } k,\, A_i = a\}},
#'   \qquad
#'   \hat F^a(k) \;=\; 1 - \prod_{j=1}^{k}(1 - \hat\lambda^a_j).}
#' Under `ipcw = FALSE` the weights `W_i` are uniform (1), yielding the
#' standard discrete-time KM. Under `ipcw = TRUE` the weights are the
#' inverse-probability-of-censoring weights `w_cens` from a fitted
#' C-hazard model (`cond_indep_cens ~ A + k + I(k^2) + I(k^3)` by default).
#'
#' KM is marginal in baseline `L_0`: no propensity model is fit and any
#' covariates carried on `pt_data` are ignored. The C-hazard model under
#' `ipcw = TRUE` is also marginal in `L_0` (treatment + time only).
#' Identifies the marginal counterfactual under marginal exchangeability
#' (e.g., randomization) and either independent censoring
#' (`ipcw = FALSE`) or no-unmeasured-confounders-of-censoring beyond `A`
#' (`ipcw = TRUE`).
#'
#' Requires baseline treatment. v1 restricts the package to a binary
#' baseline `A`.
#'
#' @keywords internal
fit_km <- function(pt_data, id_col, treatment_col, cut_times,
                   formulas, ipcw, stabilize, truncate) {

  # 1. Optionally fit IPCW machinery (C-hazard + optional stabilized
  #    numerator) -> w_cens. Skipped when ipcw = FALSE (pure KM).
  if (ipcw) {
    haz_fit <- fit_hazard_models(
      pt_data        = pt_data,
      treatment      = treatment_col,
      covariates     = character(),
      active_methods = "ipw",
      formulas       = formulas,
      ipcw           = TRUE
    )
    model_c <- haz_fit$models$model_c
    check_c <- haz_fit$checks$c

    model_c_num <- NULL
    check_c_num <- NULL
    if (identical(stabilize, "marginal")) {
      cnum_fml  <- as.formula(paste("cond_indep_cens ~", treatment_col))
      cnum_rows <- pt_data$indep_cens == 0
      cnum_fit  <- fit_logistic(
        cnum_fml, pt_data[cnum_rows, , drop = FALSE],
        "C-hazard (numerator)"
      )
      model_c_num <- cnum_fit$model
      check_c_num <- cnum_fit$check
    }

    w_cens_raw         <- ipw_cens(model_c, pt_data, id_col,
                                   model_num = model_c_num)
    pt_data$w_cens_raw <- w_cens_raw
    pt_data$w_cens     <- w_cens_raw
    trunc_out          <- apply_weight_truncation(pt_data, id_col, truncate)
    pt_data            <- trunc_out$pt_data
    flagged_ids        <- trunc_out$flagged_ids
  } else {
    model_c     <- NULL
    check_c     <- NULL
    model_c_num <- NULL
    check_c_num <- NULL
    flagged_ids <- NULL
  }

  # 2. Per-arm Hajek pooled hazard. Weights: w_cens when ipcw, else 1.
  cif_by_arm <- lapply(c(0, 1), function(a) {
    arm_rows <- pt_data[
      pt_data[[treatment_col]] == a &
        pt_data$cond_indep_cens   == 0 &
        pt_data$indep_cens == 0,
      , drop = FALSE
    ]
    weights <- if (ipcw) arm_rows$w_cens else rep(1, nrow(arm_rows))
    cum_inc_from_weighted(
      y_event   = arm_rows$y_event,
      k         = arm_rows$k,
      weights   = weights,
      cut_times = cut_times
    )
  })

  # 3. Long-format estimates.
  estimates <- make_estimates_long(cif_by_arm, cut_times)

  list(
    estimates    = estimates,
    models       = list(y = NULL, c = model_c, A = NULL,
                        A_num = NULL, c_num = model_c_num),
    model_checks = list(y = NULL, c = check_c, A = NULL,
                        A_num = NULL, c_num = check_c_num),
    weights      = if (ipcw) {
      list(
        pt_data_weighted = pt_data,
        weight_summary   = summarize_weights(pt_data),
        truncated_ids    = flagged_ids,
        truncate         = truncate
      )
    } else NULL
  )
}
