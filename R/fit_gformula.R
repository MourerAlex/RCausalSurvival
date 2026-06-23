# ----------------------------------------------------------------------------
# G-formula worker for causal_survival(method = "gformula").
#
# Exposes:
#   - fit_gformula()   internal worker dispatched from R/causal_survival.R
#
# Uses make_clone() from R/utils.R (shared with fit_ipw_msm() in R/fit_ipw.R).
# ----------------------------------------------------------------------------


#' G-Formula Cumulative Incidence Worker
#'
#' Estimand: \eqn{E[Y_k^{a, c=0}]} under conditional exchangeability given
#' baseline `L_0`. Identification: discrete-time g-formula recursion over
#' the Y-hazard.
#'
#' Single-event parametric g-formula. Fits an unweighted Y-hazard pooled
#' logistic model, then for each arm `a in {0, 1}`:
#'
#' 1. Clone baseline covariates across all `cut_times`, setting treatment to `a`.
#' 2. Predict the counterfactual hazard at every (subject, k).
#' 3. Compute per-subject survival `S_i(k) = prod_{j <= k}(1 - h(j | a, L_i))`.
#' 4. Marginalize over L by averaging across subjects: `S(k) = mean_i S_i(k)`.
#' 5. CIF = 1 - S.
#'
#' Equivalent unrolled per-arm form (kept here for explanatory reference;
#' the implementation uses an `lapply` over both arms to avoid duplication):
#'
#' ```
#' clone        <- make_clone(baseline, cut_times, treatment_col, a)
#' haz_a1       <- predict_counterfactual_hazard(model_y, clone,
#'                                               treatment_col, 1, "Y-haz a=1")
#' haz_a0       <- predict_counterfactual_hazard(model_y, clone,
#'                                               treatment_col, 0, "Y-haz a=0")
#' S_a1         <- cumprod_survival(haz_a1, clone[[id_col]])
#' S_a0         <- cumprod_survival(haz_a0, clone[[id_col]])
#' surv_a1_by_k <- tapply(S_a1, clone$k, mean)
#' surv_a0_by_k <- tapply(S_a0, clone$k, mean)
#' ```
#'
#' @keywords internal
fit_gformula <- function(pt_data, id_col, treatment_col, covariates_vec,
                         cut_times, formulas) {

  # 1. Fit unweighted Y-hazard
  fit <- fit_hazard_models(
    pt_data        = pt_data,
    treatment      = treatment_col,
    covariates     = covariates_vec,
    active_methods = "gformula",
    formulas       = formulas,
    ipcw           = FALSE
  )
  model_y <- fit$models$model_y

  # 2. Baseline per subject (k = 1, first analyzable interval under
  # the LOCKED (0, t_1] convention; left-truncation rejected upstream)
  baseline <- pt_data[pt_data$k == 1, , drop = FALSE]

  # 3. Per-arm CIF: clone -> predict -> cumprod -> mean -> 1 - S
  cif_by_arm <- lapply(c(0, 1), function(a) {
    clone <- make_clone(baseline, cut_times, treatment_col, a)
    haz   <- predict_counterfactual_hazard(
      model_y, clone, treatment_col, a,
      paste0("Y-hazard a=", a)
    )
    if (any(is.na(haz))) {
      warning(
        "G-formula: Y-hazard predictions contain ", sum(is.na(haz)),
        " NA value(s). CIF estimates will be biased.",
        call. = FALSE
      )
    }
    S_i <- cumprod_survival(haz, clone[[id_col]])
    S_k <- as.numeric(tapply(S_i, clone$k, mean))
    1 - S_k
  })

  # 4. Long-format estimates (k = integer index 1..K_max per spec §3.0.2;
  # time = cut_times[k] for human-facing display).
  estimates <- make_estimates_long(cif_by_arm, cut_times)

  list(
    estimates    = estimates,
    models       = list(y = model_y, c = NULL, A = NULL, A_num = NULL),
    model_checks = fit$checks,
    weights      = NULL
  )
}
