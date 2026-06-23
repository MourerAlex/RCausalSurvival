# ----------------------------------------------------------------------------
# IPW family for causal_survival(method = "ipw").
#
# Exposes:
#   - fit_ipw()           dispatch shim (km / msm), called from R/causal_survival.R
#   - fit_ipw_weights()   shared weight builder (propensity x IPCW)
#   - fit_ipw_msm()       weighted pooled-logistic Y-MSM survival estimator
#   - fit_ipw_km()        weighted pooled-hazard Kaplan-Meier estimator
#
# fit_ipw_msm() reuses make_clone() from R/fit_gformula.R via the package
# namespace (no explicit import needed).
# ----------------------------------------------------------------------------


#' Dispatch the IPW survival-curve estimator (`km` or `msm`)
#'
#' Thin shim that resolves `ipw_estimator` to the underlying worker
#' ([fit_ipw_km()] or [fit_ipw_msm()]) and forwards the call. Lets the
#' [causal_survival()] body show one named verb per branch (`fit_gformula`
#' vs `fit_ipw`) instead of a nested switch.
#'
#' @keywords internal
fit_ipw <- function(pt_data, id_col, treatment_col, covariates_vec,
                    cut_times, formulas, ipcw, stabilize, truncate,
                    ipw_estimator = c("km", "msm")) {
  ipw_estimator <- match.arg(ipw_estimator)
  worker <- switch(ipw_estimator, km = fit_ipw_km, msm = fit_ipw_msm)
  worker(pt_data, id_col, treatment_col, covariates_vec, cut_times,
         formulas, ipcw, stabilize, truncate)
}


#' Build IPW Weights
#'
#' Compute the inverse-probability weight for each person-time row,
#' \deqn{W_i(k) \;=\; \frac{1}{\pi(A_i \mid L_i)} \,
#'                   \prod_{j=1}^{k}\frac{1}{1 - h_C(j;\, A_i, L_i)},}
#' replacing each factor by its stabilized ratio when
#' `stabilize = "marginal"`:
#'
#' | `stabilize`  | propensity numerator | censoring numerator       |
#' |--------------|----------------------|---------------------------|
#' | `"marginal"` | `A ~ 1`              | `cond_indep_cens ~ A` (when ipcw)|
#' | NULL         | none (unstabilized)  | none (unstabilized)       |
#'
#' Truncate at the requested percentile and attach `w_a`, `w_cens`,
#' `w_combined` (and their raw counterparts) to `pt_data`. Shared by
#' [fit_ipw_msm()] and [fit_ipw_km()]; the engines differ only in the
#' downstream survival-estimation step. The censoring numerator is
#' fixed to `cond_indep_cens ~ A` (H&R Technical Point 12.2).
#'
#' @return List with
#'   - `pt_data` — augmented with weight columns.
#'   - `models`  — list of fitted glms (`A`, `A_num`, `c`, `c_num`).
#'   - `checks`  — list of per-model diagnostics.
#'   - `flagged_ids` — ids of subjects with weights flagged by
#'     truncation (or `NULL` when no truncation requested).
#' @keywords internal
fit_ipw_weights <- function(pt_data, id_col, treatment_col,
                            covariates_vec, formulas, ipcw, stabilize,
                            truncate) {

  do_stabilize <- identical(stabilize, "marginal")

  # ---------- 1. Propensity model(s) on the k = 1 row ----------
  prop_fit <- fit_propensity(
    pt_data       = pt_data,
    treatment     = treatment_col,
    covariates    = covariates_vec,
    stabilize     = do_stabilize,
    formula_full  = formulas$A,
    formula_num   = formulas$A_num
  )
  model_a     <- prop_fit$model_a
  model_a_num <- prop_fit$model_a_num

  # ---------- 2. Censoring model(s) for IPCW ----------
  model_c     <- NULL
  model_c_num <- NULL
  check_c     <- NULL
  check_c_num <- NULL
  if (ipcw) {
    haz_fit <- fit_hazard_models(
      pt_data        = pt_data,
      treatment      = treatment_col,
      covariates     = covariates_vec,
      active_methods = "ipw",
      formulas       = formulas,
      ipcw           = TRUE
    )
    model_c <- haz_fit$models$model_c
    check_c <- haz_fit$checks$c

    if (do_stabilize) {
      cnum_fml  <- as.formula(
        paste("cond_indep_cens ~", treatment_col)
      )
      cnum_rows <- pt_data$indep_cens == 0
      cnum_fit  <- fit_logistic(
        cnum_fml, pt_data[cnum_rows, , drop = FALSE],
        "C-hazard (numerator)"
      )
      model_c_num <- cnum_fit$model
      check_c_num <- cnum_fit$check
    }
  }

  # ---------- 3. Build raw weights ----------
  w_a_raw <- ipw_static_trt(
    model_full    = model_a,
    pt_data       = pt_data,
    treatment_col = treatment_col,
    id_col        = id_col,
    model_num     = model_a_num
  )
  w_cens_raw <- if (ipcw) {
    ipw_cens(model_c, pt_data, id_col, model_num = model_c_num)
  } else {
    NULL
  }

  # `*_raw` columns survive truncation so reweight() can re-apply
  # truncation without refitting upstream models.
  pt_data$w_a_raw <- w_a_raw
  pt_data$w_a     <- w_a_raw
  if (ipcw) {
    pt_data$w_cens_raw <- w_cens_raw
    pt_data$w_cens     <- w_cens_raw
  }

  # ---------- 4. Truncation ----------
  trunc_out <- apply_weight_truncation(
    pt_data  = pt_data,
    id_col   = id_col,
    truncate = truncate
  )
  pt_data <- trunc_out$pt_data

  # ---------- 5. Combined per-row weight ----------
  pt_data$w_combined <- if (ipcw) {
    pt_data$w_a * pt_data$w_cens
  } else {
    pt_data$w_a
  }

  list(
    pt_data     = pt_data,
    models      = list(A = model_a, A_num = model_a_num,
                       c = model_c, c_num = model_c_num),
    checks      = list(A = prop_fit$check_a,
                       A_num = prop_fit$check_a_num,
                       c = check_c, c_num = check_c_num),
    flagged_ids = trunc_out$flagged_ids
  )
}


#' IPW Cumulative Incidence Worker (MSM engine)
#'
#' Build IPW weights via [fit_ipw_weights()], then fit a weighted
#' pooled-logistic Y-MSM and marginalize per arm by
#' clone-predict-marginalize. The MSM is parametric in `k`
#' (`y_event ~ A + k + I(k^2) + I(k^3)` by default).
#'
#' Equivalent unrolled per-arm form (kept here for explanatory reference;
#' the implementation uses an `lapply` over both arms to avoid
#' duplication):
#'
#' ```
#' # (after fit_ipw_weights() has attached the per-row IPW weight
#' #  `w_combined` to pt_data)
#'
#' # 1. Weighted pooled-logistic Y-MSM fit. The IPW weight enters HERE
#' #    via the `weights =` argument and gets baked into model_y's
#' #    coefficients.
#' y_rows  <- pt_data$indep_cens == 0 & pt_data$cond_indep_cens == 0
#' fml_y   <- y_event ~ A + k + I(k^2) + I(k^3)
#' model_y <- fit_logistic(fml_y, pt_data[y_rows, , drop = FALSE],
#'                         label   = "Y-MSM (IPW)",
#'                         weights = pt_data$w_combined[y_rows])$model
#'
#' # 2. Per-arm marginalization. Weights no longer appear here — they
#' #    already shaped model_y's coefficients above.
#' baseline     <- pt_data[pt_data$k == 1, , drop = FALSE]
#' clone        <- make_clone(baseline, cut_times, treatment_col, a)
#' haz_a1       <- predict_counterfactual_hazard(model_y, clone,
#'                                               treatment_col, 1, "Y-MSM a=1")
#' haz_a0       <- predict_counterfactual_hazard(model_y, clone,
#'                                               treatment_col, 0, "Y-MSM a=0")
#' S_a1         <- cumprod_survival(haz_a1, clone[[id_col]])
#' S_a0         <- cumprod_survival(haz_a0, clone[[id_col]])
#' surv_a1_by_k <- tapply(S_a1, clone$k, mean)
#' surv_a0_by_k <- tapply(S_a0, clone$k, mean)
#' ```
#'
#' @keywords internal
fit_ipw_msm <- function(pt_data, id_col, treatment_col, covariates_vec,
                        cut_times, formulas, ipcw, stabilize, truncate) {

  w_out   <- fit_ipw_weights(
    pt_data        = pt_data,
    id_col         = id_col,
    treatment_col  = treatment_col,
    covariates_vec = covariates_vec,
    formulas       = formulas,
    ipcw           = ipcw,
    stabilize      = stabilize,
    truncate       = truncate
  )
  pt_data <- w_out$pt_data

  # ---------- 6. Weighted Y-MSM fit ----------
  # Default Y-MSM is marginal in covariates: weights handle confounding,
  # so no covariate adjustment in the outcome model. Users can supply
  # `formulas$y` for a covariate-conditional MSM.
  # Fit population: rows with indep_cens == 0 & cond_indep_cens == 0.
  time_terms <- "k + I(k^2) + I(k^3)"
  fml_y <- formulas$y %||% as.formula(
    paste("y_event ~",
          paste(c(treatment_col, time_terms), collapse = " + "))
  )
  y_rows  <- pt_data$indep_cens == 0 & pt_data$cond_indep_cens == 0
  msm_fit <- fit_logistic(
    formula = fml_y,
    data    = pt_data[y_rows, , drop = FALSE],
    label   = "Y-MSM (IPW)",
    weights = pt_data$w_combined[y_rows]
  )
  model_y <- msm_fit$model
  check_y <- msm_fit$check

  # ---------- 7. Per-arm CIF: clone -> predict -> cumprod -> mean ----------
  baseline   <- pt_data[pt_data$k == 1, , drop = FALSE]
  cif_by_arm <- lapply(c(0, 1), function(a) {
    clone <- make_clone(baseline, cut_times, treatment_col, a)
    haz   <- predict_counterfactual_hazard(
      model_y, clone, treatment_col, a,
      paste0("Y-MSM a=", a)
    )
    if (any(is.na(haz))) {
      warning(
        "IPW (MSM): Y-MSM predictions contain ", sum(is.na(haz)),
        " NA value(s). CIF estimates will be biased.",
        call. = FALSE
      )
    }
    S_i <- cumprod_survival(haz, clone[[id_col]])
    S_k <- as.numeric(tapply(S_i, clone$k, mean))
    1 - S_k
  })

  estimates <- make_estimates_long(cif_by_arm, cut_times)

  list(
    estimates    = estimates,
    models       = list(
      y     = model_y,
      c     = w_out$models$c,
      A     = w_out$models$A,
      A_num = w_out$models$A_num,
      c_num = w_out$models$c_num
    ),
    model_checks = c(list(y = check_y), w_out$checks),
    weights      = list(
      pt_data_weighted = pt_data,
      weight_summary   = summarize_weights(pt_data),
      truncated_ids    = w_out$flagged_ids,
      truncate         = truncate
    )
  )
}


#' IPW Cumulative Incidence — Weighted KM Engine
#'
#' Estimate the counterfactual cumulative incidence under each arm by
#' a weighted pooled-hazard Kaplan-Meier estimator:
#' \deqn{\hat\lambda^a_k \;=\;
#'   \frac{\sum_i W_i\, 1\{Y_{ik}=1,\, A_i = a\}}
#'        {\sum_i W_i\, 1\{\text{at risk at } k,\, A_i = a\}},
#'   \qquad
#'   \hat F^a(k) \;=\; 1 - \prod_{j=1}^{k}(1 - \hat\lambda^a_j).}
#' Weights `W_i` come from [fit_ipw_weights()]. The hazard is
#' nonparametric in `k`; no outcome model is fit.
#'
#' Requires baseline treatment. The arm-specific risk set is undefined
#' under time-varying `A`; use [fit_ipw_msm()] in that case. v1
#' restricts the package to baseline `A`.
#'
#' Equivalent unrolled per-arm form (kept here for explanatory reference;
#' the implementation uses an `lapply` over both arms to avoid
#' duplication):
#'
#' ```
#' # (after fit_ipw_weights() has attached w_combined to pt_data)
#' rows_a1 <- pt_data[pt_data[[treatment_col]] == 1 &
#'                    pt_data$cond_indep_cens   == 0 &
#'                    pt_data$indep_cens == 0, , drop = FALSE]
#' rows_a0 <- pt_data[pt_data[[treatment_col]] == 0 &
#'                    pt_data$cond_indep_cens   == 0 &
#'                    pt_data$indep_cens == 0, , drop = FALSE]
#' cif_a1  <- cum_inc_from_weighted(
#'              y_event = rows_a1$y_event, k = rows_a1$k,
#'              weights = rows_a1$w_combined, cut_times = cut_times)
#' cif_a0  <- cum_inc_from_weighted(
#'              y_event = rows_a0$y_event, k = rows_a0$k,
#'              weights = rows_a0$w_combined, cut_times = cut_times)
#' ```
#'
#' @keywords internal
fit_ipw_km <- function(pt_data, id_col, treatment_col, covariates_vec,
                       cut_times, formulas, ipcw, stabilize, truncate) {

  w_out   <- fit_ipw_weights(
    pt_data        = pt_data,
    id_col         = id_col,
    treatment_col  = treatment_col,
    covariates_vec = covariates_vec,
    formulas       = formulas,
    ipcw           = ipcw,
    stabilize      = stabilize,
    truncate       = truncate
  )
  pt_data <- w_out$pt_data

  # ---------- 6. Weighted pooled hazard per arm (Kaplan-Meier-style) ----------
  # Fit population: rows with indep_cens == 0 & cond_indep_cens == 0 (KM
  # symmetry with the Y-MSM fit).
  cif_by_arm <- lapply(c(0, 1), function(a) {
    arm_rows <- pt_data[
      pt_data[[treatment_col]] == a &
        pt_data$cond_indep_cens   == 0 &
        pt_data$indep_cens == 0,
      , drop = FALSE
    ]
    cum_inc_from_weighted(
      y_event   = arm_rows$y_event,
      k         = arm_rows$k,
      weights   = arm_rows$w_combined,
      cut_times = cut_times
    )
  })

  estimates <- make_estimates_long(cif_by_arm, cut_times)

  list(
    estimates    = estimates,
    models       = list(
      y     = NULL,                       # no Y outcome model under KM
      c     = w_out$models$c,
      A     = w_out$models$A,
      A_num = w_out$models$A_num,
      c_num = w_out$models$c_num
    ),
    model_checks = c(list(y = NULL), w_out$checks),
    weights      = list(
      pt_data_weighted = pt_data,
      weight_summary   = summarize_weights(pt_data),
      truncated_ids    = w_out$flagged_ids,
      truncate         = truncate
    )
  )
}
