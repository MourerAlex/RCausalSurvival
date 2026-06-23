#' Propensity Model: P(A | L)
#'
#' Fits and predicts the treatment propensity used by the IPW estimator.
#' The propensity weight `1 / pi(A | L)` extrapolates from the observed
#' treatment mixture to the counterfactual full-population estimand
#' (Hernán & Robins, *Causal Inference: What If*, ch. 12).
#'
#' @keywords internal
#' @name propensity
NULL


#' Fit Propensity Score Models
#'
#' Estimate the propensity score \eqn{\pi(a \mid l) = \Pr(A = a \mid L = l)}
#' by logistic regression. With `stabilize = TRUE`, also estimate the
#' marginal model \eqn{\pi(a) = \Pr(A = a)} for stabilized weights
#' (Robins & Hernán).
#'
#' Both fits use one row per subject: the row at `k = 1`. Each row
#' `k` carries \eqn{L_{k-1}}, \eqn{A_{k-1}}, and the event indicator
#' \eqn{Y_k} (left-edge convention: covariates and treatment apply
#' from the start of an interval, events occur at the end). In v1 the
#' treatment is point — \eqn{A_{k-1} = A_0} for all `k` — so the
#' propensity fit reduces to one row per subject at `k = 1`, carrying
#' \eqn{A_0} and \eqn{L_0} measured at `k = 0`.
#'
#' @param pt_data Person-time data frame.
#' @param treatment Character. Treatment column name.
#' @param covariates Character vector. Baseline covariate column names.
#' @param stabilize Logical. If TRUE (default) also fit the marginal
#'   numerator model A ~ 1.
#' @param formula_full Optional formula override for the conditional
#'   denominator model (default: `A ~ L1 + L2 + ...`).
#' @param formula_num Optional formula override for the numerator model.
#'   Default `A ~ 1` (marginal stabilization). Conditional stabilization
#'   (e.g. `A ~ V_baseline`) is supported by the helper but not yet
#'   exposed by the public API — see `dev/TODO.md`.
#'
#' @return List with elements:
#'   - `model_a` — fitted full propensity glm
#'   - `model_a_num` — fitted numerator glm, or NULL if `stabilize = FALSE`
#'   - `check_a` — diagnostics from [check_fitted_positivity()]
#'   - `check_a_num` — diagnostics for the numerator model, or NULL
#' @export
#' @keywords internal
fit_propensity <- function(pt_data, treatment, covariates,
                            stabilize = TRUE,
                            formula_full = NULL, formula_num = NULL) {

  # 1. Baseline subset: one row per subject at k = 1, the time point
  #    where treatment is assigned under the v1 point-treatment scope.
  baseline <- pt_data[pt_data$k == 1, ]

  # 2. Fit the full conditioning denominator model `A ~ L_1 + ... + L_p`
  #    (or `A ~ 1` when covariates is empty / user-overridden formula).
  if (is.null(formula_full)) {
    rhs <- if (length(covariates) > 0) {
      paste(covariates, collapse = " + ")
    } else {
      "1"
    }
    formula_full <- as.formula(paste(treatment, "~", rhs))
  }
  full <- fit_logistic(formula_full, baseline, "Propensity (full)")

  # 3. Optional numerator for stabilization (marginal `A ~ 1` by default;
  #    Robins/Hernán form). Skipped when stabilize = FALSE.
  num <- NULL
  if (stabilize) {
    if (is.null(formula_num)) {
      formula_num <- as.formula(paste(treatment, "~ 1"))
    }
    num <- fit_logistic(formula_num, baseline, "Propensity (numerator)")
  }

  list(
    model_a     = full$model,
    model_a_num = if (is.null(num)) NULL else num$model,
    check_a     = full$check,
    check_a_num = if (is.null(num)) NULL else num$check
  )
}
