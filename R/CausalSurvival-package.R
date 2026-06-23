#' CausalSurvival: Causal Inference for Single-Event Survival Outcomes
#'
#' Discrete-time pooled logistic regression estimators for causal survival
#' analysis with a single failure event and right censoring. Parametric
#' g-formula (plug-in) and inverse probability weighting (IPW) for cumulative
#' incidence under static, baseline-only treatment regimes. Bootstrap
#' percentile confidence intervals, identifying-assumption accessors, and
#' weight diagnostics.
#'
#' @section Design philosophy:
#' Simplicity and auditability above all. Pure base R, hardcoded simple
#' defaults, small public surface, explicit error messages. See the package
#' specification (`dev/CAUSAL_SURVIVAL_SPEC.md` in the CausalCompetingRisks
#' repository) for the full design.
#'
#' @section Two-package ecosystem:
#' CausalSurvival is the foundation. CausalCompetingRisks (CCR) imports
#' CausalSurvival's primitives and adds the separable-effects framework for
#' competing events.
#'
#' @keywords internal
"_PACKAGE"
