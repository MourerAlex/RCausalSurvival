#' Identifying assumptions for a causal_survival_fit
#'
#' Returns the hardcoded baseline identification block for the
#' counterfactual risk \eqn{E[Y_k^{a, c = 0}]} targeted by
#' [causal_survival()]. Each entry is a named list with fields
#' `name`, `statement`, `status`, and `pointer` (spec §4.6). The
#' `pointer` field directs the user to the diagnostic accessor that
#' can supply a necessary-but-not-sufficient check, when one exists.
#'
#' Untestable conditions are reported as `status = "untestable"`;
#' testable ones as `status = "testable"`. No data-driven flag is
#' set inside this accessor — the corresponding diagnostic must be
#' invoked separately.
#'
#' @param fit A `"causal_survival_fit"` object from
#'   [causal_survival()].
#'
#' @return An S3 object of class `"causal_survival_assumptions"`
#'   with one list element per assumption.
#'
#' @seealso [causal_survival()]
#' @family accessors
#' @export
causal_assumptions <- function(fit) {
  stopifnot(inherits(fit, "causal_survival_fit"))
  if (fit$method == "km") {
    return(.causal_assumptions_km(fit))
  }
  structure(
    list(
      list(
        name      = "Consistency",
        statement = paste0(
          "If A = a and C_bar_k = 0, the counterfactual equals the ",
          "observed value. No hidden version of the treatment."
        ),
        formula   = "A = a,\\ \\bar{C}_k = 0 \\;\\Rightarrow\\; Y^{a,\\bar{c}=0}_k = Y_k",
        status    = "untestable",
        pointer   = NA_character_,
        citation  = "Hernan & Robins (2020), What If, Chapter 3."
      ),
      list(
        name      = "Exchangeability",
        statement = paste0(
          "What would happen under no censoring is independent of which ",
          "treatment was received, given baseline covariates L_0. Holds ",
          "by randomisation; in observational data it is the ",
          "no-unmeasured-confounding condition."
        ),
        formula   = "Y^{a,\\bar{c}=0} \\perp A \\mid L_0",
        status    = "untestable",
        pointer   = NA_character_,
        citation  = "Hernan & Robins (2020), What If, Chapter 3."
      ),
      list(
        name      = "Positivity",
        statement = paste0(
          "For every covariate pattern with positive density, both ",
          "treatments occur with positive probability. Else the ",
          "identifying formulas divide by zero."
        ),
        formula   = "P(L_0 = l) > 0 \\;\\Rightarrow\\; P(A = a \\mid L_0 = l) > 0,\\ a \\in \\{0,1\\}",
        status    = "testable",
        pointer   = "fit$weights$weight_summary (ipw); fit$model_checks (propensity)",
        citation  = "Hernan & Robins (2020), What If, Chapter 3."
      ),
      list(
        name      = "No interference",
        statement = paste0(
          "One subject's treatment does not affect another subject's ",
          "potential outcomes."
        ),
        formula   = NA_character_,
        status    = "untestable",
        pointer   = NA_character_,
        citation  = "Hernan & Robins (2020), What If, Chapter 1."
      ),
      list(
        name      = "Correct model specification",
        statement = paste0(
          "The fitted discrete-time hazard models for Y and (under ipcw) ",
          "C, and the propensity model for A, are correctly specified."
        ),
        formula   = NA_character_,
        status    = "untestable",
        pointer   = "fit$model_checks",
        citation  = "Hernan & Robins (2020), What If, Chapter 18."
      ),
      list(
        name      = "Independent censoring (E2)",
        statement = paste0(
          "What would happen under no censoring is independent of whether ",
          "censoring occurs, given observed history. Ignorable censoring."
        ),
        formula   = "Y^{a,\\bar{c}=0}_{k+1} \\perp C_{k+1} \\mid Y_k = \\bar{C}_k = 0,\\, L_0,\\, A",
        status    = "untestable",
        pointer   = NA_character_,
        citation  = "Hernan & Robins (2020), What If, Chapter 12."
      )
    ),
    class = "causal_survival_assumptions"
  )
}


#' Identification block under method = "km"
#'
#' KM identifies the marginal counterfactual under MARGINAL
#' exchangeability (e.g., randomization) rather than the conditional
#' exchangeability that g-formula / IPW invoke. No propensity model is
#' fit, so positivity reduces to the trivial marginal version. The
#' censoring block depends on `fit$ipcw`: when `FALSE`, the user is
#' asserting independent censoring (no correction applied); when
#' `TRUE`, a conditionally-independent-censoring assumption +
#' the C-hazard specification check kick in.
#'
#' @keywords internal
.causal_assumptions_km <- function(fit) {
  base_assumptions <- list(
    list(
      name      = "Consistency",
      statement = paste0(
        "If A = a and C_bar_k = 0, the counterfactual equals the ",
        "observed value. No hidden version of the treatment."
      ),
      formula   = "A = a,\\ \\bar{C}_k = 0 \\;\\Rightarrow\\; Y^{a,\\bar{c}=0}_k = Y_k",
      status    = "untestable",
      pointer   = NA_character_,
      citation  = "Hernan & Robins (2020), What If, Chapter 3."
    ),
    list(
      name      = "Marginal exchangeability",
      statement = paste0(
        "Y^a is independent of A (unconditional). For confounder ",
        "adjustment use method = \"gformula\" or method = \"ipw\"."
      ),
      formula   = "Y^{a,\\bar{c}=0} \\perp A",
      status    = "untestable",
      pointer   = NA_character_,
      citation  = "Hernan & Robins (2020), What If, Chapter 2."
    ),
    list(
      name      = "No interference",
      statement = paste0(
        "One subject's treatment does not affect another subject's ",
        "potential outcomes."
      ),
      formula   = NA_character_,
      status    = "untestable",
      pointer   = NA_character_,
      citation  = "Hernan & Robins (2020), What If, Chapter 1."
    )
  )
  censoring_assumptions <- if (isTRUE(fit$ipcw)) {
    list(
      list(
        name      = "Conditionally-independent censoring",
        statement = paste0(
          "Censoring is independent of the outcome given A and k. ",
          "The IPCW correction targets this conditionally-independent ",
          "mechanism (cond_indep_cens)."
        ),
        formula   = NA_character_,
        status    = "untestable",
        pointer   = NA_character_,
        citation  = "Hernan & Robins (2020), What If, Chapter 12."
      ),
      list(
        name      = "C-positivity",
        statement = paste0(
          "P(cond_indep_cens = 0 | A, k) > 0 across the observed support — ",
          "no (arm, interval) with near-0 P(cond_indep_cens = 0 | A, k). ",
          "Violations inflate the inverse-probability-of-censoring ",
          "weights and the resulting survival bands. If structural, ",
          "they point to an estimand that is not possible to target ",
          "given the data."
        ),
        formula   = "P(\\bar{C}_k = 0 \\mid A, k) > 0",
        status    = "testable",
        pointer   = "fit$model_checks$c (min_fitted, max_fitted)",
        citation  = "Hernan & Robins (2020), What If, Chapter 3."
      ),
      list(
        name      = "Correct C-hazard specification",
        statement = paste0(
          "The fitted discrete-time hazard model for C is correctly ",
          "specified (the only model fit under method = \"km\" + ",
          "ipcw = TRUE)."
        ),
        formula   = NA_character_,
        status    = "untestable",
        pointer   = "fit$model_checks$c",
        citation  = "Hernan & Robins (2020), What If, Chapter 18."
      )
    )
  } else {
    list(
      list(
        name      = "Independent censoring",
        statement = paste0(
          "Censoring is independent of T^a — the censoring mechanism ",
          "does not depend on prognosis. No correction is applied ",
          "(ipcw = FALSE). For dependent censoring set ipcw = TRUE."
        ),
        formula   = NA_character_,
        status    = "untestable",
        pointer   = NA_character_,
        citation  = "Hernan & Robins (2020), What If, Chapter 12."
      )
    )
  }
  structure(c(base_assumptions, censoring_assumptions),
            class = "causal_survival_assumptions")
}


#' Print a causal_survival_assumptions object
#'
#' Renders the baseline identification block as a numbered list with
#' the status tag and (when applicable) the diagnostic pointer.
#'
#' @param x A `"causal_survival_assumptions"` object.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.causal_survival_assumptions <- function(x, ...) {
  cat("Identifying assumptions (causal_survival_assumptions)\n")
  cat("-----------------------------------------------------\n")
  for (i in seq_along(x)) {
    a <- x[[i]]
    tag <- if (a$status == "testable") "[testable]" else "[untestable]"
    cat(sprintf("%d. %s  %s\n", i, a$name, tag))
    cat("   ", a$statement, "\n", sep = "")
    if (!is.null(a$formula) && !is.na(a$formula)) {
      cat("   Formula: ", a$formula, "\n", sep = "")
    }
    if (!is.na(a$pointer)) {
      cat("   See: ", a$pointer, "\n", sep = "")
    }
    if (!is.null(a$citation) && !is.na(a$citation)) {
      cat("   Citation: ", a$citation, "\n", sep = "")
    }
    if (i < length(x)) cat("\n")
  }
  invisible(x)
}
