#' Counterfactual risk under a binary point treatment
#'
#' Estimate the counterfactual risk
#' \eqn{E[Y_k^{a, c = 0}] = \Pr(Y_k^{a, c = 0} = 1)} that the event of
#' interest has occurred by the end of discrete interval `k` under
#' arm `a` of a binary treatment assigned at baseline (`t = 0`) and
#' held fixed thereafter, with the treatment-dependent censoring
#' component held at zero. Equivalently, the cumulative incidence
#' \eqn{F^a(t_k) = \Pr(T^a \le t_k, C^{d,a} = 0)}; the survival
#' \eqn{S^a(t_k) = 1 - F^a(t_k)} follows. Time is discretized into a
#' finite grid of intervals `k = 1, …, K_max` under the half-open
#' convention `(t_{k-1}, t_k]`.
#'
#' Identification rests on consistency, exchangeability conditional
#' on the recorded baseline covariates `L_0`, and positivity in both
#' the treatment and the censoring mechanisms. Censoring is split
#' into a treatment-dependent component `cond_indep_cens` (the component
#' `c` set to zero in the estimand) and an administrative or
#' otherwise independent component `indep_cens` (handled by the
#' at-risk set, not by intervention).
#'
#' Two estimators are exposed:
#'
#' - `method = "gformula"`: discrete-time parametric g-formula. Fit
#'   the discrete hazards of `Y` (and, when needed, of `C`) on the
#'   pooled person-time data, then simulate the counterfactual risk
#'   `E[Y_k^{a, c = 0}]` in each arm by recursive substitution.
#' - `method = "ipw"`: inverse-probability weighting. Build the
#'   product weight (treatment and, when applicable, censoring) and
#'   estimate the counterfactual hazard in each arm by a weighted
#'   pooled-hazard Kaplan-Meier estimator (default, nonparametric in `k`) or by a
#'   weighted pooled-logistic marginal structural model
#'   (`.ipw_estimator = "msm"`); the risk is recovered from the
#'   hazards.
#' - `method = "km"`: unadjusted Kaplan-Meier per arm. No propensity
#'   model is fit (any covariates on `pt_data` are silently ignored).
#'   Optionally accepts `ipcw = TRUE` to correct for treatment-
#'   conditional dependent censoring via inverse-probability weights.
#'   Identifies the marginal counterfactual under MARGINAL
#'   exchangeability (e.g., randomization) and independent censoring
#'   (or conditional-on-treatment independent censoring when
#'   `ipcw = TRUE`).
#'
#' Standard errors are not produced by this function; obtain them by
#' bootstrap on the subject-level input.
#'
#' @param pt_data A `person_time` object returned by
#'   [to_person_time()].
#' @param method Estimator. One of `"gformula"`, `"ipw"`, or `"km"`.
#' @param formulas Optional named list of model formula overrides.
#'   Keys: `y` (Y-hazard), `c` (C-hazard / IPCW denominator), `A`
#'   (propensity denominator), `A_num` (propensity numerator). Any
#'   absent key falls back to the default linear formula.
#' @param truncate `NULL` or length-2 numeric `c(lower, upper)`
#'   percentile bounds for IPW weight truncation. `NULL` leaves the
#'   weights untruncated.
#' @param ipcw `NULL` or logical. `NULL` resolves to the
#'   method-conditional default (`TRUE` under `"ipw"`, `FALSE` under
#'   `"gformula"`).
#' @param stabilize One of `"marginal"` or `NULL`. Drives both the
#'   treatment and the censoring numerators; v1 supports marginal
#'   stabilization only.
#' @param verbose Logical.
#' @param keep_data Logical. When `TRUE`, the `pt_data` (and the
#'   subject-level input it was built from) are retained on the
#'   returned fit.
#' @param .ipw_estimator Internal. Survival-curve estimator under
#'   `method = "ipw"`: `"km"` (default, weighted pooled-hazard Kaplan-Meier,
#'   nonparametric in `k`) or `"msm"` (weighted pooled logistic with
#'   the cubic-in-`k` default). The leading dot flags this as a
#'   developer-facing knob.
#'
#' @return An S3 object of class `"causal_survival_fit"`.
#' @export
causal_survival <- function(pt_data,
                            method         = "gformula",
                            formulas       = NULL,
                            truncate       = NULL,
                            ipcw           = NULL,
                            stabilize      = "marginal",
                            verbose        = FALSE,
                            keep_data      = TRUE,
                            .ipw_estimator = "km") {

  # Freeze argument VALUES into the stored call so bootstrap()'s
  # eval(fit$call) is self-contained: an argument passed as a variable
  # (e.g. method = chosen_method inside a wrapper function) must not need
  # to still be in scope at bootstrap time. pt_data is left as its original
  # symbol — bootstrap overwrites it per replicate.
  cl <- match.call()
  pf <- parent.frame()
  for (nm in setdiff(names(cl)[-1L], "pt_data")) {
    cl[[nm]] <- eval(cl[[nm]], pf)
  }

  # 1. Validate / canonicalize args (also unpacks pt_data attrs onto args)
  args <- .validate_causal_survival_args(
            pt_data, method, formulas, truncate, ipcw, stabilize,
            verbose, keep_data, .ipw_estimator, cl = cl)

  if (args$verbose) {
    message("causal_survival(): fitting method = '", args$method, "'")
  }

  # 2. Fit (worker dispatch — fit_gformula, fit_km, or fit_ipw;
  #    glm warnings collected)
  worker_out <- .with_collected_warnings(
    if (args$method == "gformula") {
      fit_gformula(
        pt_data        = args$pt_data,
        id_col         = args$id_col,
        treatment_col  = args$treatment_col,
        covariates_vec = args$covariates_vec,
        cut_times      = args$cut_times,
        formulas       = args$formulas
      )
    } else if (args$method == "km") {
      fit_km(
        pt_data       = args$pt_data,
        id_col        = args$id_col,
        treatment_col = args$treatment_col,
        cut_times     = args$cut_times,
        formulas      = args$formulas,
        ipcw          = args$ipcw,
        stabilize     = args$stabilize,
        truncate      = args$truncate
      )
    } else {
      fit_ipw(
        pt_data        = args$pt_data,
        id_col         = args$id_col,
        treatment_col  = args$treatment_col,
        covariates_vec = args$covariates_vec,
        cut_times      = args$cut_times,
        formulas       = args$formulas,
        ipcw           = args$ipcw,
        stabilize      = args$stabilize,
        truncate       = args$truncate,
        ipw_estimator  = args$.ipw_estimator
      )
    }
  )

  # 3. Assemble S3 fit
  .assemble_causal_survival_fit(worker_out, args)
}


# ----------------------------------------------------------------------------
# Internal helpers for causal_survival(): validate, warning-collect, IPW
# dispatch shim, output assembly. All plumbing lives here so the public API
# body above stays a clean methodological pipeline (one verb per tile).
# ----------------------------------------------------------------------------

#' Validate and canonicalize `causal_survival()` arguments
#'
#' All input-shape checks surface here. Resolves `ipcw = NULL` to its
#' method-conditional default (`TRUE` under `"ipw"`, `FALSE` under
#' `"gformula"`). Also unpacks the `person_time` attributes (`cut_times`,
#' `id_col`, `treatment_col`, `covariates`, `treatment_levels`) onto the
#' returned list so downstream workers and the assembler can read them as
#' plain fields. Returns a canonical, normalized argument list.
#'
#' @keywords internal
.validate_causal_survival_args <- function(pt_data, method, formulas,
                                           truncate, ipcw, stabilize,
                                           verbose, keep_data,
                                           .ipw_estimator, cl) {

  # pt_data class check (stricter than CCR: classed only)
  if (!inherits(pt_data, "person_time")) {
    stop("pt_data must inherit 'person_time'. ",
         "Run to_person_time() on subject-level data first.",
         call. = FALSE)
  }

  # method (single value, not vector)
  valid_methods <- c("gformula", "ipw", "km")
  if (length(method) != 1L || !method %in% valid_methods) {
    stop("method must be one of: ",
         paste(shQuote(valid_methods), collapse = ", "),
         ". Got: ", paste(method, collapse = ", "), call. = FALSE)
  }

  # .ipw_estimator: only validated and stored when method == "ipw". For
  # gformula / km the IPW engine choice is meaningless, so we coerce to
  # NULL without complaint (avoids storing a phantom value on the fit).
  if (method == "ipw") {
    valid_ipw_estimators <- c("km", "msm")
    if (length(.ipw_estimator) != 1L ||
        !.ipw_estimator %in% valid_ipw_estimators) {
      stop(".ipw_estimator must be one of: ",
           paste(shQuote(valid_ipw_estimators), collapse = ", "),
           ". Got: ", paste(.ipw_estimator, collapse = ", "),
           call. = FALSE)
    }
  } else {
    .ipw_estimator <- NULL
  }

  # ipcw: NULL → method-conditional default. "ipw" defaults to TRUE;
  # "gformula" and "km" default to FALSE. Under method = "km" with
  # ipcw = TRUE we still fit a C-hazard model (used for inverse-
  # probability-of-censoring weights), so TRUE is a valid user choice.
  if (is.null(ipcw)) {
    ipcw <- (method == "ipw")
  } else if (!is.logical(ipcw) || length(ipcw) != 1L || is.na(ipcw)) {
    stop("ipcw must be NULL or a single TRUE/FALSE.", call. = FALSE)
  }

  # stabilize (v1: NULL or "marginal"). Silent-coerce to NULL when no
  # weights are produced (i.e., method = "gformula", or method = "km"
  # with ipcw = FALSE). Avoids storing a misleading "marginal" on a fit
  # that did not actually stabilize anything.
  if (!is.null(stabilize) && !identical(stabilize, "marginal")) {
    stop("stabilize must be NULL (no stabilization) or \"marginal\" ",
         "(v1 only allows these two).", call. = FALSE)
  }
  if (method == "gformula" || (method == "km" && !ipcw)) {
    stabilize <- NULL
  }

  # formulas keys
  valid_formula_keys <- c("y", "c", "A", "A_num")
  if (!is.null(formulas)) {
    if (!is.list(formulas) || is.null(names(formulas)) ||
        any(names(formulas) == "")) {
      stop("`formulas` must be a named list (keys: ",
           paste(valid_formula_keys, collapse = ", "), ").", call. = FALSE)
    }
    bad <- setdiff(names(formulas), valid_formula_keys)
    if (length(bad) > 0L) {
      stop("Unknown formula key(s): ",
           paste(shQuote(bad), collapse = ", "),
           ". Valid keys: ",
           paste(shQuote(valid_formula_keys), collapse = ", "), ".",
           call. = FALSE)
    }
    # method = "km" only fits a C-hazard model (when ipcw = TRUE);
    # never a Y-hazard or propensity model. Reject formula overrides
    # that have no model to attach to.
    if (method == "km") {
      bad_km <- intersect(names(formulas), c("y", "A", "A_num"))
      if (length(bad_km) > 0L) {
        stop("method = \"km\" does not fit Y-, A-, or A-numerator models. ",
             "Unsupported formula key(s): ",
             paste(shQuote(bad_km), collapse = ", "), ". ",
             "Only `formulas$c` (used when ipcw = TRUE) is meaningful here.",
             call. = FALSE)
      }
    }
  }

  # truncate (length-2 percentiles in [0,1] with lower < upper)
  if (!is.null(truncate)) {
    if (!is.numeric(truncate) || length(truncate) != 2L ||
        any(is.na(truncate)) ||
        truncate[1] < 0 || truncate[2] > 1 ||
        truncate[1] >= truncate[2]) {
      stop("truncate must be NULL or c(lower, upper) percentiles in [0, 1] ",
           "with lower < upper.", call. = FALSE)
    }
  }

  # Person-time covariates: method-conditional consumption. fit_km
  # ignores them by design (KM is marginal in L_0). The validator just
  # carries them through; downstream workers decide what to do.
  covariates_vec <- attr(pt_data, "covariates")

  list(pt_data = pt_data, method = method, formulas = formulas,
       truncate = truncate, ipcw = ipcw, stabilize = stabilize,
       verbose = verbose, keep_data = keep_data,
       .ipw_estimator = .ipw_estimator, call = cl,
       # Person-time metadata pulled off pt_data attrs (once, for downstream)
       cut_times        = attr(pt_data, "cut_times"),
       id_col           = attr(pt_data, "id_col"),
       treatment_col    = attr(pt_data, "treatment_col"),
       covariates_vec   = covariates_vec,
       treatment_levels = attr(pt_data, "treatment_levels"))
}


#' Collect glm warnings from a fit expression
#'
#' Evaluates `expr` under a calling handler that muffles `warning()` calls
#' into a vector and re-emits a single grouped notice at the end. Returns
#' `expr`'s value with `$warnings` attached. Used by [causal_survival()] so
#' inner-fitter warnings surface as `fit$warnings` rather than streaming
#' mid-run.
#'
#' Works because R promises evaluate `expr` lazily — the expression is
#' forced inside the `withCallingHandlers` frame, so handlers are active
#' when the workers run.
#'
#' @keywords internal
.with_collected_warnings <- function(expr, context = "causal_survival()") {
  collected <- character()
  out <- withCallingHandlers(
    expr,
    warning = function(w) {
      collected <<- c(collected, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (length(collected) > 0L) {
    warning(context, ": ", length(collected),
            " warning(s) collected during fit. See fit$warnings for details.",
            call. = FALSE)
  }
  out$warnings <- collected
  out
}


#' Assemble the `causal_survival_fit` S3 object
#'
#' Packs worker output and canonical args (including person-time metadata
#' previously unpacked by the validator) into the fit list and stamps the
#' S3 class.
#'
#' @keywords internal
.assemble_causal_survival_fit <- function(worker_out, args) {
  ci_list <- list(gformula = NULL, ipw = NULL, km = NULL)
  ci_list[[args$method]] <- worker_out$estimates

  fit <- list(
    call                 = args$call,
    method               = args$method,
    ipw_estimator        = if (args$method == "ipw") args$.ipw_estimator else NULL,
    cumulative_incidence = ci_list,
    weights              = worker_out$weights,
    models               = worker_out$models,
    model_checks         = worker_out$model_checks,
    model_diagnostics    = NULL,
    warnings             = worker_out$warnings,
    pt_data              = if (args$keep_data) args$pt_data else NULL,
    cut_times            = args$cut_times,
    treatment_levels     = args$treatment_levels,
    id_col               = args$id_col,
    treatment_col        = args$treatment_col,
    covariates           = args$covariates_vec,
    stabilize            = args$stabilize,
    ipcw                 = args$ipcw,
    truncate             = args$truncate
  )

  class(fit) <- "causal_survival_fit"
  fit
}
