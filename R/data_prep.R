#' Standardize Treatment to `{0, 1}`
#'
#' Coerces a 2-level treatment column to integer `{0, 1}` and returns
#' both the new vector and the original level mapping. The mapping is
#' stashed on the person-time output so display layers (`print`,
#' `summary`, `plot`) can relabel arms back to the user's original
#' codes.
#'
#' Mapping rule by input type:
#' - factor: `levels(.)[1]` -> 0, `levels(.)[2]` -> 1
#' - character: sorted unique values; lower -> 0, higher -> 1
#' - logical: `FALSE` -> 0, `TRUE` -> 1
#' - numeric: must already be `{0, 1}`; passed through unchanged
#'
#' @param x A 2-level treatment vector.
#' @return A list with two elements:
#'   - `values`: integer vector in `{0, 1}` of the same length as `x`.
#'   - `levels`: length-2 character vector. `levels[1]` is the original
#'     label for `A = 0`; `levels[2]` for `A = 1`.
#' @keywords internal
standardize_treatment <- function(x) {

  if (is.factor(x)) {
    if (length(levels(x)) != 2) {
      stop("Treatment factor must have exactly 2 levels. Got: ",
           length(levels(x)), call. = FALSE)
    }
    list(
      values = as.integer(x) - 1L,
      levels = levels(x)
    )

  } else if (is.character(x)) {
    lv <- sort(unique(x))
    if (length(lv) != 2) {
      stop("Treatment character vector must have exactly 2 unique ",
           "values. Got: ", length(lv), call. = FALSE)
    }
    list(
      values = as.integer(match(x, lv)) - 1L,
      levels = lv
    )

  } else if (is.logical(x)) {
    list(
      values = as.integer(x),
      levels = c("FALSE", "TRUE")
    )

  } else if (is.numeric(x)) {
    u <- unique(x)
    if (!setequal(u, c(0, 1))) {
      stop("Numeric treatment must be coded as {0, 1}. Found: ",
           paste(u, collapse = ", "), call. = FALSE)
    }
    list(
      values = as.integer(x),
      levels = c("0", "1")
    )

  } else {
    stop("Unsupported treatment type: ", class(x)[1],
         ". Use factor, character, numeric, or logical.",
         call. = FALSE)
  }
}


#' Prepare Subject-Level Data for Causal Survival Analysis
#'
#' Converts subject-level data (one row per subject) into discrete-time
#' person-time format using [survSplit()]. Derives the
#' three-way censoring split (`y_event`, `cond_indep_cens`, `indep_cens`) on
#' terminal rows, standardizes treatment to `{0, 1}`, materializes the
#' admin-truncation convention, and attaches metadata for downstream
#' estimation.
#'
#' Schema and structural ordering are LOCKED — see
#' `dev/CAUSAL_SURVIVAL_SPEC.md` §3.0.
#'
#' @param data A subject-level data.frame (one row per subject).
#' @param id Character. Subject identifier column.
#' @param time Character. Event/censoring time column. Must be strictly
#'   positive — `time = 0` events are a hard error (no home interval
#'   under the `(0, t_1]` convention).
#' @param status Character. Binary event indicator column
#'   (`1` = event, `0` = censored).
#' @param ipcw Logical, scalar OR length-`nrow(data)` vector.
#'   `TRUE` = the subject's censoring is conditionally-independent —
#'   independent of the outcome only after conditioning on covariates
#'   `L` and treatment `A` (LTFU / treatment-switch / etc.). It
#'   contributes to the c-hazard fit and gets IPCW-weighted. `FALSE` =
#'   independent without conditioning (administrative end-of-study is
#'   the usual example) — excluded from the c-hazard fit and assigned
#'   weight 1. Label by the independence assumption you are willing to
#'   make, not by the administrative reason. This is the user's
#'   a-priori labeling of each subject's censoring mechanism; never
#'   inferred from data. Honored only when `status = 0` AND
#'   `time <= T_max`; otherwise ignored. Default `TRUE`.
#' @param T_max Numeric scalar, end of the analyzable time grid. The
#'   reporting grid is `(0, T_max]` partitioned into `K_max` intervals
#'   (see `cut_points`). Subjects with `time > T_max` are
#'   administratively truncated (no exit row, contribute at-risk rows
#'   up to `k = K_max`). Defaults to `max(data[[time]])`. Hard error
#'   if `T_max > max(data[[time]])` (would induce empty trailing
#'   intervals).
#' @param treatment Character. 2-level treatment column (factor,
#'   character, logical, or numeric `{0, 1}`).
#' @param covariates Character vector. Baseline covariate columns.
#' @param cut_points Time grid specification over `(0, T_max]`:
#'   - `NULL` (default): 12 equi-spaced intervals.
#'   - Single positive integer: that many equi-spaced intervals.
#'   - Numeric vector of length >= 2: explicit interior cut points,
#'     strictly within `(0, T_max)`.
#' @param time_varying Reserved for future use. Must be `NULL` in v1.
#' @param ... Unused. Reserved to catch misdirected
#'   CausalCompetingRisks-style calls with a clear redirect error when
#'   both packages are attached and this function masks the other.
#'
#' @return A data.frame of class `c("person_time", "data.frame")` with
#'   columns `id, k, A, <covariates>, y_event, cond_indep_cens, indep_cens`
#'   and attributes `cut_times`, `T_max`, `K_max`, `treatment_levels`,
#'   `id_col`, `treatment_col`, `covariates`. Per row, at most one of
#'   `{y_event, cond_indep_cens, indep_cens}` is `1` (exit row); otherwise all
#'   three are `0` (at-risk row, including the final row of admin-
#'   truncated subjects).
#'
#' @details
#' ## Time grid (LOCKED §3.0.2)
#' Continuous time over `(0, T_max]` is partitioned into `K_max`
#' analyzable intervals. Intervals are left-open right-closed: a
#' subject's event time `t` lands in interval `k` iff
#' `t_{k-1} < t <= t_k`. `k = 0` is pre-baseline and not in the data.
#' Estimates are reported at `k = 1, ..., K_max`.
#'
#' Boundary handling at `t = T_max`: `time = T_max` is inside the last
#' interval `(t_{K_max-1}, T_max]` (events fire normally at
#' `k = K_max`). `time > T_max` is outside all intervals — those
#' subjects are admin-truncated.
#'
#' ## Structural ordering (LOCKED §3.0.1)
#' Within each interval: `C_admin -> C_dep -> Y`. Admin censoring is
#' placed first within the interval — subjects with `time > T_max`
#' have no exit row regardless of `status`.
#'
#' ## Row encoding (LOCKED §3.0.4)
#'
#' | input                                              | indicator        |
#' |----------------------------------------------------|------------------|
#' | `status = 1, time <= T_max`                        | `y_event = 1`    |
#' | `status = 1, time > T_max`                         | (no exit row)    |
#' | `status = 0, time <= T_max, ipcw[i] = TRUE`        | `cond_indep_cens = 1`   |
#' | `status = 0, time <= T_max, ipcw[i] = FALSE`       | `indep_cens = 1` |
#' | `status = 0, time > T_max`                         | (no exit row)    |
#' | at-risk (incl. admin-truncated subjects)           | all `0`          |
#'
#' ## Diagnostics
#' - **Hard error** if `T_max > max(data[[time]])` — empty trailing
#'   intervals.
#' - **Hard error** if any `status = 1` row has `time = 0` — no home
#'   interval. Lists affected subject ids.
#' - **Warning** if `mean(admin_reach) < 0.5` after encoding — "hazard
#'   at K_max from thin risk set; CIF at K_max unreliable". Fit
#'   proceeds.
#'
#' ## Treatment standardization
#' The treatment column is coerced to integer `{0, 1}` via
#' [standardize_treatment()]. The original level mapping is stashed as
#' `attr(<output>, "treatment_levels")` for later display by accessors
#' (`print`, `summary`, `plot`).
#'
#' ## Pre-split mode
#' Dropped from v1 (see spec §3.0.9). Users with already-classified
#' censoring must run their classification through `status` + `ipcw`
#' (per-subject vector). The function always discretizes from
#' subject-level input.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   id = 1:100, time = rexp(100), status = rbinom(100, 1, 0.6),
#'   A = sample(c("ctrl", "trt"), 100, TRUE),
#'   age = rnorm(100, 60, 10)
#' )
#' pt <- to_person_time(
#'   df, id = "id", time = "time", status = "status",
#'   treatment = "A", covariates = "age", cut_points = 10
#' )
#' attr(pt, "treatment_levels")  # c("ctrl", "trt")
#'
#' # Per-subject censoring classification:
#' ipcw_vec <- df$reason %in% c("ltfu", "switch")  # TRUE = cond_indep_cens
#' pt2 <- to_person_time(df, id = "id", time = "time", status = "status",
#'                       treatment = "A", ipcw = ipcw_vec)
#' }
#'
#' @seealso [validate_subject_level()] for the input contract.
#' @export
to_person_time <- function(data,
                           id = "id",
                           time = "time",
                           status = "status",
                           ipcw = TRUE,
                           T_max = NULL,
                           treatment = "A",
                           covariates = character(),
                           cut_points = NULL,
                           time_varying = NULL,
                           ...) {

  # 0. Masking guard: with CausalCompetingRisks also attached, the LAST
  #    attached package masks the other's to_person_time(). Catch the
  #    sibling's arguments -> clear redirect instead of R's bare
  #    "unused argument".
  if (...length() > 0) {
    dots     <- names(list(...))
    ccr_args <- intersect(dots, c("event", "event_y", "event_d", "event_c"))
    if (length(ccr_args) > 0) {
      stop("Argument(s) ", paste(ccr_args, collapse = ", "),
           " belong to CausalCompetingRisks::to_person_time(). This call ",
           "reached CausalSurvival's to_person_time() (attached last, so ",
           "it masks the other). Call ",
           "CausalCompetingRisks::to_person_time() explicitly.",
           call. = FALSE)
    }
    stop("unused argument(s): ",
         paste(if (is.null(dots)) "<unnamed>" else dots, collapse = ", "),
         call. = FALSE)
  }

  # 1. Input shape: non-NULL, data.frame, > 0 rows.
  validate_input_shape(data, "data")

  # 2. Subject-level checks: required columns, id uniqueness, time /
  #    status / treatment / cut_points / covariate quality. Hard errors
  #    on structural problems that would corrupt the long-format output.
  validate_subject_level(
    data, id = id, time = time, status = status, treatment = treatment,
    covariates = covariates, cut_points = cut_points,
    time_varying = time_varying
  )

  # 3. Canonicalize args: resolve T_max default, reject time = 0 events,
  #    coerce ipcw to a length-nrow logical vector, resolve cut_points
  #    into cut_times over (0, T_max], and standardize treatment to
  #    integer {0, 1}. Returns a canonical args list with T_max,
  #    cut_times, K_max, ipcw_vec, and the trt mapping (values + original
  #    level labels).
  args <- .canonicalize_to_person_time_args(
            data, id, time, status, ipcw, T_max, treatment,
            covariates, cut_points, time_varying)

  # 4. Discretize subject-level data into person-time. Caps event times
  #    at T_max (admin-truncated subjects keep their at-risk rows up to
  #    k = K_max), runs survSplit() on cut_times, and tags
  #    each row with the integer interval index k under the
  #    (0, t_1], ..., (t_{K_max-1}, T_max] convention.
  pt <- .discretize_person_time(args)

  # 5. Assemble the three-way exit flag triple (y_event, cond_indep_cens,
  #    indep_cens) on terminal rows. Admin-truncated subjects' terminal
  #    row is forced back to at-risk (all three flags 0); otherwise the
  #    row splits by status + ipcw per spec §3.0.4.
  pt <- .assemble_exit_flags(pt, args)

  # 6. Drop survSplit scaffolding (tstart / tstop / event_indicator),
  #    the original time / status columns, and the broadcast helper
  #    columns; reorder to the canonical (id, k, treatment, covariates,
  #    y_event, cond_indep_cens, indep_cens) layout.
  pt <- .tidy_person_time_columns(pt, args)

  # 7. Diagnostic: warn when more than 50% of individuals are censored
  #    before reaching end of follow-up or having an event — i.e., the
  #    union of (subjects reaching k = K_max) and (subjects with an
  #    event at any k) covers less than 50% of the sample. Fit proceeds
  #    either way.
  .check_admin_reach(pt, args)

  # 8. Stamp metadata attributes (cut_times, T_max, K_max,
  #    treatment_levels, id_col, treatment_col, covariates) and the S3
  #    class c("person_time", "data.frame") onto the output.
  .stamp_person_time_class(pt, args)
}


# ----------------------------------------------------------------------------
# Internal helpers for to_person_time(). All plumbing (arg canonicalization,
# survSplit prep, flag assembly, column tidy, diagnostics, S3 stamping) lives
# here so the public API body above stays a clean methodological pipeline
# (one named verb per tile).
# ----------------------------------------------------------------------------

#' Canonicalize `to_person_time()` arguments
#'
#' Resolves `T_max` from `data[[time]]` when `NULL`, rejects any
#' `time = 0` rows (no home interval under the `(0, t_1]` convention),
#' coerces `ipcw` to a length-`nrow(data)` logical vector, resolves
#' `cut_points` to a sorted `cut_times` vector over `(0, T_max]`, and
#' standardizes the treatment column to integer `{0, 1}` (capturing the
#' original level mapping). Returns a canonical, normalized argument
#' list consumed by all downstream helpers in this file.
#'
#' Shape-level validators (`validate_input_shape()`,
#' `validate_subject_level()`) are called upstream in
#' [to_person_time()] before this helper runs.
#'
#' @keywords internal
.canonicalize_to_person_time_args <- function(data, id, time, status, ipcw,
                                              T_max, treatment, covariates,
                                              cut_points, time_varying) {

  # Resolve T_max
  max_time <- max(data[[time]])
  if (is.null(T_max)) {
    T_max <- max_time
  } else if (T_max > max_time) {
    stop("T_max (", T_max, ") exceeds max(data[[time]]) (", max_time,
         "). Would create empty trailing intervals.", call. = FALSE)
  }

  # Hard error: any time = 0 (no home interval under (0, t_1])
  zero_time <- which(data[[time]] == 0)
  if (length(zero_time) > 0) {
    ids <- data[[id]][zero_time]
    stop("time = 0 not supported (no home interval under (0, t_1] ",
         "convention; would be silently dropped by survSplit). ",
         "Affected ", id, ": ",
         paste(ids, collapse = ", "), ". ",
         "Drop or recode these subjects upstream.", call. = FALSE)
  }

  # Resolve ipcw shape
  if (length(ipcw) == 1L) {
    ipcw_vec <- rep(as.logical(ipcw), nrow(data))
  } else if (length(ipcw) == nrow(data)) {
    ipcw_vec <- as.logical(ipcw)
  } else {
    stop("ipcw must be a scalar logical or length-", nrow(data),
         " logical vector. Got length ", length(ipcw), ".",
         call. = FALSE)
  }
  if (any(is.na(ipcw_vec))) {
    stop("ipcw contains NA values. Must be TRUE/FALSE per subject.",
         call. = FALSE)
  }

  # Resolve cut_points -> cut_times over (0, T_max]
  if (is.null(cut_points)) {
    cut_times <- seq(0, T_max, length.out = 12L + 1L)[-1L]
  } else if (length(cut_points) == 1L) {
    cut_times <- seq(0, T_max, length.out = cut_points + 1L)[-1L]
  } else {
    cut_times <- sort(cut_points)
    if (cut_times[1] <= 0 || cut_times[length(cut_times)] >= T_max) {
      stop("Explicit cut_points must lie strictly within (0, T_max). ",
           "T_max = ", T_max, ".", call. = FALSE)
    }
    cut_times <- c(cut_times, T_max)
  }
  K_max <- length(cut_times)

  # Standardize treatment (raw -> integer {0, 1} + original level mapping)
  trt <- standardize_treatment(data[[treatment]])

  list(
    data       = data,
    id         = id,
    time       = time,
    status     = status,
    treatment  = treatment,
    covariates = covariates,
    T_max      = T_max,
    cut_times  = cut_times,
    K_max      = K_max,
    ipcw_vec   = ipcw_vec,
    trt        = trt
  )
}


#' Discretize subject-level data into person-time via `survSplit`
#'
#' Caps event times at `T_max` so admin-truncated subjects retain
#' at-risk rows up to `k = K_max` (their terminal `event_indicator = 1`
#' is suppressed later by [.assemble_exit_flags()]). Broadcasts the
#' per-subject `ipcw`, admin-truncation flag, and original status onto
#' helper columns so flag assembly can split rows by reason. Tags each
#' row with the integer interval index `k` under the
#' `(0, t_1], ..., (t_{K_max-1}, T_max]` convention.
#'
#' @keywords internal
.discretize_person_time <- function(args) {

  # Admin-truncated: survived past the window (time > T_max), OR censored
  # at exactly T_max. The boundary case has no later interval, so it must
  # stay at-risk for the final event instead of being excluded (which
  # collapses the last denominator). Reclassify and warn.
  past_window <- args$data[[args$time]] >  args$T_max
  at_boundary <- args$data[[args$time]] == args$T_max & args$data[[args$status]] == 0
  admin_trunc <- past_window | at_boundary

  if (any(at_boundary)) {
    warning(sum(at_boundary), " subjects censored at T_max reclassified as ",
            "administrative censoring; set T_max explicitly to override.",
            call. = FALSE)
  }

  # Prepare survSplit input: capped time, mapped treatment, scaffolding
  # columns, broadcast subject-level helpers.
  df <- args$data
  df[[args$treatment]]    <- args$trt$values
  df[[args$time]]         <- pmin(df[[args$time]], args$T_max)
  df$event_indicator      <- 1L
  df$tstart               <- 0
  df$tstop                <- df[[args$time]]
  df$.ipcw_subject        <- args$ipcw_vec
  df$.admin_trunc_subject <- admin_trunc
  df$.status_subject      <- df[[args$status]]

  pt <- survSplit(
    data  = df,
    cut   = args$cut_times,
    start = "tstart",
    end   = "tstop",
    event = "event_indicator"
  )

  pt$k <- findInterval(pt$tstop, c(0, args$cut_times), left.open = TRUE)
  pt
}


#' Build the three-way exit flag triple (`y_event`, `cond_indep_cens`, `indep_cens`)
#'
#' `event_indicator == 1` marks each subject's terminal row. For admin-
#' truncated subjects the terminal row is forced back to at-risk (all
#' three flags 0). Otherwise the row splits by `status` + `ipcw` per
#' spec §3.0.4: `status = 1` -> `y_event`; `status = 0` & `ipcw = TRUE`
#' -> `cond_indep_cens`; `status = 0` & `ipcw = FALSE` -> `indep_cens`.
#'
#' @keywords internal
.assemble_exit_flags <- function(pt, args) {
  is_terminal <- pt$event_indicator == 1L & !pt$.admin_trunc_subject
  is_event    <- is_terminal & pt$.status_subject == 1
  is_dep      <- is_terminal & pt$.status_subject == 0 &  pt$.ipcw_subject
  is_indep    <- is_terminal & pt$.status_subject == 0 & !pt$.ipcw_subject

  pt$y_event    <- as.integer(is_event)
  pt$cond_indep_cens   <- as.integer(is_dep)
  pt$indep_cens <- as.integer(is_indep)
  pt
}


#' Drop survSplit scaffolding + reorder to the canonical person-time layout
#'
#' Removes `tstart`, `tstop`, `event_indicator`, the original `time` and
#' `status` columns, and the three broadcast `.*_subject` helpers used
#' by [.assemble_exit_flags()]. Reorders the surviving columns to
#' `(id, k, treatment, covariates..., y_event, cond_indep_cens, indep_cens)`.
#'
#' @keywords internal
.tidy_person_time_columns <- function(pt, args) {
  pt$tstart                <- NULL
  pt$tstop                 <- NULL
  pt$event_indicator       <- NULL
  pt[[args$time]]          <- NULL
  pt[[args$status]]        <- NULL
  pt$.ipcw_subject         <- NULL
  pt$.admin_trunc_subject  <- NULL
  pt$.status_subject       <- NULL

  ordered_cols <- c(args$id, "k", args$treatment, args$covariates,
                    "y_event", "cond_indep_cens", "indep_cens")
  pt[, ordered_cols, drop = FALSE]
}


#' Warn when most individuals are censored before event or end of follow-up
#'
#' Computes the fraction of subjects who either (a) reach `k = K_max`
#' (with or without admin truncation) or (b) have an event
#' (`y_event == 1`) at any `k`. The complement of this union is the
#' fraction of subjects censored before reaching end of follow-up or
#' having an event — those contributing only partial information at
#' `k = K_max`. When the union covers less than 50% of subjects, the
#' function emits a warning; fit proceeds either way.
#'
#' @keywords internal
.check_admin_reach <- function(pt, args) {
  n_subjects     <- length(unique(args$data[[args$id]]))
  ids_at_kmax    <- unique(pt[[args$id]][pt$k == args$K_max])
  ids_with_event <- unique(pt[[args$id]][pt$y_event == 1])
  covered_frac   <- length(union(ids_at_kmax, ids_with_event)) / n_subjects

  if (covered_frac < 0.5) {
    warning(
      "More than 50% of individuals are censored before reaching ",
      "end of follow-up or having an event ",
      "(covered fraction = ",
      formatC(covered_frac, digits = 3, format = "f"), ").",
      call. = FALSE
    )
  }
  invisible(NULL)
}


#' Stamp `person_time` metadata attributes and S3 class
#'
#' Attaches `cut_times`, `T_max`, `K_max`, `treatment_levels`, `id_col`,
#' `treatment_col`, and `covariates` as attributes on `pt`, and prepends
#' `"person_time"` to its S3 class vector so downstream code can dispatch
#' on it via [inherits()].
#'
#' @keywords internal
.stamp_person_time_class <- function(pt, args) {
  attr(pt, "cut_times")        <- args$cut_times
  attr(pt, "T_max")            <- args$T_max
  attr(pt, "K_max")            <- args$K_max
  attr(pt, "treatment_levels") <- args$trt$levels
  attr(pt, "id_col")           <- args$id
  attr(pt, "treatment_col")    <- args$treatment
  attr(pt, "covariates")       <- args$covariates
  class(pt) <- c("person_time", class(pt))
  pt
}
