#' Validate Basic Input Data Frame Shape
#'
#' Checks that an input object is a non-empty data.frame.
#'
#' @param x The object to check.
#' @param name Character. Name of the argument (for error messages).
#'
#' @return Invisibly TRUE if all checks pass; otherwise throws a helpful
#'   error.
#' @export
#' @keywords internal
validate_input_shape <- function(x, name) {
  if (is.null(x)) {
    stop(name, " is NULL.", call. = FALSE)
  }
  if (!is.data.frame(x)) {
    stop(
      name, " must be a data.frame. Got class: ",
      paste(class(x), collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(x) == 0) {
    stop(name, " has 0 rows.", call. = FALSE)
  }
  invisible(TRUE)
}


#' Shared Covariate Checks
#'
#' Hard errors on (a) NAs, (b) unsupported types (anything other than
#' numeric / integer / factor / character / logical — these break glm
#' downstream). Warnings on (c) constant columns, (d) high-cardinality
#' factors.
#'
#' @param data A data.frame.
#' @param covariates Character vector of covariate column names.
#' @return Invisibly NULL.
#' @export
#' @keywords internal
check_covariate_quality <- function(data, covariates) {

  errors        <- character()
  warnings_msgs <- character()

  for (cov in covariates) {
    col <- data[[cov]]
    col_class <- class(col)[1]

    # --- Type (ERROR) ---
    if (!col_class %in% c("numeric", "integer", "factor",
                          "character", "logical")) {
      errors <- c(errors,
        paste0("  - '", cov, "': unsupported type '", col_class,
               "' (expected numeric, integer, factor, character, or logical)")
      )
      next
    }

    # --- NAs (ERROR) ---
    if (any(is.na(col))) {
      errors <- c(errors,
        paste0("  - '", cov, "': contains ", sum(is.na(col)), " NA value(s)")
      )
    }

    # --- Constant (WARNING) ---
    n_unique_nona <- length(unique(col[!is.na(col)]))
    if (n_unique_nona < 2) {
      warnings_msgs <- c(warnings_msgs,
        paste0("  - '", cov, "': fewer than 2 unique values")
      )
    }

    # --- Cardinality, categorical only (WARNING) ---
    if (col_class %in% c("factor", "character")) {
      n_levels <- if (col_class == "factor") length(levels(col))
                  else length(unique(col))
      if (n_levels > 20) {
        warnings_msgs <- c(warnings_msgs,
          paste0("  - '", cov, "': ", n_levels, " levels ",
                 "(high cardinality - consider grouping)")
        )
      }
    }
  }

  if (length(errors) > 0) {
    stop(
      "Covariate validation failed:\n",
      paste(errors, collapse = "\n"),
      call. = FALSE
    )
  }
  if (length(warnings_msgs) > 0) {
    warning(
      "Covariate quality issues:\n",
      paste(warnings_msgs, collapse = "\n"),
      call. = FALSE
    )
  }

  invisible(NULL)
}


#' Validate Subject-Level Input
#'
#' Checks subject-level data (one row per subject) before person-time
#' construction in [to_person_time()]. Hard-errors on structural problems
#' that would corrupt the long-format output or downstream estimation.
#'
#' @param data A data.frame in subject-level format.
#' @param id,time,status,treatment Character. Column names.
#' @param covariates Character vector of baseline covariate column names.
#' @param cut_points `NULL`, a length-1 positive integer (number of
#'   equi-spaced intervals), or a numeric vector of length >= 2 (explicit
#'   interior cut points strictly within `(0, max(time))`).
#' @param time_varying Reserved for future use. Must be `NULL` in v1.
#'
#' @return Invisibly `TRUE` on success.
#'
#' @details
#' ## Hard errors
#' - Missing required columns
#' - `id`: NAs or duplicates (one row per subject)
#' - `time`: not numeric, negative, or NAs
#' - `status`: not in `{0, 1}`, or NAs
#' - `treatment`: NAs, fewer than 2 levels, or more than 2 levels
#' - covariates: NAs or unsupported types (via [check_covariate_quality()])
#' - `cut_points`: invalid form, non-integer scalar, non-increasing
#'   vector, or values outside `(0, max(time))`
#' - `time_varying != NULL`
#'
#' @keywords internal
validate_subject_level <- function(data, id, time, status, treatment,
                                   covariates, cut_points, time_varying) {

  # --- v1 block ---
  if (!is.null(time_varying)) {
    stop("time_varying covariates are not supported in v1.", call. = FALSE)
  }

  # --- Required columns present ---
  required <- c(id, time, status, treatment, covariates)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop("Subject-level data is missing required column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  # --- id: no NAs, no duplicates ---
  if (any(is.na(data[[id]]))) {
    stop("Column '", id, "' contains NA values.", call. = FALSE)
  }
  if (anyDuplicated(data[[id]]) > 0) {
    stop("Column '", id, "' contains duplicate values. ",
         "Subject-level data must have one row per subject.",
         call. = FALSE)
  }

  # --- time: numeric, non-negative, no NAs ---
  t_vals <- data[[time]]
  if (!is.numeric(t_vals)) {
    stop("Column '", time, "' must be numeric. Got: ",
         class(t_vals)[1], call. = FALSE)
  }
  if (any(is.na(t_vals))) {
    stop("Column '", time, "' contains NA values.", call. = FALSE)
  }
  if (any(t_vals < 0)) {
    stop("Column '", time, "' contains negative values.", call. = FALSE)
  }

  # --- status: {0, 1}, no NAs ---
  s_vals <- data[[status]]
  if (any(is.na(s_vals))) {
    stop("Column '", status, "' contains NA values.", call. = FALSE)
  }
  if (!all(s_vals %in% c(0, 1))) {
    stop("Column '", status, "' must contain only 0 and 1. ",
         "Found: ", paste(unique(s_vals), collapse = ", "),
         call. = FALSE)
  }

  # --- treatment: NAs + 2 unique levels ---
  if (any(is.na(data[[treatment]]))) {
    stop("Column '", treatment, "' contains NA values.", call. = FALSE)
  }
  trt_vals <- unique(data[[treatment]])
  if (length(trt_vals) > 2) {
    stop("Column '", treatment, "' has more than 2 unique values. ",
         "Multi-arm treatment is not supported in v1. Found: ",
         paste(trt_vals, collapse = ", "), call. = FALSE)
  }
  if (length(trt_vals) < 2) {
    stop("Column '", treatment, "' has fewer than 2 unique values. ",
         "Cannot estimate counterfactuals without both arms.",
         call. = FALSE)
  }

  # --- cut_points: form + range ---
  if (!is.null(cut_points)) {
    if (length(cut_points) == 1) {
      if (!is.numeric(cut_points) || cut_points <= 0 ||
          cut_points != as.integer(cut_points)) {
        stop("cut_points scalar must be a positive integer (number of ",
             "equi-spaced intervals). Got: ", cut_points,
             call. = FALSE)
      }
    } else {
      if (!is.numeric(cut_points)) {
        stop("cut_points vector must be numeric.", call. = FALSE)
      }
      if (any(is.na(cut_points))) {
        stop("cut_points contains NA values.", call. = FALSE)
      }
      max_t <- max(t_vals)
      if (any(cut_points <= 0 | cut_points > max_t)) {
        stop("cut_points must lie within (0, max(", time, ")]. ",
             "max(", time, ") = ", max_t, "; ",
             "found: ", paste(cut_points, collapse = ", "),
             call. = FALSE)
      }
      if (any(diff(cut_points) <= 0)) {
        stop("cut_points must be strictly increasing.", call. = FALSE)
      }
    }
  }

  # --- covariate quality ---
  if (length(covariates) > 0) {
    check_covariate_quality(data, covariates)
  }

  invisible(TRUE)
}
