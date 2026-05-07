# ================================
# PARAMS CONSTRUCTOR
# ================================

#' Create an Empty LMM Params Object
#'
#' Inspects an `lmm_spec` and returns a params object with all required slots
#' set to NULL. Prints a summary of what needs to be filled via `patch()`.
#'
#' @param spec An `lmm_spec` object.
#' @return An object of class `lmm_params`.
#' @export
lmm_params <- function(spec) {
  if (!inherits(spec, "lmm_spec")) stop("spec must be an lmm_spec object")
  if (is.null(spec$fixed$formula)) stop("spec has no fixed effects; call set_fixed() first")

  # beta slots: named by term ("1" -> "intercept")
  beta_names <- ifelse(spec$fixed$terms == "1", "intercept", spec$fixed$terms)
  beta <- stats::setNames(vector("list", length(beta_names)), beta_names)

  # dispersion: fixed for binomial/poisson, NULL otherwise
  dispersion <- switch(spec$response$dist,
    gaussian = NULL,
    gamma    = NULL,
    poisson  = 1L,
    binomial = 1L
  )

  # random slots: sd per term, corr NULL (optional)
  random <- lapply(spec$random, function(g) {
    sd_names <- ifelse(g$terms == "1", "intercept", g$terms)
    list(
      sd   = stats::setNames(vector("list", length(sd_names)), sd_names),
      corr = NULL
    )
  })

  params <- structure(
    list(spec = spec, beta = beta, dispersion = dispersion, random = random),
    class = "lmm_params"
  )

  invisible(params)
}


# ================================
# PATCH
# ================================

#' Patch Parameter Values
#'
#' Sets or updates values in an `lmm_params` object. Validates keys and values
#' immediately. Can be used both for initial fill and updates.
#'
#' Accepted keys in the update list:
#' - `beta`: named list of fixed effect values, e.g. `list(intercept = 0, x1 = 0.5)`
#' - `dispersion`: single positive numeric
#' - `random_sd`: named list per group, e.g. `list(id = list(intercept = 1, x1 = 0.3))`
#' - `random_corr`: named list per group with either:
#'     - `list(structure = "cs", r = 0.3)` for shorthand structures
#'     - `list(matrix = m)` for a raw correlation matrix
#'
#' @param params An `lmm_params` object.
#' @param values A named list of updates (see Details).
#' @return An updated `lmm_params` object.
#' @export
patch <- function(params, values) {
  if (!inherits(params, "lmm_params")) stop("params must be an lmm_params object")
  if (!is.list(values) || is.null(names(values)) || any(names(values) == "")) {
    stop("values must be a named list")
  }

  unknown <- setdiff(names(values), c("beta", "dispersion", "random_sd", "random_corr"))
  if (length(unknown) > 0) {
    stop(
      "unknown key(s): ", paste(unknown, collapse = ", "),
      ". Expected: beta, dispersion, random_sd, random_corr"
    )
  }

  if (!is.null(values$beta)) params <- .patch_beta(params, values$beta)
  if (!is.null(values$dispersion)) params <- .patch_dispersion(params, values$dispersion)
  if (!is.null(values$random_sd)) params <- .patch_random_sd(params, values$random_sd)
  if (!is.null(values$random_corr)) params <- .patch_random_corr(params, values$random_corr)

  params
}


.patch_beta <- function(params, beta) {
  if (!is.list(beta) || is.null(names(beta)) || any(names(beta) == "")) {
    stop("beta must be a named list")
  }

  unknown <- setdiff(names(beta), names(params$beta))
  if (length(unknown) > 0) {
    stop(
      "unknown beta term(s): ", paste(unknown, collapse = ", "),
      ". Available: ", paste(names(params$beta), collapse = ", ")
    )
  }

  for (nm in names(beta)) {
    if (!is.numeric(beta[[nm]]) || length(beta[[nm]]) != 1) {
      stop("beta$", nm, " must be a single numeric value")
    }
    params$beta[[nm]] <- beta[[nm]]
  }
  params
}


.patch_dispersion <- function(params, value) {
  if (params$spec$response$dist %in% c("binomial", "poisson")) {
    stop("dispersion is fixed at 1 for '", params$spec$response$dist, "'")
  }
  if (!is.numeric(value) || length(value) != 1 || value <= 0) {
    stop("dispersion must be a single positive numeric value")
  }
  params$dispersion <- value
  params
}


.patch_random_sd <- function(params, random_sd) {
  if (!is.list(random_sd) || is.null(names(random_sd))) {
    stop("random_sd must be a named list")
  }

  unknown_groups <- setdiff(names(random_sd), names(params$random))
  if (length(unknown_groups) > 0) {
    stop(
      "unknown random group(s): ", paste(unknown_groups, collapse = ", "),
      ". Available: ", paste(names(params$random), collapse = ", ")
    )
  }

  for (group in names(random_sd)) {
    sd_vals <- random_sd[[group]]
    if (!is.list(sd_vals) || is.null(names(sd_vals))) {
      stop("random_sd$", group, " must be a named list")
    }

    valid_terms <- names(params$random[[group]]$sd)
    unknown_terms <- setdiff(names(sd_vals), valid_terms)
    if (length(unknown_terms) > 0) {
      stop(
        "unknown term(s) in random_sd$", group, ": ",
        paste(unknown_terms, collapse = ", "),
        ". Available: ", paste(valid_terms, collapse = ", ")
      )
    }

    for (term in names(sd_vals)) {
      val <- sd_vals[[term]]
      if (!is.numeric(val) || length(val) != 1 || val < 0) {
        stop("random_sd$", group, "$", term, " must be a single non-negative numeric")
      }
      params$random[[group]]$sd[[term]] <- val
    }
  }
  params
}


.patch_random_corr <- function(params, random_corr) {
  if (!is.list(random_corr) || is.null(names(random_corr))) {
    stop("random_corr must be a named list")
  }

  unknown_groups <- setdiff(names(random_corr), names(params$random))
  if (length(unknown_groups) > 0) {
    stop(
      "unknown random group(s): ", paste(unknown_groups, collapse = ", "),
      ". Available: ", paste(names(params$random), collapse = ", ")
    )
  }

  for (group in names(random_corr)) {
    spec_g <- params$spec$random[[group]]
    term_names <- ifelse(spec_g$terms == "1", "intercept", spec_g$terms)

    if (length(term_names) < 2) {
      stop("random_corr for group '", group, "' requires at least 2 random terms")
    }

    corr_spec <- random_corr[[group]]
    m <- .make_corr_matrix(term_names, corr_spec)
    params$random[[group]]$corr <- m
  }
  params
}


# ================================
# CHECK_PARAMS
# ================================

#' Check Parameter Completeness
#'
#' Prints which parameters are set and which are missing.
#' Call explicitly before running simulations to confirm everything is filled.
#'
#' @param params An `lmm_params` object.
#' @param quiet Logical; if TRUE returns missing list silently (default FALSE).
#' @return Invisibly returns a character vector of missing param names, or NULL if complete.
#' @export
check_params <- function(params, quiet = FALSE) {
  if (!inherits(params, "lmm_params")) stop("params must be an lmm_params object")

  result <- list()

  # Check beta parameters
  missing_beta <- character()
  for (nm in names(params$beta)) {
    if (is.null(params$beta[[nm]])) missing_beta <- c(missing_beta, nm)
  }
  result$beta <- if (length(missing_beta) > 0) missing_beta else NULL

  # Check dispersion
  if (params$spec$response$dist %in% c("binomial", "poisson")) {
    result$dispersion <- NULL
  } else {
    result$dispersion <- if (is.null(params$dispersion)) TRUE else NULL
  }

  # Check random parameters
  result$random <- lapply(names(params$random), function(group) {
    missing_sd <- character()
    for (term in names(params$random[[group]]$sd)) {
      if (is.null(params$random[[group]]$sd[[term]])) {
        missing_sd <- c(missing_sd, term)
      }
    }
    spec_g <- params$spec$random[[group]]
    n_terms <- length(spec_g$terms)
    missing_corr <- if (n_terms >= 2 && is.null(params$random[[group]]$corr)) {
      TRUE
    } else {
      NULL
    }
    list(
      sd = if (length(missing_sd) > 0) missing_sd else NULL,
      corr = missing_corr
    )
  })
  names(result$random) <- names(params$random)

  has_missing <- any(c(
    length(result$beta) > 0,
    isTRUE(result$dispersion),
    any(vapply(result$random, function(g) length(g$sd) > 0 || isTRUE(g$corr), logical(1)))
  ))

  if (!quiet) {
    if (!has_missing) {
      cli::cli_alert_success("All parameters set.")
    } else {
      cli::cli_alert_warning("Missing parameters:")
      if (length(result$beta) > 0) {
        cli::cli_ul("Fixed effects (beta):")
        for (m in result$beta) cli::cli_li("{m}")
        cli::cli_end()
      }
      if (isTRUE(result$dispersion)) cli::cli_li("Dispersion parameter")
      for (group in names(result$random)) {
        g <- result$random[[group]]
        if (length(g$sd) > 0 || isTRUE(g$corr)) {
          cli::cli_ul("Random effects ({group}):")
          if (length(g$sd) > 0) {
            for (m in g$sd) cli::cli_li("SD: {m}")
          }
          if (isTRUE(g$corr)) cli::cli_li("Correlation structure")
          cli::cli_end()
        }
      }
    }
  }

  invisible(if (!has_missing) NULL else result)
}
