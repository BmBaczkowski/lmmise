# ================================
# AS_LME4
# ================================

#' Translate Spec to lme4 Arguments
#'
#' Returns a list with `formula` and `family` ready to pass to `lme4::lmer()`
#' or `lme4::glmer()`. Does not require params — structure only.
#'
#' @param spec An `lmm_spec` object.
#' @param include_random Logical; include random effects in formula (default TRUE).
#' @return A list with elements `formula` and `family`.
#' @export
as_lme4 <- function(spec, include_random = TRUE) {
  if (!inherits(spec, "lmm_spec")) stop("spec must be an lmm_spec object")
  if (is.null(spec$fixed$formula)) stop("spec has no fixed effects; call set_fixed() first")

  fixed_str <- paste(spec$fixed$terms, collapse = " + ")
  random_str <- if (include_random) .build_random_str(spec$random) else ""

  formula <- stats::as.formula(paste(spec$response$name, "~", fixed_str, random_str))
  family <- do.call(spec$response$dist, list(link = spec$response$link))

  list(formula = formula, family = family)
}


.build_random_str <- function(random) {
  if (length(random) == 0) {
    return("")
  }

  parts <- mapply(function(group, g) {
    term_names <- g$terms
    if (g$model_corr) {
      paste0("(", paste(term_names, collapse = " + "), " | ", group, ")")
    } else if (length(term_names) == 1) {
      paste0("(", term_names, " | ", group, ")")
    } else {
      re <- c()
      if ("1" %in% term_names) re <- c(re, paste0("(1 | ", group, ")"))
      slopes <- setdiff(term_names, "1")
      if (length(slopes) > 0) re <- c(re, paste0("(0 + ", slopes, " | ", group, ")"))
      paste(re, collapse = " + ")
    }
  }, names(random), random, SIMPLIFY = TRUE)

  paste("+", paste(parts, collapse = " + "))
}

# ================================
# AS_LMLIST
# ================================

#' Translate Spec to lmList Arguments
#'
#' Returns a list with `formula` and `family` ready to pass to `lme4::lmList()`.
#' The formula uses the `response ~ fixed | cluster` syntax expected by lmList.
#' Does not require params — structure only.
#'
#' @param spec An `lmm_spec` object.
#' @param cluster A single character string naming the grouping factor (must
#'   appear in `spec$random`).
#' @return A list with elements `formula` and `family`.
#' @export
as_lmList <- function(spec, cluster) {
  if (!inherits(spec, "lmm_spec")) stop("spec must be an lmm_spec object")
  if (is.null(spec$fixed$formula)) stop("spec has no fixed effects; call set_fixed() first")
  if (missing(cluster)) stop("cluster must be provided")
  if (!is.character(cluster) || length(cluster) != 1) {
    stop("cluster must be a single character string")
  }
  if (!cluster %in% names(spec$random)) {
    stop("cluster '", cluster, "' not found in spec random groups")
  }

  fixed_str <- paste(spec$fixed$terms, collapse = " + ")
  formula <- stats::as.formula(paste(spec$response$name, "~", fixed_str, "|", cluster))
  family <- do.call(spec$response$dist, list(link = spec$response$link))

  list(formula = formula, family = family)
}

# ================================
# AS_SIMSTUDY
# ================================

#' Translate Spec + Params to simstudy Arguments
#'
#' @param spec An `lmm_spec` object.
#' @param params An `lmm_params` object (fully set).
#' @param type One of `"rand_effects"` (args for `addCorData()`) or
#'   `"outcome"` (args for `defDataAdd()`).
#' @return A list of arguments for the corresponding simstudy function.
#' @export
as_simstudy <- function(spec, params, type = "rand_effects") {
  if (!inherits(spec, "lmm_spec")) stop("spec must be an lmm_spec object")
  if (!inherits(params, "lmm_params")) stop("params must be an lmm_params object")
  if (!type %in% c("rand_effects", "outcome")) {
    stop("type must be 'rand_effects' or 'outcome'")
  }

  .check_spec_params_compat(spec, params)

  missing <- check_params(params, quiet = TRUE)
  if (!is.null(missing)) {
    stop("params has unset values; run check_params() for details")
  }

  switch(type,
    rand_effects = .as_simstudy_rand_effects(spec, params),
    outcome      = .as_simstudy_outcome(spec, params)
  )
}


.as_simstudy_rand_effects <- function(spec, params) {
  # Returns a named list (one entry per random group) of argument lists
  # for simstudy::addCorData().
  #
  # Usage:
  #   args <- as_simstudy(spec, params, type = "rand_effects")
  #   do.call(addCorData, c(list(dtOld = dt, idname = "id"), args$id))

  lapply(names(spec$random), function(group) {
    g <- spec$random[[group]]
    term_names <- g$terms
    u_names <- g$u_names

    sd_names <- ifelse(term_names == "1", "intercept", term_names)
    sigma <- sapply(sd_names, function(t) params$random[[group]]$sd[[t]])

    corr <- params$random[[group]]$corr

    corr_args <- if (is.null(corr)) {
      list(rho = 0, corstr = "ind")
    } else {
      # always stored as a named matrix internally
      list(corMatrix = corr)
    }

    c(
      list(mu = rep(0, length(term_names)), sigma = unname(sigma), cnames = unname(u_names)),
      corr_args
    )
  }) |> stats::setNames(names(spec$random))
}


.as_simstudy_outcome <- function(spec, params) {
  # Returns an argument list for simstudy::defDataAdd().
  #
  # Usage:
  #   args <- as_simstudy(spec, params, type = "outcome")
  #   do.call(defDataAdd, args)

  # build term -> u_names lookup across all groups
  u_map <- list()
  for (group in names(spec$random)) {
    for (term in names(spec$random[[group]]$u_names)) {
      u_map[[term]] <- c(u_map[[term]], spec$random[[group]]$u_names[[term]])
    }
  }

  parts <- sapply(spec$fixed$terms, function(term) {
    beta_name <- if (term == "1") "intercept" else term
    beta <- params$beta[[beta_name]]
    u_names <- u_map[[term]]

    if (term == "1") {
      paste(c(beta, u_names), collapse = " + ")
    } else {
      coef_str <- if (!is.null(u_names)) {
        paste0("(", paste(c(beta, u_names), collapse = " + "), ")")
      } else {
        as.character(beta)
      }
      paste0(coef_str, " * ", term)
    }
  })

  dist <- switch(spec$response$dist,
    gaussian = "normal",
    spec$response$dist
  )

  list(
    varname  = spec$response$name,
    formula  = paste(parts, collapse = " + "),
    dist     = dist,
    link     = spec$response$link,
    variance = params$dispersion
  )
}


# ================================
# COMPATIBILITY CHECK
# ================================

.check_spec_params_compat <- function(spec, params) {
  expected_beta <- sort(ifelse(spec$fixed$terms == "1", "intercept", spec$fixed$terms))
  actual_beta <- sort(names(params$beta))
  if (!identical(expected_beta, actual_beta)) {
    stop(
      "params beta terms (", paste(actual_beta, collapse = ", "),
      ") do not match spec (", paste(expected_beta, collapse = ", "), ")"
    )
  }

  expected_groups <- sort(names(spec$random))
  actual_groups <- sort(names(params$random))
  if (!identical(expected_groups, actual_groups)) {
    stop(
      "params random groups (", paste(actual_groups, collapse = ", "),
      ") do not match spec (", paste(expected_groups, collapse = ", "), ")"
    )
  }

  for (group in names(spec$random)) {
    expected_terms <- sort(ifelse(spec$random[[group]]$terms == "1", "intercept", spec$random[[group]]$terms))
    actual_terms <- sort(names(params$random[[group]]$sd))
    if (!identical(expected_terms, actual_terms)) {
      stop("random terms for group '", group, "' differ between spec and params")
    }
  }
}
