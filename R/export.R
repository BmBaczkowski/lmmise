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
  if (length(random) == 0L) {
    return("")
  }

  parts <- mapply(
    function(group, g) {
      term_names <- g$terms

      if (g$model_corr) {
        if ("1" %in% term_names) {
          rhs <- paste(term_names, collapse = " + ")
        } else {
          rhs <- paste(c("0", term_names), collapse = " + ")
        }

        paste0("(", rhs, " | ", group, ")")
      } else {
        re <- character()

        if ("1" %in% term_names) {
          re <- c(re, paste0("(1 | ", group, ")"))
        }

        slopes <- setdiff(term_names, "1")

        if (length(slopes) > 0L) {
          re <- c(
            re,
            paste0("(0 + ", slopes, " | ", group, ")")
          )
        }

        paste(re, collapse = " + ")
      }
    },
    names(random),
    random,
    SIMPLIFY = TRUE
  )

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
#' @param exclude_terms A character vector of variable names to exclude from the
#'   formula. Interaction terms containing any excluded variable will also be
#'   removed. Comma-separated strings are automatically split. For example,
#'   use \code{c("x1", "x2")} or \code{"x1, x2"} to exclude \code{x1} and
#'   \code{x2} (and any interactions involving them).
#' @return A list with elements `formula` and `family`.
#' @export
as_lmList <- function(spec, cluster, exclude_terms = character(0)) {
  if (!inherits(spec, "lmm_spec")) stop("spec must be an lmm_spec object")
  if (is.null(spec$fixed$formula)) stop("spec has no fixed effects; call set_fixed() first")
  if (missing(cluster)) stop("cluster must be provided")
  if (!is.character(cluster) || length(cluster) != 1) {
    stop("cluster must be a single character string")
  }
  if (!cluster %in% names(spec$random)) {
    stop("cluster '", cluster, "' not found in spec random groups")
  }

  if (length(exclude_terms) > 0) {
    # Allow comma-separated input for convenience
    exclude_terms <- unique(trimws(unlist(strsplit(exclude_terms, ","))))
    spec$fixed$terms <- spec$fixed$terms[sapply(spec$fixed$terms, function(term) {
      components <- unlist(strsplit(term, ":"))
      !any(components %in% exclude_terms)
    })]
  }
  if (length(spec$fixed$terms) == 0) {
    stop("No fixed terms remaining after exclusions")
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

    term_clean <- gsub(":", "*", term) # Replace : with * for interactions

    if (term == "1") {
      paste(c(beta, u_names), collapse = " + ")
    } else {
      coef_str <- if (!is.null(u_names)) {
        paste0("(", paste(c(beta, u_names), collapse = " + "), ")")
      } else {
        as.character(beta)
      }
      paste0(coef_str, " * ", term_clean)
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

#' Translate Spec to metafor Arguments
#'
#' Returns model components suitable for constructing a multilevel model
#' with metafor. Variables appearing in `mods` are treated as level-2
#' moderators and are excluded from the level-1 formula.
#'
#' @param spec An `lmm_spec` object.
#' @param mods A one-sided formula specifying level-2 moderators, e.g.
#'   \code{~ x2}, \code{~ x2 + x3}, or \code{~ x2 * treatment}.
#'   Defaults to \code{~ 1}.
#'
#' @return A list with elements `lme4_like`, `formula`, `cluster`, `mods`,
#'   and `family`.
#' @export
as_metafor <- function(spec, mods = ~ 1) {
  if (!inherits(spec, "lmm_spec")) {
    stop("spec must be an lmm_spec object")
  }

  if (is.null(spec$fixed$formula)) {
    stop("spec has no fixed effects; call set_fixed() first")
  }

  # Validate moderator formula
  if (!inherits(mods, "formula")) {
    stop("mods must be a formula, e.g. ~ x2 or ~ x2 + x3")
  }

  if (length(mods) != 2L) {
    stop("mods must be a one-sided formula, e.g. ~ x2")
  }

  mod_vars <- all.vars(mods)

  # Moderator variables must appear in the fixed-effects specification
  fixed_vars <- unique(unlist(
    lapply(spec$fixed$terms, function(term) {
      if (term == "1") {
        return(character(0))
      }

      all.vars(stats::as.formula(paste("~", term)))
    })
  ))

  missing_mods <- setdiff(mod_vars, fixed_vars)

  if (length(missing_mods) > 0L) {
    stop(
      "Moderator variable",
      if (length(missing_mods) > 1L) "s " else " ",
      paste0("'", missing_mods, "'", collapse = ", "),
      if (length(missing_mods) > 1L) " are" else " is",
      " not present in the fixed-effects formula"
    )
  }

  # Infer cluster from spec$random
  random_groups <- names(spec$random)

  if (length(random_groups) == 0L) {
    stop("spec has no random-effects grouping factor")
  }

  if (length(random_groups) > 1L) {
    stop(
      "as_metafor() currently requires exactly one random-effects grouping factor; ",
      "found: ",
      paste(random_groups, collapse = ", ")
    )
  }

  cluster <- random_groups[[1L]]

  # Full lme4-style formula
  lme4_like <- as_lme4(spec, include_random = TRUE)$formula

  # Remove terms involving level-2 moderators
  level1_terms <- spec$fixed$terms[
    vapply(spec$fixed$terms, function(term) {
      if (term == "1") {
        return(TRUE)
      }

      term_vars <- all.vars(
        stats::as.formula(paste("~", term))
      )

      !any(term_vars %in% mod_vars)
    }, logical(1))
  ]

  if (length(level1_terms) == 0L) {
    stop("No level-1 fixed terms remaining after removing moderators")
  }

  formula <- stats::as.formula(
    paste(
      spec$response$name,
      "~",
      paste(level1_terms, collapse = " + ")
    )
  )

  family <- do.call(
    spec$response$dist,
    list(link = spec$response$link)
  )

  list(
    lme4_like = lme4_like,
    formula = formula,
    cluster = cluster,
    mods = mods,
    family = family
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
