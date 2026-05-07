# ================================
# SPEC CONSTRUCTOR
# ================================

#' Create an LMM Spec
#'
#' Defines the structure of a linear mixed model: response variable,
#' distribution, and link function. Use `set_fixed()` and `set_random()`
#' to add terms via pipe.
#'
#' @param response Character; name of the response variable (default "y").
#' @param dist Character; one of "gaussian", "binomial", "poisson", "gamma".
#' @param link Character; link function. Defaults to canonical link for dist.
#' @return An object of class `lmm_spec`.
#' @export
lmm_spec <- function(response = "y", dist = "gaussian", link = NULL) {
  supported_dists <- c("gaussian", "binomial", "poisson", "gamma")

  if (!is.character(response) || length(response) != 1) {
    stop("response must be a single character string")
  }
  if (!dist %in% supported_dists) {
    stop("dist must be one of: ", paste(supported_dists, collapse = ", "))
  }

  if (is.null(link)) {
    link <- switch(dist,
      gaussian = "identity",
      binomial = "logit",
      poisson  = "log",
      gamma    = "inverse"
    )
  }

  structure(
    list(
      response = list(name = response, dist = dist, link = link),
      fixed    = list(formula = NULL, terms = character()),
      random   = list()
    ),
    class = "lmm_spec"
  )
}


# ================================
# SET_FIXED
# ================================

#' Set Fixed Effects
#'
#' @param spec An `lmm_spec` object.
#' @param formula A formula, e.g. `y ~ 1 + x1`. The response name overrides
#'   the one set in `lmm_spec()`.
#' @return An updated `lmm_spec`.
#' @export
set_fixed <- function(spec, formula) {
  if (!inherits(spec, "lmm_spec")) stop("spec must be an lmm_spec object")
  if (!inherits(formula, "formula")) stop("formula must be of class 'formula'")
  if (grepl("\\|", deparse(formula))) stop("use set_random() for random effects")
  if (!is.null(spec$fixed$formula)) stop("fixed effects already set; create a new spec to change structure")

  response_name <- if (length(formula) == 3) as.character(formula[[2]]) else spec$response$name
  terms_obj <- stats::terms(formula)

  spec$response$name <- response_name
  spec$fixed$formula <- formula
  spec$fixed$terms <- c(
    if (attr(terms_obj, "intercept")) "1",
    attr(terms_obj, "term.labels")
  )

  spec
}


# ================================
# SET_RANDOM
# ================================

#' Set Random Effects for a Group
#'
#' Call once per random group, chained via pipe.
#'
#' @param spec An `lmm_spec` object.
#' @param formula A formula of the form `group ~ terms`, e.g. `id ~ 1 + x1`.
#' @param model_corr Logical; correlated random effects in the *fitted model*
#'   (default FALSE). TRUE → `(1 + x1 | id)`; FALSE → `(1 | id) + (0 + x1 | id)`.
#'   Independent of the simulation correlation set in `patch()`.
#' @return An updated `lmm_spec`.
#' @export
set_random <- function(spec, formula, model_corr = FALSE) {
  if (!inherits(spec, "lmm_spec")) stop("spec must be an lmm_spec object")
  if (!inherits(formula, "formula")) stop("formula must be of class 'formula'")
  if (length(formula) != 3) stop("formula must be of the form 'group ~ terms'")
  if (!is.logical(model_corr)) stop("model_corr must be TRUE or FALSE")

  group <- as.character(formula[[2]])
  if (group %in% names(spec$random)) stop("random group '", group, "' already defined")

  rhs <- formula[[3]]
  term_labels <- attr(stats::terms(formula[-2]), "term.labels")
  term_names <- c(if (grepl("(^|\\+)\\s*1", deparse(rhs))) "1", term_labels)

  if (length(term_names) == 0) stop("formula must have at least one term")

  counter <- sum(vapply(spec$random, function(g) length(g$terms), integer(1)))
  u_names <- stats::setNames(paste0("u", counter + seq_along(term_names) - 1), term_names)

  spec$random[[group]] <- list(
    terms      = term_names,
    u_names    = u_names,
    model_corr = model_corr
  )

  spec
}
