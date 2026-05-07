# ================================
# VALIDATION HELPERS
# ================================

.validate_named_numeric <- function(values, fn_name) {
  if (length(values) == 0) stop(fn_name, " requires at least one named argument")
  if (is.null(names(values)) || any(names(values) == "")) stop("all arguments to ", fn_name, " must be named")
  if (!all(sapply(values, is.numeric))) stop("all values in ", fn_name, " must be numeric")
}

.validate_corr_params <- function(structure, r) {
  if (structure %in% c("cs", "ar1")) {
    if (is.null(r) || !is.numeric(r) || length(r) != 1) stop("structure '", structure, "' requires a single numeric r")
    if (r < -1 || r > 1) stop("r must be between -1 and 1")
  }
  if (structure == "toep") {
    if (is.null(r) || !is.numeric(r)) stop("structure 'toep' requires a numeric vector r")
    if (any(r < -1 | r > 1)) stop("all values of r must be between -1 and 1")
  }
}

.validate_corr_matrix <- function(m) {
  if (!is.matrix(m)) stop("matrix must be a matrix object")
  if (nrow(m) != ncol(m)) stop("correlation matrix must be square")
  if (!all(diag(m) == 1)) stop("diagonal of correlation matrix must be all 1s")
  if (!all(m >= -1 & m <= 1)) stop("all values must be between -1 and 1")
  if (!all(abs(m - t(m)) < .Machine$double.eps^0.5)) stop("correlation matrix must be symmetric")
}


# ================================
# CORRELATION MATRIX BUILDERS
# ================================

.make_corr_matrix <- function(term_names, spec) {
  n <- length(term_names)

  if (!is.null(spec$matrix)) {
    m <- spec$matrix
    if (!is.null(rownames(m)) && !identical(rownames(m), term_names)) {
      stop(
        "matrix rownames do not match terms: expected ",
        paste(term_names, collapse = ", ")
      )
    }
    rownames(m) <- term_names
    colnames(m) <- term_names
    .validate_corr_matrix(m)
    return(m)
  }

  structure <- spec$structure
  r <- spec$r

  valid_structures <- c("cs", "ind", "ar1", "toep")
  if (!structure %in% valid_structures) {
    stop("structure must be one of: ", paste(valid_structures, collapse = ", "))
  }

  .validate_corr_params(structure, r)

  m <- if (structure == "ind") {
    diag(n)
  } else if (structure == "cs") {
    mm <- matrix(r, nrow = n, ncol = n)
    diag(mm) <- 1
    mm
  } else if (structure == "ar1") {
    outer(1:n, 1:n, function(i, j) r^abs(i - j))
  } else if (structure == "toep") {
    if (length(r) != n - 1) stop("toep requires ", n - 1, " value(s) in r")
    mm <- diag(n)
    for (lag in seq_along(r)) {
      for (i in 1:(n - lag)) {
        mm[i, i + lag] <- r[lag]
        mm[i + lag, i] <- r[lag]
      }
    }
    mm
  }

  rownames(m) <- term_names
  colnames(m) <- term_names
  m
}
