# ================================
# PRINT LMM_SPEC
# ================================

#' @export
print.lmm_spec <- function(x, ...) {
  cat("== lmm_spec ==\n")
  cat("Response : ", x$response$name,
    " (", x$response$dist, " / ", x$response$link, ")\n",
    sep = ""
  )

  if (!is.null(x$fixed$formula)) {
    cat("Fixed    : ", paste(x$fixed$terms, collapse = " + "), "\n", sep = "")
  } else {
    cat("Fixed    : <not set>\n")
  }

  if (length(x$random) > 0) {
    for (group in names(x$random)) {
      g <- x$random[[group]]
      terms_str <- paste(ifelse(g$terms == "1", "1", g$terms), collapse = " + ")
      corr_tag <- if (length(g$terms) > 1) {
        if (g$model_corr) " [correlated]" else " [independent]"
      } else {
        ""
      }
      cat("Random   : [", group, "] ", terms_str, corr_tag, "\n", sep = "")
    }
  } else {
    cat("Random   : <not set>\n")
  }

  invisible(x)
}


# ================================
# PRINT LMM_PARAMS
# ================================

#' @export
print.lmm_params <- function(x, ...) {
  cat("== lmm_params ==\n\n")

  # --- Beta ---
  cat("Beta:\n")
  for (nm in names(x$beta)) {
    val <- x$beta[[nm]]
    cat("  ", .pad(nm), " = ", if (is.null(val)) "<unset>" else val, "\n", sep = "")
  }
  cat("\n")

  # --- Dispersion ---
  dist <- x$spec$response$dist
  if (dist %in% c("binomial", "poisson")) {
    cat("Dispersion: 1 (fixed for", dist, ")\n")
  } else {
    val <- x$dispersion
    cat("Dispersion:", if (is.null(val)) "<unset>" else val, "\n")
  }
  cat("\n")

  # --- Random ---
  groups <- names(x$random)
  for (i in seq_along(groups)) {
    group <- groups[[i]]
    g <- x$random[[group]]
    term_names <- names(g$sd)

    .print_divider(paste("Random:", group))

    cat("  SD:\n")
    for (term in term_names) {
      val <- g$sd[[term]]
      cat("    ", .pad(term), " = ", if (is.null(val)) "<unset>" else val, "\n", sep = "")
    }

    if (length(term_names) >= 2) {
      cat("\n  Corr:\n")
      if (is.null(g$corr)) {
        cat("  (not set, defaults to independent)\n")
        .print_corr_matrix(.make_identity(term_names), term_names, indent = 4)
      } else {
        .print_corr_matrix(g$corr, term_names, indent = 4)
      }
    } else {
      cat("\n  (single term, no correlation)\n")
    }

    .print_divider()
    if (i < length(groups)) cat("\n")
  }

  invisible(x)
}


# ================================
# PRINT HELPERS
# ================================

.pad <- function(x, width = 12) formatC(x, width = -width, flag = "-")

.print_divider <- function(label = NULL, width = 48) {
  if (!is.null(label)) {
    prefix <- "\u2500\u2500 "
    suffix_len <- width - nchar(label) - nchar(prefix) - 1
    cat(prefix, label, " ", strrep("\u2500", max(suffix_len, 2)), "\n", sep = "")
  } else {
    cat(strrep("\u2500", width), "\n", sep = "")
  }
}

.make_identity <- function(term_names) {
  m <- diag(length(term_names))
  rownames(m) <- colnames(m) <- term_names
  m
}

.print_corr_matrix <- function(m, term_names, indent = 0) {

  col_width <- max(nchar(term_names)) + 2

  pad <- strrep(" ", indent)
  cat(pad, sep = "")
  for (t in term_names) cat(.pad(t, col_width))
  cat("\n")

  for (i in seq_along(term_names)) {
    cat(pad, .pad(term_names[i], col_width), sep = "")
    for (j in seq_along(term_names)) {
      cat(.pad(formatC(m[i, j], format = "f", digits = 2), col_width))
    }
    cat("\n")
  }
}
