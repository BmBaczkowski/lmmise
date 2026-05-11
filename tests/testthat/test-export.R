describe("as_lme4", {
  test_that("returns formula and family", {
    spec <- set_random(
      set_fixed(lmm_spec(), y ~ 1 + x1),
      id ~ 1 + x1,
      model_corr = FALSE
    )

    result <- as_lme4(spec)
    expect_true(inherits(result$formula, "formula"))
    expect_equal(deparse(result$formula), "y ~ 1 + x1 + (1 | id) + (0 + x1 | id)")
  })

  test_that("correlated random effects use correct lme4 syntax", {
    spec <- set_random(
      set_fixed(lmm_spec(), y ~ 1 + x1),
      id ~ 1 + x1,
      model_corr = TRUE
    )

    result <- as_lme4(spec)
    expect_equal(deparse(result$formula), "y ~ 1 + x1 + (1 + x1 | id)")
  })

  test_that("include_random = FALSE drops random effects", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1 + x1), id ~ 1 + x1)

    result <- as_lme4(spec, include_random = FALSE)
    expect_equal(deparse(result$formula), "y ~ 1 + x1")
  })

  test_that("returns correct family for gaussian", {
    spec <- set_fixed(lmm_spec(), y ~ 1 + x1)
    result <- as_lme4(spec)
    expect_equal(result$family$family, "gaussian")
    expect_equal(result$family$link, "identity")
  })

  test_that("returns correct family for binomial", {
    spec <- set_fixed(lmm_spec(dist = "binomial"), y ~ 1 + x1)
    result <- as_lme4(spec)
    expect_equal(result$family$family, "binomial")
    expect_equal(result$family$link, "logit")
  })

  test_that("multiple random groups produce correct formula", {
    spec <- set_random(
      set_random(
        set_fixed(lmm_spec(), y ~ 1 + x1),
        id ~ 1,
        model_corr = FALSE
      ),
      site ~ 1,
      model_corr = FALSE
    )

    result <- as_lme4(spec)
    expect_equal(deparse(result$formula), "y ~ 1 + x1 + (1 | id) + (1 | site)")
  })

  test_that("expand main terms to include interaction", {
    spec <- set_fixed(lmm_spec(), y ~ 1 + x1 * x2)

    result <- as_lme4(spec)
    expect_equal(deparse(result$formula), "y ~ 1 + x1 + x2 + x1:x2")
  })
})

describe("as_simstudy rand_effects", {
  test_that("returns correct mu and sigma", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1 + x1), id ~ 1 + x1)
    params <- lmm_params(spec)
    params <- patch(params, list(
      beta = list(intercept = 1.0, x1 = 0.5),
      dispersion = 0.5,
      random_sd = list(id = list(intercept = 1.0, x1 = 0.3)),
      random_corr = list(id = list(structure = "cs", r = 0.3))
    ))

    result <- as_simstudy(spec, params, type = "rand_effects")
    expect_equal(result$id$mu, c(0, 0))
    expect_equal(result$id$sigma, c(1.0, 0.3))
  })

  test_that("returns correct cnames", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1 + x1), id ~ 1 + x1)
    params <- lmm_params(spec)
    params <- patch(params, list(
      beta = list(intercept = 1.0, x1 = 0.5),
      dispersion = 0.5,
      random_sd = list(id = list(intercept = 1.0, x1 = 0.3)),
      random_corr = list(id = list(structure = "cs", r = 0.3))
    ))

    result <- as_simstudy(spec, params, type = "rand_effects")
    expect_equal(result$id$cnames, c("u0", "u1"))
  })

  test_that("returns correct corstr for cs", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1 + x1), id ~ 1 + x1)
    params <- lmm_params(spec)
    params <- patch(params, list(
      beta = list(intercept = 1.0, x1 = 0.5),
      dispersion = 0.5,
      random_sd = list(id = list(intercept = 1.0, x1 = 0.3)),
      random_corr = list(id = list(structure = "cs", r = 0.3))
    ))

    result <- as_simstudy(spec, params, type = "rand_effects")

    m <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
    idnames <- c("intercept", "x1")
    rownames(m) <- idnames
    colnames(m) <- idnames
    expect_equal(result$id$corMatrix, m)
  })

  test_that("returns corMatrix for unstructured", {
    m <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1 + x1), id ~ 1 + x1)
    params <- lmm_params(spec)
    params <- patch(params, list(
      beta = list(intercept = 1.0, x1 = 0.5),
      dispersion = 0.5,
      random_sd = list(id = list(intercept = 1.0, x1 = 0.3)),
      random_corr = list(id = list(matrix = m))
    ))

    idnames <- c("intercept", "x1")
    rownames(m) <- idnames
    colnames(m) <- idnames
    result <- as_simstudy(spec, params, type = "rand_effects")
    expect_equal(result$id$corMatrix, m)
  })

  test_that("u names are globally unique across multiple groups", {
    spec <- set_random(set_random(set_fixed(lmm_spec(), y ~ 1 + x1), id ~ 1 + x1), site ~ 1)
    params <- lmm_params(spec)
    params <- patch(params, list(
      beta = list(intercept = 1.0, x1 = 0.5),
      dispersion = 0.5,
      random_sd = list(
        id = list(intercept = 1.0, x1 = 0.3),
        site = list(intercept = 0.5)
      ),
      random_corr = list(id = list(structure = "cs", r = 0.3))
    ))

    result <- as_simstudy(spec, params, type = "rand_effects")
    expect_equal(result$id$cnames, c("u0", "u1"))
    expect_equal(result$site$cnames, c("u2"))
  })
})


describe("as_lmList", {
  test_that("returns formula and family", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1 + x1), id ~ 1 + x1)

    result <- as_lmList(spec, cluster = "id")
    expect_true(inherits(result$formula, "formula"))
    expect_equal(deparse(result$formula), "y ~ 1 + x1 | id")
  })

  test_that("returns correct family for gaussian", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1 + x1), id ~ 1 + x1)

    result <- as_lmList(spec, cluster = "id")
    expect_equal(result$family$family, "gaussian")
    expect_equal(result$family$link, "identity")
  })

  test_that("returns correct family for binomial", {
    spec <- set_random(set_fixed(lmm_spec(dist = "binomial"), y ~ 1 + x1), id ~ 1 + x1)

    result <- as_lmList(spec, cluster = "id")
    expect_equal(result$family$family, "binomial")
    expect_equal(result$family$link, "logit")
  })

  test_that("works with intercept-only fixed effects", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1), id ~ 1)

    result <- as_lmList(spec, cluster = "id")
    expect_equal(deparse(result$formula), "y ~ 1 | id")
  })

  test_that("selects correct cluster from multiple random groups", {
    spec <- set_random(set_random(set_fixed(lmm_spec(), y ~ 1 + x1), id ~ 1), site ~ 1)

    result <- as_lmList(spec, cluster = "site")
    expect_equal(deparse(result$formula), "y ~ 1 + x1 | site")
  })

  test_that("errors when cluster is missing", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1 + x1), id ~ 1)

    expect_error(as_lmList(spec), "cluster must be provided")
  })

  test_that("errors when cluster is not in spec random groups", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1 + x1), id ~ 1)

    expect_error(as_lmList(spec, cluster = "site"), "not found in spec random groups")
  })

  test_that("errors when cluster is not a single string", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1 + x1), id ~ 1)

    expect_error(
      as_lmList(spec, cluster = c("id", "site")),
      "cluster must be a single character string"
    )
    expect_error(as_lmList(spec, cluster = 1L), "cluster must be a single character string")
  })

  test_that("errors when spec has no fixed effects", {
    spec <- lmm_spec()

    expect_error(as_lmList(spec, cluster = "id"), "call set_fixed\\(\\) first")
  })

  test_that("exclude_terms removes specified main effects", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1 + x1 + x2), id ~ 1)

    result <- as_lmList(spec, cluster = "id", exclude_terms = "x2")
    expect_equal(deparse(result$formula), "y ~ 1 + x1 | id")
  })

  test_that("exclude_terms removes interactions containing excluded variables", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1 + x1 * x2), id ~ 1)

    result <- as_lmList(spec, cluster = "id", exclude_terms = "x2")
    expect_equal(deparse(result$formula), "y ~ 1 + x1 | id")
  })

  test_that("exclude_terms with multiple exclusions", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1 + x1 + x2 + x3), id ~ 1)

    result <- as_lmList(spec, cluster = "id", exclude_terms = c("x1", "x3"))
    expect_equal(deparse(result$formula), "y ~ 1 + x2 | id")
  })

  test_that("exclude_terms removes x1 but keeps intercept", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1 + x1), id ~ 1)

    result <- as_lmList(spec, cluster = "id", exclude_terms = "x1")
    expect_equal(deparse(result$formula), "y ~ 1 | id")
  })

  test_that("exclude_terms errors when all terms are excluded", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 0 + x1), id ~ 1)

    expect_error(
      as_lmList(spec, cluster = "id", exclude_terms = "x1"),
      "No fixed terms remaining"
    )
  })

  test_that("exclude_terms with empty vector does nothing", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1 + x1), id ~ 1)

    result1 <- as_lmList(spec, cluster = "id")
    result2 <- as_lmList(spec, cluster = "id", exclude_terms = character(0))

    expect_equal(result1$formula, result2$formula)
  })

  test_that("exclude_terms handles comma-separated input", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1 + x1 * x2), id ~ 1)

    result <- as_lmList(spec, cluster = "id", exclude_terms = "x1, x2")
    expect_equal(deparse(result$formula), "y ~ 1 | id")
  })
})


describe("as_simstudy outcome", {
  test_that("returns correct formula string with no random effects on slope", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1 + x1), id ~ 1)
    params <- lmm_params(spec)
    params <- patch(params, list(
      beta = list(intercept = 1.0, x1 = 0.5),
      dispersion = 0.5,
      random_sd = list(id = list(intercept = 1.0))
    ))

    result <- as_simstudy(spec, params, type = "outcome")
    expect_equal(result$formula, "1 + u0 + 0.5 * x1")
  })

  test_that("returns correct formula string with random slope", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1 + x1), id ~ 1 + x1)
    params <- lmm_params(spec)
    params <- patch(params, list(
      beta = list(intercept = 1.0, x1 = 0.5),
      dispersion = 0.5,
      random_sd = list(id = list(intercept = 1.0, x1 = 0.3)),
      random_corr = list(id = list(structure = "cs", r = 0.3))
    ))

    result <- as_simstudy(spec, params, type = "outcome")
    expect_equal(result$formula, "1 + u0 + (0.5 + u1) * x1")
  })

  test_that("returns correct varname, dist, link, variance", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1 + x1), id ~ 1)
    params <- lmm_params(spec)
    params <- patch(params, list(
      beta = list(intercept = 1.0, x1 = 0.5),
      dispersion = 0.5,
      random_sd = list(id = list(intercept = 1.0))
    ))

    result <- as_simstudy(spec, params, type = "outcome")
    expect_equal(result$varname, "y")
    expect_equal(result$dist, "normal")
    expect_equal(result$link, "identity")
    expect_equal(result$variance, 0.5)
  })

  test_that("returns correct formula string with interactions", {
    spec <- set_random(
      set_fixed(lmm_spec(), y ~ 1 + x1 * x2),
      id ~ 1 + x1 * x2
    )
    params <- lmm_params(spec)
    params <- patch(params, list(
      beta = list(intercept = 1.0, x1 = 0.5, x2 = -1.0, "x1:x2" = 0.2),
      dispersion = 0.5,
      random_sd = list(id = list(
        intercept = 2, x1 = .5, x2 = 3, "x1:x2" = 2
      )),
      random_corr = list(id = list(structure = "ind"))
    ))

    result <- as_simstudy(spec, params, type = "outcome")
    expect_equal(
      result$formula,
      "1 + u0 + (0.5 + u1) * x1 + (-1 + u2) * x2 + (0.2 + u3) * x1*x2"
    )
  })
})
