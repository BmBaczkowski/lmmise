# Test patch beta
describe("patch beta", {
  test_that("sets beta params", {
    spec <- set_fixed(lmm_spec(), ~ 1 + x1)
    params <- lmm_params(spec)
    params <- patch(params, list(beta = list(intercept = 2, x1 = 1.5)))

    expect_equal(params$beta$intercept, 2)
    expect_equal(params$beta$x1, 1.5)
  })
})


# Test patch dispersion
describe("patch dispersion", {
  test_that("sets dispersion param", {
    spec <- set_fixed(lmm_spec(), ~ 1 + x1)
    params <- lmm_params(spec)
    params <- patch(params, list(dispersion = 1))
    expect_equal(params$dispersion, 1)
  })

  test_that("sets dispersion param for gamma", {
    spec <- set_fixed(lmm_spec(dist = "gamma"), ~ 1 + x1)
    params <- lmm_params(spec)
    params <- patch(params, list(dispersion = 2))
    expect_equal(params$dispersion, 2)
  })
})


# Test patch random_sd
describe("patch random_sd", {
  test_that("sets random sd params", {
    spec <- set_fixed(lmm_spec(), ~ 1 + x1) |>
      set_random(id ~ 1 + x1 + x2)
    params <- lmm_params(spec)
    params <- patch(
      params,
      list(random_sd = list(id = list(intercept = 1, x1 = 0.5, x2 = 0.8)))
    )

    expect_equal(params$random$id$sd[["intercept"]], 1)
    expect_equal(params$random$id$sd$x1, 0.5)
    expect_equal(params$random$id$sd$x2, 0.8)
  })
})


# Test patch random_corr
describe("patch random_corr", {
  test_that("sets cs correlation structure", {
    spec <- lmm_spec() |>
      set_fixed(~1) |>
      set_random(id ~ 1 + x1 + x2)
    params <- lmm_params(spec)
    params <- patch(params, list(random_corr = list(id = list(structure = "cs", r = 0.3))))

    expect_equal(nrow(params$random$id$corr), 3)
    expect_equal(params$random$id$corr[1, 2], 0.3)
  })

  test_that("sets ar1 correlation structure", {
    spec <- lmm_spec() |>
      set_fixed(~1) |>
      set_random(id ~ 1 + x1 + x2)
    params <- lmm_params(spec)
    params <- patch(params, list(random_corr = list(id = list(structure = "ar1", r = 0.5))))

    expect_equal(params$random$id$corr[1, 2], 0.5)
    expect_equal(params$random$id$corr[1, 3], 0.5^2)
  })

  test_that("sets ind correlation structure as identity matrix", {
    spec <- lmm_spec() |>
      set_fixed(~1) |>
      set_random(id ~ 1 + x1 + x2)
    params <- lmm_params(spec)
    params <- patch(params, list(random_corr = list(id = list(structure = "ind"))))

    m <- diag(3)
    idnames <- c("intercept", "x1", "x2")
    rownames(m) <- idnames
    colnames(m) <- idnames
    expect_equal(params$random$id$corr, m)
  })

  test_that("sets unstructured correlation via matrix", {
    m <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
    spec <- lmm_spec() |>
      set_fixed(~1) |>
      set_random(id ~ 1 + x1)
    params <- lmm_params(spec)
    params <- patch(params, list(random_corr = list(id = list(matrix = m))))

    idnames <- c("intercept", "x1")
    rownames(m) <- idnames
    colnames(m) <- idnames
    expect_equal(params$random$id$corr, m)
  })
})
