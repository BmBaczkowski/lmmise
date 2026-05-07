describe("check_params", {
  test_that("returns NULL and messages when all params are set", {
    spec <- set_random(set_fixed(lmm_spec(), y ~ 1 + x1), id ~ 1 + x1)
    params <- lmm_params(spec)
    params <- patch(params, list(
      beta = list(intercept = 1.0, x1 = 0.5),
      dispersion = 0.5,
      random_sd = list(id = list(intercept = 1.0, x1 = 0.5)),
      random_corr = list(id = list(structure = "cs", r = 0.3))
    ))

    expect_null(check_params(params))
  })

  test_that("detects missing beta params", {
    spec <- set_fixed(lmm_spec(), y ~ 1 + x1 + x2)
    params <- lmm_params(spec)

    result <- check_params(params)
    expect_equal(sort(result$beta), sort(c("intercept", "x1", "x2")))
  })

  test_that("detects partially missing beta params", {
    spec <- set_fixed(lmm_spec(), y ~ 1 + x1 + x2)
    params <- lmm_params(spec)
    params <- patch(params, list(beta = list(intercept = 1.0)))

    result <- check_params(params)
    expect_equal(sort(result$beta), sort(c("x1", "x2")))
  })

  test_that("detects missing dispersion for gaussian", {
    spec <- set_fixed(lmm_spec(), y ~ 1)
    params <- lmm_params(spec)
    result <- check_params(params)
    expect_true(result$dispersion)
  })

  test_that("does not flag dispersion for poisson", {
    spec <- set_fixed(lmm_spec(dist = "poisson"), y ~ 1)
    params <- lmm_params(spec)
    result <- check_params(params)
    expect_null(result$dispersion)
  })

  test_that("detects missing random sd", {
    spec <- set_fixed(lmm_spec(), y ~ 1) |> set_random(id ~ 1 + x1 + x2)
    params <- lmm_params(spec)

    result <- check_params(params)
    expect_equal(sort(result$random$id$sd), sort(c("intercept", "x1", "x2")))
  })

  test_that("detects partially missing random sd", {
    spec <- set_fixed(lmm_spec(), y ~ 1) |> set_random(id ~ 1 + x1 + x2)
    params <- lmm_params(spec)
    params <- patch(params, list(random_sd = list(id = list(intercept = 1.0))))

    result <- check_params(params)
    expect_equal(sort(result$random$id$sd), sort(c("x1", "x2")))
  })

  test_that("detects missing random corr when group has multiple terms", {
    spec <- set_fixed(lmm_spec(), y ~ 1) |> set_random(id ~ 1 + x1)
    params <- lmm_params(spec)

    result <- check_params(params)
    expect_true(result$random$id$corr)
  })

  test_that("does not flag missing corr when group has single term", {
    spec <- set_fixed(lmm_spec(), y ~ 1) |> set_random(id ~ 1)
    params <- lmm_params(spec)

    result <- check_params(params)
    expect_null(result$random$id$corr)
  })

  test_that("handles multiple random groups", {
    spec <- set_fixed(lmm_spec(), y ~ 1) |>
      set_random(id ~ 1 + x1) |>
      set_random(site ~ 1)
    params <- lmm_params(spec)
    params <- patch(params, list(random_sd = list(id = list(intercept = 1.0))))

    result <- check_params(params)
    expect_equal(result$random$id$sd, "x1")
    expect_true(result$random$id$corr)
    expect_equal(result$random$site$sd, "intercept")
    expect_null(result$random$site$corr)
  })

  test_that("errors on non-params input", {
    expect_error(check_params(list()), "params must be an lmm_params object")
  })
})
