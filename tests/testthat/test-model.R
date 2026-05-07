# =========================
# Tests for LMM Spec
# =========================

# Test lmm_spec
describe("lmm_spec", {
  test_that("creates a spec object", {
    spec <- lmm_spec()
    expect_s3_class(spec, "lmm_spec")
    expect_type(spec$fixed, "list")
    expect_type(spec$random, "list")
    expect_type(spec$response, "list")
    expect_equal(spec$response$name, "y")
    expect_equal(spec$response$dist, "gaussian")
    expect_equal(spec$response$link, "identity")
  })
  test_that("creates a spec with custom response and distribution", {
    spec <- lmm_spec(response = "outcome", dist = "poisson")
    expect_equal(spec$response$name, "outcome")
    expect_equal(spec$response$dist, "poisson")
    expect_equal(spec$response$link, "log")
  })
})

# Test set_fixed
describe("set_fixed", {
  test_that("adds fixed effects", {
    spec <- set_fixed(lmm_spec(), ~ 1 + x1 + x2)
    expect_length(spec$fixed$terms, 3)
    expect_true("1" %in% spec$fixed$terms)
    expect_true("x1" %in% spec$fixed$terms)
    expect_true("x2" %in% spec$fixed$terms)
    expect_equal(spec$response$name, "y")
    expect_equal(length(spec$random), 0)
  })
  test_that("adds fixed effects without intercept", {
    spec <- set_fixed(lmm_spec(), ~ 0 + x1 + x2)
    expect_length(spec$fixed$terms, 2)
    expect_false("1" %in% spec$fixed$terms)
    expect_true("x1" %in% spec$fixed$terms)
    expect_true("x2" %in% spec$fixed$terms)
    expect_equal(spec$response$name, "y")
    expect_equal(length(spec$random), 0)
  })
  test_that("adds fixed effects with custom response", {
    spec <- set_fixed(lmm_spec(), outcome ~ 1 + x1)
    expect_equal(spec$response$name, "outcome")
    expect_true("1" %in% spec$fixed$terms)
    expect_true("x1" %in% spec$fixed$terms)
  })
})

# Test set_random
describe("set_random", {
  test_that("adds only random effects with intercept but without correlation", {
    spec <- set_random(lmm_spec(), id ~ 1)
    expect_true("id" %in% names(spec$random))
    expect_true("1" %in% spec$random$id$terms)
    expect_false(spec$random$id$model_corr)
    expect_equal(length(spec$fixed$terms), 0)
  })
  test_that("adds only random effects without intercept", {
    spec <- set_random(lmm_spec(), id ~ x2, model_corr = TRUE)
    expect_true("id" %in% names(spec$random))
    expect_true("x2" %in% spec$random$id$terms)
    expect_true(spec$random$id$model_corr)
    expect_false("1" %in% spec$random$id$terms)
  })
  test_that("adds only random effects without intercept (-1)", {
    spec <- set_random(lmm_spec(), id ~ -1 + x2)
    expect_true("id" %in% names(spec$random))
    expect_false("1" %in% spec$random$id$terms)
  })
  test_that("adds only random effects without intercept (0)", {
    spec <- set_random(lmm_spec(), id ~ 0 + x2)
    expect_true("id" %in% names(spec$random))
    expect_false("1" %in% spec$random$id$terms)
  })
  test_that("adds multiple random effects", {
    spec <- set_random(set_random(lmm_spec(), id ~ 0 + x1), item ~ 1)
    expect_true("id" %in% names(spec$random))
    expect_true("item" %in% names(spec$random))
    expect_false("1" %in% spec$random$id$terms)
    expect_true("x1" %in% spec$random$id$terms)
    expect_true("1" %in% spec$random$item$terms)
  })
})

# Test assertions
describe("Assertions", {
  test_that("lmm_spec rejects invalid response", {
    expect_error(lmm_spec(response = c("a", "b")), "response must be")
  })
  test_that("lmm_spec rejects invalid distribution", {
    expect_error(lmm_spec(dist = "foo"), "dist must be one of")
  })
  test_that("set_fixed rejects invalid formula", {
    expect_error(set_fixed(lmm_spec(), "not a formula"), "formula must be")
  })
  test_that("set_fixed rejects duplicate fixed", {
    spec <- set_fixed(lmm_spec(), ~1)
    expect_error(set_fixed(spec, ~x1), "fixed effects already set")
  })
})
