test_that("Oja signs: basic output structure", {
  X <- matrix(rnorm(10), nrow=5, ncol=2)

  signs <- ojaSign(X, center="colMean", p=1)
  expect_equal(dim(signs), dim(X))
  expect_true(is.numeric(signs))
})

test_that("Oja signs: single point query", {
  X <- matrix(rnorm(10), nrow=5, ncol=2)

  sgn <- ojaSign(X, x=c(0,0), center="colMean", p=1)
  expect_length(sgn, ncol(X))
  expect_true(is.numeric(sgn))
})

test_that("Oja signs: different center options work", {
  X <- matrix(rnorm(20), nrow=10, ncol=2)

  sgn_median  <- ojaSign(X, center="ojaMedian", p=1)
  sgn_colmean <- ojaSign(X, center="colMean", p=1)
  sgn_vec     <- ojaSign(X, center=c(0,0), p=1)

  expect_equal(dim(sgn_median), dim(X))
  expect_equal(dim(sgn_colmean), dim(X))
  expect_equal(dim(sgn_vec), dim(X))
})

test_that("Oja signs: output dimensions preserved under affine transform", {
  X <- matrix(rnorm(8), nrow=4, ncol=2)
  A <- matrix(c(2,0, 0,3), nrow=2)
  b <- c(10, 20)

  signs_X <- ojaSign(X, center="colMean", p=1)
  signs_AXb <- ojaSign(t(apply(X, 1, function(x) A %*% x + b)), center="colMean", p=1)

  expect_equal(dim(signs_X), dim(X))
  expect_equal(dim(signs_AXb), dim(X))
})

test_that("Oja signs: rejects invalid input", {
  X <- matrix(1:6, nrow=3, ncol=2)

  expect_error(ojaSign(data.frame(a=letters[1:3], b=1:3)),
               "must be numeric")
  expect_error(ojaSign(X, x=c(1,2,3)),
               "must have the same dimension")
})

test_that("Oja signs: works with biochem data", {
  data(biochem)
  X <- as.matrix(biochem[, c("comp.1", "comp.2")])

  signs <- ojaSign(X, center="colMean", p=1)
  expect_equal(nrow(signs), nrow(X))
  expect_equal(ncol(signs), 2)
})
