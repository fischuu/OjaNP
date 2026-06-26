test_that("Oja SCM: basic output structure", {
  X <- matrix(rnorm(20), nrow=10, ncol=2)

  scm <- ojaSCM(X, p=1)
  expect_equal(dim(scm), c(2, 2))
  expect_true(isSymmetric(scm))
  expect_true(is.numeric(scm))
})

test_that("Oja SCM: with colMean centering", {
  X <- matrix(rnorm(20), nrow=10, ncol=2)

  scm <- ojaSCM(X, center="colMean", p=1)
  expect_equal(dim(scm), c(2, 2))
  expect_true(isSymmetric(scm))
})

test_that("Oja RCM: basic output structure", {
  X <- matrix(rnorm(20), nrow=10, ncol=2)

  rcm <- ojaRCM(X, p=1)
  expect_equal(dim(rcm), c(2, 2))
  expect_true(isSymmetric(rcm))
  expect_true(is.numeric(rcm))
})

test_that("Oja RCM: 3D data", {
  X <- matrix(rnorm(30), nrow=10, ncol=3)

  rcm <- ojaRCM(X, p=1)
  expect_equal(dim(rcm), c(3, 3))
  expect_true(isSymmetric(rcm))
})
