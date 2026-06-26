test_that("Oja median: evolutionary algorithm with ICS", {
  set.seed(123)
  X <- matrix(rnorm(20), nrow=10, ncol=2)

  med <- ojaMedian(X, alg="evolutionary",
                   control=ojaMedianControl(iter=2000))
  expect_length(med, ncol(X))
  expect_true(is.numeric(med))
})

test_that("Oja median: grid algorithm with ICS", {
  set.seed(123)
  X <- matrix(rnorm(20), nrow=10, ncol=2)

  med <- ojaMedian(X, alg="grid",
                   control=ojaMedianControl(samples=5, eps=0.5))
  expect_length(med, ncol(X))
  expect_true(is.numeric(med))
})

test_that("Oja median: multiple runs averaged", {
  set.seed(123)
  X <- matrix(rnorm(20), nrow=10, ncol=2)

  med <- ojaMedian(X, alg="evolutionary", sp=2,
                   control=ojaMedianControl(iter=500))
  expect_length(med, ncol(X))
  expect_true(is.numeric(med))
})

test_that("Oja median function value: returns scalar", {
  X <- matrix(c(1,2, 3,4, 5,6, 7,8), nrow=4, ncol=2, byrow=TRUE)

  fval <- ojaMedianFn(X, x=c(3,4))
  expect_length(fval, 1)
  expect_true(is.numeric(fval))
})

test_that("ojaMedianEvo barebones works directly", {
  set.seed(123)
  X <- matrix(rnorm(20), nrow=10, ncol=2)

  med <- ojaMedianEvo(X, control=ojaMedianControl(iter=1000))
  expect_length(med, ncol(X))
})

test_that("ojaMedianGrid barebones works directly", {
  set.seed(123)
  X <- matrix(rnorm(20), nrow=10, ncol=2)

  med <- ojaMedianGrid(X, control=ojaMedianControl(samples=5))
  expect_length(med, ncol(X))
})
