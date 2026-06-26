test_that("Oja ranks: basic output structure", {
  X <- matrix(c(1,2, 3,4, 5,6), nrow=3, ncol=2, byrow=TRUE)

  ranks <- ojaRank(X, p=1)
  expect_equal(dim(ranks), dim(X))
  expect_true(is.numeric(ranks))
  expect_equal(rownames(ranks), rownames(X))
})

test_that("Oja ranks: single point query", {
  X <- matrix(c(1,2, 3,4, 5,6), nrow=3, ncol=2, byrow=TRUE)

  rnk <- ojaRank(X, x=c(2,3), p=1)
  expect_length(rnk, ncol(X))
  expect_true(is.numeric(rnk))
})

test_that("Oja ranks: output dimensions preserved under location shift", {
  X <- matrix(rnorm(12), nrow=6, ncol=2)
  b <- c(100, -50)

  ranks_X    <- ojaRank(X, p=1)
  X_shift <- sweep(X, 2, b, "+")
  ranks_X_shift <- ojaRank(X_shift, p=1)

  expect_equal(dim(ranks_X), dim(X))
  expect_equal(dim(ranks_X_shift), dim(X))
})

test_that("Oja ranks: rejects invalid input", {
  X <- matrix(1:6, nrow=3, ncol=2)

  expect_error(ojaRank(data.frame(a=letters[1:3], b=1:3)),
               "must be numeric")
  expect_error(ojaRank(X, x=c(1,2,3)),
               "must have the same dimension")
})

test_that("Oja ranks: works with biochem data", {
  data(biochem)
  X <- as.matrix(biochem[, c("comp.1", "comp.2")])

  ranks <- ojaRank(X, p=1)
  expect_equal(nrow(ranks), nrow(X))
  expect_equal(ncol(ranks), 2)
})
