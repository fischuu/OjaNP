test_that("Oja signed ranks: basic output structure", {
  X <- matrix(c(1,2, 3,4, 5,6), nrow=3, ncol=2, byrow=TRUE)

  sr <- ojaSignedRank(X, p=1)
  expect_equal(dim(sr), dim(X))
  expect_true(is.numeric(sr))
  expect_equal(rownames(sr), rownames(X))
})

test_that("Oja signed ranks: single point query", {
  X <- matrix(c(1,2, 3,4, 5,6), nrow=3, ncol=2, byrow=TRUE)

  rnk <- ojaSignedRank(X, x=c(2,3), p=1)
  expect_length(rnk, ncol(X))
  expect_true(is.numeric(rnk))
})

test_that("Oja signed ranks: 3D data", {
  X <- matrix(c(1,2,3, 4,5,6, 7,8,9, 10,11,12), nrow=4, ncol=3, byrow=TRUE)

  sr <- ojaSignedRank(X, p=1)
  expect_equal(nrow(sr), nrow(X))
  expect_equal(ncol(sr), ncol(X))
})

test_that("Oja signed ranks: rejects invalid input", {
  X <- matrix(1:6, nrow=3, ncol=2)

  expect_error(ojaSignedRank(data.frame(a=letters[1:3], b=1:3)),
               "must be numeric")
  expect_error(ojaSignedRank(X, x=c(1,2,3)),
               "must have the same dimension")
})

test_that("Oja signed ranks: works with biochem data", {
  data(biochem)
  X <- as.matrix(biochem[, c("comp.1", "comp.2")])

  sr <- ojaSignedRank(X, p=1)
  expect_equal(nrow(sr), nrow(X))
  expect_equal(ncol(sr), 2)
})
