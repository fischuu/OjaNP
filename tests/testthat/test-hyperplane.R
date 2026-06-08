test_that("hyperplane: produces correct output structure (2D)", {
  X <- matrix(c(1,1, 2,3), nrow=2, ncol=2, byrow=TRUE)

  hp <- hyperplane(X)
  expect_length(hp, 3)
  expect_true(is.numeric(hp))
})

test_that("hyperplane: 3D hyperplane through 3 points", {
  X <- matrix(c(0,0,0, 1,0,0, 0,1,0), nrow=3, ncol=3, byrow=TRUE)

  hp <- hyperplane(X)
  expect_length(hp, 4)
  expect_true(is.numeric(hp))
})

test_that("hyperplane: rejects non-square matrix", {
  X <- matrix(1:6, nrow=3, ncol=2)

  expect_error(hyperplane(X), "must be a k x k matrix")
})
