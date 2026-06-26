test_that("oja1sampleTest: sign scores with approximation", {
  set.seed(123)
  X <- matrix(rnorm(20), nrow=10, ncol=2)

  res <- oja1sampleTest(X, mu=c(0,0), scores="sign", method="approximation")
  expect_s3_class(res, "htest")
  expect_true(is.numeric(res$statistic))
  expect_true(is.numeric(res$p.value))
  expect_equal(res$alternative, "two.sided")
})

test_that("oja1sampleTest: rank scores with approximation", {
  set.seed(123)
  X <- matrix(rnorm(20), nrow=10, ncol=2)

  res <- oja1sampleTest(X, mu=c(0,0), scores="rank", method="approximation")
  expect_s3_class(res, "htest")
  expect_true(is.numeric(res$statistic))
  expect_true(is.numeric(res$p.value))
})

test_that("oja1sampleTest: sign scores with permutation", {
  set.seed(123)
  X <- matrix(rnorm(12), nrow=6, ncol=2)

  res <- oja1sampleTest(X, mu=c(0,0), scores="sign", method="permutation",
                        n.simu=100)
  expect_s3_class(res, "htest")
  expect_true(is.numeric(res$statistic))
  expect_true(is.numeric(res$p.value))
})

test_that("oja1sampleTest: uses default mu=origin", {
  set.seed(123)
  X <- matrix(rnorm(20), nrow=10, ncol=2)

  res <- oja1sampleTest(X, scores="sign", method="approximation")
  expect_s3_class(res, "htest")
})

test_that("oja1sampleTest: rejects invalid scores", {
  X <- matrix(rnorm(12), nrow=6, ncol=2)

  expect_error(oja1sampleTest(X, scores="invalid"))
})
