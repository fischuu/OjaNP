test_that("ojaCsampleTest: two-sample default method sign scores", {
  set.seed(123)
  X <- matrix(rnorm(14), nrow=7, ncol=2)
  Y <- matrix(rnorm(14, mean=1), nrow=7, ncol=2)

  res <- ojaCsampleTest(X, Y, scores="sign", method="approximation")
  expect_s3_class(res, "htest")
  expect_true(is.numeric(res$statistic))
  expect_true(is.numeric(res$p.value))
  expect_equal(res$alternative, "two.sided")
})

test_that("ojaCsampleTest: two-sample rank scores", {
  set.seed(123)
  X <- matrix(rnorm(16), nrow=8, ncol=2)
  Y <- matrix(rnorm(16, mean=0.5), nrow=8, ncol=2)

  res <- ojaCsampleTest(X, Y, scores="rank", method="approximation")
  expect_s3_class(res, "htest")
  expect_true(is.numeric(res$statistic))
  expect_true(is.numeric(res$p.value))
})

test_that("ojaCsampleTest: two-sample permutation", {
  set.seed(123)
  X <- matrix(rnorm(8), nrow=4, ncol=2)
  Y <- matrix(rnorm(8, mean=1), nrow=4, ncol=2)

  res <- ojaCsampleTest(X, Y, scores="sign", method="permutation",
                        n.simu=100)
  expect_s3_class(res, "htest")
  expect_true(is.numeric(res$p.value))
})

test_that("ojaCsampleTest: formula method with biochem data", {
  data(biochem)

  res <- ojaCsampleTest(cbind(comp.1, comp.2) ~ group, data=biochem,
                         scores="sign", method="approximation")
  expect_s3_class(res, "htest")
})
