test_that("Paper Section 4.2: Oja median evolutionary matches paper", {
  data(biochem, package = "OjaNP")
  set.seed(1)
  X <- as.matrix(biochem[, 1:2])

  OMev <- ojaMedian(X, alg = "evolutionary",
                    control = ojaMedianControl(iter = 1e6))

  expect_equal(round(OMev, 3), c(comp.1 = 1.150, comp.2 = 0.425),
               tolerance = 0.01)
})

test_that("Paper Section 4.2: Oja signs with compMedian center match", {
  data(biochem, package = "OjaNP")
  X <- as.matrix(biochem[, 1:2])

  signs <- ojaSign(X, center = "compMedian")

  expect_equal(nrow(signs), 22)
  expect_equal(ncol(signs), 2)
  expect_equal(round(signs[1, 1], 4), 0.0107, tolerance = 0.001)
  expect_equal(round(signs[1, 2], 4), 0.1145, tolerance = 0.001)
  expect_equal(round(signs[2, 1], 4), -0.0634, tolerance = 0.001)
  expect_equal(round(signs[2, 2], 4), 0.0195, tolerance = 0.001)
})

test_that("Paper Section 4.2: Oja SCM is symmetric positive", {
  data(biochem, package = "OjaNP")
  X <- as.matrix(biochem[, 1:2])

  scm <- ojaSCM(X)

  expect_equal(dim(scm), c(2, 2))
  expect_true(isSymmetric(scm))
  expect_true(all(diag(scm) > 0))
})

test_that("Paper Section 4.2: 1-sample sign test matches paper", {
  data(biochem, package = "OjaNP")
  X <- as.matrix(biochem[, 1:2])

  res <- oja1sampleTest(X[1:10, ], mu = c(1, 0.5))

  expect_s3_class(res, "htest")
  expect_equal(round(res$statistic, 4), c("Q.S" = 3.3745))
  expect_equal(res$parameter, c("df" = 2))
  expect_true(res$p.value > 0.1)
})

test_that("Paper Section 4.2: 1-sample test permutation", {
  data(biochem, package = "OjaNP")
  X <- as.matrix(biochem[, 1:2])

  res <- oja1sampleTest(X[1:10, ], mu = c(1, 0.5),
                        method = "permutation")

  expect_s3_class(res, "htest")
  expect_equal(round(res$statistic, 4), c("Q.S" = 3.3745))
  expect_true(res$p.value > 0.1)
})

test_that("Paper Section 4.2: C-sample rank test matches paper", {
  data(biochem, package = "OjaNP")
  X <- as.matrix(biochem[, 1:2])
  GROUP <- biochem$group

  res <- ojaCsampleTest(X ~ GROUP, scores = "rank",
                        method = "permutation")

  expect_s3_class(res, "htest")
  expect_equal(round(res$statistic, 2), c("Q.R" = 15.17))
})

test_that("Paper Section 4.2: C-sample sign test with formula", {
  data(biochem)

  res <- ojaCsampleTest(cbind(comp.1, comp.2) ~ group, data = biochem,
                         scores = "sign", method = "approximation")

  expect_s3_class(res, "htest")
  expect_true(res$p.value < 0.01)
})
