################################################################################
# Reproduction of Section 4.2 from:
# Fischer D, Mosler K, Mottonen J, Nordhausen K, Pokotylo O, Vogel D (2020).
# "Computing the Oja Median in R: The Package OjaNP"
# Journal of Statistical Software, 92(8), 1-36.
# doi: 10.18637/jss.v092.i08
################################################################################

library(OjaNP)
data(biochem)

cat("===========================================================================\n")
cat("Section 4.2: Illustration on biochemical data\n")
cat("Journal of Statistical Software, 92(8), 2020\n")
cat("===========================================================================\n\n")

X <- as.matrix(biochem[, 1:2])

cat("Data: biochem (Brown & Hettmansperger, 1987)\n")
cat("  n =", nrow(X), "observations\n")
cat("  k =", ncol(X), "variables\n")
cat("  Groups: Control (n=10), Treat (n=12)\n\n")
print(biochem)
library("OjaNP")
data("biochem", package = "OjaNP")
set.seed(1)
X <- as.matrix(biochem[, 1:2])
GROUP <- biochem$group
GRlabel <- as.numeric(GROUP)

# ---------------------------------------------------------------------------
# 4.2.1  Oja median
# ---------------------------------------------------------------------------
cat("\n-----------------------------------------------------------------------\n")
cat("4.2.1  Oja median\n")
cat("-----------------------------------------------------------------------\n\n")

set.seed(1)
OMev <- ojaMedian(X)
OMev

OMex <- ojaMedian(X, alg = "exact")
OMex
cat("\nBounded exact algorithm:\n")
set.seed(1)
OM_bounded <- ojaMedian(X, alg = "bounded_exact")
print(OM_bounded)

cat("\nGrid algorithm:\n")
OM_grid <- ojaMedian(X, alg = "grid",
                     control = ojaMedianControl(samples = 100))
print(OM_grid)

# ---------------------------------------------------------------------------
# 4.2.2  Oja signs
# ---------------------------------------------------------------------------
cat("\n-----------------------------------------------------------------------\n")
cat("4.2.2  Oja signs\n")
cat("-----------------------------------------------------------------------\n\n")
set.seed(1)
cat("Oja signs with componentwise median center:\n")
signs_comp <- head(ojaSign(X,center = "ojaMedian"))
print(signs_comp)

cat("\nOja signs centered at Oja median (compMedian):\n")
signs_om <- ojaSign(X, center = "compMedian")

print(signs_om)

# ---------------------------------------------------------------------------
# 4.2.3  Oja SCM (Sign Covariance Matrix)
# ---------------------------------------------------------------------------
cat("\n-----------------------------------------------------------------------\n")
cat("4.2.3  Oja SCM (Sign Covariance Matrix)\n")
cat("-----------------------------------------------------------------------\n\n")

cat("Oja SCM with Oja median center (exact):\n")
scm <- ojaSCM(X, center = "ojaMedian", alg = "exact")
print(scm)



cat("\nSCM equals (n-1)/n * cov(signs) when signs are centered:\n")
signs <- ojaSign(X, center = "ojaMedian", alg = "exact")
scm_from_cov <- (nrow(X) - 1) / nrow(X) * cov(signs)
cat("max|scm - (n-1)/n*cov(signs)| =", max(abs(scm - scm_from_cov)), "\n")

cat("\nOja SCM with column mean center:\n")
scm_colmean <- ojaSCM(X, center = "colMean")
print(scm_colmean)

# ---------------------------------------------------------------------------
# 4.2.4  Oja RCM (Rank Covariance Matrix)
# ---------------------------------------------------------------------------
cat("\n-----------------------------------------------------------------------\n")
cat("4.2.4  Oja RCM (Rank Covariance Matrix)\n")
cat("-----------------------------------------------------------------------\n\n")

cat("Oja RCM:\n")
rcm <- ojaRCM(X, p = 1)
print(rcm, digits = 5)

cat("\nRCM is symmetric:", isSymmetric(rcm), "\n")

# ---------------------------------------------------------------------------
# 4.2.5  1-sample location test
# ---------------------------------------------------------------------------
cat("\n-----------------------------------------------------------------------\n")
cat("4.2.5  1-sample location test\n")
cat("-----------------------------------------------------------------------\n\n")

cat("1-sample sign test, H0: mu = c(1, 0.5), first 10 obs:\n")
cat("  Approximation method:\n")
res_1s_approx <- oja1sampleTest(X[1:10, ], mu = c(1, 0.5))
print(res_1s_approx)

cat("\n  Permutation method (1000 permutations, set.seed = 1):\n")
set.seed(1)
res_1s_perm <- oja1sampleTest(X[1:10, ], mu = c(1, 0.5),
                              method = "permutation")
print(res_1s_perm)

cat("\n  Signed rank test (scores = 'rank'):\n")
res_1s_sr <- oja1sampleTest(X[1:10, ], mu = c(1, 0.5), scores = "rank")
print(res_1s_sr)

# ---------------------------------------------------------------------------
# 4.2.6  C-sample location test
# ---------------------------------------------------------------------------
cat("\n-----------------------------------------------------------------------\n")
cat("4.2.6  C-sample location test\n")
cat("-----------------------------------------------------------------------\n\n")

GROUP <- biochem$group

cat("C-sample rank test (permutation):\n")
set.seed(1)
res_cs_rank <- ojaCsampleTest(X ~ GROUP, scores = "rank",
                               method = "permutation")
print(res_cs_rank)

cat("\nC-sample sign test (formula interface, approximation):\n")
set.seed(1)
res_cs_sign <- ojaCsampleTest(cbind(comp.1, comp.2) ~ group,
                               data = biochem,
                               scores = "sign",
                               method = "approximation")
print(res_cs_sign)

cat("\n===========================================================================\n")
cat("Reproduction complete.\n")
cat("===========================================================================\n")
