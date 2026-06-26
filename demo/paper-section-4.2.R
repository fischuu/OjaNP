################################################################################
# Reproduction of Section 4.2 from:
# Fischer D, Mosler K, Mottonen J, Nordhausen K, Pokotylo O, Vogel D (2020).
# "Computing the Oja Median in R: The Package OjaNP"
# Journal of Statistical Software, 92(8), 1-36.
# doi: 10.18637/jss.v092.i08
################################################################################

## NOTE (fixes_j branch): removed extra blank line between header and code


cat("===========================================================================\n")
cat("Section 4.2: Illustration on biochemical data\n")
cat("Journal of Statistical Software, 92(8), 2020\n")
cat("===========================================================================\n\n")


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


# 4.2.1  Oja median
cat("\n-----------------------------------------------------------------------\n")
cat("4.2.1  Oja median\n")
cat("-----------------------------------------------------------------------\n\n")
## NOTE (fixes_j): added label; removed redundant set.seed(1); removed Grid algorithm section
cat("Default configuration:\n")
OMev <- ojaMedian(X)
OMev

cat("\nExact algorithm:\n")

OMex <- ojaMedian(X, alg = "exact")
OMex

cat("\nBounded exact algorithm:\n")
OM_bounded <- ojaMedian(X, alg = "bounded_exact")
print(OM_bounded)

# ---------------------------------------------------------------------------
# 4.2.2  Oja signs
# ---------------------------------------------------------------------------
cat("\n-----------------------------------------------------------------------\n")
cat("4.2.2  Oja signs\n")
cat("-----------------------------------------------------------------------\n\n")
## NOTE (fixes_j): simplified section removed set.seed(1); call now uses default center="ojaMedian"
cat("Oja signs default:\n")
signs_comp <- head(ojaSign(X))
print(signs_comp)

## NOTE (fixes_j): wrapped in head() for consistent output length
cat("\nOja signs centered at Oja median (compMedian):\n")
signs_om <- head(ojaSign(X, center = "compMedian"))

print(signs_om)

# ---------------------------------------------------------------------------
# 4.2.3  Oja SCM (Sign Covariance Matrix)
# ---------------------------------------------------------------------------
cat("\n-----------------------------------------------------------------------\n")
cat("4.2.3  Oja SCM (Sign Covariance Matrix)\n")
cat("-----------------------------------------------------------------------\n\n")

## NOTE (fixes_j): removed comparison sections (SCM vs cov(signs), colMean SCM)
cat("Oja SCM with Oja median center (exact):\n")
scm <- ojaSCM(X)
print(scm)


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

## NOTE (fixes_j): removed redundant set.seed(1); removed signed rank test section
cat("1-sample sign test, H0: mu = c(1, 0.5), first 10 obs:\n")

res_1s_approx <- oja1sampleTest(X[1:10, ], mu = c(1, 0.5))
print(res_1s_approx)

cat("\n  Permutation method (1000 permutations):\n")

res_1s_perm <- oja1sampleTest(X[1:10, ], mu = c(1, 0.5),
                              method = "permutation")
print(res_1s_perm)

# ---------------------------------------------------------------------------
# 4.2.6  C-sample location test
# ---------------------------------------------------------------------------
cat("\n-----------------------------------------------------------------------\n")
cat("4.2.6  C-sample location test\n")
cat("-----------------------------------------------------------------------\n")

cat("C-sample rank test (permutation):\n")
res_cs_rank <- ojaCsampleTest(X ~ GROUP, scores = "rank",
                              method = "permutation")
print(res_cs_rank)