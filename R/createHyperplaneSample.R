`createHyperplaneSample` <- function(N, pN, n, k) {
  
  FailureMessage <- "Can't allocate enough memory: matrix 'X' is too large."
  
  if (pN < N) {
    
    # Draw pN distinct numbers from 1,...,N
    Sample <- suppressWarnings(
      try(sample(N, size = pN, replace = FALSE), silent = TRUE)
    )
    
    if (inherits(Sample, "try-error")) {
      Sample <- suppressWarnings(
        try(floor((N * runif(n = pN)) + 0.5), silent = TRUE)
      )
    }
    
    if (inherits(Sample, "try-error")) {
      stop(FailureMessage)
    }
    
    CombinationsSample <- suppressWarnings(
      try(
        t(combn(x = n, m = k)[, Sample, drop = FALSE]),
        silent = TRUE
      )
    )
    
    if (inherits(CombinationsSample, "try-error")) {
      CombinationsSample <- suppressWarnings(
        try(
          combinationsSample(n = n, k = k, vec = Sample),
          silent = TRUE
        )
      )
    }
    
  } else {
    
    CombinationsSample <- suppressWarnings(
      try(
        t(combn(x = n, m = k)),
        silent = TRUE
      )
    )
    
  }
  
  if (inherits(CombinationsSample, "try-error")) {
    stop(FailureMessage)
  }
  
  CombinationsSample
}