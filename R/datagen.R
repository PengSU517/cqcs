#' Data generation function for survival analysis
#' This function generates survival data based on specified parameters.
#' @param n Number of samples to generate.
#' @param p Number of features.
#' @param tau Quantile level for the error term.
#' @param distri Distribution of the error term, either "normal" or "cauchy".
#' @param heter Heterogeneity parameter for the error term.

datagen = function(n, p, tau, distri, heter) {

  # n = 100
  # p = 10
  # tau = 0.5
  # distri = "normal" # "cauchy"
  # heter = 1

  beta = c(1, 0, -1, 0, 1)
  lambda = c(1,0,1,0,1)
  mean = rep(0, p)
  # sigmax = diag(rep(1, p))
  # r = 0.5
  # for (i in 1:p) {
  #   for (j in 1:p) {
  #     if (i != j)
  #       sigmax[i, j] = sqrt(sigmax[i, i] * sigmax[j, j]) * r^abs(i - j)
  #   }
  # }

  ## generate datasets
  x = matrix(rnorm(n * p), nrow = n, ncol = p)
  # x = MASS::mvrnorm(n, mean, sigmax)
  if (distri == "cauchy") {
    disfunc = rcauchy
    quantfunc = qcauchy
  } else {
    disfunc = rnorm
    quantfunc = qnorm
  }

  eps = (disfunc(n, 0, 1) - quantfunc(tau)) *  exp(heter * (x[,1:5] %*% lambda)) +  quantfunc(tau)
  t = x[,1:5] %*% beta + eps
  c = runif(n, 0, 3) ## generate censoring time
  y = pmin(t, c)
  d = as.numeric(t < c)
  # sum(d)

  return(list(
    y = y,
    d =  d,
    x = x,
    t = t,
    c = c,
    rate = 1 - sum(d) / n,
    beta = beta,
    lambda = lambda
  ))
}
