#' @title Quantile Correlation
#' @description Computes the quantile correlation between a response variable and a set of predictors.
#' @param y A numeric vector representing the response variable.
#' @param x A numeric matrix or data frame representing the predictors.
#' @param tau A numeric value representing the quantile level (between 0 and 1).
#' @return A numeric vector of quantile correlations, ordered by their absolute values in descending order.
qcor <- function(y, x, tau) {
  n = length(y)
  sum((tau - (y < quantile(y, tau))) * (x - mean(x))) /
    sqrt(n * (tau - tau^2) * sum((x - mean(x))^2))
}



#' @title Censored Quantile Correlation
#' @description Computes the censored quantile correlation between a response variable, a censoring indicator, and a set of predictors.
#' @param y A numeric vector representing the response variable.
#' @param d A numeric vector representing the censoring indicator (1 for uncensored, 0 for censored).
#' @param x A numeric matrix or data frame representing the predictors.
#' @param tau A numeric value representing the quantile level (between 0 and 1).
#' @return A numeric vector of censored quantile correlations, ordered by their absolute values in descending order.
cqcor <- function(y, d, x, tau) {
  n <- length(y)
  y <- as.numeric(y)
  x <- as.numeric(x)
  Te <- as.data.frame(cbind(y, d, x, P = NA))[order(y), ]
  lambda <- ((n - c(1:n)) / (n - c(1:n) + 1))^Te$d
  S <- cumprod(lambda)
  for (i in 1:n)  {
    if ((Te$d[i] == 0) && ((1 - S[i]) < tau)) {
      Te$P[i] <- (tau - (1 - S[i])) / S[i]
    } else {
      Te$P[i] <- (Te$y[i] < quantile(Te$y, tau))
    }

  }
  cqcov <- sum((tau - Te$P) * Te$x) / n
  covx <- sum((Te$x - mean(Te$x))^2) / n
  cqcov / sqrt((tau - tau^2) * covx)
}


#' @title Quantile Correlation Screening
#' @param y A numeric vector representing the response variable.
#' @param x A numeric matrix or data frame representing the predictors.
#' @param tau A numeric value representing the quantile level (between 0 and 1).
#' @return A numeric vector of quantile correlations, ordered by their absolute values in descending order.
qcs <- function(y, x, tau, len) {
  n <- length(y)
  cpcvec <- numeric(ncol(x))

  for (k in seq_len(ncol(x))) {
    cpcvec[k] <- qcor(y, x[, k], tau)
  }

  orders <- order(abs(cpcvec), decreasing = TRUE)

  if (is.null(len)) {
    return(orders[1:floor(n / log(n))])
  } else {
    return(orders[1:len])
  }
}

#' @title Censored Quantile Correlation Screening
#' @param y A numeric vector representing the response variable.
#' @param d A numeric vector representing the censoring indicator (1 for uncensored, 0 for censored).
#' @param x A numeric matrix or data frame representing the predictors.
#' @param tau A numeric value representing the quantile level (between 0 and 1).
#' @return A numeric vector of censored quantile correlations, ordered by their absolute values in descending order.
#'
cqcs <- function(y, d, x, tau, len) {
  n <- length(y)
  cpcvec <- numeric(ncol(x))

  for (k in seq_len(ncol(x))) {
    cpcvec[k] <- cqcor(y, d, x[, k], tau)
  }

  orders <- order(abs(cpcvec), decreasing = TRUE)

  if (is.null(len)) {
    return(orders[1:floor(n / log(n))])
  } else {
    return(orders[1:len])
  }
}

#' @title Censored Quantile Correlation Screening for Complete Data
#' @param y A numeric vector representing the response variable.
#' @param d A numeric vector representing the censoring indicator (1 for uncensored, 0 for censored).
#' @param x A numeric matrix or data frame representing the predictors.
#' @param tau A numeric value representing the quantile level (between 0 and 1).
#' @return A numeric vector of censored quantile correlations for complete data, ordered by their absolute values in descending order.
#'
coqcs <- function(y, d, x, tau, len) {
  n <- length(y)
  cpcvec <- numeric(ncol(x))

  for (k in seq_len(ncol(x))) {
    cpcvec[k] <- qcor(y[d == 1], x[d == 1, k], tau)
  }

  orders <- order(abs(cpcvec), decreasing = TRUE)

  if (is.null(len)) {
    return(orders[1:floor(n / log(n))])
  } else {
    return(orders[1:len])
  }
}
