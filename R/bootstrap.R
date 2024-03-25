#' Generate bootstrap sample of censoring time
#'
#' @param n An integer of sample size.
#' @param rev.km A \code{\link[stats]{stepfun}} of the Kaplan-Meier estimate
#' for the censoring distribution computed from \code{\link{estimate.km}}.
#'
#' @return A numeric vector of bootstrap sample of censoring time.
#' @keywords internal
resample.censor <- function(n, rev.km) {
  extract <- extract.km(rev.km)
  x.selected <- extract$x
  p.selected <- extract$p
  cens.star <- sample(x.selected, size = n, replace = TRUE, prob = p.selected)
  return(cens.star)
}

#' Generate bootstrap sample of event time
#'
#' @param n An integer of sample size.
#' @param skm.df A \code{data.frame} of the Smoothed LCM of KME
#' computed from \code{\link{smooth.kme}}.
#' @param tau.c A numeric of the end of the study time.
#'
#' @return A numeric vector of bootstrap sample of event time.
#' @importFrom stats runif
#' @keywords internal
resample.survival <- function(n, skm.df, tau.c) {
  t.star <- numeric(n)
  xx <- skm.df[["x"]]
  skm <- skm.df[["skm"]]
  u.star <- runif(n)
  for (l in 1:n) {
    t.star[l] <- xx[which(skm >= u.star[l])[1]]
  }
  t.star[which(is.na(t.star))] <- 100000
  return(t.star)
}

#' Generate bootstrap sample of smooth Grenander density estimate
#'
#'
#' @param rev.km A \code{\link[stats]{stepfun}} of the Kaplan-Meier estimate
#' for the censoring distribution computed from \code{\link{estimate.km}}.
#' @param gcmlcm.fit A list of LCM of KME returned from \code{\link[fdrtool]{gcmlcm}}.
#' @param skm.df A \code{data.frame} of the Smoothed LCM of KME
#' computed from \code{\link{smooth.kme}}.
#' @param tau.c A numeric of the end of the study time.
#' @param n An integer of sample size.
#'
#' @return A list containing bootstrap samples, including \eqn{\hat{f}_{nh}^{SG^\ast}(\tau_{G})}.
#'
#' @importFrom fdrtool gcmlcm
#' @importFrom survival survfit Surv
#' @keywords internal
smooth.grenander.bootstrap <- function(rev.km, gcmlcm.fit, skm.df, tau.c, n) {
  c.star <- resample.censor(n, rev.km)
  t.star <- resample.survival(n, skm.df, tau.c)
  y.star <- pmin(t.star, c.star)
  delta.star <- as.numeric(c.star >= t.star)

  # KME
  km.star.fit <- survfit(Surv(y.star, delta.star) ~ 1)

  # Cure fraction from the KME
  p.hat.star <- min(km.star.fit$surv)

  # Construct the LCM of 1 - KME
  km.star.dist <- 1 - km.star.fit$surv
  km.star.times <- km.star.fit$time
  if (min(km.star.times > 0)) {
    km.star.times <- c(0, km.star.times)
    km.star.dist <- c(0, km.star.dist)
  }

  gcmlcm.star.fit <- gcmlcm(km.star.times, km.star.dist, type = 'lcm')

  gren.star.times <- gcmlcm.star.fit$x.knots
  gren.star.val <- c(gcmlcm.star.fit$slope.knots,
                     gcmlcm.star.fit$slope.knots[length(gcmlcm.star.fit$slope.knots)])

  # Grenander estimator
  gren.star.fit <- function(u) f.hat(u, gren.star.times, gren.star.val)

  # Bandwidth for smoothing
  bw.star <- min(tau.c*(n^(-1/5)), tau.c / 2)
  # Smoothed Grenander estimator with boundary correction
  grid.size <- 100L
  x.grid <- seq(0, tau.c, length = grid.size)
  sg.star <- boundary.smoother(x.grid,
                               jumptimes = gcmlcm.star.fit$x.knots,
                               jumpheights = gcmlcm.star.fit$slope.knots,
                               bandwidth = bw.star, a = 0, b = tau.c)

  sg.star.tau.c <- sg.star[grid.size]
  return(list('sg.star.tau.c' = sg.star.tau.c,
              'sg.star' = data.frame('x' = x.grid, 'sg.star' = sg.star),
              'gren.star' = gren.star.fit, 'p.hat.star' = p.hat.star,
              'km.star' = km.star.fit))
}

