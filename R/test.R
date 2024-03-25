#' Testing sufficient follow-up using Smoothed Grenander Estimator \eqn{\hat{f}_{nh}^{SG}}
#'
#' Testing the hypotheses:
#' \deqn{
#'    \tilde{H}_{0}: q_{1 - \epsilon} \geq \tau_{G}
#'    \quad\text{versus}\quad
#'    \tilde{H}_{a}: q_{1 - \epsilon} < \tau_{G},
#' }
#' using the smoothed Grenander estimator, \eqn{\hat{f}_{nh}^{SG}}.
#'
#' @param y A numeric vector of the observed survival times, \eqn{Y_i}.
#' @param delta A numeric vector of the censoring indicators, \eqn{\Delta_i}.
#' @param n.boot An integer of number of bootstrap samples.
#' @param tau.c A numeric of the end of the study time.
#' \code{NULL} indicates the maximum observed survival time \code{max(y)} is used.
#' @param skm.bw.multi A numeric of scaling to the bandwidth, \eqn{h_0}, for smoothing
#' the LCM of the KME.
#'
#' @return A list containing the results including the estimate
#' \eqn{\hat{f}_{nh}^{SG}(\tau_{G})} and
#' the bootstrap samples \eqn{\hat{f}_{nh}^{SG^\ast}(\tau_{G})}.
#' @importFrom fdrtool gcmlcm
#' @importFrom survival survfit Surv
#' @keywords internal
smooth.grenander.test <- function(y, delta, n.boot = 1000L, tau.c = NULL, skm.bw.multi = 0.7) {
  if (is.null(tau.c)) {
    tau.c <- max(y)
  }
  n <- length(y)

  # KME
  km.fit <- survfit(Surv(y, delta) ~ 1)

  # Cure fraction from the KME
  p.hat <- min(km.fit$surv)

  # Construct the LCM of 1 - KME
  km.dist <- 1 - km.fit$surv
  km.times <- km.fit$time

  gcmlcm.fit <- gcmlcm(c(0, km.times), c(0, km.dist), type = 'lcm')


  gren.times <- gcmlcm.fit$x.knots
  gren.val <- c(gcmlcm.fit$slope.knots,
                gcmlcm.fit$slope.knots[length(gcmlcm.fit$slope.knots)])

  # Grenander estimator
  gren.fit <- function(u) f.hat(u, gren.times, gren.val)

  # Bandwidth for smoothing
  bw <- min(tau.c*(n^(-1/5)), tau.c / 2)
  # Smoothed Grenander estimator with boundary correction
  grid.size <- 100L
  sg <- boundary.smoother(seq(0, tau.c, length = grid.size),
                          jumptimes = gcmlcm.fit$x.knots,
                          jumpheights = gcmlcm.fit$slope.knots,
                          bandwidth = bw, a = 0, b = tau.c)

  sg.tau.c <- sg[grid.size]

  # Smooth KME
  skm.df <- smooth.kme(gcmlcm.fit, tau.c, n, skm.bw.multi)
  skm.grid.size <- nrow(skm.df)
  f.star.tau.c <-
    (skm.df[["skm"]][skm.grid.size] - skm.df[["skm"]][skm.grid.size - 1]) /
    (skm.df[["x"]][skm.grid.size] - skm.df[["x"]][skm.grid.size - 1])

  # Reverse KME for estimating the censoring distribution
  rev.km <- estimate.km(y, 1L - delta)

  b.res <- lapply(seq(n.boot), function(j) {
    if (j %% 100 == 0) cat("\r", sprintf("Bootstrap iteration: %d / %d", j, n.boot))
    seed.out <- .Random.seed
    boot.res <- smooth.grenander.bootstrap(rev.km, gcmlcm.fit, skm.df, tau.c, n)
    return(list("bs.iter" = j, "bs.seed" = seed.out,
                "sg.star.tau.c" = boot.res$sg.star.tau.c))
  })
  cat("\n")

  sg.tau.c.star <- sapply(b.res, function(b) {
    return(b$sg.star.tau.c)
  })

  output.list <- list("sg.tau.c.fit" = sg.tau.c,
                      "sg.tau.c.star" = sg.tau.c.star,
                      "km.fit" = km.fit, "rev.km.fit" = rev.km,
                      "skm.df" = skm.df, "f.star.tau.c" = f.star.tau.c,
                      "boot.res" = b.res, "gren.fit" = gren.fit,
                      "n" = n, "tau.c" = tau.c)
  return(output.list)
}


#' Lower bound of the density at \eqn{q_{1-\epsilon}}
#'
#' The lower bound of \eqn{f(q_{1-\epsilon})}:
#' \deqn{{f(q_{1-\epsilon}) \geq
#'        \frac{(\epsilon - \eta)p}{\tau - \tau_{G}} \approx
#'        \frac{\epsilon p}{\tau - \tau_{G}}}.}
#'
#' @param p.hat A numeric of the estimated cure fraction \eqn{p} from KME.
#' @param tau A numeric of the parameter \eqn{\tau} for the test.
#' @param tau.c A numeric of the end of the study time.
#' @param eps A numeric of the \eqn{\epsilon} for the test.
#'
#' @return A numeric of the lower bound
#' @keywords internal
tail.approx <- function(p.hat, tau, tau.c, eps = 0.05) {
  return(eps * (1 - p.hat) / (tau - tau.c))
}


#' Compute critical value from the bootstrap samples for the test using \eqn{\hat{f}_{nh}^{SG}}
#'
#' Compute critical value from the bootstrap samples:
#' \deqn{
#'  \frac{\epsilon \hat{F}_{n}(\tau_{G})}{\tau - \tau_{G}} + n^{-2/5}Q_{\alpha}^{SG},
#' }
#' where \eqn{\hat{F}_{n}} is the KME,
#' \eqn{Q_{\alpha}^{SG}} is the \eqn{\alpha}-quantile of the bootstrap estimates
#' \eqn{n^{2/5}\{\hat{f}_{nh}^{SG^{\ast}}(\tau_{G})-\tilde{f}_{nh_{0}}(\tau_{G})\}},
#' with \eqn{\tilde{f}_{nh_{0}}} being the derivative of the smoothed LCM of KME from
#' \code{\link{smooth.kme}}.
#'
#' @param sg A list from \code{\link{smooth.grenander.test}}.
#' @param tau A numeric of the parameter \eqn{\tau} for the test.
#' @param tau.c A numeric of the end of the study time.
#' @param eps A numeric of the \eqn{\epsilon} for the test.
#' @param alpha A numeric of the significance level.
#'
#' @importFrom stats quantile
#' @return A numeric of the computed critical value
#' @keywords internal
smooth.grenander.crit.val <- function(sg, tau, tau.c, eps = 0.01, alpha = 0.05) {
  p.hat <- min(sg[["km.fit"]]$surv)
  q.alpha <- quantile(sg[["sg.tau.c.star"]] - sg[["f.star.tau.c"]], alpha)
  crit.val <- tail.approx(p.hat, tau, tau.c, eps) + q.alpha
  return(crit.val)
}


#' Compute \eqn{p}-value from the bootstrap samples for the test \eqn{\hat{f}_{nh}^{SG}}
#'
#' @param sg A list from \code{\link{smooth.grenander.test}}.
#' @param tau A numeric of the parameter \eqn{\tau} for the test.
#' @param tau.c A numeric of the end of the study time.
#' @param eps A numeric of the \eqn{\epsilon} for the test.
#'
#' @return A numeric of the computed \eqn{p}-value
#' @keywords internal
smooth.grenander.p.val <- function(sg, tau, tau.c, eps = 0.01) {
  p.hat <- min(sg[["km.fit"]]$surv)
  sg.star <- sg[["sg.tau.c.star"]] - sg[["f.star.tau.c"]]
  t.val <- sg[["sg.tau.c.fit"]] - tail.approx(p.hat, tau, tau.c, eps)
  p.val <- mean(sg.star < t.val)
  return(p.val)
}


#' Testing sufficient follow-up using Grenander Estimator \eqn{\hat{f}^{G}}
#'
#' Testing the hypotheses:
#' \deqn{
#'    \tilde{H}_{0}: q_{1 - \epsilon} \geq \tau_{G}
#'    \quad\text{versus}\quad
#'    \tilde{H}_{a}: q_{1 - \epsilon} < \tau_{G},
#' }
#' using the Grenander estimator, \eqn{\hat{f}^{G}}.
#'
#' @param y A numeric vector of the observed survival times, \eqn{Y_i}.
#' @param delta A numeric vector of the censoring indicators, \eqn{\Delta_i}.
#' @param tau.c A numeric of the end of the study time.
#' \code{NULL} indicates the maximum observed survival time \code{max(y)} is used.
#'
#' @return A list containing the results including the estimate
#' \eqn{\hat{f}^{G}(\tau_{G})}.
#' @importFrom fdrtool gcmlcm
#' @importFrom survival survfit Surv
#' @keywords internal
grenander.test <- function(y, delta, tau.c = NULL) {
  if (is.null(tau.c)) {
    # Estimate tau.c by max{y}
    tau.c <- max(y)
  }
  n <- length(y)

  # KME
  km.fit <- survfit(Surv(y, delta) ~ 1)

  # Cure fraction from the KME
  p.hat <- min(km.fit$surv)

  # Construct the LCM of 1 - KME
  km.dist <- 1 - km.fit$surv
  km.times <- km.fit$time

  gcmlcm.fit <- gcmlcm(c(0, km.times), c(0, km.dist), type = 'lcm')


  gren.times <- gcmlcm.fit$x.knots
  gren.val <- c(gcmlcm.fit$slope.knots,
                gcmlcm.fit$slope.knots[length(gcmlcm.fit$slope.knots)])

  # Grenander estimator
  gren.fit <- function(u) f.hat(u, gren.times, gren.val)

  # Modified reverse KME for estimating the censoring distribution
  rev.km <- estimate.km(y, 1L - delta, TRUE)

  output.list <- list("gren.fit" = gren.fit,
                      "km.fit" = km.fit, "rev.km" = rev.km,
                      "n" = n, "tau.c" = tau.c)
  return(output.list)
}


#' Compute the test statistic for the test using \eqn{\hat{f}^{G}}
#'
#' Compute the test statistic for the test using \eqn{\hat{f}^{G}}:
#' \deqn{\hat{f}_{n}^{G}(\tau_{G} - cn^{-a})}.
#'
#' @param gren A list from \code{\link{grenander.test}}.
#' @param a A numeric of the constant \eqn{a\in(1/3,1)} controls
#' the rate of convergence of the estimator.
#' @param c.shift A positive number of the constant \eqn{c}.
#' @param eps A numeric of the \eqn{\epsilon} for the test.
#'
#'
#' @return A numeric of the computed test statistic.
#' @keywords internal
grenander.stat <- function(gren, a, c.shift = 1, eps = 0.01) {
  n <- gren[["n"]]
  tau.c <- gren[["tau.c"]]
  gren.fit <- gren[["gren.fit"]]
  x <- tau.c - c.shift / (n ^ a)
  gren.fit.tau.c <- gren.fit(x)
  gren.fit.tau.c <- ifelse(abs(gren.fit.tau.c) < .Machine$double.eps, 0, gren.fit.tau.c)
  return(gren.fit.tau.c)
}


#' Compute critical value for the test using \eqn{\hat{f}^{G}}
#'
#' Compute critical value:
#' \deqn{
#'  \frac{\epsilon\hat{F}_{n}(\tau_{G})}{\tau - \tau_{G}} - A_{1}^{-1}n^{-(1-a)/2}Q_{1 - \alpha}^{G}
#' }
#' where \eqn{\hat{F}_{n}} is the KME,
#' \eqn{A_1=\sqrt{c[1-G(\tau_G-)]/f(\tau_G)}}, and
#' \eqn{Q_{\alpha}^{G}} is the \eqn{\alpha}-quantile of the distribution \eqn{D_R[W(t)](1)},
#' with \eqn{W(t)} denotes a Brownian motion and
#' \eqn{D_R[Z(t)](b)} is the right derivative of the LCM on
#' \eqn{[0,\infty)} of the process \eqn{Z(t)} at the point \eqn{t=b}.
#'
#' @param gren A list from \code{\link{grenander.test}}.
#' @param tau A numeric of the parameter \eqn{\tau} for the test.
#' @param c.shift A positive number of the constant \eqn{c}.
#' @param eps A numeric of the \eqn{\epsilon} for the test.
#' @param alpha A numeric of the significance level.
#'
#' @importFrom stats quantile
#' @return A numeric of the computed critical value
#' @keywords internal
grenander.crit.val <- function(gren, a, tau, c.shift = 1, eps = 0.01, alpha = 0.05) {
  n <- gren[["n"]]
  tau.c <- gren[["tau.c"]]
  p.hat <- min(gren[["km.fit"]]$surv)
  rev.kme <- gren[["rev.km"]]
  rev.kme.tau.c.minus <- rev.kme(tau.c - sqrt(.Machine$double.eps))

  gren.fit.tau.c <- grenander.stat(gren, a, c.shift, eps)
  q <- quantile(dw, 1 - alpha)
  q <- q * sqrt(gren.fit.tau.c) / (n ^ (0.5 * (1 - a)) *
                                     sqrt(c.shift * (1 - rev.kme.tau.c.minus)))
  crit.val <- tail.approx(p.hat, tau, tau.c, eps) - q
  return(crit.val)
}


#' Compute \eqn{p}-value for the test using \eqn{\hat{f}^{G}}
#'
#'
#' @param gren A list from \code{\link{grenander.test}}.
#' @param tau A numeric of the parameter \eqn{\tau} for the test.
#' @param c.shift A positive number of the constant \eqn{c}.
#' @param eps A numeric of the \eqn{\epsilon} for the test.
#'
#' @importFrom stats quantile
#' @return A numeric of the computed \eqn{p}-value
#' @keywords internal
grenander.p.val <- function(gren, a, tau, c.shift = 1, eps = 0.01) {
  n <- gren[["n"]]
  tau.c <- gren[["tau.c"]]
  p.hat <- min(gren[["km.fit"]]$surv)
  rev.kme <- gren[["rev.km"]]
  rev.kme.tau.c.minus <- rev.kme(tau.c - sqrt(.Machine$double.eps))
  norm.const <- n ^ (0.5 * (1 - a)) * sqrt(c.shift * (1 - rev.kme.tau.c.minus))

  gren.fit.tau.c <- grenander.stat(gren, a, c.shift, eps)
  t.val <- tail.approx(p.hat, tau, tau.c, eps) - gren.fit.tau.c
  p.val <- mean(dw * sqrt(gren.fit.tau.c) / norm.const > t.val)

  return(p.val)
}
