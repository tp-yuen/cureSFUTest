#' Test for sufficient follow-up in survival data with a cure fraction
#' using shape constrained density estimator
#'
#' Test for sufficient follow-up in survival data with a cure fraction
#' using the smoothed Grenander estimator \eqn{\hat{f}_{nh}^{SG}} or
#' the Grenander estimator \eqn{\hat{f}_{n}^{G}}.
#' Let \eqn{F_{u}} be the distribution of the survival time for the uncured subject and
#' \eqn{G} be the distribution of the random right censoring time.
#' The test hypothesis is:
#' \deqn{
#'    \tilde{H}_{0}: q_{1 - \epsilon} \geq \tau_{G}
#'    \quad\text{versus}\quad
#'    \tilde{H}_{a}: q_{1 - \epsilon} < \tau_{G},
#' }
#' where \eqn{q_{1 - \epsilon} = \inf\lbrace{t \geq 0}:{F_{u}(t) \geq 1 - \epsilon}\rbrace} and
#' \eqn{\tau_{G}=\sup\lbrace{t\geq0}:{G(t)<1}\rbrace}.
#'
#' @param y a numeric vector of the observed survival times, \eqn{Y_i}.
#' @param delta a numeric vector of the censoring indicators, \eqn{\Delta_i}.
#' @param tau a numeric of the parameter \eqn{\tau} for the test.
#'  It must be greater than \code{tau.g}. See 'Details' below.
#' @param tau.g a numeric of the end of the study time, \eqn{\tau_{G}}.
#'  See 'Details' for the meaning of \code{NULL} (the default).
#' @param eps a numeric in \eqn{(0, 1)} of the \eqn{\epsilon} for the test.
#'  Default is \code{0.01}.
#' @param method a character string indicating which shape constrained density
#'  estimator is to be used for the test.
#'  One of "g" (indicating Grenander estimator \eqn{\hat{f}_{n}^{G}}), or
#'  "sg" (smoothed Grenander estimator \eqn{\hat{f}_{nh}^{SG}}).
#' @param alpha A numeric in \eqn{(0, 1)} of the significance level for the test.
#'  Default is \code{0.05}.
#' @param n.boot An integer of number of bootstrap samples when
#'  the smoothed Grenander estimator \eqn{\hat{f}_{nh}^{SG}} is used, i.e. \code{method = "sg"}.
#'  Default is \code{1000L}.
#'
#'
#' @details
#'  If \code{tau.g} is \code{NULL},
#'  the maximum observed survival time \code{max(y)} is used for the test.
#'  The choice of \eqn{\tau} can be based on some prior knowledge, for example,
#'  one can take \eqn{\tau > \tau_{G}} such that it is almost impossible
#'  for the event to happen after \eqn{\tau}.
#'
#' @return A list with class \code{sfu} containing the following components:
#' \describe{
#' \item{\code{statistic}}{the value of the test statistic.}
#' \item{\code{p.value}}{the \eqn{p}-value of the test.}
#' \item{\code{crit.val}}{the critical value of the test using the significance level \code{alpha}.}
#' \item{\code{method}}{a character string indicating how the density was estimated.}
#' \item{\code{alpha}}{the value of the significance level used.}
#' \item{\code{tau.g}}{the value of \code{tau.g} used.}
#' \item{\code{tau}}{the value of \code{tau} used.}
#' \item{\code{eps}}{the value of \code{eps} used.}
#' \item{\code{data.name}}{a character string giving the names of the data.}
#' }
#' @export
sfu.test <- function(y, delta, tau, tau.g = NULL, eps = 0.01,
                     method = c("g", "sg"), alpha = 0.05, n.boot = 1000L) {
  method <- match.arg(method)
  DNAME <- paste0("(", deparse(substitute(y)), ", ", deparse(substitute(delta)), ")")

  if(length(y) != length(delta))
    stop("'y' and 'delta' must have the same length")
  if(!is.numeric(y)) stop("'y' must be a numeric vector")
  if(!is.numeric(delta)) stop("'delta' must be a numeric vector")
  if(is.null(tau.g)) tau.g <- max(y)
  if(tau <= tau.g) stop("'tau' must be greater than 'tau.g'")
  if(eps <= 0 || eps >= 1) stop("'eps' must be in (0, 1)")
  if(alpha <= 0 || alpha >= 1) stop("'alpha' must be in (0, 1)")

  switch (
    method,
    "sg" = {
      res <- smooth.grenander.test(y = y, delta = delta, n.boot = n.boot, tau.c = tau.g)
      test.stat <- res[["sg.tau.c.fit"]]
      crit.val <- smooth.grenander.crit.val(res, tau = tau, tau.c = tau.g, eps = eps, alpha = alpha)
      p.val <- smooth.grenander.p.val(res, tau = tau, tau.c = tau.g, eps = eps)
    },
    "g" = {
      gren.a <- 0.34
      res <- grenander.test(y = y, delta = delta, tau.c = tau.g)
      test.stat <- grenander.stat(res, a = gren.a, c.shift = tau.g, eps = eps)
      crit.val <- grenander.crit.val(res, gren.a, tau, c.shift = tau.g, eps = eps, alpha = alpha)
      p.val <- grenander.p.val(res, gren.a, tau, c.shift = tau.g, eps = eps)
    }
  )
  sfu.res <- list(
    "statistic" = test.stat,
    "p.value" = p.val,
    "method" = method,
    "crit.val" = crit.val,
    "alpha" = alpha,
    "tau.g" = tau.g,
    "tau" = tau,
    "eps" = eps,
    "data.name" = DNAME
  )
  class(sfu.res) <- "sfu"
  return(sfu.res)
}


#' Print a sfu object
#'
#' Print a summary of the sufficient follow-up test result.
#'
#' @details
#' The value of the test statistic, p-value of the test, and the setting are printed.
#'
#' @param x An \code{sfu} object.
#' @param \dots Additional print arguments.
#' @return A printout mentioned in Details.
#' @method print sfu
#' @export
print.sfu <- function (x, ...)
{
  cat("\t", sprintf("Test for sufficient follow-up using %sGrenander estimator",
                    ifelse(x[["method"]] == "sg", "Smoothed ", "")), "\n\n")
  cat("data:\t", x[["data.name"]], "\n")
  cat(sprintf("test stat. = %0.4f, p-value = %0.4f", x[["statistic"]], x[["p.value"]]), "\n")
  cat(sprintf("alternative hypothesis: q(%g) < tau.g (sufficient follow-up)",
              1 - x[["eps"]]), "\n")
  cat("setting:\n")
  print(c("eps" = x[["eps"]], "tau.g" = x[["tau.g"]], "tau" = x[["tau"]]))
}


#' Test for sufficient follow-up in survival data with a cure fraction
#' and categorical covariates
#' using shape constrained density estimator
#'
#' Test for sufficient follow-up in survival data with a cure fraction
#' and categorical covariates \eqn{X}
#' using the smoothed Grenander estimator \eqn{\hat{f}_{nh}^{SG}}.
#' Let \eqn{F_{u}(\cdot|x)} be the conditional distribution of
#' the survival time for the uncured subject given \eqn{X=x} and
#' \eqn{G(\cdot|x)} be the conditional distribution of
#' the random right censoring time given \eqn{X=x}.
#' The test hypothesis is:
#' \deqn{
#'    {H}_{0}: q_{1 - \epsilon}(x) \geq \tau_{G}(x)\text{ for some }x
#'    \quad\text{versus}\quad
#'    {H}_{1}: q_{1 - \epsilon}(x) < \tau_{G}(x)\text{ for all }x,
#' }
#' where \eqn{q_{1 - \epsilon}(x) = \inf\lbrace{t \geq 0}:{F_{u}(t|x) \geq 1 - \epsilon}\rbrace} and
#' \eqn{\tau_{G}(x)=\sup\lbrace{t\geq0}:{G(t|x)<1}\rbrace}.
#' The test hypothesis for fixed covariate value \eqn{x} is:
#' \deqn{
#'    {H}_{0x}: q_{1 - \epsilon}(x) \geq \tau_{G}(x)
#'    \quad\text{versus}\quad
#'    {H}_{1x}: q_{1 - \epsilon}(x) < \tau_{G}(x),
#' }
#'
#' @param y a numeric vector of the observed survival times, \eqn{Y}.
#' @param delta a numeric vector of the censoring indicators, \eqn{\Delta}.
#' @param X a numeric vector of the categorical covariates, \eqn{X}.
#' @param tau a numeric of the parameter \eqn{\tau} for the test.
#'  It must be greater than the observed survival time \code{max(y)}. See 'Details' below.
#' @param eps a numeric in \eqn{(0, 1)} of the \eqn{\epsilon} for the test.
#'  Default is \code{0.01}.
#' @param alpha A numeric in \eqn{(0, 1)} of the significance level for the test.
#'  Default is \code{0.05}.
#' @param gamma A numeric in \eqn{(0, 1)} of the parameter used to estimate \eqn{x_\ast}.
#' @param n.boot An integer of number of bootstrap samples. Default is \code{1000L}.
#'
#'
#' @details
#'  The choice of \eqn{\tau} can be based on some prior knowledge, for example,
#'  one can take \eqn{\tau > \tau_{G}(x)} such that it is almost impossible
#'  for the event to happen after \eqn{\tau} given \eqn{X = x} for all \eqn{x}.
#'
#' @return A list with class \code{sfu.cov} containing the following components:
#' \describe{
#' \item{\code{statistic}}{the value of the test statistic for \eqn{H_{0x}}.}
#' \item{\code{p.value}}{the \eqn{p}-value of the test for \eqn{H_{0x}}.}
#' \item{\code{crit.val}}{the critical value of the test for \eqn{H_{0x}} using the significance level \code{alpha}.}
#' \item{\code{alpha}}{the value of the significance level used.}
#' \item{\code{y.max.x}}{the value of \eqn{Y_{(n)}^x}.}
#' \item{\code{x.ast.hat}}{the value of \eqn{\hat{x}_{\ast,n}}.}
#' \item{\code{tau}}{the value of \code{tau} used.}
#' \item{\code{eps}}{the value of \code{eps} used.}
#' \item{\code{x.unique}}{the unique value of \code{X}.}
#' \item{\code{data.name}}{a character string giving the names of the data.}
#' }
#' @export
sfu.cov.test <- function(y, delta, X, tau, eps = 0.01,
                         alpha = 0.05, gamma = 0.025, n.boot = 1000L) {
  DNAME <- paste0("(", deparse(substitute(y)), ", ",
                  deparse(substitute(delta)), ", ", deparse(substitute(X)), ")")

  if(length(y) != length(delta))
    stop("'y' and 'delta' must have the same length")
  if(!is.numeric(y)) stop("'y' must be a numeric vector")
  if(!is.numeric(delta)) stop("'delta' must be a numeric vector")
  if(!is.numeric(X)) stop("'X' must be a numeric matrix")
  if(length(y) != length(X)) stop("The number of rows 'X' must equal to the length of 'y'")
  if(tau <= max(y)) stop("'tau' must be greater than 'max(y)'")
  if(eps <= 0 || eps >= 1) stop("'eps' must be in (0, 1)")
  if(alpha <= 0 || alpha >= 1) stop("'alpha' must be in (0, 1)")

  d <- sort(unique(X))
  ord <- order(y, 1 - delta)
  y.ord <- y[ord]
  status.ord <- delta[ord]
  X.ord <- X[ord]
  res <- lapply(d, function(xx) {
    X.ord.xx.idx <- which(X.ord == xx)
    y.max.xx <- max(y.ord[X.ord.xx.idx])
    sg.xx.test <- smooth.grenander.test.undersmooth(
      y.ord[X.ord.xx.idx], status.ord[X.ord.xx.idx], n.boot = n.boot,
      tau.c = y.max.xx, skm.bw.multi = 1.0)
    tail.xx <- tail.approx(min(sg.xx.test[["km.fit"]]$surv), tau, y.max.xx, eps)
    test.stat.xx <- sg.xx.test[["sg.tau.c.fit"]] - tail.xx
    crit.val.xx <- smooth.grenander.crit.val(sg.xx.test, tau = tau,
                                             tau.c = y.max.xx, eps = eps,
                                             alpha = alpha) - tail.xx
    p.val.xx <- smooth.grenander.p.val(sg.xx.test, tau = tau, tau.c = y.max.xx,
                                       eps = eps)
    q.ast <- quantile(
      sg.xx.test[["sg.tau.c.star"]] -
        eps * (1 - sg.xx.test[['p.hat.star']]) / (tau - sg.xx.test[['tau.c.star']]),
      1 - gamma)
    return(list('x' = xx, 'res.x' = sg.xx.test, 'y.max.x' = y.max.xx,
                'test.stat.x' = test.stat.xx, 'crit.val.x' = crit.val.xx,
                'p.val.x' = p.val.xx, 'q.ast' = q.ast))
  })

  test.stats.x <- sapply(res, function(r) {
    return(r[['test.stat.x']])
  })
  p.val.x <- sapply(res, function(r) {
    return(r[['p.val.x']])
  })
  crit.val.x <- sapply(res, function(r) {
    return(r[['crit.val.x']])
  })
  y.max.x <- sapply(res, function(r) {
    return(r[['y.max.x']])
  })
  q.ast <- sapply(res, function(r) {
    return(r[['q.ast']])
  })
  x.ast.hat <- which.max(q.ast)

  sfu.res <- list(
    "statistic" = test.stats.x,
    "p.value" = p.val.x,
    "crit.val" = crit.val.x,
    "alpha" = alpha,
    "y.max.x" = y.max.x,
    "x.ast.hat" = x.ast.hat,
    "tau" = tau,
    "eps" = eps,
    "x.unique" = d,
    "data.name" = DNAME
  )
  class(sfu.res) <- "sfu.cov"
  return(sfu.res)
}

#' Print a sfu.cov object
#'
#' Print a summary of the sufficient follow-up with categorical covariates test result.
#' Let \eqn{F_{u}(\cdot|x)} be the conditional distribution of
#' the survival time for the uncured subject given \eqn{X=x} and
#' \eqn{G(\cdot|x)} be the conditional distribution of
#' the random right censoring time given \eqn{X=x}.
#' The test hypothesis is:
#' \deqn{
#'    {H}_{0}: q_{1 - \epsilon}(x) \geq \tau_{G}(x)\text{ for some }x
#'    \quad\text{versus}\quad
#'    {H}_{1}: q_{1 - \epsilon}(x) < \tau_{G}(x)\text{ for all }x,
#' }
#' where \eqn{q_{1 - \epsilon}(x) = \inf\lbrace{t \geq 0}:{F_{u}(t|x) \geq 1 - \epsilon}\rbrace} and
#' \eqn{\tau_{G}(x)=\sup\lbrace{t\geq0}:{G(t|x)<1}\rbrace}.
#' The test hypothesis for fixed covariate value \eqn{x} is:
#' \deqn{
#'    {H}_{0x}: q_{1 - \epsilon}(x) \geq \tau_{G}(x)
#'    \quad\text{versus}\quad
#'    {H}_{1x}: q_{1 - \epsilon}(x) < \tau_{G}(x),
#' }
#'
#' @details
#' The value of the test statistic, p-value of the test \eqn{H_{0x}},
#' the estimate \eqn{\hat{x}_{\ast,n}} and the setting are printed.
#'
#' @param x An \code{sfu.cov} object.
#' @param \dots Additional print arguments.
#' @return A printout mentioned in Details.
#' @method print sfu.cov
#' @export
print.sfu.cov <- function (x, ...)
{
  cat("\t", "Test for sufficient follow-up using Smoothed Grenander estimator", "\n\n")
  cat("data:\t", x[["data.name"]], "\n")
  test.H0x.res.df <- data.frame(
    "x" = x[['x.unique']],
    "test stat for H0x" = x[['statistic']],
    "p-value for H0x" = x[['p.value']],
    "x.ast.hat" = rep("", length(x[['x.unique']]))
  )
  colnames(test.H0x.res.df) <- c("x", "test stat for H0x",
                                 "p-value for H0x", "x.ast.hat")
  test.H0x.res.df[x[['x.ast.hat']], "x.ast.hat"] <- "*"
  print(test.H0x.res.df, row.names = FALSE)
  cat("\n")
  cat(sprintf("alternative hypothesis: q(%g)(x) < tau.g(x) for all x (sufficient follow-up)",
              1 - x[["eps"]]), "\n")
  cat("setting:\n")
  print(c("eps" = x[["eps"]], "tau" = x[["tau"]]))
}
