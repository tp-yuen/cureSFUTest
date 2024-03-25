#' Tri-weight kernel function
#'
#' Tri-weight kernel function:
#' \eqn{k(u) = \frac{35}{32}(1-u^2)^3I_{[-1,1]}(u)}.
#'
#' @param x A numeric vector or matrix of points for which the value of
#' tri-weight kernel is computed.
#'
#' @return A numeric vector or matrix of the computed values.
#' @keywords internal
kern <- function(x) {
  K <- rep(0, length(x))
  is.nz <- (x > -1 & x < 1)
  K[is.nz] <- (35/32)*((1-x[is.nz]^2)^3)
  if (is.matrix(x)) {
    K <- matrix(K, nrow(x), ncol(x))
  }
  return(K)
}


#' Integrated tri-weight kernel function
#'
#' Integrated tri-weight kernel function:
#' \deqn{\int_{-\infty}^x k(u)du,}
#' where \eqn{k(u) = \frac{35}{32}(1-u^2)^3I_{[-1,1]}(u)}.
#'
#' @param x A numeric vector or matrix of points for which the value of
#' \code{integrated.kern} is computed.
#'
#' @return A numeric vector or matrix of the computed values.
#' @keywords internal
integrated.kern <- function(x){
  K <- rep(0, length(x))
  is.nz <- (x >= -1 & x <= 1)
  x.sel <- x[is.nz]
  K[is.nz] <- ((35/32)*(x.sel-x.sel^3+3*x.sel^5/5-x.sel^7/7)+0.5)
  K <- K + (x > 1)
  if (is.matrix(x)) {
    K <- matrix(K, nrow(x), ncol(x))
  }
  return(K)
}


#' Coefficients for the boundary kernel
#'
#' Compute the following:
#' \deqn{\int_{-\infty}^s k(u)du,}
#' where \eqn{k(u) = \frac{35}{32}(1-u^2)^3I_{[-1,1]}(u)}.
#'
#' @param s A numeric vector or matrix of points for which the value of
#'          the above is computed.
#'
#' @return A numeric vector or matrix of the computed values.
#' @keywords internal
i0 <- function(s) {
  i0 <- integrated.kern(s)
  return(i0)
}


#' Coefficients for the boundary kernel
#'
#' Compute the following:
#' \deqn{\int_{-\infty}^s uk(u)du,}
#' where \eqn{k(u) = \frac{35}{32}(1-u^2)^3I_{[-1,1]}(u)}.
#'
#' @param s A numeric vector or matrix of points for which the value of
#'          the above is computed.
#'
#' @return A numeric vector or matrix of the computed values.
#' @keywords internal
i1 <- function(s){
  K <- rep(0, length(s))
  is.nz <- (s >= -1 & s <= 1)
  s.sel <- s[is.nz]
  K[is.nz] <- ((35/32)*(s.sel^2/2-3*s.sel^4/4+s.sel^6/2-s.sel^8/8-1/8))
  if (is.matrix(s)) {
    K <- matrix(K, nrow(s), ncol(s))
  }
  return(K)
}


#' Coefficients for the boundary kernel
#'
#' Compute the following:
#' \deqn{\int_{-\infty}^s u^2k(u)du,}
#' where \eqn{k(u) = \frac{35}{32}(1-u^2)^3I_{[-1,1]}(u)}.
#'
#' @param s A numeric vector or matrix of points for which the value of
#'          the above is computed.
#'
#' @return A numeric vector or matrix of the computed values.
#' @keywords internal
i2 <- function(s) {
  K <- rep(0, length(s))
  is.nz <- (s >= -1 & s <= 1)
  s.sel <- s[is.nz]
  K[is.nz] <- (35/32)*((s.sel^3/3-3*s.sel^5/5+3*s.sel^7/7-s.sel^9/9)+16/(35*9))
  K <- K + (s > 1) / 9
  if (is.matrix(s)) {
    K <- matrix(K, nrow(s), ncol(s))
  }
  return(K)
}


#' Coefficients for the boundary kernel
#'
#' Compute the coefficients \eqn{\phi(x)} for \eqn{k(u)}.
#'
#' @param x A numeric vector or matrix of points for which the value of
#'          the above is computed.
#'
#' @return A numeric vector or matrix of the computed values.
#' @keywords internal
phi <- function(x) {
  return(i2(x)/(i0(x)*i2(x)-(i1(x))^2))
}


#' Coefficients for the boundary kernel
#'
#' Compute the coefficients \eqn{\psi(x)} for \eqn{uk(u)}.
#'
#' @param x A numeric vector or matrix of points for which the value of
#'          the above is computed.
#'
#' @return A numeric vector or matrix of the computed values.
#' @keywords internal
psi <- function(x) {
  return(-i1(x)/(i0(x)*i2(x)-(i1(x))^2))
}


#' Kernel smoothing with boundary correction using tri-weight kernel
#'
#' @param x A numeric vector of points at which the smoothing is computed.
#' @param jumptimes A numeric vector of points at which there is a jump.
#' @param jumpheights A numeric vector of jump heights.
#' @param bandwidth A positive number of the bandwidth parameter.
#' @param a A numeric of left endpoint of the support.
#' @param b A numeric of right endpoint of the support.
#'
#' @return A numeric vector of the computed values at \code{x}.
#' @keywords internal
boundary.smoother <- function(x, jumptimes, jumpheights, bandwidth, a, b) {
  x.smoothed <- rep(NA, length(x))
  is.lb <- (x < a + bandwidth)
  is.not.lb.ub <- !is.lb
  is.ub <- TRUE
  if (!is.na(b)) {
    is.ub <- (x > b - bandwidth)
    is.not.lb.ub <- (is.not.lb.ub & !is.ub)
  }

  # No correction
  xt <- outer(x, jumptimes, "-") / bandwidth
  if (any(is.not.lb.ub)) {
    k.xt <- integrated.kern(xt[is.not.lb.ub, ])
    if (sum(is.not.lb.ub) > 1) {
      x.smoothed[is.not.lb.ub] <- -1 * drop(jumpheights %*% diff(t(k.xt)))
    } else {
      x.smoothed[is.not.lb.ub] <- -1 * sum(jumpheights * diff(k.xt))
    }
  }

  # Left correction
  if (any(is.lb)) {
    phi.lb <- phi((x[is.lb] - a) / bandwidth)
    i0.lb <- i0(xt[is.lb, ])
    psi.lb <- psi((x[is.lb] - a) / bandwidth)
    i1.lb <- i1(xt[is.lb, ])
    if (sum(is.lb) > 1) {
      x.smoothed[is.lb] <- -phi.lb * drop(jumpheights %*% diff(t(i0.lb))) - psi.lb * drop(jumpheights %*% diff(t(i1.lb)))
    } else {
      x.smoothed[is.lb] <- -phi.lb * sum(jumpheights * diff(i0.lb)) - psi.lb * sum(jumpheights * diff(i1.lb))
    }
  }

  # Right correction
  if (!is.na(b)) {
    if (any(is.ub)) {
      phi.ub <- phi((b - x[is.ub]) / bandwidth)
      i0.ub <- i0(xt[is.ub, ])
      psi.ub <- psi((b - x[is.ub]) / bandwidth)
      i1.ub <- i1(xt[is.ub, ])
      if (sum(is.ub) > 1) {
        x.smoothed[is.ub] <- -phi.ub * drop(jumpheights %*% diff(t(i0.ub))) + psi.ub * drop(jumpheights %*% diff(t(i1.ub)))
      } else {
        x.smoothed[is.ub] <- -phi.ub * sum(jumpheights * diff(i0.ub)) + psi.ub * sum(jumpheights * diff(i1.ub))
      }
    }
  }
  return(x.smoothed)
}
