#' Evaluate the Grenander density estimate at \code{x}
#'
#'
#' @param x A numeric vector for which the value of Grenander density estimate is computed.
#' @param gren.times A numeric vector of points at which there is a jump.
#' @param gren.val A numeric vector of Grenander density estimate corresponds to `gren.times`.
#'
#' @return A numeric vector of the computed values
#' @keywords internal
f.hat <- function(x, gren.times, gren.val) {
  f.hat.val <- x
  int.idx <- findInterval(x, gren.times, left.open = FALSE)
  has.val <- (int.idx != 0)
  f.hat.val[has.val] <- gren.val[int.idx[has.val]]
  f.hat.val[!has.val] <- gren.val[1]
  at.rendpt <- (int.idx == length(gren.times))
  f.hat.val[at.rendpt] <- gren.val[length(gren.val)]
  return(f.hat.val)
}


# #' Obtain the Least Concave Majorant (LCM)
# #'
# #' Obtain the Least Concave Majorant from \code{\link[fdrtool]{gcmlcm}}.
# #'
# #' @param gcmlcm.fit A list of LCM of KME returned from \code{\link[fdrtool]{gcmlcm}}.
# #'
# #' @return A function from \code{\link[stats]{approxfun}} of the LCM
# #' @importFrom stats approxfun
# #' @keywords internal
# lcm.F.hat <- function(gcmlcm.fit) {
#   return(approxfun(gcmlcm.fit$x.knots, gcmlcm.fit$y.knots,
#                    yleft = 0, yright = max(gcmlcm.fit$y.knots)))
# }


#' Kenrnel smoothing of the LCM of the KME (near left endpoint)
#'
#'
#' @param gcmlcm.fit A list of LCM of KME returned from \code{\link[fdrtool]{gcmlcm}}.
#' @param tau.c A numeric of the right endpoint of the support.
#' @param x.grid.left A numeric vector of the points for which the value of Smoothed LCM of KME is computed.
#' @param bw A numeric of bandwidth
#'
#' @return A numeric vector of the computed value at \code{x.grid.left}.
#' @keywords internal
smooth.kme.left <- function(gcmlcm.fit, tau.c, x.grid.left, bw) {
  ### Left side x \in [0, bw]
  ## Step 1
  gcmlcm.slope.knots <- gcmlcm.fit$slope.knots
  gcmlcm.t.knots <- gcmlcm.fit$x.knots
  gcmlcm.intercept.knots <- gcmlcm.fit$y.knots[-1] - gcmlcm.slope.knots * gcmlcm.t.knots[-1]
  gcmlcm.slope.knots <- c(gcmlcm.slope.knots, gcmlcm.slope.knots[length(gcmlcm.slope.knots)])
  gcmlcm.intercept.knots <- c(gcmlcm.intercept.knots, gcmlcm.intercept.knots[length(gcmlcm.intercept.knots)])

  ## Step 2
  phi.x.grid.left <- phi(x.grid.left / bw)
  psi.x.grid.left <- psi(x.grid.left / bw)
  x.grid.left.size <- length(x.grid.left)

  ## Step 3
  t.knots <- c(gcmlcm.t.knots, sqrt(.Machine$double.xmax))
  t.minus.x.grid <- -1.0 * outer(t.knots, x.grid.left, "-") / bw
  knots.size <- length(t.knots)

  ## Step 4
  x.b <- matrix(rep(x.grid.left / bw, knots.size),
                knots.size, x.grid.left.size, byrow = TRUE)
  integral.left.0.u <- i0(pmin(x.b, t.minus.x.grid))
  integral.left.1.u <- i1(pmin(x.b, t.minus.x.grid))
  integral.left.2.u <- i2(pmin(x.b, t.minus.x.grid))

  neg.one <- matrix(-1, nrow(t.minus.x.grid), ncol(t.minus.x.grid))
  integral.left.0.l <- i0(pmax(neg.one, t.minus.x.grid))
  integral.left.1.l <- i1(pmax(neg.one, t.minus.x.grid))
  integral.left.2.l <- i2(pmax(neg.one, t.minus.x.grid))

  ## Step 5
  integral.left.A <- integral.left.0.u[-knots.size, ] -
    integral.left.0.l[-1, ]
  integral.left.B <- integral.left.1.u[-knots.size, ] -
    integral.left.1.l[-1, ]
  integral.left.C <- integral.left.2.u[-knots.size, ] -
    integral.left.2.l[-1, ]

  ## Step 6
  ax.plus.b <- outer(gcmlcm.slope.knots, x.grid.left) +
    matrix(rep(gcmlcm.intercept.knots, x.grid.left.size),
           nrow = knots.size - 1L, ncol = x.grid.left.size, byrow = FALSE)
  a.bw <- matrix(rep(bw * gcmlcm.slope.knots, x.grid.left.size),
                 nrow = knots.size - 1L, ncol = x.grid.left.size, byrow = FALSE)
  integral.left.I <- phi.x.grid.left * drop(colSums(ax.plus.b * integral.left.A - a.bw * integral.left.B))
  integral.left.II <- psi.x.grid.left * drop(colSums(ax.plus.b * integral.left.B - a.bw * integral.left.C))
  integral.left <- integral.left.I + integral.left.II
  return(integral.left)
}


#' Kenrnel smoothing of the LCM of the KME (no boundary correction)
#'
#'
#' @param gcmlcm.fit A list of LCM of KME returned from \code{\link[fdrtool]{gcmlcm}}.
#' @param tau.c A numeric of the right endpoint of the support.
#' @param x.grid.middle A numeric vector of the points for which the value of Smoothed LCM of KME is computed.
#' @param bw A numeric of bandwidth
#'
#' @return A numeric vector of the computed value at \code{x.grid.middle}.
#' @keywords internal
smooth.kme.middle <- function(gcmlcm.fit, tau.c, x.grid.middle, bw) {
  ### Middle x \in (bw, tau.c - bw)
  ## Step 1
  gcmlcm.slope.knots <- gcmlcm.fit$slope.knots
  gcmlcm.t.knots <- gcmlcm.fit$x.knots
  gcmlcm.intercept.knots <- gcmlcm.fit$y.knots[-1] - gcmlcm.slope.knots * gcmlcm.t.knots[-1]
  gcmlcm.slope.knots <- c(gcmlcm.slope.knots, gcmlcm.slope.knots[length(gcmlcm.slope.knots)])
  gcmlcm.intercept.knots <- c(gcmlcm.intercept.knots, gcmlcm.intercept.knots[length(gcmlcm.intercept.knots)])

  #### Step 2
  x.grid.middle.size <- length(x.grid.middle)
  t.knots <- c(gcmlcm.t.knots, sqrt(.Machine$double.xmax))
  t.minus.x.grid <- -1.0 * outer(t.knots, x.grid.middle, "-") / bw
  knots.size <- length(t.knots)

  #### Step 3
  one <- matrix(1, nrow(t.minus.x.grid), ncol(t.minus.x.grid))
  integral.middle.0.u <- i0(pmin(one, t.minus.x.grid))
  integral.middle.1.u <- i1(pmin(one, t.minus.x.grid))

  neg.one <- matrix(-1, nrow(t.minus.x.grid), ncol(t.minus.x.grid))
  integral.middle.0.l <- i0(pmax(neg.one, t.minus.x.grid))
  integral.middle.1.l <- i1(pmax(neg.one, t.minus.x.grid))

  #### Step 4
  integral.middle.A <- integral.middle.0.u[-knots.size, ] -
    integral.middle.0.l[-1, ]
  integral.middle.B <- integral.middle.1.u[-knots.size, ] -
    integral.middle.1.l[-1, ]

  #### Step 5
  ax.plus.b <- outer(gcmlcm.slope.knots, x.grid.middle) +
    matrix(rep(gcmlcm.intercept.knots, x.grid.middle.size),
           nrow = knots.size - 1L, ncol = x.grid.middle.size, byrow = FALSE)
  a.bw <- matrix(rep(bw * gcmlcm.slope.knots, x.grid.middle.size),
                 nrow = knots.size - 1L, ncol = x.grid.middle.size, byrow = FALSE)
  integral.middle <- drop(colSums(ax.plus.b * integral.middle.A - a.bw * integral.middle.B))
  return(integral.middle)
}


#' Kenrnel smoothing of the LCM of the KME (near right endpoint)
#'
#'
#' @param gcmlcm.fit A list of LCM of KME returned from \code{\link[fdrtool]{gcmlcm}}.
#' @param tau.c A numeric of the right endpoint of the support.
#' @param x.grid.right A numeric vector of the points for which the value of Smoothed LCM of KME is computed.
#' @param bw A numeric of bandwidth
#'
#' @return A numeric vector of the computed value at \code{x.grid.right}.
#' @keywords internal
smooth.kme.right <- function(gcmlcm.fit, tau.c, x.grid.right, bw) {
  ### Right x \in [tau.c - bw, tau.c]
  ## Step 1
  gcmlcm.slope.knots <- gcmlcm.fit$slope.knots
  gcmlcm.t.knots <- gcmlcm.fit$x.knots
  gcmlcm.intercept.knots <- gcmlcm.fit$y.knots[-1] - gcmlcm.slope.knots * gcmlcm.t.knots[-1]
  gcmlcm.slope.knots <- c(gcmlcm.slope.knots, gcmlcm.slope.knots[length(gcmlcm.slope.knots)])
  gcmlcm.intercept.knots <- c(gcmlcm.intercept.knots, gcmlcm.intercept.knots[length(gcmlcm.intercept.knots)])

  #### Step 2
  x.grid.right.size <- length(x.grid.right)
  t.knots <- c(gcmlcm.t.knots, sqrt(.Machine$double.xmax))
  t.minus.x.grid <- -1.0 * outer(t.knots, x.grid.right, "-") / bw
  knots.size <- length(t.knots)

  #### Step 2
  phi.x.grid.right <- phi((tau.c - x.grid.right) / bw)
  psi.x.grid.right <- psi((tau.c - x.grid.right) / bw)

  #### Step 4
  one <- matrix(1, nrow(t.minus.x.grid), ncol(t.minus.x.grid))
  integral.right.0.u <- i0(pmin(one, t.minus.x.grid))
  integral.right.1.u <- i1(pmin(one, t.minus.x.grid))
  integral.right.2.u <- i2(pmin(one, t.minus.x.grid))

  x.tau.b <- matrix(rep((x.grid.right - tau.c) / bw, knots.size),
                    knots.size, x.grid.right.size, byrow = TRUE)
  integral.right.0.l <- i0(pmax(x.tau.b, t.minus.x.grid))
  integral.right.1.l <- i1(pmax(x.tau.b, t.minus.x.grid))
  integral.right.2.l <- i2(pmax(x.tau.b, t.minus.x.grid))

  #### Step 5
  integral.right.A <- integral.right.0.u[-knots.size, ] -
    integral.right.0.l[-1, ]
  integral.right.B <- integral.right.1.u[-knots.size, ] -
    integral.right.1.l[-1, ]
  integral.right.C <- integral.right.2.u[-knots.size, ] -
    integral.right.2.l[-1, ]

  #### Step 6
  ax.plus.b <- outer(gcmlcm.slope.knots, x.grid.right) +
    matrix(rep(gcmlcm.intercept.knots, x.grid.right.size),
           nrow = knots.size - 1L, ncol = x.grid.right.size, byrow = FALSE)
  a.bw <- matrix(rep(bw * gcmlcm.slope.knots, x.grid.right.size),
                 nrow = knots.size - 1L, ncol = x.grid.right.size, byrow = FALSE)
  integral.right.I <- phi.x.grid.right * drop(colSums(ax.plus.b * integral.right.A - a.bw * integral.right.B))
  integral.right.II <- psi.x.grid.right * drop(colSums(ax.plus.b * integral.right.B - a.bw * integral.right.C))
  integral.right <- integral.right.I - integral.right.II
  return(integral.right)
}


#' Kenrnel smoothing of the LCM of the KME
#'
#' Kenrnel smoothing of the LCM of the KME:
#' \deqn{
#'  \int\frac{1}{h_{0}}k_{B,t}({\frac{t-u}{h_{0}}})\hat{F}_{n}^{G}(u)du,
#' }
#' where \eqn{k_{B,t}} is the boundary kernel, \eqn{h_{0}} is the bandwidth,
#' \eqn{\hat{F}_{n}^{G}} is the LCM of KME.
#'
#'
#'
#' @param gcmlcm.fit A list of LCM of KME returned from \code{\link[fdrtool]{gcmlcm}}.
#' @param tau.c A numeric of the end of the study time.
#' @param n A numeric of the sample size to determine the bandwidth.
#' @param skm.bw.multi A numeric of scaling to the bandwidth for smoothing
#' the LCM of the KME.
#' @param grid.size An integer of the grid size to generate a grid of points
#' for which the value of Smoothed LCM of KME is computed.
#'
#' @return A \code{data.frame}.
#' @keywords internal
smooth.kme <- function(gcmlcm.fit, tau.c, n,
                       skm.bw.multi = 0.7, grid.size = 5000L) {
  bw.skm <- min(skm.bw.multi*tau.c*(n^(-1/9)),tau.c/2)

  x.grid <- seq(0, tau.c, length = grid.size)
  smooth.km <- rep(NA, grid.size)
  x.left <- (x.grid <= bw.skm)
  x.middle <- (x.grid > bw.skm & x.grid < tau.c - bw.skm)
  x.right <- (x.grid >= tau.c - bw.skm & x.grid <= tau.c)
  smooth.km[x.left] <- smooth.kme.left(gcmlcm.fit, tau.c, x.grid[x.left], bw.skm)
  smooth.km[x.middle] <- smooth.kme.middle(gcmlcm.fit, tau.c, x.grid[x.middle], bw.skm)
  smooth.km[x.right] <- smooth.kme.right(gcmlcm.fit, tau.c, x.grid[x.right], bw.skm)
  return(data.frame("x" = x.grid, "skm" = smooth.km))
}
