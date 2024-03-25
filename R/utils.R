#' Kaplan-Meier estimate of the survivor function
#'
#' @param y A numeric vector (length \code{n}) of the observed event time.
#' @param delta  A numeric vector (length \code{n}) of the censoring indicator.
#' @param adjust A logical to indicate if imposing an adjustment in computing
#'               the hazard components to avoid the KME being zero when
#'               time greater than the largest observation time.
#'
#' @return A \code{\link[stats]{stepfun}} of the estimate.
#' @importFrom stats stepfun
#' @keywords internal
estimate.km <- function(y, delta, adjust = FALSE) {
  # Distinct observed event times in ascending order
  y.star <- sort(unique(y))
  # Risk set at t_{(j)} at the j-th row (n x length(y.star))
  risk.set <- sapply(y.star, function(x) {
    c(y >= x, # Risk set
      sum(y == x) - sum(delta[y == x] == 0) # No. of events excluding censored
    )
  })
  # No. of events at t_{(j)} at the j-th entry
  d.star <- drop(risk.set[nrow(risk.set), ])
  risk.set <- risk.set[-nrow(risk.set), ]

  product.component <- if (!adjust) {
    # Hazard components
    hazard <- d.star / drop(colSums(risk.set))
    cumprod(c(1, 1 - hazard))
  } else {
    risk.sum <- drop(colSums(risk.set))
    cumprod(c(1, (risk.sum - d.star + 1) / (risk.sum + 1)))
  }

  km.stepfun <- stepfun(y.star, product.component, f = 0, right = FALSE)
  return(km.stepfun)
}


#' Extract probability masses from Kaplan-Meier estimate
#'
#' @param km A \code{\link[stats]{stepfun}} of the Kaplan-Meier estimate
#' computed from \code{\link{estimate.km}}.
#'
#' @return A list of the extract masses and the corresponding time
#' @keywords internal
extract.km <- function(km) {
  km.x <- get("x", environment(km))
  km.y <- get("y", environment(km))
  km.y.diff <- diff(c(1, km.y))
  is.jump <- which(km.y.diff != 0)
  x.selected <- km.x[is.jump]
  p.selected <- -1 * km.y.diff[is.jump]
  if (any(p.selected <= 0 | p.selected >= 1)) {
    stop("Error from retrieving probability from the Kaplan-Meier estiamte")
  }
  return(list(x = x.selected, p = p.selected))
}
