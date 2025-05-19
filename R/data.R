#' Simulated dataset
#'
#' A simulated dataset using the setting with two binary covariates described
#' in Yuen, Musta and van Keilegom (2025) with
#' \eqn{\tau_G(0,0)=q_{0.975}(0,0)},
#' \eqn{\tau_G(0,1)=q_{0.975}(0,1)},
#' \eqn{\tau_G(1,0)=q_{0.975}(1,0)},
#' \eqn{\tau_G(1,1)=q_{0.95}(1,1)},
#' where \eqn{q_{1 - \epsilon}(x) = \inf\lbrace{t \geq 0}:{F_{u}(t|x) \geq 1 - \epsilon}\rbrace}.
#'
#' @format
#' \describe{
#'   \item{\eqn{Y}}{follow-up time}
#'   \item{\eqn{delta}}{censoring indicator}
#'   \item{\eqn{X}}{categorical covariate, 1, 2, 3 and 4 corresponds to
#'                  (0,0), (0,1), (1,0) and (1,1), respectively}
#' }
#'
#' @references Yuen, T. P., Musta, E., and Van Keilegom, I. (2025)
#' \emph{Testing for sufficient follow-up in survival data with categorical covariates}
"simdata"

