#' Evaluate the log-density of a stationary Gaussian AR(1) process.
#'
#' Evaluate the log-density of a stationary Gaussian AR(1) process, observed at
#' times \code{times} taking values \code{x}.
#' @param x A vector of observed values.
#' @param mu A vector of expected values.
#' @param times A vector of the time points of observation.
#' @param rho A real number strictly less than 1 in absolute value.
#' @param sigma A positive real number.
#' @return A scalar, the log density.
#' @export
#' @examples
#' x <- rnorm(5) + 1:5
#' t <- c(1, 3, 5:6, 10)
#' rho <- 0.5
#' sigma <- 1
#' # zero mean
#' ar1_lpdf(x, t, rho, sigma)
#' # means equal times
#' mu <- t
#' ar1_lpdf(x + mu, t, rho, sigma, mu)
ar1_lpdf <- function(x, times, rho, sigma, mu = 0) {
  if (!(length(mu) %in% c(1, length(x)))) {
    stop("mu must be of length 1 or length(x).")
  }
  if (length(x) != length(times)) {
    stop("x must be the same length as the argument times.")
  }
  if (length(sigma) != 1 || length(rho) != 1) {
    stop("sigma and rho must be scalars.")
  }
  if (abs(rho) >= 1) stop("rho must be less than 1 in magnitude.")
  if (sigma < 0) stop("sigma must be positive")
  if (length(mu) == 1) {
    mu <- rep(mu, length(x))
  }

  return(ar1_lpdf_cpp(x, mu, times, rho, sigma))
}