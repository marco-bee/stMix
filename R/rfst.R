#' Random number simulation from a skewed t distribution
#'
#' This function simulates random numbers from a skewed t distribution
#' @param n positive integer: number of simulated random numbers.
#' @param gammap positive real: skewness parameter.
#' @param mu real: location parameter.
#' @param sigma positive real: scale parameter.
#' @param nu positive real: number of degrees of freedom.
#' @return n iid random numbers from the skewed t distribution.
#' @export
#' @examples
#' ySim <- rfst(100,1.2,1.3,0.5,3)

rfst <- function(n,gammap,mu,sigma,nu)
{
  u <- runif(n)
  x <- qfst(u,gammap,mu,sigma,nu)
  return(x)
}
