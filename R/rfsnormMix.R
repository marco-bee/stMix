#' Random number simulation from a two-population mixture of skewed-normal
#' distributions
#'
#' This function simulates random numbers from a two-population
#' mixture of skewed-normal distributions
#' @param n positive integer: number of simulated random numbers.
#' @param p real 0<=p<=1: mixing weight.
#' @param gamma1 positive real: skewness parameter of the first population.
#' @param mu1 real: location parameter of the first population.
#' @param sigma1 positive real: scale parameter of the first population.
#' @param gamma2 positive real: skewness parameter of the second population.
#' @param mu2 real: location parameter of the second population.
#' @param sigma2 positive real: scale parameter of the second population.
#' @return A vector containing n iid random numbers from a mixture of skewed-
#' normal distributions.
#' @export
#' @examples
#' ySim <- rfsnormMix(100,1.2,1.3,0.5,.8,.5,1)
rfsnormMix <- function(n,p,gamma1,mu1,sigma1,gamma2,mu2,sigma2)
{
  u <- runif(n)
  x <- qfsnormMix(u,p,gamma1,mu1,sigma1,gamma2,mu2,sigma2)
  return(x)
}
