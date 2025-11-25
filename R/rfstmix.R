#' Random number simulation from a two-population mixture of skewed t
#' distributions
#'
#' This function simulates random numbers from a two-population
#' mixture of skewed t distributions
#' @param n positive integer: number of simulated random numbers.
#' @param p real 0<=p<=1: mixing weight.
#' @param gamma1 positive real: skewness parameter of the first population.
#' @param mu1 real: location parameter of the first population.
#' @param sigma1 positive real: scale parameter of the first population.
#' @param nu1 positive real: number of degrees of freedomr of the first
#' population.
#' @param gamma2 positive real: skewness parameter of the second population.
#' @param mu2 real: location parameter of the second population.
#' @param sigma2 positive real: scale parameter of the second population.
#' @param nu2 positive real: number of degrees of freedom of the second
#' population.
#' @export
#' @examples
#' ySim <- rfstMix(100,1.2,1.3,0.5,3,.8,.5,1,8)
rfstMix <- function(n,p,gamma1,mu1,sigma1,nu1,gamma2,mu2,sigma2,nu2)
{
  u <- runif(n)
  x <- qfstMix(u,p,gamma1,mu1,sigma1,nu1,gamma2,mu2,sigma2,nu2)
  return(x)
}
