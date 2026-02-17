#' Cumulative distribution function of a two-population mixture of
#' skewed-normal distributions
#'
#' This function evaluates the cdf of a two-population mixture of
#' skewed-normal distributions for a vectorial input
#' @param x real vector: values where the cdf has to be evaluated.
#' @param p real 0<=p<=1: mixing weight.
#' @param gamma1 positive real: skewness parameter of the first population.
#' @param mu1 real: location parameter of the first population.
#' @param sigma1 positive real: scale parameter of the first population.
#' @param gamma2 positive real: skewness parameter of the second population.
#' @param mu2 real: location parameter of the second population.
#' @param sigma2 positive real: scale parameter of the second population.
#' @return Value of the cdf of a two-population mixture of
#' skewed-normal distributions.
#' @export
#' @examples
#' yd <- pfsnormMix(1,0.5,1.3,-1,0.5,.7,2,1)
#' @importFrom Rdpack reprompt
pfsnormMix <- function(x,p,gamma1,mu1,sigma1,gamma2,mu2,sigma2)
{
  f <- p * pfsnorm(x,gamma1,mu1,sigma1) + (1-p) * pfsnorm(x,gamma2,mu2,sigma2)
  return(f)
}
