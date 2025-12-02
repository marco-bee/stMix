#' Cumulative distribution function of a two-population mixture of
#' skewed t distributions
#'
#' This function evaluates the cdf of a two-population mixture of
#' skewed t distributions for a vectorial input
#' @param x real vector: values where the cdf has to be evaluated.
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
#' @return value of the cdf of a two-population mixture of
#' skewed t distributions.
#' @export
#' @examples
#' yd <- pfst(1,1.3,1,0.5,3)
#' @importFrom Rdpack reprompt
pfstMix <- function(x,p,gamma1,mu1,sigma1,nu1,gamma2,mu2,sigma2,nu2)
{
  f <- p * pfst(x,gamma1,mu1,sigma1,nu1) + (1-p) * pfst(x,gamma2,mu2,sigma2,nu2)
  return(f)
}
