#' Density function of the a mixture of two skewed-normal distributions
#'
#' This function computes the density of a mixture of two skewed-normal
#' distributions.
#' @param x real: values where the density has to be evaluated.
#' @param gamma1 positive real: skewness parameter of the first population.
#' @param mu1 real: location parameter of the first population.
#' @param sigma1 positive real: scale parameter of the first population.
#' @param gamma2 positive real: skewness parameter of the second population.
#' @param mu2 real: location parameter of the second population.
#' @param sigma2 positive real: scale parameter of the second population.
#' @return Value of the density function of a mixture of two skewed normal
#' distributions.
#' @export
#' @examples
#' yd <- dfsnormMix(1,1.3,1,0.5,.8,3,.7)
#' 
dfsnormMix <- function(x,p,gamma1,mu1,sigma1,gamma2,mu2,sigma2)
{
  if (gamma1!=1 & gamma2!=1)
    f <- p * dsnorm(x,gamma1,mu1,sigma1) + (1-p) * dsnorm(x,gamma2,mu2,sigma2)
  if (gamma1==1 & gamma2!=1)
    f <- p * dnorm(x,mu1,sigma1) + (1-p) * dsnorm(x,gamma2,mu2,sigma2)
  if (gamma1!=1 & gamma2==1)
    f <- p * dsnorm(x,gamma1,mu1,sigma1) + (1-p) * dnorm(x,mu2,sigma2)
  if (gamma1==1 & gamma2==1)
    f <- p * dnorm(x,mu1,sigma1) + (1-p) * dnorm(x,mu2,sigma2)
  return(f)
}
