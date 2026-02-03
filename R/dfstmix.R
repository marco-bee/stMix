#' Density function of the a mixture of two skewed t distributions
#'
#' This function computes the density of a mixture of two skewed t
#' distributions.
#' @param x real: values where the density has to be evaluated.
#' @param gamma1 positive real: skewness parameter of the first population.
#' @param mu1 real: location parameter of the first population.
#' @param sigma1 positive real: scale parameter of the first population.
#' @param nu1 real: number of degrees of freedom of the first population.
#' @param gamma2 positive real: skewness parameter of the second population.
#' @param mu2 real: location parameter of the second population.
#' @param sigma2 positive real: scale parameter of the second population.
#' @param nu2 real: number of degrees of freedom of the second population.
#' @return Value of the density function of the skewed t distribution. Notice
#' that, if \eqn{\nu_i=+\infty} (\eqn{i=1,2}), the corresponding component
#' becomes a skewed normal density.
#' @export
#' @examples
#' yd <- dfstMix(1,1.3,1,0.5,3,.8,3,.7,5)
#' 
dfstMix <- function(x,p,gamma1,mu1,sigma1,nu1,gamma2,mu2,sigma2,nu2)
{
  if (gamma1!=1 & gamma2!=1)
    f <- p * dfst(x,gamma1,mu1,sigma1,nu1) + (1-p) * dfst(x,gamma2,mu2,sigma2,nu2)
  if (gamma1==1 & gamma2!=1)
    f <- p * extraDistr::dlst(x,nu1,mu1,sigma1) + (1-p) * dfst(x,gamma2,mu2,sigma2,nu2)
  if (gamma1!=1 & gamma2==1)
    f <- p * dfst(x,gamma1,mu1,sigma1,nu1) + (1-p) * extraDistr::dlst(x,nu2,mu2,sigma2)
  if (gamma1==1 & gamma2==1)
    f <- p * extraDistr::dlst(x,nu1,mu1,sigma1) + (1-p) * extraDistr::dlst(x,nu2,mu2,sigma2)
  return(f)
}
