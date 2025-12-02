#' Quantile function of a mixture of skewed Student t distributions
#'
#' This function evaluates the Quantile function of a mixture of skewed
#' t distributions
#' @param w real, 0<=w<=1: quantile levels.
#' @param p real 0<=p<=1: mixing weight.
#' @param gamma1 positive real: skewness parameter of the first population.
#' @param mu1 real: location parameter of the first population.
#' @param sigma1 positive real: scale parameter of the first population.
#' @param nu1 positive real: number of degrees of freedom of the first
#' population.
#' @param gamma2 positive real: skewness parameter of the second population.
#' @param mu2 real: location parameter of the second population.
#' @param sigma2 positive real: scale parameter of the second population.
#' @param nu2 positive real: number of degrees of freedom of the second
#' population.
#' @return value of the quantile function of a mixture of skewed t
#' distributions. Notice that, if \eqn{\nu_i=+\infty} (\eqn{i=1,2}), the
#' corresponding component becomes a skewed normal density.
#' @export
#' @examples
#' yd <- qfst(.1,1.2,1,0.5,3)
#' @importFrom Rdpack reprompt

qfstMix <- function(w,p,gamma1,mu1,sigma1,nu1,gamma2,mu2,sigma2,nu2)
{
  nw <- length(w)
  x <- rep(NA,nw)
  minmax <- range(qfst(w,gamma1,mu1,sigma1,nu1),qfst(w,gamma2,mu2,sigma2,nu2)) 
  for (i in 1:nw)
  {
    x[i] <- uniroot(function(x) pfstMix(x,p,gamma1,mu1,sigma1,nu1,gamma2,mu2,sigma2,nu2)-w[i],
                    interval = minmax,
                    tol = 10^{-16})$root  
  }
  return(x)
}
