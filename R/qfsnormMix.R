#' Quantile function of a mixture of skewed normal t distributions
#'
#' This function evaluates the Quantile function of a mixture of skewed
#' normal distributions
#' @param w real, 0<=w<=1: quantile levels.
#' @param p real 0<=p<=1: mixing weight.
#' @param gamma1 positive real: skewness parameter of the first population.
#' @param mu1 real: location parameter of the first population.
#' @param sigma1 positive real: scale parameter of the first population.
#' @param gamma2 positive real: skewness parameter of the second population.
#' @param mu2 real: location parameter of the second population.
#' @param sigma2 positive real: scale parameter of the second population.
#' @return value of the quantile function of a mixture of skewed-normal
#' distributions. 
#' @export
#' @examples
#' yd <- qfsnormMix(.9,.1,.8,-1,0.5,1.3,2,1)
#' @importFrom Rdpack reprompt

qfsnormMix <- function(w,p,gamma1,mu1,sigma1,gamma2,mu2,sigma2)
{
  nw <- length(w)
  x <- rep(NA,nw)
  minmax <- range(qfsnorm(w,gamma1,mu1,sigma1),qfsnorm(w,gamma2,mu2,sigma2)) 
  for (i in 1:nw)
  {
    x[i] <- uniroot(function(x) pfsnormMix(x,p,gamma1,mu1,sigma1,gamma2,mu2,sigma2)-w[i],
                    interval = minmax,
                    tol = 10^{-16})$root  
  }
  return(x)
}
