#' Cumulative distribution function of the skewed Student t distribution
#'
#' This function evaluates the cdf of the skewed Student t distribution
#' for a vectorial input
#' @param x real vector: values where the cdf has to be evaluated.
#' @param gammap positive real: skewness parameter.
#' @param mu real: location parameter.
#' @param sigma positive real: scale parameter.
#' @param nu positive real: number of degrees of freedom.
#' @return Value of the cdf of the skewed t distribution.
#' @export
#' @examples
#' yd <- pfst(1,1.3,1,0.5,3)
#' @importFrom Rdpack reprompt

pfst <- function(x,gammap,mu,sigma,nu)
{
  F <- rep(0,length(x))
  xst <- (x-mu)/sigma
  indici1 <- xst>=0
  x1 <- xst[indici1]
  indici2 <- xst<0
  x2 <- xst[indici2]
  F[indici2] <- (2/(1+gammap^2)) * pt(x2*gammap,nu)
  F[indici1] <- (1/(1+gammap^2)) +
    (2*gammap/(gammap+1/gammap)) * (pt(x1/gammap,nu) - 1/2)
  return(F)
}

