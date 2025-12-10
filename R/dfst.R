#' Density function of the skewed t distribution
#'
#' This function computes the density of a skewed t distribution.
#' @param x real: values where the density has to be evaluated.
#' @param gammap positive real: skewness parameter.
#' @param mu real: location parameter.
#' @param sigma positive real: scale parameter.
#' @param nu real: number of degrees of freedom.
#' @return Value of the density function of the skewed t distribution.
#' @export
#' @examples
#' yd <- dfst(1,1.3,1,0.5,3)

dfst <- function(x,gammap,mu,sigma,nu)
{
  if (gammap < 0) 
    stop("invalid skewness")
  if (nu < 2) 
    stop("invalid DoF")
  if (nu == Inf)
    nu = 1.7e+307
  xst <- (x-mu)/(sigma*sqrt(nu/(nu-2)))
  indici1 <- xst>=0
  x1 <- xst[indici1]
  indici2 <- xst<0
  x2 <- xst[indici2]
  f1 <- rep(0,length(x1))
  f2 <- rep(0,length(x2))
  f1[indici1] <- (1/sigma) * dt(x1/gammap,nu)
  f1[!indici1] = 0
  f2[indici2] <- (1/sigma) * dt(x2*gammap,nu)
  f2[!indici2] = 0
  f <- pmax((2/(gammap+(1/gammap))) * (f1+f2),1e-323)
  return(f)
}
