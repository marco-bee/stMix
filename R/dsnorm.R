#' Density function of the skewed normal distribution
#'
#' This function computes the density of a skewed normal distribution.
#' @param x real: values where the density has to be evaluated.
#' @param gammap real: skewness parameter.
#' @param mu real: location parameter.
#' @param sigma positive real: scale parameter.
#' @return value of the density function of the skewed normal distribution.
#' @export
#' @examples
#' yd <- dsnorm(1,1.3,1,0.5)

dsnorm <- function(x,gammap,mu,sigma)
{
  if (gammap < 0) 
    stop("invalid skewness")
  xst <- (x-mu)/sigma
  indici1 <- xst>=0
  x1 <- xst[indici1]
  indici2 <- xst<0
  x2 <- xst[indici2]
  f1 <- rep(0,length(x1))
  f2 <- rep(0,length(x2))
  f1[indici1] <- dnorm(x1/gammap)/sigma
  f1[!indici1] = 0
  f2[indici2] <- dnorm(x2*gammap)/sigma
  f2[!indici2] = 0
  f <- (2/(gammap+(1/gammap))) * (f1+f2)
  return(f)
}
