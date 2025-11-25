#' Quantile function of the skewed Student t distribution
#'
#' This function evaluates the Quantile function of the skewed
#' Student t distribution
#' @param p real, 0<=x<=1: quantile levels.
#' @param gammap positive real: skewness parameter.
#' @param mu real: location parameter.
#' @param sigma positive real: scale parameter.
#' @param nu positive real: number of degrees of freedom.
#' @return value of the quantile function of the skewed t distribution.
#' @export
#' @examples
#' yd <- qfst(.1,1.2,1,0.5,3)
#' @importFrom Rdpack reprompt

qfst <- function(p,gammap,mu,sigma,nu)
{
  x <- rep(0,length(p))
  alphaStar <- pfst(mu,gammap,mu,sigma,nu)
  indici1 <- p < alphaStar
  p1 <- p[indici1]
  indici2 <- p >= alphaStar
  p2 <- p[indici2]
  x[indici1] <- sigma*qt(p1*(gammap^2+1)/2,nu)/gammap + mu
  x[indici2] <- mu + sigma*gammap * 
    qt((1/(2*gammap^2)) * (p2*(gammap^2+1)-1)+1/2,nu)
  return(x)
}
