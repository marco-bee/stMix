#' Quantile function of the skewed Student t distribution
#'
#' This function evaluates the quantile function of the skewed
#' Student t distribution
#' @param w real, 0<=w<=1: quantile levels.
#' @param gammap positive real: skewness parameter.
#' @param mu real: location parameter.
#' @param sigma positive real: scale parameter.
#' @param nu positive real: number of degrees of freedom.
#' @return Value of the quantile function of the skewed t distribution.
#' @export
#' @examples
#' yd <- qfst(.1,1.2,1,0.5,3)
#' @importFrom Rdpack reprompt

qfst <- function(w,gammap,mu,sigma,nu)
{
  x <- rep(0,length(w))
  alphaStar <- pfst(mu,gammap,mu,sigma,nu)
  indici1 <- w < alphaStar
  w1 <- w[indici1]
  indici2 <- w >= alphaStar
  w2 <- w[indici2]
  x[indici1] <- sigma*qt(w1*(gammap^2+1)/2,nu)/gammap + mu
  x[indici2] <- mu + sigma*gammap * 
    qt((1/(2*gammap^2)) * (w2*(gammap^2+1)-1)+1/2,nu)
  return(x)
}
