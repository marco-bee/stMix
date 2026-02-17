#' Quantile function of the skewed-normal distribution
#'
#' This function evaluates the quantile function of the skewed
#' Student t distribution
#' @param w real, 0<=w<=1: quantile levels.
#' @param gammap positive real: skewness parameter.
#' @param mu real: location parameter.
#' @param sigma positive real: scale parameter.
#' @return Value of the quantile function of the skewed t distribution.
#' @export
#' @examples
#' yd <- qfsnorm(.1,1.2,1,0.5)
#' @importFrom Rdpack reprompt

qfsnorm <- function(w,gammap,mu,sigma)
{
  x <- rep(0,length(w))
  alphaStar <- pfsnorm(mu,gammap,mu,sigma)
  indici1 <- w < alphaStar
  w1 <- w[indici1]
  indici2 <- w >= alphaStar
  w2 <- w[indici2]
  x[indici1] <- sigma*qnorm(w1*(gammap^2+1)/2)/gammap + mu
  x[indici2] <- mu + sigma*gammap * 
    qnorm((1/(2*gammap^2)) * (w2*(gammap^2+1)-1)+1/2)
  return(x)
}
