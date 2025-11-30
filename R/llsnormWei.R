#' Weighted log-likelihood of a skewed normal distribution
#'
#' This function evaluates the log-likelihood of a skewed normal distribution
#'  with weighted observations.
#' @param x numerical vector (3x1): values of the parameters \eqn{\gamma},
#' \eqn{\mu} and \eqn{sigma}. 
#' @param y numerical vector (nx1): observed data.
#' @param post numerical vector (nx1) with elements in (0,1): weights
#' of the observations (in the EM algorithm, posterior probabilities).
#' @return llik real: numerical value of the log-likelihood
#' function
#' @export
#' @examples
#' y <- rfst(100,1.2,1.3,0.5,300)
#' x <- c(1.2,1.3,0.5,3)
#' post <- runif(100)
#' llik <- llsnormWei(x,y,post)
#'
#' @importFrom Rdpack reprompt
llsnormWei <- function(x,y,post)
{
  gammap <- exp(x[1])
  mu <- x[2]
  sigma <- exp(x[3])
  ll <- dsnorm(y,gammap,mu,sigma)
  llik <- sum(log(ll)*post)
  return(llik)
}
