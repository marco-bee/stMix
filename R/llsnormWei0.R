#' Weighted log-likelihood of a skewed normal distribution
#'
#' This function evaluates the log-likelihood of a skewed normal distribution
#'  with weighted observations.
#' @param x numerical vector (3x1): values of the parameters \eqn{\gamma},
#' \eqn{\mu} and \eqn{\sigma}. 
#' @param y numerical vector (nx1): observed data.
#' @param post numerical vector (nx1) with elements in (0,1): weights
#' of the observations (in the EM algorithm, posterior probabilities).
#' @return Numerical value of the log-likelihood
#' function
#' @export
#' @examples
#' y <- rfst(100,1.2,1.3,0.5,300)
#' x <- c(1.2,1.3,0.5,3)
#' post <- runif(100)
#' llik <- llsnormWei(x,y,post)
#'
#' @importFrom Rdpack reprompt
llsnormWei0 <- function(x,y,post)
{
  mu <- x[1]
  sigma <- exp(x[2])
  ll <- dsnorm(y,1,mu,sigma)
  llik <- sum(log(ll)*post)
  return(llik)
}
