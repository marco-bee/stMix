#' Weighted log-likelihood of a skewed t distribution
#'
#' This function evaluates the log-likelihood of a skewed t distribution
#'  with weighted observations.
#' @param x numerical vector (4x1): values of the parameters \eqn{\gamma},
#' \eqn{\mu}, \eqn{\sigma} and \eqn{\nu}. 
#' @param y numerical vector (nx1): observed data.
#' @param post numerical vector (nx1) with elements in (0,1): weights
#' of the observations (in the EM algorithm, posterior probabilities).
#' @return Numerical value of the log-likelihood function
#' @export
#' @examples
#' y <- rfst(100,1,1.3,0.5,3)
#' x <- c(1.3,0.5,3)
#' post <- runif(100)
#' llik <- llstWei0(x,y,post)
#'
#' @importFrom Rdpack reprompt
llstWei0 <- function(x,y,post)
{
  mu <- x[1]
  sigma <- exp(x[2])
  nu <- exp(x[3])
  ll <- dfst(y,1,mu,sigma,nu)
  llik <- sum(log(ll)*post)
  return(llik)
}
