#' Log-likelihood of a skewed t distribution
#'
#' This function evaluates the log-likelihood of a skewed t distribution.
#' @param x numerical vector (4x1): values of the parameters \eqn{\gamma},
#' \eqn{\mu}, \eqn{\sigma} and \eqn{\nu}. 
#' @param y numerical vector (nx1): observed data.
#' @return Numerical value of the log-likelihood
#' function
#' @export
#' @examples
#' y <- rfst(100,1.2,1.3,0.5,3)
#' x <- c(1.2,1.3,0.5,3)
#' llik <- llfst(x,y)
#'
#' @importFrom Rdpack reprompt

llfst <- function(x,y)
{
  gammap <- exp(x[1])
  mu <- x[2]
  sigma <- exp(x[3])
  tau <- exp(x[4])
  ll <- sum(log(dfst(y,gammap,mu,sigma,tau)))
  return(ll)
}
