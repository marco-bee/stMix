#' Mixture estimation via EM
#'
#' This function estimates a of a two-population mixture of
#' skewed t distributions by means of the EM algorithm.
#' @param x0 numerical vector (9x1): initial values of the parameters \emph{p},
#' \eqn{\mu_1,} \eqn{\sigma_1}, \eqn{\nu_1}, 
#' \eqn{\gamma_2}, \eqn{\mu_2}, \eqn{\sigma_2}, \eqn{\nu_2}. 
#' @param y vector: observed data.
#' @param maxiter positive integer: maximum number of iterations of the EM algorithm.
#' @return A list with the following elements is returned:
#' "p" = estimated value of p,
#' "post" = posterior probabilities of all observations,
#' "mu1" = estimated value of \eqn{\mu_1},
#' "sigma1 " = estimated value of \eqn{\sigma_1},
#' "nu1" = estimated value of \eqn{\nu_1},
#' "gamma2" = estimated value of \eqn{\gamma_2},
#' "mu2" = estimated value of \eqn{\mu_2},
#' "sigma2 " = estimated value of \eqn{\sigma_2},
#' "nu2" = estimated value of \eqn{\nu_2},
#' "nit" = number of iterations.
#' @export
#' @examples
#' y <- rfstMix(100,.5,1,1.3,0.5,3,.8,1,1,7)
#' x <- c(.6,1.3,1,0.8,3,1.2,.7,.9,9)
#' res <- stMixEM(y,.6,0.8,1,3,1.7,.9,9,1000)
#'
#' @importFrom Rdpack reprompt

stMixEM01 <- function(Y,p1,mu1,sigma1,nu1,gamma2,mu2,sigma2,nu2,maxiter)
{
  p2 <- 1 - p1
  n <- length(Y)
  nit <- 1
  epsilon <- 1.e-10
  change <- 500

  # start iterations

  while (2>1)
  {
    parold <- c(p1,mu1,sigma1,gamma2,mu2,sigma2)

    # evaluate component densities and mixture density

    f1 <- dfst(Y,1,mu1,sigma1,nu1)
    f2 <- dfst(Y,gamma2,mu2,sigma2,nu2)
    f <- p1 * f1 + p2 * f2
    loglik <- sum(log(f))

    # update posterior probabilities: E step

    post1 <- p1*f1 / f
    post2 <- p2*f2 / f
    post <- cbind(post1, post2)

    # update prior probabilities: M step

    prior <- colMeans(post)
    p1 <- as.double(prior[1])
    p2 <- as.double(prior[2])

    # update means: M-step for means

    res1 <- optim(c(mu1,log(sigma1),pmin(log(nu1),1.7e+307)),
                  llstWei0, gr = NULL, control = list(fnscale = -1), Y, post1)
    mu1 <- res1$par[1]
    sigma1 <- exp(res1$par[2])
    nu1 <- exp(res1$par[3])
    res2 <- optim(c(log(gamma2),mu2,log(sigma2),pmin(log(nu2),1.7e+307)),
          llstWei, gr = NULL, control = list(fnscale = -1), Y, post2)
    gamma2 <- exp(res2$par[1])
    mu2 <- res2$par[2]
    sigma2 <- exp(res2$par[3])
    nu2 <- exp(res2$par[4])
    loglik <- sum(log(dfstMix(Y,p1,1,mu1,sigma1,nu1,gamma2,mu2,sigma2,nu2)))
    
    # check convergence

    diffpar <- parold - c(p1,mu1,sigma1,gamma2,mu2,sigma2)
    change <- max(abs(diffpar))
    if (nit > maxiter || change < epsilon)
      break
    nit <- nit + 1
    estPars <- c(p1,mu1,sigma1,nu1,gamma2,mu2,sigma2,nu2)
    results <- list(estPars=estPars,post=post1,nit=nit,loglik=loglik)
  }
  return(results)
}
