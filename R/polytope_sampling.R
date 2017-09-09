# Here is the script where bilateral data will be Bayesian-processed. This is
# where we call JAGS, where we sample into a polytope and where our custom C++
# JAGS module is used.

#' Get path to jags files used within bimark
#'
#' Bimark uses internal jags files. Use this function to retrieve their location
#' on your system.
#'
#' @param name string the name of the model
#'
#' @return a path to the correspindong jags file
#'
#' @examples
#' GetBayesianModel('test')
#'
#' @export

GetBayesianModel <- function(name) { # {{{
  filename <- paste0(name, '.jags')
  jagsFile <- system.file('jags', filename, package='bimark')
  if (!file.exists(jagsFile))
    throw(paste("Cannot find", filename, "model file."))
  return(jagsFile)
}
# }}}

#' Run a dummy Bayesian MCMC procedure
#'
#' If this runs well, it means that the whole R + JAGS procedure is okay:
#' model files are found, \code{JAGS} is reachable, \code{rjags} package is here
#' and works fine.
#'
#' The dummy procedure draws N random data from an univariate Gaussian, then
#' estimates the mean and variance of this Gaussian with one MCMC. Estimators
#' are searched with uniform priors within a short range around the actual
#' values.
#'
#' @param N integer, the number of data to draw
#' @param n.iter integer, the number of MCMC iterations
#' @param mu real, dummy mean for the Gaussian
#' @param sigma positive real, dummy standard deviance for the Gaussian
#' @param priors range of uniform priors around mu and sigma: \code{c(lowerMu,
#' upperMu, lowerSigma, upperSigma)}
#'
#' @examples
#' DummyJags()
#'
#' @export

DummyJags <- function(N=1e3, n.iter=1e3, mu=15., sigma=.3, # {{{
                      priors=c(-20., 20., 0., 10.)) {

  # retrieve the dummy model
  jagsFile <- GetBayesianModel('test')

  # JAGS's dnorm does not use sigma but tau -_-
  tau <- 1. / sigma ^ 2
  data <- rnorm(N, mu, sigma)
  jm <- rjags::jags.model(jagsFile,
                          data=list(a=data, N=N, priors=priors), quiet=TRUE)
  mc <- rjags::coda.samples(jm,
                            variable.names=c('mu', 'tau'),
                            n.iter=n.iter,
                            n.chains=1)[[1]] # only one chain
  mus <- as.numeric(mc[1:n.iter, 'mu'])
  taus <- as.numeric(mc[1:n.iter, 'tau'])
  sigmas <- 1 / sqrt(taus)
  cat(paste0(paste0("mu estimate: ", mean(mus),
                    " vs real: ", mu, '\n')))
  print(paste0(paste0("sigma estimate: ", mean(sigmas),
                      " vs real: ", sigma, '\n')))
  par(mfrow=c(2, 1))
  plot(mus, type='l')
  plot(taus, type='l')

  cat("Okay, everything went well.\n")

}
# }}}

#' Estimate latent counts with bayesian sampling inside the polytope
#'
#' For now, it is just a sandbox to make JAGS work again from this repo. Start
#' with naive sampling algorithm where we sample inside the bounding box and
#' reject points outside of the polytope.
#'
#' @param model a BimarkModel object
#' @param method string the jags model to use, for now 'test' is the only one
#' available
#'
#' @return an updated BimarkModel object with the estimators, the chains traces
#' and their properties.
#'
#' @examples
#' m <- BimarkSimulationModel(N=20, T=5)
#' EstimateLatentCounts(m)
#'
#' @export

EstimateLatentCounts <- function(model, method='boundingBox', # {{{
                                 priors=list(h.P=c(1, 1),
                                             h.delta=rep(1, nb.capture.events)))
{

  # temp for building
  model <- BimarkSimulationModel(N=50, T=5)
  method <- 'boundingBox'

  # retrieve the desired model
  jagsFile <- GetBayesianModel(method)

  # Prepare values for all required fixed nodes:
  # basic needed values:
  LS <- model$LS
  LL <- model$LL
  LR <- model$LR
  LU <- LL + LR
  LB <- LL * LR
  LM <- LU + LB
  L <- LS + LL + LR + LB
  # data counts
  F <- model$F
  # special LS = 0 case to handle
  if (LS > 0) {
    is.Sf.empty <- 0
    jagsLS <- LS
    jagsS.f <- F[1:LS]
    L.f <- F[(LS+1):(LS+LL)]
    R.f <- F[(LS+LL+1):(LS+LL+LR)]
  } else {
    is.Sf.empty <- 1
    jagsLS <- 1 # fill it with one dummy piece of data
    jagsS.f <- 0 # dummy data read nowhere except in JAGS: n <- sum(X)
    L.f <- F[1:LL]
    R.f <- F[(LL+1):(LL+LR)]
  }
  # maximum values for latent counts
  # or there would be too many left- or right-histories
  box.ceilings <- vapply(1:LB, function(i) # watch canonical order of B-hists!
                     min(L.f[((i - 1) %/% LR) + 1],
                         R.f[((i - 1) %% LR) + 1]), 1)

  # see `boundingBox.jags` for data node descriptions
  data.list <- list(T = model$T,
                    nbCaptureEvents = nb.capture.events,
                    jagsLS = jagsLS,
                    LL = LL,
                    LR = LR,
                    LU = LU,
                    L  =  L,
                    is.Sf.empty = is.Sf.empty,
                    jagsS.f = jagsS.f,
                    L.f = L.f,
                    R.f = R.f,
                    Omega = GetOmega(model),
                    box.ceilings = box.ceilings,
                    B = as.matrix(ComputeB(LL, LR)))

  jm <- rjags::jags.model(jagsFile,
                          data=c(data.list, priors), quiet=TRUE)

  return(model)

}
# }}}

