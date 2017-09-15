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

  # TODO: better (degenerated-)multi-case organisation: one common processing
  # then each part should add its small contribution to the building of
  # `data.list`, `inits`, `monitor`, etc.

  # Prepare values for all required fixed nodes:
  # basic needed values:
  LS <- model$LS
  LL <- model$LL
  LR <- model$LR
  LU <- LL + LR
  LB <- LL * LR
  LM <- LU + LB
  L <- LS + LM

  # in case the model is degenerated, switch to the more adapted one:
  degenerated <- LB == 0
  if (degenerated) {
    method <- 'degenerated'
    cat(paste0('No unobservable histories in the model. Switching to ',
               method, ' method.\n'))
  }

  # retrieve the desired model
  jagsFile <- GetBayesianModel(method)

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
  # degenerated cases:
  if (LL == 0) L.f <- integer(0)
  if (LR == 0) R.f <- integer(0)
  x0 <- c(L.f, R.f, rep(0, LB))
  # maximum values for latent counts
  # or there would be too many left- or right-histories
  box.ceilings <- vapply(1:LB, function(i) # watch canonical order of B-hists!
                     min(L.f[((i - 1) %/% LR) + 1],
                         R.f[((i - 1) %% LR) + 1]), 1)

  # Fill data node (see `ins/jags/*.jags` for node descriptions)
  # common nodes
  data.list <- list(T = model$T,
                    nbCaptureEvents = nb.capture.events,
                    jagsLS = jagsLS,
                    LL = LL,
                    LR = LR,
                    LU = LU,
                    LM = LM,
                    LB = LB,
                    L  =  L,
                    is.Sf.empty = is.Sf.empty,
                    jagsS.f = jagsS.f,
                    Omega = GetOmega(model),
                    PoissonC = 1e5,
                    PoissonTrick = 1
                    )
  # nodes to monitor/track during MCMC process (many for debugging)
  monitor <- c('ln.Likelihood', 'jagsS.f', 'n', 'N', 'G', 'P', 'p', 'delta')
  # initial nodes
  inits <- list()

  # case-specific nodes
  if (degenerated) {
    data.list <- c(data.list, list(fixed.X = c(jagsS.f, L.f, R.f)))
  } else {
    data.list <- c(data.list, list(
                                   L.f = L.f,
                                   R.f = R.f,
                                   x0 = x0,
                                   box.ceilings = box.ceilings,
                                   B = as.matrix(ComputeB(LL, LR)),
                                   BernoulliTrick = 1
                                   ))
    monitor <- c(monitor, c('X', 'B', 'x', 'd', 'x0', 'allxPos',
                            'BernoulliTrick', 'BernoulliTrickP'))
    # Initial values must be consistent with BernoulliTrick at least: all zeroes
    # for the polytope, say.
    inits <- c(inits,
               list(unifs=rep(0, LB))) # => d=rep(0, LB), but it is "fixed node"
  }

  # TODO: mak'em parameters
  RNG <- list(.RNG.name="base::Marsaglia-Multicarry", .RNG.seed=8)

  # build
  jm <- rjags::jags.model(jagsFile,
                          data=c(data.list, priors),
                          # seed RNGs
                          inits=c(inits, RNG), # TODO: parameters
                          quiet=TRUE,
                          )

  # Run!
  # length of the chain
  n.iter <- 1e3 # TODO: make it a parameter
  # only ONE chain will be processed. Launch in parallel if you need several.
  n.chains <- 1
  mc <- rjags::coda.samples(jm,
                            variable.names=monitor,
                            n.iter=n.iter,
                            n.chains=n.chains)[[1]] # only one chain

  return(mc)

}
# set.seed(6)  # got another bug here: 'lgamma' output + R error
# set.seed(11) # got another bug there: degenerated case: solved
# m <- BimarkSimulationModel(N=20, T=5)
# mc <- EstimateLatentCounts(m)
# }}}

