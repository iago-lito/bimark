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
#' getBayesianModel('test')
#'
#' @export

getBayesianModel <- function(name) {
  filename <- paste0(name, '.jags')
  jagsFile <- system.file('jags', filename, package='bimark')
  if (!file.exists(jagsFile))
    throw(paste("Cannot find", filename, "model file."))
  return(jagsFile)
}

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
#' @export

estimateLatentCounts <- function(model, method='test') {

  # retrieve the right model
  jagsFile <- getBayesianModel('test')

  return(model)

}

