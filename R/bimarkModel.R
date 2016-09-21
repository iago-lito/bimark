# Here is R script to mock object organization of the bimark package

# convenience string for the scripts:
bmclass <- "BimarkModel"

#' Create the data object that'll endure all this package processing
#'
#' This object will contain every relevant data. And they will be *consistent*
#' at any time.. provided one do not try to set the slots by oneself: nothing is
#' protected!
#'
#' R's S4 object opportunities sound nice but.. terrible to use. Furthermore, we
#' do not need such a sophisticated object architecture. As an alternative, we
#' are falling back on this plain named-and-typed list which will gather
#' information as the analysis goes by. Items are simili-slots. Their values are
#' set to NULL if they are not defined.
#'
#' @slot L raw latent histories matrix if available (available if simulated)
#' @slot N number of latent histories: actual number of individuals!
#' @slot T number of capture occasions
#' @slot M raw observation matrix
#'
#' @seealso print.BimarkModel
#'
#' @examples
#' m <- createBimarkModel()
#'
#' @export

createBimarkModel <- function() { # {{{
  model <- list(
                L=NULL,
                T=NULL,
                N=NULL,
                M=NULL
                )
  class(model) <- bmclass
  return(model)
} # }}}

#' Print BimarkModel pseudo-object to the console
#'
#' @seealso createBimarkModel
#'
#' @examples
#' m <- createBimarkModel()
#' m
#' print(m)
#'
#' @export

print.BimarkModel <- function(model, ...) { # {{{

  # Consider it as empty if there are no raw histories matrices:
  if (is.null(model$L) && is.null(model$M)) {
    cat("empty", bmclass, '\n')
    return()
  }

  cat(bmclass, "object:\n")

  # About Simulated data:
  if (!is.null(model$L)) {
    tab0 <- 2
    cat('\n', strrep(' ', 2), "Simulated data ($L): \n", sep='')
    tab1 <- 4
    l <- list(N="individuals",
              T="capture occasions",
              P="capture probabilities",
              delta="capture events probabilities"
              )
    maxLabel <- max(vapply(l, nchar, 1)) + tab1
    maxName <- max(vapply(names(l), nchar, 1))
    for (name in names(l)) {
      label <- l[[name]]
      nl <- nchar(label)
      nn <- nchar(name)
      beforeLabel <- strrep(' ', tab1 + (maxLabel - tab1 - nl))
      beforeName <- strrep(' ', maxName - nn)
      cat(beforeLabel, label, " : ", beforeName, name, " = ",
          do.call(paste, as.list(model[[name]])), "\n", sep='')
    }
  }

} # }}}

#' Populate model object with simulated latent data
#'
#' @seealso generateLatentHistories
#'
#' @param model a BimarkModel object
#' @param ... further arguments to generateLatentHistories
#'
#' @examples
#' m <- createBimarkModel()
#' m <- addLatentHistories(m)
#' addLatentHistories(m)
#'
#' @export

addLatentHistories <- function(model, N     = 500, # {{{
                                      T     = 5,
                                      P     = rep(.1, T),
                                      delta = c(0., 4., 4., 2., 3.)) {

  # TODO: better handle the case where they already exist. This'll need to check
  # for further dependent data and erase them.
  if (!is.null(model$L) | !is.null(model$M))
    throw("This model object is not empty. Better use a fresh one.")

  L <- generateLatentHistories(N, T, P, delta)

  model$T <- T
  model$P <- P
  model$delta <- delta
  model$L <- L
  model$N <- nrow(L)
  model$T <- ncol(L)

  return(model)

} # }}}

