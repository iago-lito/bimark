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
#' @slot n number of observed histories (side and/or individuals)
#' @slot F frequencies of observed histories in canonical order
#' @slot LS number of S-histories observed
#' @slot LL number of L-histories observed
#' @slot LR number of R-histories observed
#' @slot iOmega latent histories matrix (stored as ids in canonical order)
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
                M=NULL,
                n=NULL,
                F=NULL,
                LS=NULL,
                LL=NULL,
                LR=NULL,
                iOmega=NULL
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
  if (is.null(model$M)) {
    cat("empty", bmclass, '\n')
    return()
  }

  cat(bmclass, "object:\n")

  # common script to the two sections?
  arranger <- function(title, tab0=2, tab1=4, l) {
    cat('\n', strrep(' ', tab0), title, '\n', sep='')
    maxLabel <- max(vapply(l, nchar, 1)) + tab1
    maxName <- max(vapply(names(l), nchar, 1))
    for (name in names(l)) {
      label <- l[[name]]
      nl <- nchar(label)
      nn <- nchar(name)
      beforeLabel <- strrep(' ', tab1 + (maxLabel - tab1 - nl))
      beforeName <- strrep(' ', maxName - nn)
      value <- model[[name]]
      if (is.null(value)) value <- '??'
      cat(beforeLabel, label, " : ", beforeName, name, " = ",
          do.call(paste, as.list(value)), "\n", sep='')
    }
  }

  # About Simulated data:
  arranger("Latent data ($L):", 2, 4,
           list(N="individuals",
                T="capture occasions",
                P="capture probabilities",
                delta="capture events probabilities"
                ))

  # About observed data:
  arranger("Observed data ($M):", 2, 4,
           list(T="capture occasions",
                n="number of records",
                LS="number of simulaneous histories",
                LL="number of left-histories",
                LR="number of right-histories"
                ))

} # }}}

#' Populate model object with simulated latent data
#'
#' @seealso generateLatentHistories
#'
#' @param model a virgin BimarkModel object
#' @param ... further arguments to generateLatentHistories
#'
#' @examples
#' m <- createBimarkModel()
#' m <- addLatentHistories(m)
#'
#' @return the model object updated
#'
#' @export

addLatentHistories <- function(model, N     = 500, # {{{
                                      T     = 5,
                                      P     = rep(.1, T),
                                      delta = c(0., 4., 4., 2., 3.)) {

  L <- generateLatentHistories(N, T, P, delta)

  model$T <- T
  model$P <- P
  model$delta <- delta
  model$L <- L
  model$N <- nrow(L)
  model$T <- ncol(L)

  # then simulate the observation process of these latent histories:
  M <- observeHist(L)
  model <- addObservedHistories(model, M)

  return(model)

} # }}}

#' Populate model with observed raw histories matrix M
#'
#' This will also compute every little data summarizing M.. and preparing the
#' polytope processing.
#'
#' @param model a virgin BimarkModel object
#' @param M a raw observation matrix
#'
#' @return the model object updated
#'
#' @examples
#' m <- createBimarkModel()
#' m <- addObservedHistories(m, example.M)
#'
#' @seealso Hist2Id getFrequencies getOmega.B getA getB orderHists
#'
#' @export

addObservedHistories <- function(model, M) { # {{{

  T <- ncol(M)
  counts <- getFrequencies(M)
  hists <- ID2Hist(counts$id, T) # unique hists
  # reorganize the canonical way:
  o <- orderHists(hists)
  LS <- length(o$S)
  LL <- length(o$L)
  LR <- length(o$R)
  canonicalOrder <- c(o$S, o$L, o$R)
  iOmega.SLR <- counts$id[canonicalOrder]
  F <- counts$F[canonicalOrder]

  # Compute Omega.B and nullspace informations:
  Omega.B <- getOmega.B(hists[o$L,], hists[o$R,])
  iOmega <- c(iOmega.SLR, Hist2ID(Omega.B))

  # That's it!
  model$T <- ncol(M)
  model$n <- nrow(M)
  model$M <- M
  model$F <- F
  model$LS <- LS
  model$LL <- LL
  model$LR <- LR
  model$iOmega <- iOmega
  return(model)

} # }}}

