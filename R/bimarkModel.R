# Here is R script to mock object organization of the bimark package

# convenience string for the scripts:
bmclass <- "BimarkModel"

#' Manipulate a simple \code{BimarkModel} object to generate and process
#' bilateral data.
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
#' @seealso bimarkSimulationModel bimarkObservationModel print.BimarkModel
#'
#' @docType package
#' @name bimark
NULL

# Utility method to create an empty model object
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
#' @examples
#' m <- bimarkSimulationModel()
#' print(m)
#'
#' @export

print.BimarkModel <- function(model, ...) { # {{{

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
#' @seealso \code{\link{generateLatentHistories}}
#'
#' @inheritParams generateLatentHistories
#'
#' @examples
#' set.seed(12)
#' m <- bimarkSimulationModel()
#'
#' @return the model object updated
#'
#' @export

bimarkSimulationModel <- function(N     = 500, # {{{
                                  T     = 5,
                                  P     = rep(.1, T),
                                  delta = c(0., 4., 4., 2., 3.)) {

  model <- createBimarkModel()

  L <- generateLatentHistories(N, T, P, delta)

  model$T <- T
  model$P <- P
  model$delta <- delta
  model$L <- L
  model$N <- nrow(L)
  model$T <- ncol(L)

  # then simulate the observation process of these latent histories:
  M <- observeHist(L)
  model <- addObservationToModel(model, M)

  return(model)

} # }}}

#' Populate model with observed raw histories matrix M
#'
#' This will also compute every little data summarizing M.. and preparing the
#' polytope processing.
#'
#' @param M a raw observation matrix as given by \code{\link{observeHist}}
#'
#' @return the model object updated
#'
#' @examples
#' m <- bimarkObservationModel(example.M)
#'
#' @seealso Hist2Id compute.Frequencies compute.Omega.B compute.A compute.B
#' orderHists
#'
#' @export

bimarkObservationModel <- function(M) { # {{{

  model <- createBimarkModel()

  return(addObservationToModel(model, M))

}

# convenience to be called from both `bimarkSimulationModel` and
# `bimarkObservationModel`
addObservationToModel <- function(model, M) {

  T <- ncol(M)
  counts <- compute.Frequencies(M)
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
  Omega.B <- compute.Omega.B(hists[o$L,], hists[o$R,])
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

#' Protect BimarkModel pseudo-members for encapsulation
#'
#' By overloading these functions, we expect the user not to be able to write to
#' the typed lists \code{BimarkModel}.
#'
#' If truly needed, change the type of the object with
#' \code{class(model) <- "HighjackModel"} then do whatever you want. However
#' there is no guarantee then that the package will keep functional nor
#' consistent.
#'
#' @examples
#' m <- bimarkObservationModel(example.M)
#' # read ok
#' m$LR
#' m[['LR']]
#' m['LR']
#' # no write
#' m$LR <- 82
#' m[['LR']] <- 82
#' m['LR'] <- list(LR=82)
#'
#' @name BimarkModelEncapsulation
NULL

#' @rdname BimarkModelEncapsulation
#' @usage model$attribute <- value # forbidden `$<-.BimarkModel`
#' @export
`$<-.BimarkModel` <- function(o, ...) {  # {{{
  # This should not forbid use by the package. So we need to check the
  # context here: if the operator is called from a regular `bimark:*` method,
  # then let it run as expected. Otherwise it is that the user has asked for
  # internal members edition: raise an error.
  if (isCalledByAFriend()) {
    # then we've been called by a friend, do as usual:
    return(nativeWriteResult(`$<-`, o, ...))
  }
  # Or raise an error and explain to user.. sorry :)
  warnUserEncapsulationViolation('$')
  return(o)
}

#' @rdname BimarkModelEncapsulation
#' @usage model[[attribute]] <- value # forbidden `[[<-.BimarkModel`
#' @export
`[[<-.BimarkModel` <- function(o, ...) {
  # This should not forbid use by the package. So we need to check the
  # context here: if the operator is called from a regular `bimark:*` method,
  # then let it run as expected. Otherwise it is that the user has asked for
  # internal members edition: raise an error.
  if (isCalledByAFriend()) {
    # then we've been called by a friend, do as usual:
    return(nativeWriteResult(`[[<-`, o, ...))
  }
  # Or raise an error and explain to user.. sorry :)
  warnUserEncapsulationViolation('[[')
  return(o)
}

#' @rdname BimarkModelEncapsulation
#' @usage model[attribute] <- value # forbidden `[<-.BimarkModel`
#' @export
`[<-.BimarkModel` <- function(o, ...) {
  # This should not forbid use by the package. So we need to check the
  # context here: if the operator is called from a regular `bimark:*` method,
  # then let it run as expected. Otherwise it is that the user has asked for
  # internal members edition: raise an error.
  if (isCalledByAFriend()) {
    # then we've been called by a friend, do as usual:
    return(nativeWriteResult(`[<-`, o, ...))
  }
  # Or raise an error and explain to user.. sorry :)
  warnUserEncapsulationViolation('[')
  return(o)
}

# Isolate checking caller function: is it allowed to edit model object?
isCalledByAFriend <- function() {
  level <- -3 # the caller -> `*<-.BimarkModel` -> `isCalledByAFriend`
  # First: list authorized methods: methods in this package:
  bimarkEnv <- environment(generateLatentHistories)
  friends <- objects(bimarkEnv) # here are their names
  # Second, get the calling function name:
  callerName <- as.character(sys.call(level)[1])
  if (length(callerName) == 0) callerName <- "" # dirty R adjustment.. 'hate R.
  # Third: if the names match, check that they are the actual same method:
  if (callerName %in% friends) {
    # get the actual caller function, the actual friend function, and compare
    caller <- sys.function(level)
    friend <- get(callerName, bimarkEnv)
    return(identical(caller, friend)) # that's it
  }
  return(FALSE)
}

# Isolate getting the regular, native `*<-` result:
# @param method name of the native method to call
# @param object, ... arguments received by the first call to `*<-`
nativeWriteResult <- function(method, object, ...) {
    # This is R fiddle-faddle, I've found no other way :(
    # change the class temporarily so the right `$<-` is called
    cl <- class(object)                # save the class
    class(object) <- NULL              # clear the class
    args <- c(list(object), list(...)) # gather into one arguments list
    res <- do.call(method, args)   # get native result
    class(res) <- cl              # restore the class
    return(res)                   # that's it
}

# Isolated warning the user of denied access
# @param operator char the operator user has tried using
warnUserEncapsulationViolation <- function(operator) {
  cat("The ", bmclass, " object is supposed to be encapsulated. ",
      "Write access is denied to the user.\n",
      "Please do not try editing its elements with `", operator,"`. ",
      "See `help(BimarkModelEncapsulation)`.\n", sep='')
} # }}}

#' Get polytope matrices from a BimarkModel object
#'
#' See \code{\link{compute.A}}, \code{\link{compute.B}}
#'
#' @param model a valid BimarkModel object
#'
#' @seealso \code{\link{compute.A}} \code{\link{compute.B}}
#' @export

get.A <- function(model) { # {{{
  return(compute.A(model$LL, model$LR))
}

#' @rdname get.A
#' @export
get.B <- function(model) {
  return(compute.B(model$LL, model$LR))
} # }}}

#' Get various parts of the Omega matrix
#'
#' See \code{\link{compute.Omega.B}}.
#'
#' @inheritParams get.A
#'
#' @return the corresponding Omega matrix
#'
#' @seealso \code{\link{compute.Omega.B}}
#' @export

get.Omega <- function(model) { # {{{
  return(ID2Hist(model$iOmega, model$T))
}

#' @rdname get.Omega
#' @export
get.Omega.S <- function(model) {
  LS <- model$LS
  iOmega.S <- model$iOmega[1:LS]
  return(ID2Hist(iOmega.S, model$T))
}

#' @rdname get.Omega
#' @export
get.Omega.L <- function(model) {
  LS <- model$LS
  LL <- model$LL
  iOmega.L <- model$iOmega[LS + (1:LL)]
  return(ID2Hist(iOmega.L, model$T))
}

#' @rdname get.Omega
#' @export
get.Omega.R <- function(model) {
  LS <- model$LS
  LL <- model$LL
  LR <- model$LR
  iOmega.R <- model$iOmega[LS + LL + (1:LR)]
  return(ID2Hist(iOmega.R, model$T))
}

#' @rdname get.Omega
#' @export
get.Omega.B <- function(model) {
  LS <- model$LS
  LL <- model$LL
  LR <- model$LR
  LB <- LL * LR
  iOmega.B <- model$iOmega[LS + LL + LR + (1:LB)]
  return(ID2Hist(iOmega.B, model$T))
}

#' @rdname get.Omega
#' @export
get.Omega.SLR <- function(model) {
  LS <- model$LS
  LL <- model$LL
  LR <- model$LR
  iOmega.SLR <- model$iOmega[1:(LS + LL + LR)]
  return(ID2Hist(iOmega.SLR, model$T))
}

#' @rdname get.Omega
#' @export
get.Omega.LR <- function(model) {
  LS <- model$LS
  LL <- model$LL
  LR <- model$LR
  iOmega.LR <- model$iOmega[LS + 1:(LL + LR)]
  return(ID2Hist(iOmega.LR, model$T))
} # }}}

