# Here is R script to mock object organization of the bimark package

# convenience string for the scripts:
bmclass <- "BimarkModel"

#' Manipulate a simple \code{BimarkModel} object to generate and process
#' bilateral data.
#'
#' This object will contain every relevant data. And they will be consistent
#' at any time.
#'
#' R's S4 object opportunities sound nice but.. terrible to use. Furthermore, we
#' do not need such a sophisticated object architecture yet. As an alternative,
#' we are falling back on this plain named-and-typed list which will gather
#' information as the analysis goes by. Items are used as pseudo-slots. Their
#' values are set to NULL if they are not defined.
#'
#' @slot L raw latent histories matrix if available (available if simulated)
#' @slot N number of latent histories: actual number of individuals!
#' @slot T number of capture occasions
#' @slot M raw observation matrix (see \code{\link{example.M}})
#' @slot n number of observed histories (side and/or individuals)
#' @slot F frequencies of observed histories in canonical order (see
#' \code{\link{ComputeFrequencies}})
#' @slot LS number of S-histories observed
#' @slot LL number of L-histories observed
#' @slot LR number of R-histories observed
#' @slot iOmega latent histories matrix (stored as ids in canonical order) (see
#' \code{\link{GetOmega.B}})
#'
#' @seealso \code{\link{BimarkSimulationModel}}
#' \code{\link{BimarkObservationModel}} \code{\link{print.BimarkModel}}
#' \code{\link{BimarkModelEncapsulation}}
#'
#' @examples
#' model <- BimarkSimulationModel()
#'
#' class(model)
#' print(model)
#' model$n
#' model$F
#' model$iOmega
#'
#' model$n <- 12 # do not try editing these members
#'
#' @docType package
#' @name BimarkModel
NULL

# Utility method to create an empty model object
CreateBimarkModel <- function() { # {{{
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
#' @seealso \code{\link{BimarkModel}}
#'
#' @examples
#' m <- BimarkSimulationModel()
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

#' Create a bimark model object based on simulated data
#'
#' This model will be full because latent histories are known. About the actual
#' returned R object, see \code{help(BimarkModel)}.
#'
#' @seealso \code{\link{GenerateLatentHistories}} \code{\link{BimarkModel}}
#'
#' @inheritParams GenerateLatentHistories
#'
#' @return a \code{\link{BimarkModel}} object
#'
#' @examples
#' set.seed(12)
#' m <- BimarkSimulationModel()
#'
#' @export

BimarkSimulationModel <- function(N     = 500, # {{{
                                  T     = 5,
                                  P     = rep(.1, T),
                                  delta = c(0., 4., 4., 2., 3.)) {

  model <- CreateBimarkModel()

  L <- GenerateLatentHistories(N, T, P, delta)

  model$T <- T
  model$P <- P
  model$delta <- delta
  model$L <- L
  model$N <- nrow(L)
  model$T <- ncol(L)

  # then simulate the observation process of these latent histories:
  M <- ObserveHist(L)
  model <- AddObservationToModel(model, M)

  return(model)

} # }}}

#' Create a birmark model object based on observed data
#'
#' This model will not be full because latent histories are unknown. Abouth the
#' actual returned R object, see \code{help(BimarkModel)}.
#'
#' @param M a raw observation matrix formatted as \code{\link{example.M}}
#'
#' @return a \code{\link{BimarkModel}} object
#'
#' @seealso \code{\link{example.M}} \code{\link{capture.events}}
#' \code{\link{BimarkModel}}
#'
#' @examples
#' m <- BimarkObservationModel(example.M)
#'
#' @export

BimarkObservationModel <- function(M) { # {{{

  model <- CreateBimarkModel()

  return(AddObservationToModel(model, M))

}

# convenience to be called from both `BimarkSimulationModel` and
# `BimarkObservationModel`
AddObservationToModel <- function(model, M) {

  T <- ncol(M)
  counts <- ComputeFrequencies(M)
  hists <- ID2Hist(counts$id, T) # unique hists
  # reorganize the canonical way:
  o <- OrderHists(hists)
  LS <- length(o$S)
  LL <- length(o$L)
  LR <- length(o$R)
  canonicalOrder <- c(o$S, o$L, o$R)
  iOmega.SLR <- counts$id[canonicalOrder]
  F <- counts$F[canonicalOrder]

  # Compute Omega.B and nullspace informations:
  Omega.B <- ComputeOmegaB(hists[o$L,], hists[o$R,])
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
#' the typed lists \code{\link{BimarkModel}}.
#'
#' If truly needed, change the type of the object with
#' \code{class(model) <- "HighjackModel"} then do whatever you want. However
#' there is no guarantee then that the package will keep functional nor
#' consistent.
#'
#' @seealso \code{\link{BimarkModel}}
#'
#' @examples
#' m <- BimarkObservationModel(example.M)
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
  if (IsCalledByAFriend()) {
    # then we've been called by a friend, do as usual:
    return(NativeWriteResult(`$<-`, o, ...))
  }
  # Or raise an error and explain to user.. sorry :)
  WarnUserEncapsulationViolation('$')
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
  if (IsCalledByAFriend()) {
    # then we've been called by a friend, do as usual:
    return(NativeWriteResult(`[[<-`, o, ...))
  }
  # Or raise an error and explain to user.. sorry :)
  WarnUserEncapsulationViolation('[[')
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
  if (IsCalledByAFriend()) {
    # then we've been called by a friend, do as usual:
    return(NativeWriteResult(`[<-`, o, ...))
  }
  # Or raise an error and explain to user.. sorry :)
  WarnUserEncapsulationViolation('[')
  return(o)
}

# Isolate checking caller function: is it allowed to edit model object?
IsCalledByAFriend <- function() {
  level <- -3 # the caller -> `*<-.BimarkModel` -> `IsCalledByAFriend`
  # First: list authorized methods: methods in this package:
  bimarkEnv <- environment(GenerateLatentHistories)
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
NativeWriteResult <- function(method, object, ...) {
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
WarnUserEncapsulationViolation <- function(operator) {
  cat("The ", bmclass, " object is supposed to be encapsulated. ",
      "Write access is denied to the user.\n",
      "Please do not try editing its elements with `", operator,"`. ",
      "See `help(BimarkModelEncapsulation)`.\n", sep='')
} # }}}

#' Generate A observation process matrix from a BimarkModel object
#'
#' This version of the A matrix, as in Bonner2013, only takes into account the
#' histories actually observed in the data, and the only latent histories that
#' could have generated them. Each latent history (in rows) generates the
#' observed histories corresponding to the columns in which they contain a 1.
#' Thanks to polytope ordering, its generation is really simple since its
#' pattern only depends on the number of L-histories (LL) and R-histories (LR).
#' Here is the pattern:\preformatted{
#'       ______L_L_R_R_R-observable_and_ghosts___________________
#'        a L  1 0 0 0 0 L-block |
#'        b_L__0_1_0_0_0_________|
#'        1 R  0 0 1 0 0         | Identity
#'        2 R  0 0 0 1 0 R-block |
#'  A:    3_R__0_0_0_0_1_________|_______________________________
#'       a1 B  1 0 1 0 0         |
#'       a2 B  1 0 0 1 0         |
#'       a3 B  1 0 0 0 1 B-block | Combinations in polytope order
#'       b1 B  0 1 1 0 0         |
#'       b2 B  0 1 0 1 0         |
#'       b3_B__0_1_0_0_1_________|_______________________________}
#'
#' The matrix A is actually not directly useful for JAGS because the polytope is
#' better described by its kernel B.
#'
#' @param model a valid BimarkModel object
#'
#' @return the observation processs matrix A as a
#' \code{\link[Matrix]{sparseMatrix}} : a LM x LU matrix with ones and zeroes,
#' each row corresponding to a latent history sorted in polytope order, each
#' column corresponding to an observable history sorded in polytope (or
#' "canonical") order.
#'
#' @seealso \code{\link{get.B}} \code{\link{ComputeA}}
#' @export

get.A <- function(model) {
  return(ComputeA(model$LL, model$LR))
}

#' Generate B matrix, nullspace of t(A), from a BimarkModel object
#'
#' The polytope is directly described by a kernel of the observation process
#' matrix A (more exactly t(A)). Just like A, its pattern is quite easy-to-grasp
#' once the histories has been sorted in polytope order. So it is easy to
#' generate and it only depends on the number of L-histories (LL) and of
#' R-histories (LR). Here is the pattern:\preformatted{
#'       ______B_B_B_B_B_B-latent_histories______________________
#'        a L  - - - 0 0 0 L-block |
#'        b_L__0_0_0_-_-_-_________|
#'        1 R  - 0 0 - 0 0         | '-' are minus ones -1's
#'        2 R  0 - 0 0 - 0 R-block |
#'  B:    3_R__0_0_-_0_0_-_________|_______________________________
#'       a1 B  + 0 0 0 0 0         |
#'       a2 B  0 + 0 0 0 0         |
#'       a3 B  0 0 + 0 0 0 B-block | Identity : '+' are ones 1's
#'       b1 B  0 0 0 + 0 0         |
#'       b2 B  0 0 0 0 + 0         |
#'       b3_B__0_0_0_0_0_+_________|_______________________________}
#'
#
#' @inheritParams get.A
#'
#' @return the polytope matrix B as a \code{\link[Matrix]{sparseMatrix}} : a LM
#' x LB matrix with ones, zeroes and minus ones, each row corresponding to a
#' latent history sorted in polytope order, each column corresponding to an
#' unobservable B-history sorded in polytope order.
#'
#' @seealso \code{\link{get.A}} \code{\link{ComputeB}}
#' @export

get.B <- function(model) {
  return(ComputeB(model$LL, model$LR))
}

#' Get various parts of the Omega matrix
#'
#' \code{GetOmega.B} function computes all possible latent, unobservable
#' histories (aka B-histories) which may have generated the given letf- and
#' right-histories (L- and R-histories). In a nutshell, it generates Omega.B
#' from Omega.L and Omega.R. Consider it as a kind of reversed
#' \code{\link{ObserveHist}}. Omega.B is sorted with polytope order.
#'
#' Each B-history comes from the combination of one L- and one R-history. Here
#' are the calculation rules for each couple of simultaneous observed events:
#'     \itemize{
#'     \item{Rule 1: }{observed 0 and 0 together may come from a latent 0}
#'     \item{Rule 2: }{observed 0 and L together may come from a latent L}
#'     \item{Rule 3: }{observed 0 and R together may come from a latent R}
#'     \item{Rule 4: }{observed L and R together may come from a latent B}
#'     }
#'
#' Polytope order for B-histories is the nested order of canonical R-order
#' within canonical L-order, as follows: \preformatted{
#'       |          |  S - no bilateral ambiguity
#'       | Omega.S  |  S - no bilateral ambiguity
#'       |         |   L - 1 left side ambiguous history
#'       | Omega.L |   L - 2 left side ambiguous history
#'       |          |  R - 1 right side ambiguous history
#'       | Omega.R  |  R - 2 right side ambiguous history
#' Omega |          |  R - 3 right side ambiguous history
#'       |         |   B - generating L1 and R1
#'       |         |   B - generating L1 and R2
#'       | Omega.B |   B - generating L1 and R3
#'       |         |   B - generating L2 and R1
#'       |         |   B - generating L2 and R2
#'       |         |   B - generating L2 and R3}
#'
#'
#' @inheritParams get.A
#'
#' @return the corresponding Omega matrix
#'
#' @seealso \code{\link{OrderHists}} \code{\link{ComputeOmegaB}}
#' @export

GetOmega <- function(model) { # {{{
  return(ID2Hist(model$iOmega, model$T))
}

#' @rdname GetOmega
#' @export
GetOmegaS <- function(model) {
  LS <- model$LS
  iOmega.S <- model$iOmega[1:LS]
  return(ID2Hist(iOmega.S, model$T))
}

#' @rdname GetOmega
#' @export
GetOmega.L <- function(model) {
  LS <- model$LS
  LL <- model$LL
  iOmega.L <- model$iOmega[LS + (1:LL)]
  return(ID2Hist(iOmega.L, model$T))
}

#' @rdname GetOmega
#' @export
GetOmega.R <- function(model) {
  LS <- model$LS
  LL <- model$LL
  LR <- model$LR
  iOmega.R <- model$iOmega[LS + LL + (1:LR)]
  return(ID2Hist(iOmega.R, model$T))
}

#' @rdname GetOmega
#' @export
GetOmega.B <- function(model) {
  LS <- model$LS
  LL <- model$LL
  LR <- model$LR
  LB <- LL * LR
  iOmega.B <- model$iOmega[LS + LL + LR + (1:LB)]
  return(ID2Hist(iOmega.B, model$T))
}

#' @rdname GetOmega
#' @export
GetOmegaSLR <- function(model) {
  LS <- model$LS
  LL <- model$LL
  LR <- model$LR
  iOmega.SLR <- model$iOmega[1:(LS + LL + LR)]
  return(ID2Hist(iOmega.SLR, model$T))
}

#' @rdname GetOmega
#' @export
GetOmega.LR <- function(model) {
  LS <- model$LS
  LL <- model$LL
  LR <- model$LR
  iOmega.LR <- model$iOmega[LS + 1:(LL + LR)]
  return(ID2Hist(iOmega.LR, model$T))
} # }}}

