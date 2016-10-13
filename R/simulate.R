# Here is R script to simulate and manipulate bilateral data:

# Defining events of capture: # {{{

#' Capture events and their integer representations
#'
#'  \tabular{ccc}{
#'    0 : \tab no capture                     \tab : 0\cr
#'    L : \tab right capture                  \tab : 1\cr
#'    R : \tab left capture                   \tab : 2\cr
#'    B : \tab both sides capture, no match   \tab : 3\cr
#'    S : \tab simultaneous capture and match \tab : 4\cr
#'  }
#'
#' @seealso \code{\link{nbCaptureEvents}}
#'
#' @export

captureEvents <- c(
    "0"  = 0L
  , "L"  = 1L
  , "R"  = 2L
  , "B"  = 3L
  , "S"  = 4L
  )

#' Number of capture events
#'
#' This is just defined as \code{length(\link{captureEvents})}
#'
#' @seealso \code{\link{captureEvents}}
#'
#' @export

nbCaptureEvents <- length(captureEvents)

#' One small example of an observation matrix.
#'
#' The observation matrix M has n lines and T columns, with n being the number
#' of records and T the number of capture histories. Events are coded by an
#' integer accordingly to \code{\link{captureEvents}}.
#'
#' This example matrix is drawn from Bonnici's M2 report.
#'
#' @examples
#' print(example.M)
#' seeHist(example.M)
#'
#' @export

example.M <- matrix(captureEvents[
            c("0", "0", "R", "0", "0",
              "0", "0", "0", "0", "L",
              "L", "0", "B", "S", "0",
              "0", "0", "L", "0", "L",
              "0", "0", "R", "0", "0",
              "0", "R", "0", "0", "0",
              "0", "R", "0", "0", "0",
              "0", "0", "L", "0", "L",
              "0", "0", "0", "0", "L",
              "S", "0", "L", "0", "0",
              "0", "R", "0", "0", "0",
              "0", "R", "0", "0", "0",
              "0", "0", "0", "0", "L",
              "0", "0", "0", "0", "L",
              "0", "R", "0", "0", "0",
              "0", "0", "0", "0", "L",
              "0", "0", "R", "0", "0",
              "S", "0", "L", "0", "0",
              "0", "0", "L", "0", "L",
              "0", "0", "R", "0", "0",
              "0", "R", "0", "0", "0",
              "R", "0", "0", "0", "R")], ncol=5, byrow=TRUE)

# }}}

# custom error thrower for the package.
throw <- function(message){
  call <- sys.call(-1)
  callingFunctionName <- as.character(call[1])
  stop(paste0("bimark: ", callingFunctionName, ": ", message), call.=FALSE)
}

# changing default behaviour of submatricing: no more dropping of dimensions or
# I'll leave you, R!
`[` <- function(...) base::`[`(..., drop=FALSE)

#' Generate latent histories
#'
#' Simulate a virtual population where each individual undergoes a particular,
#' latent capture history.
#'
#' @param N integer number of simulated individual
#' @param T integer number of capture occasions
#' @param P \code{N}-vector: capture probabilities on each occasion (one value
#' for constant probability constant)
#' @param delta \code{\link{nbCaptureEvents}}-vector: capture probabilities for
#' each capture event (or positive weights whose sum will be normalized to 1)
#'
#' @seealso \code{\link{seeHists}} \code{\link{observeHists}}
#'
#' @return a matrix of capture events: one row per individual, one column per
#' capture occasion: each row is an individual capture history. These are the
#' latent histories.
#'
#' @examples
#' set.seed(12)
#' generateLatentHistories(N=20, T=5, P=rep(.1, 5), delta=c(0, 4, 4, 2, 3))
#'
#' @export

generateLatentHistories <- function(N     = 500, # {{{
                                    T     = 5,
                                    P     = rep(.1, T),
                                    delta = c(0., 4., 4., 2., 3.)
                                    ){

  # small checks for parameters consistency
  if (length(P) == 1)
    P <- rep(P, T)
  else if (length(P) != T)
    throw("Inconsistent number of capture probabilities.")
  if (any(P > 1.) || any(P < 0.))
    throw("Inconsistent capture probabilities.")
  if (length(delta) != nbCaptureEvents)
    throw("Inconsistent number of capture event probabilities.")
  if (any(delta < 0.))
    throw("No negative weights for capture event probabilities.")
  # normalize delta to probabilities
  delta <- delta / sum(delta)

  # First, the "structure" of capture histories (only ones and zeroes)
  L <- t(replicate(N, stats::rbinom(T, 1, P)))

  # Then, find out for each capture on which side it happened
  L[as.logical(L)] <- sample(captureEvents, sum(L), prob=delta, replace=TRUE)

  return(L)

} # }}}

#' Compute ID of a raw history
#'
#' Each history is assigned a unique identifier according to McClintock2010
#' (equation (3)). It is \code{\link{ID2Hist}} reverse function.
#'
#' The ids are integers, too big to be handled by R's `int` primitive type. So
#' we use strings to represent them.
#'
#' @param history for T capture events, either an events vector of length T or
#' a matrix of n histories \code{dim(n, T)}
#'
#' @seealso \code{\link{ID2Hist}} \code{\link{generateLatentHistories}}
#'
#' @return if \code{history} is one history, its id. If \code{history} is
#' several history, their corresponding vector of ids. They are (potentially
#' long) integers represented as strings.
#'
#' @examples
#' Hist2ID(c(0, 4, 0, 1, 3))
#' Hist2ID(matrix(c(0, 4, 0, 1, 3,
#'                  1, 1, 0, 2, 0,
#'                  0, 3, 4, 0, 0,
#'                  0, 0, 0, 4, 4), 4, 5, byrow=TRUE))
#' Hist2ID(generateLatentHistories(4, 5))
#'
#' @export

Hist2ID <- function(history){ # {{{

  # recursive call for vectorization over a matrix of histories stored as rows:
  if (!is.null(dim(history)))
    return(apply(history, 1, Hist2ID))

  # Plain arithmetics of history id: translate an integer from base
  # `nbCaptureEvents` into base 10
  T <- length(history)
  if (T == 1)
    return(as.character(history + 1))
  # use arbitrary long integers not to crush `int`'s ceiling!
  res <- gmp::as.bigz(1)
  for (i in 1:length(history))
    res <- res + history[i] * (nbCaptureEvents ^ (nbCaptureEvents - i))
  return(as.character(res))

} # }}}


#' Get a raw history from its ID
#'
#' Each positive integer can be interpreted as one capture history according to
#' McClintock2010 (equation (3)). It is \code{\link{Hist2ID}} reverse function.
#'
#' @param id ID of a history, either an string or a vector of string. An ID
#' must be a positive integer.
#' @param T total number of capture events to interpret the IDs with
#'
#' @seealso \code{\link{Hist2ID}}
#'
#' @return if \code{id} is one ID, the corresponding raw history. If \code{id}
#' is several IDs, the corresponding matrix of histories.
#'
#' @examples
#' Hist2ID(c(0, 4, 0, 1, 3))
#' set.seed(6)
#' hists <- generateLatentHistories(10)
#' Hist2ID(hists)
#'
#' @export

ID2Hist <- function(id, T){ # {{{

  # Recursive call for vectorization
  if (length(id) > 1) {
    res <- do.call(rbind, plyr::alply(id, 1, ID2Hist, T=T))
    attr(res, 'dimnames') <- NULL
    return(res)
  }

  # Plain arithmetics: translate an integer from base 10 to base
  # `nbCaptureEvents`
  # Thank you gmp for arbitrary long integer manipulation!
  res <- as.character(gmp::as.bigz(id) - 1, b=nbCaptureEvents)
  res <- captureEvents[as.integer(strsplit(res, '')[[1]]) + 1]
  # adjust the size to desired T
  lr <- length(res)
  if (lr < T)
    res <- c(rep(captureEvents['0'], T - lr), res)

  return(res)

} # }}}


#' Get canonical order for a set of histories
#'
#' Canonical order will facilitate managing and analysing the histories.
#' In particular, it will make it easy to derive the polytope associated with a
#' set of histories.
#' Order as follow into 4 blocks:\itemize{
#'    \item{0 - 0-history  : Only 0s}
#'    \item{1 - S-histories: Simultaneous histories (at least one S)}
#'    \item{2 - L-histories: Left-only histories (only 0 and R)}
#'    \item{3 - R-histories: Right-only histories (only 0 and L)}
#'    \item{4 - B-histories: Unobservable histories that generate the previous
#'    two}
#'    }
#' As far as histories within blocks are concerned, they will be ordered in
#' increasing order of their corresponding ID's.
#'
#' @note Be aware that canonical order is not well-suited for polytope analysis
#' of the B-bloc. Instead, polytope order is recommended. See
#' \code{\link{get.Omega.B}}.
#'
#' @param hists a histories matrix (see \code{\link{example.M}})
#'
#' @seealso \code{\link{isObservable}}
#'
#' @return a list of index vectors, each named item corresponding to a block of
#' histories: "0", "S", "L", "R", "B"
#'
#' @examples
#' set.seed(6)
#' hists <- generateLatentHistories(N=20)
#' orderHists(hists)
#'
#' @export

orderHists <- function(hists){ # {{{

  # original order that will be reorganized to get the result, hitchhiking along
  # the sorting process.
  hitchHike <- 1:nrow(hists)

  # Variables `hists` and `hitchHike` will be destroyed by the process

  # First sort by ID's: # {{{
  ids <- Hist2ID(hists)
  # Watch out! sort them as big integers, so they need to be padded with zeroes!
  ml <- max(vapply(ids, nchar, 1))
  ids <- vapply(ids, function(id) paste0(strrep('0', ml - nchar(id)), id), 't')
  o <- order(ids)
  hists     <- hists[o,]
  hitchHike <- hitchHike[o]
  # }}}

  # Separate the blocks: # {{{
  # Successively build boolean masks then remove the targetted histories from
  # the variables to put their index into `blocks`
  i0 <- integer(0)
  blocks <- list(null=i0, S=i0, L=i0, R=i0, B=i0)

  # Macro to execute each time a mask is built: # {{{
  # (basically, put everything the mask selects in the right block and remove it
  # from the destroyed variables)
  macro <- function(mask)
    eval(substitute(expression({
      blocks$MASK <- hitchHike[MASK.mask]
      hitchHike <- hitchHike[!MASK.mask]
      hists <- hists[!MASK.mask,]
    }), list(MASK=gsub(".mask", "", substitute(mask)), MASK.mask=mask)))
  # }}}

  # find the null histories (all zero)
  null.mask <- apply(hists, 1, function(h)
                     all(h == captureEvents["0"]))
  eval(macro(null.mask))

  # Simultaneous histories (contain one S at least)
  S.mask    <- apply(hists, 1, function(h)
                     any(h == captureEvents["S"]))
  eval(macro(S.mask))

  # Observable left-histories: (left-only)
  L.mask <- apply(hists, 1, function(h){
      collapsed <- unique(h)
      return(all(collapsed %in% captureEvents[c("0", "L")]))})
  eval(macro(L.mask))

  # Observable right-histories: (right-only)
  R.mask <- apply(hists, 1, function(h){
      collapsed <- unique(h)
      return(all(collapsed %in% captureEvents[c("0", "R")]))})
  eval(macro(R.mask))

  # All remaining histories are unobservable B-histories
  B.mask <- rep(TRUE, length(hitchHike))
  eval(macro(B.mask))
  # }}}

  # and that's it! ^ ^
  names(blocks)[1] <- '0' # was hard to use with the macro
  return(blocks)
} # }}}


#' Check whether a history is observable
#'
#' Some histories cannot be observed because they involve two different sides of
#' the individual and we have no clue to match them together. They are either
#' the null history or histories containing no S and two sides of the
#' individual. This function performs the check.
#'
#' @inheritParams Hist2ID
#'
#' @seealso \code{\link{orderHists}}
#'
#' @return boolean \code{TRUE} if the history can be observed, or a boolean mask
#' selecting the observable histories only
#'
#' @examples
#' isObservable(captureEvents[c('R', 'S', 'S', 'R', 'S')])
#' isObservable(captureEvents[c('R', 'B', 'B', 'B', '0')])
#' set.seed(12)
#' isObservable(generateLatentHistories(N=20))
#'
#' @export

isObservable <- function(history){ # {{{

  # recursive call for vectorization
  if (!is.null(dim(history)))
    return(apply(history, 1, isObservable))

  # TODO: `don't like this algorithm, improve it if it is too slow

  # Any history is observable if it contains a S event: both sides are matched
  if (any(history == captureEvents['S']))
    return(TRUE)

  # Without an S, any two-sided history is unobservable, that is either..
  if (any(history == captureEvents['B'])) # .. containing a B
    return(FALSE)
  if (all(captureEvents[c('R', 'L')] %in% history)) # .. or containing R *and* L
    return(FALSE)

  # The null history is also unobservable
  if (all(history == captureEvents['0']))
    return(FALSE)

  # But all remaining histories are observable (one-sided histories)
  return(TRUE)

}# }}}


#' Compute observed histories from latent histories
#'
#' Latent, actual histories undergone by the individuals may not be the same as
#' the ones we observe, because we may miss match informations.
#' Here are the rules, inspired from Brett2013 and Simon2010:\itemize{
#'   \item{1: every zero history is removed (unobserved)}
#'   \item{2: every history containing S is kept (two sides matched)}
#'   \item{3: every history containing only R's or only L's is kept
#'   (observable)}
#'   \item{4: every remaining history generates 2 ghosts:\itemize{
#'       \item{a L-ghost where R's become 0's and B's become L's}
#'       \item{a R-ghost where L's become 0's and B's become R's}}}
#'   \item{5: okay, rule 4 stands for rule 3 if we omit zeroes, so 3 rules are
#'   enough}}
#'
#' @param latent one raw history or a matrix of histories as given by
#' \code{\link{generateLatentHistories}}
#'
#' @seealso \code{\link{isObservable}}
#'
#' @return an observation matrix: histories one would observe if \code{latent}
#' were the latent ones.
#'
#' @examples
#' set.seed(2)
#' latent <- generateLatentHistories(N=10)
#' observeHist(latent)
#'
#' @export

observeHist <- function(latent){ # {{{

  # Recursive call to vectorize over histories:
  if (!is.null(dim(latent))) {
    return(do.call(rbind, plyr::alply(latent, 1, observeHist)))
  }

  T <- length(latent)

  # Rule 2: a latent containing a S is kept as it is
  if (captureEvents[c('S')] %in% latent)
    return(matrix(latent, nrow=1))

  # Rules 3-4: any other latent generates 2 ghosts:
  Bs <- latent == captureEvents['B'] # locate the B's
  Ls <- latent == captureEvents['L'] # locate the L's
  Rs <- latent == captureEvents['R'] # locate the R's
  L.ghosts <- replace(latent   , Bs, captureEvents['L']) # L side of B
  L.ghosts <- replace(L.ghosts , Rs, captureEvents['0']) # then no R's
  R.ghosts <- replace(latent   , Bs, captureEvents['R']) # R side of B
  R.ghosts <- replace(R.ghosts , Ls, captureEvents['0']) # then no L's

  # Rule 1: remove the null ghosts:
  if (all(L.ghosts == captureEvents['0']))
    L.ghosts <- matrix(integer(0), ncol=T)
  if (all(R.ghosts == captureEvents['0']))
    R.ghosts <- matrix(integer(0), ncol=T)

  # paste them together in the result:
  res <- matrix(c(L.ghosts, R.ghosts), ncol=T, byrow=TRUE)
  return(res)

} # }}}


#' Convenience printing of histories to the console
#'
#' @param x a matrix of raw histories or a vector of histories IDs
#' @param T number of capture events to interpret the IDs with. If not given,
#' the minimal T is choosen, or the actual length or the raw histories given.
#'
#' @return \code{NULL} only formats the histories \code{x} to the console
#'
#' @seealso \code{\link{ID2Hist}} \code{\link{Hist2ID}}
#'
#' @examples
#' seeHist(1015)
#' seeHist(c(1015, 1, 42, 2092))
#' seeHist(generateLatentHistories(N=20))
#'
#' @export

seeHist <- function(x, T=NULL) { # {{{

  # interpret the argument:
  if (!is.null(dim(x))) {# interpreted as a matrix of histories
    hist <- x
    if (is.null(T))
      T <- ncol(hist)
    else if (T < ncol(hist))
      throw(paste("Cannot interpret these histories with only", T, "events!"))
    id <- Hist2ID(hist)
  }
  else { # interpreted as a vector of IDs
    id <- x
    minimal.T <- max(gmp::sizeinbase(gmp::as.bigz(id), b=nbCaptureEvents))
    if (is.null(T))
      T <- minimal.T
    else if (T < minimal.T)
      throw(paste("Cannot interpret these histories with only", T, "events!"))
    hist <- ID2Hist(id, T)
  }

  # Translate it to actual history events and print it to the screen:
  N <- nrow(hist)
  res <- names(captureEvents)[hist + 1]
  dim(res) <- c(N, T)
  res <- cbind(res, ":", id, "\n")
  cat("", do.call(paste, as.list(c(t(res)))))

} # }}}

