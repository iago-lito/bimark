# Here is R script to simulate bilateral data:

# Defining events of capture: # {{{

#' Capture events and their integer representations
#'
#' @description
#'
#'  \tabular{ccc}{
#'    0 : \tab no capture                     \tab : 0\cr
#'    L : \tab right capture                  \tab : 1\cr
#'    R : \tab left capture                   \tab : 2\cr
#'    B : \tab both sides capture, no match   \tab : 3\cr
#'    S : \tab simultaneous capture and match \tab : 4\cr
#'  }
#'
#' @seealso nbCaptureEvents
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
#' This is just defined as length(captureEvents)
#'
#' @seealso captureEvents
#'
#' @export

nbCaptureEvents <- length(captureEvents)

# }}}

# custom error thrower:
throw <- function(message){
  call <- sys.call(-1)
  callingFunctionName <- as.character(call[1])
  stop(paste0("bimark: ", callingFunctionName, ": ", message), call.=FALSE)
}

#' Generate latent histories
#'
#' Simulate a virtual population where each individual undergoes a particular,
#' latent capture history.
#'
#' @usage generateLatentHistories(N, T, P, delta)
#'
#' @param N integer number of simulated individual
#' @param T integer number of capture occasions
#' @param P N-vector: capture probabilities on each occasion (one value for
#' constant probability constant)
#' @param delta nbCaptureEvents-vector: capture probabilities for each capture
#' event (or positive weights whose sum will be normalized to 1)
#'
#' @keywords history latent simulate generate
#'
#' @examples
#' generateLatentHistories(N=20, T=5, P=rep(.1, 5), delta=c(0, 4, 4, 2, 3))
#'
#' @return a matrix of capture events: one row per individual, one column per
#' capture occasion: each row is an individual capture history.
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
#' (equation (3)). It is ID2Hist reverse function.
#'
#' @param history for T capture events, either a vector of events of length T or
#' a matrix of n histories dim(n, T)
#'
#' @keywords history raw id
#'
#' @examples
#' Hist2ID(c(0, 4, 0, 1, 3))
#' Hist2ID(matrix(c(0, 4, 0, 1, 3,
#'                  1, 1, 0, 2, 0,
#'                  0, 3, 4, 0, 0,
#'                  0, 0, 0, 4, 4), 4, 5, byrow=TRUE))
#' Hist2ID(generateLatentHistories(4, 5))
#'
#' @seealso ID2Hist generateLatentHistories
#'
#' @return if `history` is one history, its id. If `history` is several history,
#' their corresponding vector of ids.
#'
#' @export

Hist2ID <- function(history){ # {{{

  # recursive call for vectorization over a matrix of histories stored as rows:
  if (!is.null(dim(history)))
    return(apply(history, 1, Hist2ID))

  # Plain arithmetics of history id: translate an integer from base
  # `nbCaptureEvents` into base 10
  T <- length(history)
  return(as.integer(1 + sum(history * nbCaptureEvents^((T-1):0))))

} # }}}


#' Get a raw history from its ID
#'
#' Each positive integer can be interpreted as one capture history according to
#' McClintock2010 (equation (3)). It is Hist2ID reverse function.
#'
#' @param id ID of a history, either an integer or a vector of integers. An ID
#' must be positive.
#' @param T total number of capture events to interpret the IDs with
#'
#' @keywords history raw id
#'
#' @examples
#' Hist2ID(c(0, 4, 0, 1, 3))
#' Hist2ID(matrix(c(0, 4, 0, 1, 3,
#'                  1, 1, 0, 2, 0,
#'                  0, 3, 4, 0, 0,
#'                  0, 0, 0, 4, 4), 4, 5, byrow=TRUE))
#' Hist2ID(generateLatentHistories(4, 5))
#'
#' @seealso Hist2ID
#'
#' @return if `id` is one IDs, the corresponding raw history. If `id` is several
#' IDs, the corresponding `dim(length(id), T)` matrix of histories.
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
  res <- rep(NA, T)
  id  <- id - 1
  divisor <- nbCaptureEvents^(T - 1:T)
  for (i in 1:T){
    res[i]  <- id %/% divisor[i]
    id      <- id - res[i] * divisor[i]
  }

  return(res)

} # }}}


#' Get canonical order of a set of histories
#'
#' Canonical order will facilitate managing and analysing the histories.
#' In particular, it will make it easy to derive the polytope associated with a
#' set of histories.
#' Order as follow into 4 blocks:
#'    0 - 0-history  : Only 0s
#'    1 - S-histories: Simultaneous histories (at least one S)
#'    2 - L-histories: Left-only histories (only 0 and R)
#'    3 - R-histories: Right-only histories (only 0 and L)
#'    4 - B-histories: Unobservable histories that generate the previous two
#' As far as histories within blocks are concerned, they will be ordered in
#' increasing order of their corresponding ID's.
#'
#' @note Be aware that canonical order is not well-suited for polytope analysis
#' of the B-bloc. Instead, polytope order is recommended. See `getOmega.B`.
#'
#' @param hists N raw histories stored in a dim(N, T) matrix
#'
#' @keywords histories order canonical analysis groups observable unobservable
#' blocks
#'
#' @seealso isObservable
#'
#' @return a list of index vectors, each named item corresponding to a block of
#' histories: "0", "S", "L", "R", "B"
#'
#' @export

orderHists <- function(hists){ # {{{

  # original order that will be reorganized to get the result, hitchhiking along
  # the sorting process.
  hitchHike <- 1:nrow(hists)

  # Variables `hists` and `hitchHike` will be destroyed by the process

  # First sort by ID's: # {{{
  ids <- Hist2ID(hists)
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
      # I hate you R will you STOP changing my variables types with no warning?!
      if (is.null(dim(hists))) dim(hists) <- c(1, length(hists))
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
#' @param history a raw history or a N,T-matrix of N histories as given by
#' `generateLatentHistories`
#'
#' @keywords histories observable unobservable test
#'
#' @seealso orderHists
#'
#' @return boolean TRUE is the history can be observed, or a mask selecting the
#' observable histories
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
#' Here are the rules, inspired from Brett2013 and Simon2010:
#'   1: every zero history is removed (unobserved)
#'   2: every history containing S is kept (two sides matched)
#'   3: every history containing only R's or only L's is kept (observable)
#'   4: every remaining history generates 2 ghosts:
#'       - a L-ghost where R's become 0's and B's become L's
#'       - a R-ghost where L's become 0's and B's become R's
#'   5: okay, rule 4 stands for rule 3 if we omit zeroes, so 3 rules are enough
#'
#' @param latent a raw history or a N,T-matrix of N histories as given by
#' `generateLatentHistories`
#'
#' @keywords histories latent observed observation observable unobservable
#' ghosts
#'
#' @seealso isObservable
#'
#' @return an observation matrix: histories one would observe if `latent` were
#' the latent ones.
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
#' @param x N raw histories in a N,T-matrix or a vector of histories IDs.
#' @param T number of capture events to interpret the IDs with. If not given,
#' the minimal T is choosen, or the actual length or the raw histories given.
#'
#' @keywords histories id row visualisation convenience
#'
#' @return NULL only formats the histories `x` to the console
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
    minimal.T <- max(1, ceiling(log(max(id), nbCaptureEvents)))
    if (is.null(T))
      T <- minimal.T
    else if (T < minimal.T)
      throw(paste("Cannot interpret these histories with only", T, "events!"))
    hist <- ID2Hist(id, T)
    if (is.null(dim(hist))) hist <- matrix(hist, nrow=1) # `hate you R!
  }

  # Translate it to actual history events and print it to the screen:
  N <- nrow(hist)
  res <- names(captureEvents)[hist + 1]
  dim(res) <- c(N, T)
  res <- cbind(res, ":", id, "\n")
  cat("", do.call(paste, as.list(c(t(res)))))

} # }}}

