# Here are R functions to preprocess bilateral data

#' Get frequency counts from raw histories matrix
#'
#' The raw history matrix (like the observation matrix M or the latent histories
#' matrix L) stacks every history, with one line per side or individual. The
#' frequency counts sums it up by the number of times each history has occured.
#' This is useful to get F: observed counts from the observation matrix M, or to
#' get the latent histories counts X from the latent simulated histories.
#'
#' @seealso \code{\link{observeHist}} \code{\link{Hist2ID}}
#'
#' @inheritParams observeHist
#'
#' @return a \code{data.frame} with two fields: \code{id} are observed histories
#' IDs, and \code{F} with the corresponding observed counts
#'
#' @examples
#' # on toy example
#' F <- compute.Frequencies(example.M)
#'
#' # on generated data
#' set.seed(12)
#' latent <- generateLatentHistories(N=20)
#' obs <- observeHist(latent)
#' seeHist(obs)
#' F <- compute.Frequencies(obs)
#'
#' # This also works directly with latent histories
#' X <- compute.Frequencies(latent)
#'
#' @export

compute.Frequencies <- function(M){ # {{{

  counts <- as.data.frame(table(Hist2ID(M)), stringsAsFactors=FALSE)
  counts[[1]] <- as.character(counts[[1]])
  names(counts) <- c('id', 'F')
  # Note that, here, `F` differs from Bonners' "f" since it only contains "f"
  # positive elements. Much ligther.

  return(counts)

} # }}}


#' Build B-histories from L- and R-histories
#'
#' See porcelain function \code{\link{get.Omega.B}}.
#'
#' @param Omega.L a matrix of unique, observed L-histories.
#' @param Omega.R a matrix of unique, observed R-histories.
#'
#' @seealso \code{\link{observeHist}}
#'
#' @examples
#' # raw observation data
#' obs <- compute.Frequencies(example.M)
#' Omega.LRS <- ID2Hist(obs$id, T=ncol(example.M))
#' # sort in canonical order
#' o <- orderHists(Omega.LRS)
#' Omega.S <- Omega.LRS[o$S,]
#' Omega.L <- Omega.LRS[o$L,]
#' Omega.R <- Omega.LRS[o$R,]
#' # get unobservable block in polytope order
#' Omega.B <- compute.Omega.B(Omega.L, Omega.R)
#' # and this is canonical Omega:
#' Omega <- rbind(Omega.S, Omega.L, Omega.R, Omega.B)
#' seeHist(Omega)
#'
#' @export

compute.Omega.B <- function(Omega.L, Omega.R) { # {{{

  # basic informations:
  T <- ncol(Omega.L)
  LL <- nrow(Omega.L)
  LR <- nrow(Omega.R)

  # prepare the ghost matrix
  LB <- LL * LR
  Omega.B <- matrix(0L, nrow=LB, ncol=T)
  if (LB == 0)
    return(Omega.B)

  # fill it!
  b <- 0                   # index visiting Omega.B
  for (l in 1:LL) {
    left <- Omega.L[l,]    # select a left history
    for (r in 1:LR){
      right <- Omega.R[r,] # select a right history
      b <- b + 1           # then fill the `B.hist` bloc accordingly

      Omega.B[b,] <- apply(rbind(left,right), 2, function(events){
        # `events` contains two simultaneous events in two different histories
        zeroes <- events == captureEvents["0"]      # mask for zeroes
        if (all(zeroes)) return(captureEvents["0"]) # Rule 1
        if (any(zeroes)) return(events[!zeroes])    # Rules 2 and 3
        return(captureEvents["B"])})                # Rule 4

    }
  }

  # this is it!
  return(Omega.B)

}# }}}

#' Generate A matrix from LL and LR
#'
#' See porcelain function \code{\link{get.A}}.
#'
#' @param LL int the number of observed L-histories (not the sum of their
#' frequencies but the number of *different* types of L-histories observed.
#' @param LR int the number of observed R-histories
#'
#' @seealso \code{\link{compute.Omega.B}} \code{\link{compute.B}}
#'
#' @examples
#' compute.A(LL=2, LR=3)
#'
#' @export

compute.A <- function(LL, LR) { # {{{

  # Basic needed values:
  LU <- LL + LR
  LB <- LL * LR
  LM <- LU + LB

  # Generate the matrix!
  seqLU <- seq.int(1L, LU)
  if (LB == 0) { # degenerated case
    if(LU == 0) # only S-hists
      A <- Matrix::sparseMatrix(i=integer(0), j=integer(0), dims=c(0, 0))
    else
      A <- Matrix::sparseMatrix(seqLU, seqLU, x=1L)
  } else {
    seqLB <- seq.int(1L, LB)
    seqLL <- seq.int(1L, LL)
    seqLR <- seq.int(1L, LR)
    each  <- rep.int(LR, LL)
    # rows where there are 1's
    rows <- c(seqLU, LU + rep.int(seqLB, 2L))
    # columns where there are 1's
    cols  <- c(seqLU, rep.int(seqLL, each), LL + rep.int(seqLR, LL))
    # final matrix:
    A <- Matrix::sparseMatrix(rows, cols, x=1L)
  }

  # that's it!
  A <- methods::as(A, 'dgCMatrix') # or R'll choose anything it likes
  return(A)

} # }}}

#' Generate B matrix, nullspace of t(A)
#'
#' See porcelain function \code{\link{get.B}}.
#'
#' @inheritParams compute.A
#'
#' @seealso \code{\link{compute.A}}
#'
#' @examples
#' compute.B(LL=2, LR=3)
#'
#' @export

compute.B <- function(LL, LR) { # {{{

  # Basic needed values:
  LU <- LL + LR
  LB <- LL * LR
  LM <- LU + LB

  if (LB == 0){ # degenerated case
    if(LU == 0) # only S-hists
      B <- Matrix::sparseMatrix(i=integer(0), j=integer(0), dims=c(0, 0))
    else
      B <- Matrix::sparseMatrix(i=integer(0), j=integer(0), dims=c(LM, 0))
  } else {
    seqLL <- seq.int(1L, LL)
    seqLR <- seq.int(1L, LR)
    seqLB <- seq.int(1L, LB)
    each  <- rep.int(LR, LL)
    rows <- rep.int(seqLB, 3) # those target the -1's
    cols <- c(LU + seqLB, rep.int(seqLL, each), LL + rep.int(seqLR, LL))
    B <- Matrix::t(Matrix::sparseMatrix(rows, cols, x=rep.int(c(1L, -1L),
                                                              c(LB, 2L * LB))))
  }

  # that's it!
  B <- methods::as(B, 'dgCMatrix') # or R'll choose anything it likes
  return(B)

} # }}}

