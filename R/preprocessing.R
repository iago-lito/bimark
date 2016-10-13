# Here is R script to preprocess bilateral data
# Feed it with observation matrices as given by `observeHist`

#' Get frequency counts from raws histories matrix
#'
#' The raw history matrix (like the observation matrix M or the latent histories
#' matrix) stacks every history, with one line per side or individual. The
#' frequency counts sums it up by the number of times each history has occured.
#' This is useful to get F observed counts from the observation matrix M, or to
#' get the latent histories counts X from the latent simulated histories.
#'
#' @seealso observeHist Hist2ID
#'
#' @inheritParams observeHist
#'
#' @return a data.frame with two fields: \code{id} are observed histories IDs,
#' and \code{F} with the corresponding observed counts
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
#' This function computes all possible latent, unobservable histories (aka
#' B-histories) which may have generated the given letf- and right-histories (L-
#' and R-histories). In a nutshell, it generates Omega.B from Omega.L and
#' Omega.R. Consider it as a kind of reversed \code{\link{observeHist}}. Omega.B
#' is sorted with polytope order.
#'
#' Each B-history comes from the combination of one L- and one R-history. Here
#' are the rules for each couple of simultaneous observed events:
#'     \itemize{
#'     \item{Rule 1: }{observed 0 and 0 together may come from a latent 0}
#'     \item{Rule 2: }{observed 0 and L together may come from a latent L}
#'     \item{Rule 3: }{observed 0 and R together may come from a latent R}
#'     \item{Rule 4: }{observed L and R together may come from a latent B}
#'     }
#'
#' Polytope order for B-histories is the nested order of canonical R-order
#' within canonical L-order, as follows: \preformatted{
#'       |          |  S - generating L2 and R2
#'       | Omega.S  |  S - generating L2 and R3
#'       |         |   L - 1
#'       | Omega.L |   L - 2
#'       |          |  R - 1
#'       | Omega.R  |  R - 2
#' Omega |          |  R - 3
#'       |         |   B - generating L1 and R1
#'       |         |   B - generating L1 and R2
#'       | Omega.B |   B - generating L1 and R3
#'       |         |   B - generating L2 and R1
#'       |         |   B - generating L2 and R2
#'       |         |   B - generating L2 and R3}
#'
#' @param Omega.L a matrix of unique, observed L-histories.
#' @param Omega.R a matrix of unique, observed R-histories.
#'
#' @seealso observeHist
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
#' This version of the A matrix, as in Bonner2013 only takes into account the
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
#' @param LL int the number of observed L-histories (not the sum of their
#' frequencies but the number of *different* types of L-histories observed.
#' @param LR int the number of observed R-histories
#'
#' @seealso compute.Omega.B compute.B
#'
#' @examples
#' compute.A(LL=2, LR=3)
#'
#' @return the observation matrix A as a \code{\link[Matrix]{sparseMatrix}} : a
#' LM x LU matrix with ones and zeroes, each row corresponding to a latent
#' history sorted in polytope order, each column corresponding to an observable
#' history sorded in polytope (or "canonical") order.
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
# (ugly-but-efficient as well: the pattern of B is easy to grasp and only
# depends on LL and LR)
#'
#' @seealso \code{\link{compute.A}}
#'
#' @examples
#' compute.B(LL=2, LR=3)
#'
#' @inheritParams compute.A
#'
#' @return the polytope matrix B as a \code{\link[Matrix]{sparseMatrix}} : a LM
#' x LB matrix with ones, zeroes and minus ones, each row corresponding to a
#' latent history sorted in polytope order, each column corresponding to an
#' unobservable B-history sorded in polytope order.
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

