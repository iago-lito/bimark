# Here is R script to preprocess bilateral data
# Feed it with observation matrices as given by `observeHist`

#' Get observed frequency counts from the observation matrix
#'
#' The observation matrix stacks every observed history, with one line per side
#' or individual. The frequency counts sums it up by the number of times each
#' history has been observed.
#'
#' @keywords observation frequency counts
#'
#' @seealso observeHist Hist2ID
#'
#' @param M an observation matrix with raw histories stacked in rows as given by
#' `observeHist`
#'
#' @return a data.frame with two fields: `id` are observed histories IDs, and
#' `F` with the corresponding observed counts
#'
#' @export

observedFrequencies <- function(M){ # {{{

  counts <- as.data.frame(table(Hist2ID(M)))
  counts[,1] <- as.integer(as.character(counts[,1]))
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
#' Omega.R. Consider it as a kind of reversed `observeHist`. Omega.B is sorted
#' with polytope order.
#'
#' Each B-history comes from the combination of one L- and one R-history. Here
#' are the rules for each couple of simultaneous observed events:
#'     Rule 1 : observed 0 and 0 together may come from a latent 0
#'     Rule 2 : observed 0 and L together may come from a latent L
#'     Rule 3 : observed 0 and R together may come from a latent R
#'     Rule 4 : observed L and R together may come from a latent B
#' Polytope order for B-histories is the nested order of canonical R-order
#' within canonical L-order, as follows:
#'      L - 1
#'      L - 2
#'      R - 1
#'      R - 2
#'      R - 3
#'      B - generating L1 and R1
#'      B - generating L1 and R2
#'      B - generating L1 and R3
#'      B - generating L2 and R1
#'      B - generating L2 and R2
#'      B - generating L2 and R3
#'
#' @param Omega.L a matrix of unique, observed L-histories.
#' @param Omega.R a matrix of unique, observed R-histories.
#'
#' @keywords observation latent observable ghosts
#'
#' @seealso observeHist
#'
#' @export

getOmega.B <- function(Omega.L, Omega.R) { # {{{

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

