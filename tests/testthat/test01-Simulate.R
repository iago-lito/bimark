context("Simulation procedures")

# simulation parameters:
T <- 5                         # number of capture occasions
N <- 10                        # true population size
P <- rep(.1, T)                # probabilities of capture on each occasion
delta <- c(0., 4., 4., 2., 3.) # probabilities of capture events

test_that("generateLatentHistories works", {
  set.seed(4)
  expected <- structure( # got from `dump` once it worked :P
                c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 2L, 4L, 0L, 0L, 0L, 0L, 3L, 4L,
                  0L, 0L, 0L, 0L, 4L, 0L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                  0L, 1L, 1L, 3L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                  0L, 0L, 0L, 4L, 0L), .Dim = c(10L, 5L))
  actual <- generateLatentHistories(N, T, P, delta)
  expect_equal(expected, actual)
  })

test_that("generateLatentHistories throws errors", {
  expect_error(generateLatentHistories(N, T, P[-1], delta),
               "Inconsistent number of capture probabilities.")
  expect_error(generateLatentHistories(N, T, P - 1., delta),
               "Inconsistent capture probabilities.")
  expect_error(generateLatentHistories(N, T, P, delta[-1]),
               "Inconsistent number of capture event probabilities.")
  expect_error(generateLatentHistories(N, T, P, - delta),
               "No negative weights for capture event probabilities.")
  })

test_that("Hist2ID works with single histories", {
  expect_equal(Hist2ID(c(0, 0, 0, 0, 0)), "1")
  expect_equal(Hist2ID(c(0, 4, 0, 1, 3)), "509")
  })

test_that("Hist2ID works with matrices", {
  actual <- Hist2ID(matrix(c(0, 4, 0, 1, 3,
                             1, 1, 0, 2, 0,
                             0, 3, 4, 0, 0,
                             0, 0, 0, 4, 4), 4, 5, byrow=TRUE))
  expected <- c("509", "761", "476", "25")
  expect_equal(expected, actual)
  })

test_that("ID2Hist works with single IDs", {
  actual <- ID2Hist("1", T)
  expected <- captureEvents[c('0', '0', '0', '0', '0')]
  expect_equal(actual, expected)
  actual <- ID2Hist("509", T)
  expected <- captureEvents[c('0', 'S', '0', 'L', 'B')]
  expect_equal(actual, expected)
  })

test_that("ID2Hist works with several IDs", {
  expected <- matrix(c(0, 4, 0, 1, 3,
                       1, 1, 0, 2, 0,
                       0, 3, 4, 0, 0,
                       0, 0, 0, 4, 4), 4, 5, byrow=TRUE)
  actual <- ID2Hist(c(509, 761, 476, 25), T)
  expect_equal(expected, actual)
  })

test_that("Hist2ID and ID2Hist are each other's reversion" , {
  set.seed(12)
  N <- 1e3
  # Testing hist -> ID -> hist
  # aim for a balanced set of histories:
  P <- (nbCaptureEvents - 1) / nbCaptureEvents
  delta <- c(0., 1., 1. ,1. ,1.)
  histories <- generateLatentHistories(N, T, 4./5., c(0., 1., 1., 1., 1.))
  expect_equal(histories, ID2Hist(Hist2ID(histories), T))
  # Testing ID -> hist -> ID
  # aim for a balanced set of ids:
  # (hypothesis: R integers are long enough to handle this particular test ;)
  min <- 1
  max <- as.integer(Hist2ID(rep(max(captureEvents), T)))
  ids <- as.character(floor(runif(N) * max) + min)
  expect_equal(ids, Hist2ID(ID2Hist(ids, T)))
  })

test_that("orderHists works", {
  # a tricky set of histories to order: type id
  hists <- matrix(c(1, 0, 0, 1, 0, #  1  L   631
                    1, 0, 2, 3, 0, #  2  B   691
                    0, 2, 2, 0, 2, #  3  R   303
                    0, 1, 0, 0, 1, #  4  L   127
                    0, 1, 0, 0, 1, #  5  L   127 duplicate: order conserved
                    0, 0, 0, 0, 0, #  6  0     1
                    0, 0, 4, 0, 0, #  7  S   101
                    4, 0, 2, 3, 0, #  8  S  2566
                    0, 1, 2, 2, 0, #  9  B   186
                    0, 0, 2, 2, 0, # 10  R    61
                    1, 1, 0, 0, 0  # 11  L   751
                    ), 11, 5, byrow=TRUE)
  actual <- orderHists(hists)
  expected <- list('0'= c(6),
                    S = c(7, 8),
                    L = c(4, 5, 1, 11), # order conserved
                    R = c(10, 3),
                    B = c(9, 2))
  expect_equal(actual, expected)
  })

test_that("isObservable works with single histories", {
  # null history
  expect_equal(isObservable(captureEvents[c('0', '0', '0', '0', '0')]), FALSE)
  # L-histories
  expect_equal(isObservable(captureEvents[c('0', '0', '0', '0', 'L')]), TRUE)
  expect_equal(isObservable(captureEvents[c('0', 'L', '0', '0', 'L')]), TRUE)
  # R-histories
  expect_equal(isObservable(captureEvents[c('R', '0', '0', '0', 'R')]), TRUE)
  expect_equal(isObservable(captureEvents[c('R', 'R', 'R', 'R', 'R')]), TRUE)
  # B-histories
  expect_equal(isObservable(captureEvents[c('R', '0', '0', 'L', '0')]), FALSE)
  expect_equal(isObservable(captureEvents[c('R', 'B', '0', 'L', '0')]), FALSE)
  expect_equal(isObservable(captureEvents[c('R', 'B', 'B', 'B', '0')]), FALSE)
  # S-histories
  expect_equal(isObservable(captureEvents[c('R', 'S', 'B', 'B', '0')]), TRUE)
  expect_equal(isObservable(captureEvents[c('R', 'S', '0', 'L', '0')]), TRUE)
  expect_equal(isObservable(captureEvents[c('R', 'S', 'S', 'R', 'S')]), TRUE)
  })

test_that("isObservable works with matrices", {
  set.seed(12)
  N <- 1e3
  histories <- generateLatentHistories(N, T, P, delta)
  order <- orderHists(histories)
  result <- isObservable(histories)
  # all null histories should be unobservable
  expect_equal(all(result[order$'0']), FALSE)
  # all left- and right-histories should be observable
  expect_equal(all(result[order$L]), TRUE)
  expect_equal(all(result[order$R]), TRUE)
  # all B-histories should be unobservable
  expect_equal(all(result[order$B]), FALSE)
  # all simultaneous histories should be observable
  expect_equal(all(result[order$S]), TRUE)
  })

test_that("observeHist works", {
  # .. on null histories:                                            # type
  actual.1 <- observeHist(captureEvents[c('0', '0', '0', '0', '0')]) # 0
  expected.1 <- matrix(integer(0), 0 , 5) # no observation
  expect_equal(actual.1, expected.1)
  # .. on observable histories:
  actual.2 <- observeHist(captureEvents[c('0', 'L', '0', 'L', '0')]) # L
  actual.3 <- observeHist(captureEvents[c('0', '0', 'R', 'R', 'R')]) # R
  actual.4 <- observeHist(captureEvents[c('0', 'L', 'S', 'R', '0')]) # S
  expect_equal(actual.2, observeHist(actual.2))  # unchanged
  expect_equal(actual.3, observeHist(actual.3))  # unchanged
  expect_equal(actual.4, observeHist(actual.4))  # unchanged
  # .. on ghost-generating histories
  actual.5 <- observeHist(captureEvents[c('0', 'L', 'R', 'L', '0')]) # B
  expected.5 <-    matrix(captureEvents[c('0', 'L', '0', 'L', '0',   # L-ghost
                                          '0', '0', 'R', '0', '0')], # R-ghost
                          nrow=2, byrow=TRUE)
  expect_equal(expected.5, actual.5)
  actual.6 <- observeHist(captureEvents[c('R', '0', 'B', 'B', '0')]) # B
  expected.6 <-    matrix(captureEvents[c('0', '0', 'L', 'L', '0',   # L-ghost
                                          'R', '0', 'R', 'R', '0')], # R-ghost
                          nrow=2, byrow=TRUE)
  expect_equal(expected.6, actual.6)
  # and on matrices!
  hists <- rbind(actual.1, actual.2, actual.3, actual.4, actual.5, actual.6)
  expected <- rbind(expected.1, actual.2, actual.3, actual.4,
                    expected.5, expected.6)
  actual <- observeHist(hists)
  expect_equal(actual, expected)
  })

test_that("seeHist does not throw errors", {
  set.seed(2)
  # with ids:
  expected <- paste0(" 0 0 : 1 \n",
                     " 0 L : 2 \n",
                     " 0 R : 3 \n",
                     " 0 B : 4 \n",
                     " 0 S : 5 \n",
                     " L 0 : 6 \n",
                     " L L : 7 \n",
                     " L R : 8 \n",
                     " L B : 9 \n",
                     " L S : 10 ")
  # given as integers OR characters
  expect_output(seeHist(1:10), expected)
  expect_output(seeHist(as.character(1:10)), expected)
  # with raw histories:
  hists <- generateLatentHistories(N, T, P, delta)
  expected <- paste0(" 0 0 0 0 S : 5 \n",
                     " L 0 0 0 0 : 626 \n",
                     " 0 0 0 0 0 : 1 \n",
                     " 0 S 0 0 0 : 501 \n",
                     " 0 0 0 0 0 : 1 \n",
                     " 0 0 0 L 0 : 6 \n",
                     " 0 0 0 0 0 : 1 \n",
                     " 0 0 0 0 0 : 1 \n",
                     " L 0 0 0 S : 630 \n",
                     " 0 B 0 0 0 : 376 ")
  expect_output(seeHist(hists), expected)
  })

test_that("seeHist does throw errors", {
  set.seed(12)
  hists <- generateLatentHistories(N, T, P, delta)
  expect_error(seeHist(150, T=2), "Cannot interpret")
  expect_error(seeHist(hists, T=2), "Cannot interpret")
  })


