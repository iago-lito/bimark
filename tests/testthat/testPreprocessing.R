context("Preprocessing procedures")

# playground raw observation matrix from Bonnici M2 report:
T <- 5
M <- matrix(captureEvents[
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
# Observable part of Omega in canonical order
Omega.S <- matrix(captureEvents[
                  c("L", "0", "B", "S", "0",
                    "S", "0", "L", "0", "0")], ncol=5, byrow=TRUE)
LS <- nrow(Omega.S)
Omega.L <- matrix(captureEvents[
                  c("0", "0", "0", "0", "L",
                    "0", "0", "L", "0", "L")], ncol=5, byrow=TRUE)
LL <- nrow(Omega.L)
Omega.R <- matrix(captureEvents[
                  c("0", "0", "R", "0", "0",
                    "0", "R", "0", "0", "0",
                    "R", "0", "0", "0", "R")], ncol=5, byrow=TRUE)
LR <- nrow(Omega.R)
Omega.SLR <- rbind(Omega.S, Omega.L, Omega.R)
# Associated frequency counts
F <- c(1, 2, 5, 3, 4, 6, 1)
# Unobservable part of Omega (potential latent histories) in polytope order
Omega.B <- matrix(captureEvents[
                 c("0", "0", "R", "0", "L",
                   "0", "R", "0", "0", "L",
                   "R", "0", "0", "0", "B",
                   "0", "0", "B", "0", "L",
                   "0", "R", "L", "0", "L",
                   "R", "0", "L", "0", "B")], ncol=5, byrow=TRUE)


test_that("getFrequencies works", {
  expected <- data.frame(id=Hist2ID(Omega.SLR), F=F)
  actual <- getFrequencies(M)
  # canonic order to compare
  oa <- unlist(orderHists(ID2Hist(actual$id, T)))
  oe <- unlist(orderHists(ID2Hist(expected$id, T)))
  expect_equal(actual[oa,], expected[oe,], check.attributes=FALSE)
  })

test_that("getOmega.B works", {
  o <- orderHists(Omega.SLR)
  Omega.L <- Omega.SLR[o$L,]
  Omega.R <- Omega.SLR[o$R,]
  actual <- getOmega.B(Omega.L, Omega.R)
  expected <- Omega.B # we are expecting polytope order
  expect_equal(actual, expected)
  })

test_that("getA works", {

  # regular
  actual <- getA(LL, LR)
  expected <- c(1, 0, 0, 0, 0,
                0, 1, 0, 0, 0,
                0, 0, 1, 0, 0,
                0, 0, 0, 1, 0,
                0, 0, 0, 0, 1,
                1, 0, 1, 0, 0,
                1, 0, 0, 1, 0,
                1, 0, 0, 0, 1,
                0, 1, 1, 0, 0,
                0, 1, 0, 1, 0,
                0, 1, 0, 0, 1)
  expected <- Matrix::Matrix(matrix(expected, ncol=LL + LR, byrow=TRUE))
  expect_equal(actual, expected)

  # degenerated:
  expected <- c(1, 0,
                0, 1)
  expected <- Matrix::Matrix(matrix(expected, ncol=LL + 0, byrow=TRUE),
                             sparse=TRUE)
  expected <- methods::as(expected, 'dgCMatrix')
  # left AND right:
  actual <- getA(LL, 0)
  expect_equal(actual, expected)
  actual <- getA(0, LL)
  expect_equal(actual, expected)

  # super-degenerated:
  actual <- getA(0, 0)
  expected <- Matrix::Matrix(matrix(integer(0), ncol=0, byrow=TRUE),
                             sparse=TRUE)
  expected <- methods::as(expected, 'dgCMatrix')
  expect_equal(actual, expected)

  })

test_that("getB works", {

  # regular
  actual <- getB(LL, LR)
  expected <- c(-1, -1, -1,  0,  0,  0,
                 0,  0,  0, -1, -1, -1,
                -1,  0,  0, -1,  0,  0,
                 0, -1,  0,  0, -1,  0,
                 0,  0, -1,  0,  0, -1,
                 1,  0,  0,  0,  0,  0,
                 0,  1,  0,  0,  0,  0,
                 0,  0,  1,  0,  0,  0,
                 0,  0,  0,  1,  0,  0,
                 0,  0,  0,  0,  1,  0,
                 0,  0,  0,  0,  0,  1)
  expected <- Matrix::Matrix(matrix(expected, ncol=LL * LR, byrow=TRUE))
  expect_equal(actual, expected)

  # degenerated:
  expected <- Matrix::Matrix(matrix(integer(0), ncol=0, nrow=LL, byrow=TRUE),
                             sparse=TRUE)
  expected <- methods::as(expected, 'dgCMatrix')
  # left AND right:
  actual <- getB(LL, 0)
  expect_equal(actual, expected)
  actual <- getB(0, LL)
  expect_equal(actual, expected)

  # super-degenerated:
  actual <- getB(0, 0)
  expected <- Matrix::Matrix(matrix(integer(0), ncol=0, byrow=TRUE),
                             sparse=TRUE)
  expected <- methods::as(expected, 'dgCMatrix')
  expect_equal(actual, expected)

  })

