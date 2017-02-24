context("Preprocessing procedures")

# playground raw observation matrix from Bonnici M2 report:
T <- 5
M <- matrix(capture.events[
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
Omega.S <- matrix(capture.events[
                  c("L", "0", "B", "S", "0",
                    "S", "0", "L", "0", "0")], ncol=5, byrow=TRUE)
LS <- nrow(Omega.S)
Omega.L <- matrix(capture.events[
                  c("0", "0", "0", "0", "L",
                    "0", "0", "L", "0", "L")], ncol=5, byrow=TRUE)
LL <- nrow(Omega.L)
Omega.R <- matrix(capture.events[
                  c("0", "0", "R", "0", "0",
                    "0", "R", "0", "0", "0",
                    "R", "0", "0", "0", "R")], ncol=5, byrow=TRUE)
LR <- nrow(Omega.R)
Omega.SLR <- rbind(Omega.S, Omega.L, Omega.R)
# Associated frequency counts
F <- c(1, 2, 5, 3, 4, 6, 1)
# Unobservable part of Omega (potential latent histories) in polytope order
Omega.B <- matrix(capture.events[
                  c("0", "0", "R", "0", "L",
                    "0", "R", "0", "0", "L",
                    "R", "0", "0", "0", "B",
                    "0", "0", "B", "0", "L",
                    "0", "R", "L", "0", "L",
                    "R", "0", "L", "0", "B")], ncol=5, byrow=TRUE)
Omega <- rbind(Omega.SLR, Omega.B)
# polytope matrices:
A <- c(1, 0, 0, 0, 0,
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
A <- Matrix::Matrix(matrix(A, ncol=LL + LR, byrow=TRUE))
B <- c(-1, -1, -1,  0,  0,  0,
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
B <- Matrix::Matrix(matrix(B, ncol=LL * LR, byrow=TRUE))


test_that("ComputeFrequencies works", {
  expected <- data.frame(id=Hist2ID(Omega.SLR), F=F, stringsAsFactors=FALSE)
  actual <- ComputeFrequencies(M)
  # canonic order to compare
  oa <- unlist(OrderHists(ID2Hist(actual$id, T)))
  oe <- unlist(OrderHists(ID2Hist(expected$id, T)))
  expect_equal(actual[oa,], expected[oe,], check.attributes=FALSE)
  })

test_that("ComputeOmegaB works", {
  o <- OrderHists(Omega.SLR)
  Omega.L <- Omega.SLR[o$L,]
  Omega.R <- Omega.SLR[o$R,]
  actual <- ComputeOmegaB(Omega.L, Omega.R)
  expected <- Omega.B # we are expecting polytope order
  expect_equal(actual, expected)
  # test degenerated cases
  expected <- matrix(integer(0), nrow=0, ncol=ncol(Omega.SLR))
  actual <- ComputeOmegaB(Omega.L[integer(0),], Omega.R) # empty L
  expect_equal(actual, expected)
  actual <- ComputeOmegaB(Omega.L, Omega.R[integer(0),]) # empty R
  expect_equal(actual, expected)
  actual <- ComputeOmegaB(Omega.L[integer(1),], Omega.R[integer(0),]) # both
  expect_equal(actual, expected)
  })

test_that("ComputeA works", {

  # regular
  actual <- ComputeA(LL, LR)
  expected <- A
  expect_equal(actual, expected)

  # degenerated:
  expected <- c(1, 0,
                0, 1)
  expected <- Matrix::Matrix(matrix(expected, ncol=LL + 0, byrow=TRUE),
                             sparse=TRUE)
  expected <- methods::as(expected, 'dgCMatrix')
  # left AND right:
  actual <- ComputeA(LL, 0)
  expect_equal(actual, expected)
  actual <- ComputeA(0, LL)
  expect_equal(actual, expected)

  # super-degenerated:
  actual <- ComputeA(0, 0)
  expected <- Matrix::Matrix(matrix(integer(0), ncol=0, byrow=TRUE),
                             sparse=TRUE)
  expected <- methods::as(expected, 'dgCMatrix')
  expect_equal(actual, expected)

  })

test_that("ComputeB works", {

  # regular
  actual <- ComputeB(LL, LR)
  expected <- B
  expect_equal(actual, expected)

  # degenerated:
  expected <- Matrix::Matrix(matrix(integer(0), ncol=0, nrow=LL, byrow=TRUE),
                             sparse=TRUE)
  expected <- methods::as(expected, 'dgCMatrix')
  # left AND right:
  actual <- ComputeB(LL, 0)
  expect_equal(actual, expected)
  actual <- ComputeB(0, LL)
  expect_equal(actual, expected)

  # super-degenerated:
  actual <- ComputeB(0, 0)
  expected <- Matrix::Matrix(matrix(integer(0), ncol=0, byrow=TRUE),
                             sparse=TRUE)
  expected <- methods::as(expected, 'dgCMatrix')
  expect_equal(actual, expected)

  })

test_that("get Omega from bimarkModel object", {
  model <- BimarkObservationModel(M)
  expect_equal(GetOmega(model), Omega)
  expect_equal(GetOmegaS(model), Omega.S)
  expect_equal(GetOmega.L(model), Omega.L)
  expect_equal(GetOmega.R(model), Omega.R)
  expect_equal(GetOmega.B(model), Omega.B)
  expect_equal(GetOmegaSLR(model), Omega.SLR)
  expect_equal(GetOmega.LR(model), rbind(Omega.L, Omega.R))
  })

test_that("get matrices from bimarkModel object", {
  model <- BimarkObservationModel(M)
  actual <- get.A(model)
  expected <- A
  expect_equal(actual, expected)
  actual <- get.B(model)
  expected <- B
  expect_equal(actual, expected)
  })

