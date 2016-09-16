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
Omega.L <- matrix(captureEvents[
                  c("0", "0", "0", "0", "L",
                    "0", "0", "L", "0", "L")], ncol=5, byrow=TRUE)
Omega.R <- matrix(captureEvents[
                  c("0", "0", "R", "0", "0",
                    "0", "R", "0", "0", "0",
                    "R", "0", "0", "0", "R")], ncol=5, byrow=TRUE)
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


test_that("observedFrequencies works", {
  expected <- data.frame(id=Hist2ID(Omega.SLR), F=F)
  actual <- observedFrequencies(M)
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

