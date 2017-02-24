context("Bayesian latent count estimations.")

test_that(".jags files can be accessed", {
  jagsFile <- getBayesianModel('test')
  expect(file.exists(jagsFile), "Returned file not found.")
  })

test_that("dummy jags procedure works", {
  expect_output(dummyJags(), "Okay")
  })

test_that("estimateLatentCounts works", {
  model <- bimarkSimulationModel()
  res <- estimateLatentCounts(model)
  expect_is(res, class(model))
  })

