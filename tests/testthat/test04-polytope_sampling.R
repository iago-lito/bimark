context("Bayesian latent count estimations.")

test_that(".jags files can be accessed", {
  jagsFile <- GetBayesianModel('test')
  expect(file.exists(jagsFile), "Returned file not found.")
  })

test_that("dummy jags procedure works", {
  expect_output(DummyJags(), "Okay")
  })

test_that("EstimateLatentCounts works", {
  model <- BimarkSimulationModel()
  res <- EstimateLatentCounts(model)
  expect_is(res, class(model))
  })

