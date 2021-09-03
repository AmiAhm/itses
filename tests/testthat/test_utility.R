
test_that("Soft-threshold behaves as expected",{
  expect_equal(itses:::soft.threshold.estimator(c(1,2), 1), c(0, 1))
  expect_equal(itses:::soft.threshold.estimator(c(-2, -0.5,1,2), 1), c(-1, 0, 0, 1))
})


test_that("Hard-threshold behaves as expected",{
  expect_equal(itses:::hard.threshold.estimator(c(1,2), 1), c(0, 2))
  expect_equal(itses:::hard.threshold.estimator(c(-2, -0.5,1,2), 1), c(-2, 0, 0, 2))
})


test_that("Loss functions behaves as expected",{
  expect_equal(itses:::loss.w.ht(lambda = 1, y = c(1,2,3), theta = c(0,0,0)), 13)
  expect_equal(itses:::loss.w.st(lambda = 1, y = c(1,2,3), theta = c(0,0,0)), 5)
})


test_that("Risk with soft-threshold has expected properties",{
  expect_equal(itses:::risk.st(theta = 1, lambda = 1), 2-(pnorm(0)-pnorm(-2))-(2)*dnorm(0))

  expect_length(itses:::risk.st(theta = -10:10, lambda = 3), 1)
  expect_length(itses:::risk.st(theta = -5:5, lambda = 3), 1)

  expect_equal(itses:::risk.st(theta = 1:5, lambda = 3)-itses:::risk.st(theta = -(1:5), lambda = 3), 0)
  expect_equal(itses:::risk.st(theta = 1:5, lambda = 2)-itses:::risk.st(theta = -(1:5), lambda = 2), 0)

  expect_lte(itses:::risk.st(theta = 0, lambda = 2), itses:::risk.st(theta = 0, lambda = 2))
  expect_lte(itses:::risk.st(theta = 0, lambda = 10), itses:::risk.st(theta = 0, lambda = 10))
  expect_lte(itses:::risk.st(theta = 5, lambda = 2), itses:::risk.st(theta = 5, lambda = 2)+5^2)

  expect_lte(itses:::dlambda.risk.st(theta = 0, lambda = 0.2), 0)
  expect_lte(itses:::dlambda.risk.st(theta = 0, lambda = 0.4), 0)
  expect_lte(itses:::dlambda.risk.st(theta = 0, lambda = 0.65), 0)
  expect_lte(itses:::dlambda.risk.st(theta = 0, lambda = 10), 0)
  expect_lte(itses:::dlambda.risk.st(theta = 0, lambda = 1e2), 0)
  expect_lte(itses:::dlambda.risk.st(theta = 0, lambda = 1e3), 0)

  expect_gte(itses:::dlambda.risk.st(theta = Inf, lambda = 1), 0)
  expect_gte(itses:::dlambda.risk.st(theta = Inf, lambda = 2), 0)
  expect_gte(itses:::dlambda.risk.st(theta = Inf, lambda = 10), 0)
  expect_gte(itses:::dlambda.risk.st(theta = Inf, lambda = 50), 0)

  expect_equal(itses:::dlambda.risk.st(theta = 1:5, lambda = 3)-itses:::dlambda.risk.st(theta = -(1:5), lambda = 3), 0)
  expect_equal(itses:::dlambda.risk.st(theta = 1:5, lambda = 2)-itses:::dlambda.risk.st(theta = -(1:5), lambda = 2), 0)

  expect_equal(itses:::d2lambda.risk.st(theta = 1:5, lambda = 3)-itses:::d2lambda.risk.st(theta = -(1:5), lambda = 3), 0)
  expect_equal(itses:::d2lambda.risk.st(theta = 1:5, lambda = 2)-itses:::d2lambda.risk.st(theta = -(1:5), lambda = 2), 0)
})


test_that("Risk with hard-threshold has expected properties",{
  expect_equal(itses:::risk.ht(theta = 1, lambda = 1), (pnorm(0)-pnorm(-2))+2-pnorm(0)-pnorm(2)+2*dnorm(2))

  expect_length(itses:::risk.ht(theta = -10:10, lambda = 3), 1)
  expect_length(itses:::risk.ht(theta = -5:5, lambda = 3), 1)

  expect_equal(itses:::risk.ht(theta = 1:5, lambda = 3)-itses:::risk.ht(theta = -(1:5), lambda = 3), 0)
  expect_equal(itses:::risk.ht(theta = 1:5, lambda = 2)-itses:::risk.ht(theta = -(1:5), lambda = 2), 0)

  expect_equal(itses:::dlambda.risk.ht(theta = 1:5, lambda = 3)-itses:::dlambda.risk.ht(theta = -(1:5), lambda = 3), 0)
  expect_equal(itses:::dlambda.risk.ht(theta = 1:5, lambda = 2)-itses:::dlambda.risk.ht(theta = -(1:5), lambda = 2), 0)

  expect_equal(itses:::dlambda.risk.ht(theta = 1:5, lambda = 3)-itses:::dlambda.risk.ht(theta = -(1:5), lambda = 3), 0)
  expect_equal(itses:::dlambda.risk.ht(theta = 1:5, lambda = 2)-itses:::dlambda.risk.ht(theta = -(1:5), lambda = 2), 0)
})

