
test_that("itses output as expected",{
  expect_length(itses::itses(1:10)$lambda, 1)
  expect_length(itses::itses(-10:10)$lambda, 1)
  expect_equal(itses::itses(-10:10, method = "ST")$method, "ST")
  expect_equal(itses::itses(-10:10, method = "HT")$method, "HT")
  expect_equal(itses::itses(-10:10, method = "ST", remove.zero = FALSE)$sd, itses::sparse.mad.estimator(-10:10))
  expect_equal(itses::itses(-10:10, method = "ST", remove.zero = FALSE, sd = 1)$sd, 1)
  expect_false(itses::itses(-10:10, method = "ST", remove.zero = TRUE)$sd == itses::sparse.mad.estimator(-10:10))
})


test_that("itses esimates noise correctly",{
  y <- seq(0, 1, length.out = 20)
  y.removed <- y[y!=0]
  sd <- sparse.mad.estimator(y.removed)
  expect_equal(itses(y)$sd, sd)
  expect_true(itses(y, sd = 2)$sd != sd)
  expect_equal(itses(y, sd = 2)$sd, 2)

  sd <- sparse.mad.estimator(y.removed, h = 0.2)
  expect_equal(itses(y, h = 0.2)$sd, sd)
  expect_true(itses(y)$sd != sd)

  sd <- sparse.mad.estimator(y)
  expect_equal(itses(y, remove.zero = FALSE)$sd, sd)


})

test_that("oracle is best: ", {
  theta <- c(0.0000000, 0.0000000, 0.0000000, 6.8506766, 7.8604660, 0.0000000, 0.0000000,
             0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 10.1909446,
             0.0000000, 6.5928580, 0.0000000, 0.0000000, 14.5193173, 0.0000000, 0.0000000,
             0.0000000, 0.0000000, 7.1626734, 0.0000000, 0.0000000, 0.0000000,
             0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, -1.0255067,
             0.0000000, 6.0754743, 0.0000000, 5.5105186, 0.0000000, 0.0000000,
             7.5529163, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 7.6240590, 0.0000000,
             -0.2568003, 0.0000000, 0.0000000, 0.0000000)
  
  y <- c(0.92892225, 1.85712114, -0.37883526, 4.99732520, 7.73890943, -0.04127880, -0.39256058, 1.12654543, -1.13659235,
         -0.28246305, 1.12224822, -1.30898171,
         8.42664455, 0.01983619, 6.42798258, 0.29942231, 0.29905206, 13.69868092, 0.24198747, -0.15800451, -0.80545987,
         -1.06137907, 7.43251169, -0.45389181,
         -0.87050720, 0.04672230, -1.24278365, 0.35549245, -0.35195347, -1.46368956, 0.58973000, -0.67609783, -1.36102443,
         0.75173650, 5.66741423, -0.79724988,
         5.26261203, -0.19744906, -0.55706885, 7.05735683, -0.68885342, -0.85796593, -1.32803186, 0.58041771, 8.27165800,
         0.62419179, -3.14574250, 0.90639487, -0.75056536, 0.97532792)
  itses.threshold <- itses(y, method = "HT")$lambda
  itses.risk <- itses:::risk.ht(theta, itses.threshold)
  oracle.threshold <- itses:::get.risk.oracle.threshold(theta, method = "HT", also.check = itses.threshold)$lambda
  oracle.risk <- itses:::risk.ht(theta, oracle.threshold)
  expect_gte(itses.risk, oracle.risk)

  itses.threshold <- itses(y, method = "ST")$lambda
  itses.risk <- itses:::risk.ht(theta, itses.threshold)
  oracle.threshold <- itses:::get.risk.oracle.threshold(theta, method = "ST", also.check = itses.threshold)$lambda
  oracle.risk <- itses:::risk.ht(theta, oracle.threshold)
  expect_gte(itses.risk, oracle.risk)
})