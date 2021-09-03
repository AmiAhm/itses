test_that("Oracle properties as expected",{
  expect_lt(itses::itses(1:10, method = "ST", sd = 1)$lambda,
            itses:::get.visu.threshold(1:10))

  expect_lt(itses::itses(-5:10, method = "ST", sd = 1)$lambda,
            itses:::get.visu.threshold(-5:10))

  expect_lt(itses:::get.risk.oracle.threshold(theta = -5:10, method = "ST")$lambda,
            itses:::get.visu.threshold(-5:10))

  expect_lt(itses:::get.risk.oracle.threshold(theta = -5:100, method = "ST")$lambda,
            itses:::get.visu.threshold(-5:10))

  expect_lt(itses:::get.risk.oracle.threshold(theta = -5:1000, method = "ST")$lambda,
            itses:::get.visu.threshold(-5:10))

  expect_lt(itses:::get.risk.oracle.threshold(theta = -100:-50, method = "ST")$lambda,
            itses:::get.visu.threshold(-5:10))
})

test_that("Grid function behaves as expected",{
  expect_length(itses:::get.lambda_grid(seq(from = 0.1, to = 0.88, length.out = 10), k = 15), 15)
  expect_length(itses:::get.lambda_grid(seq(from = 0.1, to = 0.88, length.out = 10), k = 8), 8)
  expect_length(itses:::get.lambda_grid(seq(from = 0.1, to = 0.88, length.out = 50), k = 8), 8)
  expect_length(itses:::get.lambda_grid(seq(from = 0.1, to = 0.88, length.out = 44), k = 20), 20)
})
