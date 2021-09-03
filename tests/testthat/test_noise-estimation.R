test_that("Nose estimation behaves as expected",{
  expect_length(sparse.mad.estimator(1:10), 1)
  expect_length(itses:::mad.estimator(1:10), 1)
  expect_gt(itses:::mad.estimator(1:10), 0)
  expect_lte(itses:::sparsity.estimator(sparse.mad.estimator(1:10)), 1)
  expect_gte(itses:::sparsity.estimator(sparse.mad.estimator(1:10)), 0)
  expect_equal(sparse.mad.estimator(1:10, h = 1.2), itses:::mad.estimator(1:10))
  expect_true(sparse.mad.estimator(1:10, h = 0) != itses:::mad.estimator(1:10))
})
test_that("Threshold project behaves as expected",{
  expect_equal(itses:::get.thresholded.lambda(lambda = 10, max.threshold = 8, min.threshold = 2), 8)
  expect_equal(itses:::get.thresholded.lambda(lambda = 0, max.threshold = 8, min.threshold = 2), 2)
  expect_equal(itses:::get.thresholded.lambda(lambda = 5, max.threshold = 8, min.threshold = 2), 5)
})