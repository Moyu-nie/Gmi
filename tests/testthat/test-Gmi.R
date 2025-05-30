

test_that("eigen_matmul(): function return the soft thresholding", {
  expect_equal(eigen_matmul(matrix(3,3),matrix(3,3)), 3)
})
