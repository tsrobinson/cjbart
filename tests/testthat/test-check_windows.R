test_that("runs model with seed warning on Windows", {

  test_df <- data.frame(Y = rbinom(50,1, prob = 0.5),
                        X1 = runif(50),
                        X2 = runif(50),
                        id = 1:50)

  if (.Platform$OS.type=='unix') {
    testthat::expect_no_warning(
      cjbart::cjbart(data = test_df, Y = "Y", id = "id", seed = 1, cores = 1,
                     ntree = 5, numcut = 10)
    )
  } else {
    testthat::expect_warning(
      cjbart::cjbart(data = test_df, Y = "Y", id = "id", seed = 1, cores = 1,
                     ntree = 5, numcut = 10)
    )
  }

})
