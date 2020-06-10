test_that("blblm works", {
  fit.blblm <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
  expect_s3_class(fit.blblm, "blblm")
  co <- coef(fit.blblm)
  expect_equal(length(co), 4)
})



