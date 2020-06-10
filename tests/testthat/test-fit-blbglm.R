test_that("blbglm works", {
  fit.blbglm <- blbglm(Species ~ Sepal.Length*Sepal.Width, data = iris[1:100,], m = 3, B = 100,family = binomial())
  expect_s3_class(fit.blbglm, "blblm")
  co <- coef(fit.blbglm)
  expect_equal(length(co), 4)
})