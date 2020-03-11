# blblm

<!-- badges: start -->
<!-- badges: end -->

## Examples

``` r
library(blblm)
fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
coef(fit)
#> (Intercept)          wt          hp       wt:hp 
#> 40.24493577 -4.69858460 -0.08122181  0.01325703
confint(fit, c("wt", "hp"))
#>          2.5%       97.5%
#> wt -6.2715264 -3.16035538
#> hp -0.1117288 -0.05037275
sigma(fit)
#> [1] 1.016356
sigma(fit, confidence = TRUE)
#>     sigma       lwr       upr 
#> 1.0163565 0.7578867 1.2286557
predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)))
#>        1        2 
#> 21.28659 19.10256
predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)), confidence = TRUE)
#>        fit      lwr      upr
#> 1 21.28659 20.68118 21.83728
#> 2 19.10256 18.46423 19.70810
```
