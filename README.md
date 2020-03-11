# blblm

<!-- badges: start -->
<!-- badges: end -->


## Example

This is a basic example which shows you how to solve a common problem:

``` r
> library(blblm)
> fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
> coef(fit)
 (Intercept)           wt           hp        wt:hp 
 56.82338917 -10.54673065  -0.16372116   0.04191719 
> confint(fit, c("wt", "hp"))
         2.5%      97.5%
wt -12.453054 -8.3635021
hp  -0.198616 -0.1283812
> sigma(fit, confidence = TRUE)
   sigma      lwr      upr 
1.394981 1.072991 1.667282 
> predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)), confidence = TRUE)
       fit      lwr      upr
1 21.61733 20.82602 22.37032
2 18.72836 17.93236 19.51310
```


## To create an R package
Steps
- create a new R package called `blblm`
- execute
```r
> library(devtools)
Loading required package: usethis
> use_readme_md()
✔ Setting active project to '/Users/Randy/Documents/blblm'
✔ Writing 'README.md'
● Modify 'README.md'
> use_r("blb.R")
● Modify 'R/blb.R'
> use_test("blb")
✔ Adding 'testthat' to Suggests field in DESCRIPTION
✔ Creating 'tests/testthat/'
✔ Writing 'tests/testthat.R'
● Call `use_test()` to initialize a basic test file and open it for editing.
✔ Increasing 'testthat' version to '>= 2.1.0' in DESCRIPTION
✔ Writing 'tests/testthat/test-blb.R'
● Modify 'tests/testthat/test-blb.R'
> use_package("purrr")
✔ Adding 'purrr' to Imports field in DESCRIPTION
● Refer to functions with `purrr::fun()`
> use_package("magrittr")
✔ Adding 'magrittr' to Imports field in DESCRIPTION
● Refer to functions with `magrittr::fun()`
> use_roxygen_md()
✔ Setting active project to '/Users/Randy/Documents/blblm'
✔ Setting Roxygen field in DESCRIPTION to 'list(markdown = TRUE)'
✔ Setting RoxygenNote field in DESCRIPTION to '7.0.2'
● Run `devtools::document()`
> use_namespace()
✔ Setting active project to '/Users/Randy/Documents/blblm'
Overwrite pre-existing file 'NAMESPACE'?

1: Yes
2: No
3: No way

Selection: 1
✔ Writing 'NAMESPACE'
```
