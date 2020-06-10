#' @import purrr
#' @import stats
#' @import parallel
#' @import readr
#' @importFrom utils capture.output
#' @importFrom magrittr %>%
#' @aliases NULL
#' @details
#' Generalized Linear Models with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#' Bag of little bootsraps for generalized linear model.
#'
#' @param formula a formula fitting a linear regression model.
#' @param data dataset for the model
#' @param m number of subsets
#' @param B number of simulations
#' @param parallel a logical operator, T for parallelization
#' @param select_names character list for file names
#' @param family family type for models, such as binomial, gaussian, gamma etc.
#'
#' @return.
#'
#' @export
blbglm <- function(formula, data, m = 10, B = 5000,parallel = FALSE,select_names = NULL,family) {
  # When user do not select specific sub files
  if (is.null(select_names)){data_list <- split_data(data, m)}
  # When user select specific sub files, provided a list of file names
  else{
    data_list <- prepare_files(data,m,select_names)
  }
  # When user set parallel = TRUE, use parallelization
  if (parallel){
    cl <- makeCluster(detectCores()) # Detect cores and make cluster
    clusterExport(cl=cl, varlist=c("glm_each_boot","glm1","blbcoef","blbsigma"),
                  envir=environment()) # Export the functions for processing
    estimates <- parLapply(cl,data_list,
                           glm_each_subsample, formula = formula, n = nrow(data), B = B,family = family)
    stopCluster(cl)}
  # When user set parallel = FALSE(default), don't use parallelization
  else{
    estimates <- map(
      data_list,
      ~ glm_each_subsample(formula = formula, data = ., n = nrow(data), B = B,family = family))}
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}


#' split data into m parts of approximated equal sizes
#'
#' @inheritParams blbglm
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}

#' Prepare the subfiles, select the specific files
#'
#' @inheritParams blbglm
prepare_files <- function(data,m,select_names){
  # Prepare the sub files
  #data <- split_data(data, m)
  #dir.create("files", showWarnings = FALSE)
  #set.seed(141)

  #for(i in names(data)){
  #  write_csv(data[[i]], file.path("files", paste0(i,".csv")))}

  ### Generated test files stored in folder called `files`
  file_names <- file.path("files", select_names)
  # Read the selected files
  data_list <- file_names %>% map(read_csv)
}


#' compute the estimates
#'
#' @param n number of rows of dataset
#' @inheritParams blbglm
glm_each_subsample <- function(formula, data, n, B,family) {
  replicate(B, glm_each_boot(formula, data, n,family), simplify = FALSE)
}


#' compute the regression estimates for a blb dataset
#'
#' @inheritParams glm_each_subsample
glm_each_boot <- function(formula, data, n,family) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  glm1(formula, data, freqs,family)
}


#' estimate the regression estimates based on given the number of repetitions
#'
#' @param freqs frequency for each bootstrapping
#' @inheritParams blbglm
glm1 <- function(formula, data, freqs,family) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- glm(formula,family = family, data = data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' compute the coefficients from fit
#'
#' @param fit fitting the bootstrapping results
blbcoef <- function(fit) {
  coef(fit)
}


#' compute sigma from fit
#'
#' @inheritParams blbcoef
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' @export
#' @method print blblm
print.blbglm <- function(x, ...) {
  cat("blbglm model:", capture.output(x$formula))
  cat("\n")
}


#' @export
#' @method sigma blblm
sigma.blbglm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' @export
#' @method coef blblm
coef.blbglm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}


#' @export
#' @method confint blblm
confint.blbglm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' @export
#' @method predict blblm
predict.blbglm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
               apply(1, mean_lwr_upr, level = level) %>%
               t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
