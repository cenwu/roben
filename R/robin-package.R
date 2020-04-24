#' @useDynLib robin, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @docType package
#' @keywords overview
#' @name robin-package
#' @title Robust Bayesian Variable Selection for Gene-Environment Interactions
#' @aliases robin-package
#' @description Gene-environment (G\eqn{\times}E) interactions have important implications to elucidate the etiology of complex diseases beyond the main genetic and environmental effects. Outliers and data contamination in disease phenotypes of G\eqn{\times}E studies have been commonly encountered, leading to the development of a broad spectrum of robust penalization methods. Nevertheless, within the Bayesian framework, the issue has not been taken care of in existing studies. We develop a robust Bayesian variable selection method for G\eqn{\times}E interaction studies. The proposed Bayesian method can effectively accommodate heavy--tailed errors and outliers in the response variable while conducting variable selection by accounting for structural sparsity. In particular, the spike--and--slab priors have been imposed on both individual and group levels to identify important main and interaction effects. An efficient Gibbs sampler has been developed to facilitate fast computation.
#'
#' @details The user friendly, integrated interface robin() allows users to flexibly choose the fitting methods they prefer. There are three arguments in robin() that control the fitting method
#' \tabular{rl}{
#' robust: \tab whether to robust methods. \cr\cr
#' sparse: \tab whether to use the spike-and-slab priors to achieve sparsity. \cr\cr
#' structure: \tab structural identification. Three choices are available: \cr \tab "sparsegroup", "group" and “individual”.
#' }
#'
#' robin() returns a robin object that contains the posterior estimates of each coefficients.
#' S3 generic functions BVSelection(), predict() and print() are implemented for robin objects.
#' BVSelection() takes a robin object and returns the variable selection results.
#' predict() takes a robin object and returns the predicted values for new observations.
#'
#' @references
#' Ren, J., Zhou, F., Li, X., Ma, S., Jiang, Y. and Wu, C. (2020). Robust Bayesian variable selection for gene-environment interactions.
#'
#' @seealso \code{\link{robin}}
NULL
