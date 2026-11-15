#' blcfa: Bayesian Model Modification in Confirmatory Factor Analysis
#'
#' Detects residual covariances (and cross-loadings) with the Bayesian
#' covariance lasso and passes the modified CFA model to Mplus.
#'
#' @keywords internal
#' @importFrom MASS mvrnorm
#' @importFrom MCMCpack rwish
#' @importFrom statmod rinvgauss
#' @importFrom stats rgamma rnorm runif pnorm sd
#' @importFrom utils getFromNamespace
"_PACKAGE"
