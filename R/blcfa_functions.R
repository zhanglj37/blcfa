#' Bayesian covariance lasso CFA with Mplus second step
#'
#' Step 1 estimates the CFA with a Bayesian covariance lasso prior on the
#' residual covariance matrix (and, for loadings coded -1 in an IDY matrix, a
#' lasso prior on the loadings) and detects the residual covariances (and
#' cross-loadings) that differ from zero. Step 2 frees these parameters in an
#' Mplus input file (ML/MLM and/or Bayes) and runs it if Mplus is available.
#'
#' @param filename Path of the ASCII data file (\code{.dat}/\code{.txt},
#'   whitespace-delimited, no header row).
#' @param varnames Names of all columns of the data file.
#' @param usevar Names of the variables used in the analysis.
#' @param myModel The CFA model, either as lavaan-style syntax, e.g.
#'   \code{"f1 =~ y1 + y2 + y3\nf2 =~ y4 + y5 + y6"} (the first indicator of
#'   each factor is fixed at 1; comments start with \code{#}), or as an IDY
#'   matrix with one row per \code{usevar} and one column per factor and
#'   entries 9 (fixed at 1), 1 (free), -1 (free with lasso shrinkage) or
#'   0 (fixed at 0).
#' @param estimation Estimator of the second step in Mplus: \code{"ml"}
#'   (default; MLM if the data are non-normal), \code{"bayes"} or \code{"both"}.
#' @param ms Missing-value code in the data file (\code{NA} by default).
#' @param MCMAX Total number of MCMC iterations.
#' @param N.burn Number of burn-in iterations (discarded).
#' @param bloutput If \code{TRUE}, detailed first-step results are written to
#'   the \code{results} folder.
#' @param interval \code{TRUE}: significance by the 95\% HPD interval;
#'   \code{FALSE}: by the posterior p-value (< .05). Ignored if
#'   \code{criterion} is given.
#' @param conver_check If \code{TRUE}, two chains are run and convergence is
#'   checked with the EPSR (< 1.2 over the last 100 burn-in iterations); if
#'   \code{FALSE}, one chain is run without the check.
#' @param criterion Optional: \code{"hpd"}, \code{"pvalue"} or
#'   \code{"threshold"} (|residual correlation| > .15, |loading| > .1).
#' @param outdir Directory for all output files (default: working directory).
#' @param mplus_path Optional path of the Mplus executable if it is not on
#'   the system \code{PATH}.
#' @return Invisibly, a list with \code{blcfa_est} (\code{ppp}, \code{ly},
#'   \code{mu}, \code{phi}, \code{psx}), \code{estimation},
#'   \code{convergence}, \code{criterion}, \code{sig_psx} (significant
#'   residual correlations), \code{sig_ly} (significant lasso loadings) and
#'   \code{outdir}.
#' @references
#' Zhang, L., Pan, J., Dube, L., & Ip, E. H. (2021). blcfa: An R package for
#' Bayesian model modification in confirmatory factor analysis.
#' \emph{Structural Equation Modeling}, 28(4), 649-658.
#'
#' Pan, J., Ip, E. H., & Dube, L. (2017). An alternative to post hoc model
#' modification in confirmatory factor analysis: The Bayesian lasso.
#' \emph{Psychological Methods}, 22(4), 687-704.
#' @examples
#' \dontrun{
#' filename <- system.file("extdata", "simu_data.txt", package = "blcfa")
#' varnames <- paste0("y", 1:10)
#' usevar <- varnames
#' myModel <- "
#' f1 =~ y1 + y2 + y3 + y4 + y5
#' f2 =~ y6 + y7 + y8 + y9 + y10
#' "
#' results <- blcfa(filename, varnames, usevar, myModel, estimation = "both",
#'                  MCMAX = 5000, N.burn = 2500, bloutput = TRUE,
#'                  interval = TRUE)
#' results
#' sum_second(results)
#' }
#' @export
blcfa <- function(filename, varnames, usevar, myModel, estimation = "ml",
                  ms = NA, MCMAX = 10000, N.burn = 5000, bloutput = FALSE,
                  interval = TRUE, conver_check = TRUE, criterion = NULL,
                  outdir = getwd(), mplus_path = NULL) {
  invisible(run_blcfa("full", filename, varnames, usevar, myModel, estimation,
                      ms, MCMAX, N.burn, bloutput, interval, conver_check,
                      criterion, outdir, mplus_path))
}

#' Bayesian lasso CFA for cross-loadings only
#'
#' Like \code{\link{blcfa}}, but the residual covariance matrix is diagonal,
#' so only loadings coded -1 in the IDY matrix are subject to lasso shrinkage
#' and selection (Chen, Guo, Zhang & Pan, 2021).
#'
#' @inheritParams blcfa
#' @param IDY0 The IDY matrix (see \code{myModel} in \code{\link{blcfa}});
#'   model syntax is accepted as well.
#' @return As \code{\link{blcfa}}; \code{blcfa_est$psx} contains residual
#'   variances only.
#' @export
blcfa_ly <- function(filename, varnames, usevar, IDY0, estimation = "ml",
                     ms = NA, MCMAX = 10000, N.burn = 5000, bloutput = FALSE,
                     interval = TRUE, conver_check = TRUE, criterion = NULL,
                     outdir = getwd(), mplus_path = NULL) {
  invisible(run_blcfa("ly", filename, varnames, usevar, IDY0, estimation,
                      ms, MCMAX, N.burn, bloutput, interval, conver_check,
                      criterion, outdir, mplus_path))
}

#' Bayesian covariance lasso CFA with a fixed loading structure
#'
#' Uses the \code{gibbs_psx_fun} sampler: loadings are estimated without
#' shrinkage (IDY codes 9/1/0 only) and all of them are kept in the Mplus
#' model; residual covariances are detected as in \code{\link{blcfa}}.
#'
#' @inheritParams blcfa
#' @param model Model syntax or IDY matrix without -1 entries.
#' @return As \code{\link{blcfa}}.
#' @export
blcfa_psx <- function(filename, varnames, usevar, model, estimation = "ml",
                      ms = NA, MCMAX = 10000, N.burn = 5000, bloutput = FALSE,
                      interval = TRUE, conver_check = TRUE, criterion = NULL,
                      outdir = getwd(), mplus_path = NULL) {
  invisible(run_blcfa("psx", filename, varnames, usevar, model, estimation,
                      ms, MCMAX, N.burn, bloutput, interval, conver_check,
                      criterion, outdir, mplus_path))
}
