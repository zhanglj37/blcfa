#' Summarise the second-step (Mplus) analysis
#'
#' Reads the Mplus output file(s) written by \code{\link{blcfa}} and returns
#' model fit summaries and (standardized) parameter estimates. Mplus itself is
#' not needed, only the \code{.out} files.
#'
#' @param blcfa_results The object returned by \code{\link{blcfa}},
#'   \code{\link{blcfa_ly}} or \code{\link{blcfa_psx}}.
#' @return For \code{estimation = "both"}: a list with \code{bayes_fit},
#'   \code{bayes_par_est}, \code{bayes_par_est_std}, \code{ml_fit},
#'   \code{ml_par_est}, \code{ml_par_est_std}. For a single estimator: a list
#'   with \code{estimates} and \code{std_estimates}.
#' @export
sum_second <- function(blcfa_results) {
  estimation <- normalize_estimation(blcfa_results$estimation)
  if (isFALSE(blcfa_results$convergence)) {
    cat("Error: the first step did not converge, so no second-step (Mplus)",
        "analysis was run.\n")
    return(invisible(NULL))
  }
  outdir <- if (is.null(blcfa_results$outdir)) getwd() else blcfa_results$outdir
  files <- switch(estimation,
                  bayes = "blcfa_bayes.out",
                  both = c("blcfa_bayes.out", "blcfa_ml.out"),
                  ml = "blcfa_ml.out")
  paths <- file.path(outdir, files)
  if (!all(file.exists(paths))) {
    cat("Error: Failed to run the Mplus software and get the results of the",
        "second-step analysis (missing: ",
        paste(paths[!file.exists(paths)], collapse = ", "), ")\n")
    return(invisible(NULL))
  }

  read_one <- function(path) {
    cat("Reading model: ", basename(path), "\n")
    m <- MplusAutomation::readModels(path, quiet = TRUE)
    std <- m$parameters$stdyx.standardized
    if (is.null(std) && length(m$parameters) >= 3) std <- m$parameters[[3]]
    list(fit = m$summaries, est = m$parameters$unstandardized, std = std)
  }

  if (length(paths) == 1) {
    r <- read_one(paths)
    list(estimates = r$est, std_estimates = r$std)
  } else {
    b <- read_one(paths[1])
    m <- read_one(paths[2])
    list(bayes_fit = b$fit, bayes_par_est = b$est, bayes_par_est_std = b$std,
         ml_fit = m$fit, ml_par_est = m$est, ml_par_est_std = m$std)
  }
}
