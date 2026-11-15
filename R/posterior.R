# Posterior summaries of the retained draws (iterations N.burn+1 .. MCMAX).
# Point estimate = posterior mean, se = posterior sd,
# p-value = posterior probability of the opposite sign
#           (share of draws < 0 if the mean is > 0, else share of draws > 0),
# HPD = 95% highest posterior density interval (coda::HPDinterval).

# column means computed exactly as mean() (same values as the original code)
col_mean <- function(x) apply(x, 2, mean)

tail_prob <- function(x) {
  # x: draws x parameters
  m <- col_mean(x)
  ifelse(m > 0, colMeans(x < 0), colMeans(x > 0))
}

#' @keywords internal
caculate_results <- function(chain, MCMAX, N.burn) {
  keep <- (N.burn + 1):MCMAX
  ELY <- chain$ELY[keep, , drop = FALSE]
  EMU <- chain$EMU[keep, , drop = FALSE]
  EPHI <- chain$EPHI[keep, , drop = FALSE]
  EPSX <- chain$EPSX[keep, , , drop = FALSE]
  NY <- dim(EPSX)[2]
  EPSX2 <- matrix(EPSX, nrow = length(keep), ncol = NY * NY)

  list(
    EmLY = col_mean(ELY), SELY = apply(ELY, 2, stats::sd), PLY = tail_prob(ELY),
    EmMU = col_mean(EMU), SEMU = apply(EMU, 2, stats::sd), PMU = tail_prob(EMU),
    EmPHI = col_mean(EPHI), SEPHI = apply(EPHI, 2, stats::sd),
    PPHI = tail_prob(EPHI),
    EmPSX = matrix(col_mean(EPSX2), NY, NY),
    SEPSX = matrix(apply(EPSX2, 2, stats::sd), NY, NY),
    PPSX = matrix(tail_prob(EPSX2), NY, NY),
    Empostp = mean(chain$Epostp)
  )
}

#' @keywords internal
hpd_fun <- function(chain, MCMAX, N.burn) {
  keep <- (N.burn + 1):MCMAX
  hpd <- function(x) {
    out <- coda::HPDinterval(coda::mcmc(x[keep, , drop = FALSE]))
    unname(out[, 1:2, drop = FALSE])
  }
  list(HPD_LY = hpd(chain$ELY), HPD_MU = hpd(chain$EMU),
       HPD_PHI = hpd(chain$EPHI),
       # rows follow the lower triangle by row: (1,1), (2,1), (2,2), (3,1), ...
       HPD_PSX = hpd(chain$chainpsx))
}

# Lower-triangle (row-wise, i >= j) index pairs, matching chain$chainpsx.
lower_tri_pairs <- function(NY) {
  i <- rep(seq_len(NY), seq_len(NY))
  j <- sequence(seq_len(NY))
  cbind(i = i, j = j)
}

# Missing data: replace each missing (standardised) value by the mean of its
# retained draws stored by the sampler (computed in run_one_chain() as
# chain$Y_missing_mean) and write data_imputed.txt.
# Returns TRUE if any value was missing.
#' @keywords internal
impute_ms <- function(Y, chain, file = "data_imputed.txt") {
  missing_ind <- chain$missing_ind
  if (sum(missing_ind) == 0) return(FALSE)
  dataset_new <- t(Y)
  idx <- which(missing_ind == 1, arr.ind = TRUE)       # rows: item, cols: person
  dataset_new[idx[, c(2, 1), drop = FALSE]] <- chain$Y_missing_mean[idx]
  utils::write.table(dataset_new, file, col.names = FALSE, row.names = FALSE,
                     sep = "\t")
  TRUE
}
