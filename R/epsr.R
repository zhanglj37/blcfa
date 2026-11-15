# Convergence check: estimated potential scale reduction (EPSR; Gelman, 1996)
# computed from two MCMC chains for every monitored parameter, for each
# iteration count n = 2, ..., N.burn + 1, using only the first n draws.
# The formula is identical to sna::potscalered.mcmc(); it is evaluated here
# with cumulative sums so all iterations are computed at once.

# Matrix of monitored draws (iterations x parameters) for one chain:
# loadings, intercepts, factor (co)variances and the lower triangle
# (including the diagonal) of the residual covariance matrix.
epsr_param_matrix <- function(chain) {
  cbind(chain$ELY, chain$EMU, chain$EPHI, chain$chainpsx)
}

# EPSR "chain" for each parameter; row k = EPSR based on draws 1..(k+1).
# x1, x2: iterations x parameters, two chains.
epsr_running <- function(x1, x2) {
  n <- seq(2, nrow(x1))
  run_mean_var <- function(x) {
    # centre on the first draw to limit floating-point cancellation
    xc <- sweep(x, 2, x[1, ], "-")
    s1 <- apply(xc, 2, cumsum)[n, , drop = FALSE]
    s2 <- apply(xc^2, 2, cumsum)[n, , drop = FALSE]
    m <- s1 / n
    v <- (s2 - n * m^2) / (n - 1)
    v[v < 0] <- 0
    list(m = sweep(m, 2, x[1, ], "+"), v = v)
  }
  a <- run_mean_var(x1)
  b <- run_mean_var(x2)
  # J = 2 chains
  mtot <- (a$m + b$m) / 2
  B <- n * ((a$m - mtot)^2 + (b$m - mtot)^2)          # n/(J-1) * sum
  W <- (a$v + b$v) / 2
  sqrt((((n - 1) / n) * W + B / n) / W)
}

#' @keywords internal
caculate_epsr <- function(N.burn, chain1, chain2, threshold = 1.2) {
  rows <- seq_len(N.burn + 1)
  x1 <- epsr_param_matrix(chain1)[rows, , drop = FALSE]
  x2 <- epsr_param_matrix(chain2)[rows, , drop = FALSE]
  epsr <- epsr_running(x1, x2)

  # converged if every parameter has EPSR < threshold over the last
  # 100 iterations of the burn-in period
  last <- epsr[(N.burn - 99):N.burn, , drop = FALSE]
  convergence <- any(!is.na(last)) && isTRUE(all(last < threshold, na.rm = TRUE))
  list(epsr = epsr, convergence = convergence)
}

#' @keywords internal
EPSR_figure <- function(epsr, N.burn, file = "EPSR.png") {
  mycolsi <- grDevices::rainbow(ncol(epsr), s = 1, v = 1, start = 0,
                                end = max(1, ncol(epsr) - 1) / ncol(epsr),
                                alpha = 1)
  repxlim <- seq_len(N.burn)
  epsr <- epsr[repxlim, , drop = FALSE]
  grDevices::png(file)
  on.exit(grDevices::dev.off())
  graphics::plot(x = repxlim, y = epsr[, 1], type = "l",
                 xlab = "iterations", ylab = "EPSR", ylim = c(0, 3),
                 col = mycolsi[1])
  for (i in seq_len(ncol(epsr))[-1]) {
    graphics::lines(x = repxlim, y = epsr[, i], col = mycolsi[i])
  }
  graphics::abline(h = 1.2, col = "black")
  invisible(file)
}
