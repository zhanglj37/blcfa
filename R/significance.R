# Tables of estimates and selection of significant parameters.
#
# Selection criteria (argument `criterion`):
#   "hpd"       : 0 is outside the 95% HPD interval              (default; interval = TRUE)
#   "pvalue"    : one-sided posterior p-value < .05 (Muthen, 2010) (interval = FALSE, as in the paper)
#   "threshold" : |residual correlation| > 0.15, |loading| > 0.1   (README "Tips3")

resolve_criterion <- function(criterion, interval) {
  if (is.null(criterion)) {
    return(if (isTRUE(interval)) "hpd" else "pvalue")
  }
  match.arg(criterion, c("hpd", "pvalue", "threshold"))
}

is_sig <- function(lower, upper, p, value, criterion, threshold) {
  switch(criterion,
         hpd = lower * upper > 0,
         pvalue = p < 0.05,
         threshold = abs(value) > threshold)
}

# Two-sided normal-approximation p-value, used by the original code for
# MU, PHI, residual variances (blcfa_ly) and loadings (blcfa_psx).
normal_p <- function(est, se) {
  z <- est / se
  ifelse(z > 0, 2 * (1 - stats::pnorm(z)), 2 * stats::pnorm(z))
}

# ---- residual (co)variances --------------------------------------------------

#' @keywords internal
sig_psx_fun <- function(res, hpd, usevar, criterion) {
  NY <- length(usevar)
  pr <- lower_tri_pairs(NY)
  EmPSX <- res$EmPSX
  d <- sqrt(diag(EmPSX))
  CORPSX <- EmPSX / outer(d, d)

  idx <- pr[, c("i", "j"), drop = FALSE]
  OUTPSX <- cbind(est = EmPSX[idx], se = res$SEPSX[idx], cor = CORPSX[idx],
                  "p-value" = res$PPSX[idx],
                  HPD_lower = hpd$HPD_PSX[, 1], HPD_upper = hpd$HPD_PSX[, 2])
  rownames(OUTPSX) <- paste(usevar[pr[, "i"]], "with", usevar[pr[, "j"]])

  off <- pr[, "i"] > pr[, "j"]
  sig <- off & is_sig(OUTPSX[, "HPD_lower"], OUTPSX[, "HPD_upper"],
                      OUTPSX[, "p-value"], OUTPSX[, "cor"], criterion, 0.15)
  SIGPSX <- OUTPSX[sig, c("cor", "p-value", "HPD_lower", "HPD_upper"),
                   drop = FALSE]
  list(OUTPSX = OUTPSX, SIGPSX = SIGPSX,
       sigloc = pr[sig, , drop = FALSE],
       sigpsxname = rownames(SIGPSX))
}

# Residual variances only (blcfa_ly: the residual covariance matrix is diagonal).
sig_psx_diag <- function(res, hpd, usevar) {
  NY <- length(usevar)
  pr <- lower_tri_pairs(NY)
  dg <- which(pr[, "i"] == pr[, "j"])
  est <- diag(res$EmPSX)
  se <- diag(res$SEPSX)
  OUTPSX <- cbind(est = est, se = se, cor = 1, "p-value" = normal_p(est, se),
                  HPD_lower = hpd$HPD_PSX[dg, 1], HPD_upper = hpd$HPD_PSX[dg, 2])
  rownames(OUTPSX) <- usevar
  empty <- OUTPSX[0, c("cor", "p-value", "HPD_lower", "HPD_upper"), drop = FALSE]
  list(OUTPSX = OUTPSX, SIGPSX = empty, sigloc = pr[0, , drop = FALSE],
       sigpsxname = character(0))
}

# ---- loadings ----------------------------------------------------------------

# keep_all = TRUE (blcfa_psx): every free loading is kept, p-values from the
# normal approximation (as in the original sig_ly_set_fun), and loadings are
# in gibbs_psx_fun's item-by-item order.
#' @keywords internal
sig_ly_fun <- function(res, hpd, IDY0, usevar, factor_names, criterion,
                       keep_all = FALSE) {
  free <- which(IDY0 != 0 & IDY0 != 9, arr.ind = TRUE)
  # order of the sampler's ELY columns: gibbs_fun / gibbs_ly_fun store
  # LY[IDY != 0] (column-major: factor by factor); gibbs_psx_fun stores the
  # free loadings item by item (row-major)
  free <- if (keep_all) free[order(free[, 1], free[, 2]), , drop = FALSE] else
    free[order(free[, 2], free[, 1]), , drop = FALSE]
  code <- IDY0[free]
  p <- if (keep_all) normal_p(res$EmLY, res$SELY) else res$PLY

  OUTLY <- cbind(est = res$EmLY, se = res$SELY, "p-value" = p,
                 HPD_lower = hpd$HPD_LY[, 1], HPD_upper = hpd$HPD_LY[, 2])
  rownames(OUTLY) <- paste(factor_names[free[, 2]], "by", usevar[free[, 1]])

  sig <- is_sig(OUTLY[, "HPD_lower"], OUTLY[, "HPD_upper"], OUTLY[, "p-value"],
                OUTLY[, "est"], criterion, 0.1)
  # loadings coded 1 are part of the confirmatory model and always kept;
  # lasso loadings (-1) are kept only if significant
  keep <- if (keep_all) rep(TRUE, nrow(OUTLY)) else (code == 1 | (code == -1 & sig))
  lasso_sig <- code == -1 & sig
  list(OUTLY = OUTLY,
       SIGLY = OUTLY[lasso_sig, c("est", "p-value", "HPD_lower", "HPD_upper"),
                     drop = FALSE],
       keep = keep, loc = free)
}

# Loadings retained in the second-step (Mplus) model, per factor:
# marker first, then retained free loadings in item order.
retained_indicators <- function(IDY0, usevar, sigly) {
  NZ <- ncol(IDY0)
  out <- lapply(seq_len(NZ), function(k) {
    marker <- usevar[IDY0[, k] == 9]
    rows <- sigly$loc[sigly$keep & sigly$loc[, 2] == k, 1]
    c(marker, usevar[sort(rows)])
  })
  one <- which(lengths(out) < 2)
  if (length(one) > 0) {
    warning("Factor(s) ", paste(one, collapse = ", "), " keep only the marker ",
            "indicator in the Mplus model (no other loading retained); the ",
            "Mplus model is not identified.", call. = FALSE)
  }
  out
}

# ---- other tables ------------------------------------------------------------

mu_table <- function(res, hpd, usevar) {
  MU <- cbind(EstMU = res$EmMU, SeMU = res$SEMU,
              "p-value" = normal_p(res$EmMU, res$SEMU),
              HPD_lower = hpd$HPD_MU[, 1], HPD_upper = hpd$HPD_MU[, 2])
  rownames(MU) <- usevar
  MU
}

phi_tables <- function(res, hpd, factor_names) {
  NZ <- length(factor_names)
  EmPHI <- matrix(res$EmPHI, NZ, NZ)
  SEPHI <- matrix(res$SEPHI, NZ, NZ)
  PPHI <- matrix(normal_p(res$EmPHI, res$SEPHI), NZ, NZ)
  d <- sqrt(diag(EmPHI))
  CORPHI <- EmPHI / outer(d, d)
  colnames(CORPHI) <- factor_names

  ij <- which(upper.tri(EmPHI, diag = TRUE), arr.ind = TRUE)
  ij <- ij[order(ij[, 1], ij[, 2]), , drop = FALSE]        # (1,1), (1,2), ..., (2,2), ...
  lin <- (ij[, 2] - 1) * NZ + ij[, 1]
  OUTPHI <- cbind(est = EmPHI[ij], se = SEPHI[ij], cor = CORPHI[ij],
                  "p-value" = PPHI[ij],
                  HPD_lower = hpd$HPD_PHI[lin, 1], HPD_upper = hpd$HPD_PHI[lin, 2])
  rownames(OUTPHI) <- paste(factor_names[ij[, 1]], "with", factor_names[ij[, 2]])
  list(OUTPHI = OUTPHI, CORPHI = CORPHI)
}
