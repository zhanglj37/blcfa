# First-step result tables; with bloutput = TRUE also written to ./results/
#   EPSR.png ppp.csv ly.csv mu.csv phi.csv phi_cormatrix.csv psx.csv psx_sig.csv

#' @keywords internal
write_results <- function(res, hpd, sigpsx, sigly, usevar, factor_names,
                          epsr, N.burn, bloutput, dir = "results") {
  MU <- mu_table(res, hpd, usevar)
  phi <- phi_tables(res, hpd, factor_names)

  if (bloutput) {
    dir.create(dir, showWarnings = FALSE)
    f <- function(x) file.path(dir, x)
    utils::write.csv(res$Empostp, file = f("ppp.csv"), row.names = FALSE)
    utils::write.csv(sigly$OUTLY, file = f("ly.csv"), row.names = TRUE)
    utils::write.csv(MU, file = f("mu.csv"), row.names = TRUE)
    utils::write.csv(phi$OUTPHI, file = f("phi.csv"), row.names = TRUE)
    utils::write.csv(phi$CORPHI, file = f("phi_cormatrix.csv"), row.names = FALSE)
    utils::write.csv(sigpsx$OUTPSX, file = f("psx.csv"), row.names = TRUE)
    if (nrow(sigpsx$SIGPSX) > 0) {
      utils::write.csv(sigpsx$SIGPSX, file = f("psx_sig.csv"), row.names = TRUE)
    } else {
      utils::write.csv("no sig residual correlation", file = f("psx_sig.csv"),
                       row.names = TRUE)
    }
    if (nrow(sigly$SIGLY) > 0) {
      utils::write.csv(sigly$SIGLY, file = f("ly_sig.csv"), row.names = TRUE)
    }
    if (!is.null(epsr)) EPSR_figure(epsr, N.burn, f("EPSR.png"))
  }

  list(ppp = res$Empostp, ly = round(sigly$OUTLY, 3), mu = round(MU, 3),
       phi = round(phi$OUTPHI, 3), psx = round(sigpsx$OUTPSX, 3))
}
