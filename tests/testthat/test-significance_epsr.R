# synthetic posterior summaries for NY = 3 items, NZ = 1 factor
fake <- function() {
  NY <- 3
  EmPSX <- matrix(c(1, 0.5, 0.05, 0.5, 1, 0.2, 0.05, 0.2, 1), 3)
  res <- list(EmPSX = EmPSX, SEPSX = matrix(0.1, 3, 3),
              PPSX = matrix(c(0, 0.001, 0.30, 0.001, 0, 0.04, 0.30, 0.04, 0), 3))
  # lower triangle by row: (1,1) (2,1) (2,2) (3,1) (3,2) (3,3)
  hpd <- list(HPD_PSX = rbind(c(.8, 1.2), c(.3, .7), c(.8, 1.2),
                              c(-.1, .2), c(-.01, .4), c(.8, 1.2)))
  list(res = res, hpd = hpd, usevar = c("a", "b", "c"))
}

test_that("residual covariance selection: hpd, pvalue, threshold", {
  f <- fake()
  s <- sig_psx_fun(f$res, f$hpd, f$usevar, "hpd")
  expect_equal(rownames(s$SIGPSX), "b with a")
  expect_equal(nrow(s$OUTPSX), 6)
  s <- sig_psx_fun(f$res, f$hpd, f$usevar, "pvalue")
  expect_equal(rownames(s$SIGPSX), c("b with a", "c with b"))
  s <- sig_psx_fun(f$res, f$hpd, f$usevar, "threshold")
  expect_equal(rownames(s$SIGPSX), c("b with a", "c with b"))
  expect_equal(colnames(s$SIGPSX), c("cor", "p-value", "HPD_lower", "HPD_upper"))
})

test_that("zero significant pairs gives an empty table, not a crash", {
  f <- fake(); f$hpd$HPD_PSX[2, ] <- c(-.1, .7)
  s <- sig_psx_fun(f$res, f$hpd, f$usevar, "hpd")
  expect_equal(nrow(s$SIGPSX), 0)
  expect_length(s$sigpsxname, 0)
})

test_that("criterion defaults follow 'interval'", {
  expect_equal(resolve_criterion(NULL, TRUE), "hpd")
  expect_equal(resolve_criterion(NULL, FALSE), "pvalue")
  expect_equal(resolve_criterion("threshold", TRUE), "threshold")
})

test_that("loadings: 1-coded always kept, -1 coded only if significant", {
  IDY0 <- cbind(c(9, 1, -1, -1, 0), c(-1, 0, 9, 1, 1))
  # free loadings in column-major order: (2,1) (3,1) (4,1) (1,2) (4,2) (5,2)
  res <- list(EmLY = c(.8, .3, .02, .5, .9, .7), SELY = rep(.1, 6),
              PLY = c(0, .001, .4, 0, 0, 0))
  hpd <- list(HPD_LY = rbind(c(.6, 1), c(.1, .5), c(-.2, .2), c(.3, .7),
                             c(.7, 1.1), c(.5, .9)))
  s <- sig_ly_fun(res, hpd, IDY0, paste0("y", 1:5), c("f1", "f2"), "hpd")
  expect_equal(rownames(s$OUTLY), c("f1 by y2", "f1 by y3", "f1 by y4",
                                    "f2 by y1", "f2 by y4", "f2 by y5"))
  expect_equal(unname(s$keep), c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE))
  expect_equal(rownames(s$SIGLY), c("f1 by y3", "f2 by y1"))
  ind <- retained_indicators(IDY0, paste0("y", 1:5), s)
  expect_equal(ind, list(c("y1", "y2", "y3"), c("y3", "y1", "y4", "y5")))
})

test_that("vectorised EPSR equals sna::potscalered.mcmc", {
  skip_if_not_installed("sna")
  set.seed(3)
  x1 <- matrix(rnorm(201 * 3), 201); x2 <- matrix(rnorm(201 * 3, .1, 1.2), 201)
  e <- epsr_running(x1, x2)
  ref <- sapply(1:3, function(p) sapply(2:201, function(b)
    sna::potscalered.mcmc(cbind(x1[1:b, p], x2[1:b, p]))))
  expect_equal(e, ref, tolerance = 1e-10)
})

test_that("EPSR monitors the full lower triangle of PSX", {
  ch <- function(s) {
    set.seed(s)
    list(ELY = matrix(rnorm(300), 150), EMU = matrix(rnorm(300), 150),
         EPHI = matrix(rnorm(150), 150), chainpsx = matrix(rnorm(150 * 6), 150))
  }
  e <- caculate_epsr(120, ch(1), ch(2))
  expect_equal(ncol(e$epsr), 2 + 2 + 1 + 6)
  expect_equal(nrow(e$epsr), 120)
  expect_type(e$convergence, "logical")
})

test_that("blcfa_psx loadings follow gibbs_psx_fun's item-by-item order", {
  IDY0 <- cbind(c(9, 0, 1, 0, 1), c(0, 9, 0, 1, 1))
  # free: (3,1) (4,2) (5,1) (5,2) item by item
  res <- list(EmLY = 1:4 / 10, SELY = rep(.01, 4), PLY = rep(0, 4))
  hpd <- list(HPD_LY = cbind(1:4 / 10 - .01, 1:4 / 10 + .01))
  s <- sig_ly_fun(res, hpd, IDY0, paste0("y", 1:5), c("f1", "f2"), "hpd",
                  keep_all = TRUE)
  expect_equal(rownames(s$OUTLY), c("f1 by y3", "f2 by y4", "f1 by y5", "f2 by y5"))
  expect_equal(unname(s$OUTLY[, "est"]), 1:4 / 10)
  s2 <- sig_ly_fun(res, hpd, IDY0, paste0("y", 1:5), c("f1", "f2"), "hpd")
  expect_equal(rownames(s2$OUTLY), c("f1 by y3", "f1 by y5", "f2 by y4", "f2 by y5"))
})

test_that("EPSR window of NaN is not 'converged'", {
  ch <- list(ELY = matrix(1, 150, 1), EMU = matrix(1, 150, 1),
             EPHI = matrix(1, 150, 1), chainpsx = matrix(1, 150, 1))
  expect_false(caculate_epsr(120, ch, ch)$convergence)
})
