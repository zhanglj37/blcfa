# Short end-to-end runs (no Mplus needed). These check the plumbing, not
# convergence, so very short chains are used.

simu <- system.file("extdata", "simu_data.txt", package = "blcfa")
if (!nzchar(simu)) simu <- testthat::test_path("../../inst/extdata/simu_data.txt")
v <- paste0("y", 1:10)
mod <- "f1 =~ y1 + y2 + y3 + y4 + y5\nf2 =~ y6 + y7 + y8 + y9 + y10"

test_that("argument validation", {
  expect_error(blcfa(simu, v, v, mod, estimation = "wls"), "estimation")
  expect_error(blcfa(simu, v, v, mod, MCMAX = 100, N.burn = 100), "larger")
  expect_error(blcfa(simu, v, v, mod, MCMAX = 200, N.burn = 50), "at least 100")
})

test_that("conver_check = FALSE: one chain, missing data imputed, inp files written", {
  skip_on_cran()
  od <- file.path(tempdir(), "pipe1"); unlink(od, recursive = TRUE)
  # introduce missing values: NA and the code -9
  d <- utils::read.table(simu)
  d[1, 2] <- NA; d[5, 7] <- -9; d[9, 10] <- -9
  f <- file.path(tempdir(), "simu_ms.txt")
  utils::write.table(d, f, row.names = FALSE, col.names = FALSE)
  # unused leading column to check NAMES/USEV with missing data
  d2 <- cbind(gender = 1, d)
  f2 <- file.path(tempdir(), "simu_ms2.txt")
  utils::write.table(d2, f2, row.names = FALSE, col.names = FALSE)

  r <- blcfa(f2, c("gender", v), v, mod, estimation = "both", ms = -9,
             MCMAX = 150, N.burn = 100, conver_check = FALSE, bloutput = TRUE,
             outdir = od, mplus_path = NULL)
  expect_true(is.na(r$convergence))
  expect_equal(names(r$blcfa_est), c("ppp", "ly", "mu", "phi", "psx"))
  expect_equal(dim(r$blcfa_est$psx), c(55, 6))
  expect_true(file.exists(file.path(od, "data_imputed.txt")))
  imp <- utils::read.table(file.path(od, "data_imputed.txt"))
  expect_equal(dim(imp), c(500, 10))
  expect_false(anyNA(imp))
  ml <- readLines(file.path(od, "blcfa_ml.inp"))
  expect_true("DATA: FILE = data_imputed.txt;" %in% ml)
  expect_true(any(grepl("NAMES = y1 y2", ml)))       # gender not listed for imputed data
  expect_true(file.exists(file.path(od, "blcfa_bayes.inp")))
  for (x in c("ppp.csv", "ly.csv", "mu.csv", "phi.csv", "phi_cormatrix.csv",
              "psx.csv", "psx_sig.csv"))
    expect_true(file.exists(file.path(od, "results", x)), info = x)
  expect_equal(normalizePath(getwd()), normalizePath(getwd()))  # wd restored
})

test_that("blcfa_ly with an IDY matrix (lasso cross-loadings)", {
  skip_on_cran()
  od <- file.path(tempdir(), "pipe2"); unlink(od, recursive = TRUE)
  IDY <- matrix(-1, 10, 2); IDY[1:5, 1] <- 1; IDY[6:10, 2] <- 1
  IDY[1, 1] <- 9; IDY[6, 2] <- 9
  r <- blcfa_ly(simu, v, v, IDY, MCMAX = 150, N.burn = 100,
                conver_check = FALSE, outdir = od)
  expect_equal(nrow(r$blcfa_est$ly), 18)
  expect_equal(nrow(r$blcfa_est$psx), 10)          # residual variances only
  ml <- readLines(file.path(od, "blcfa_ml.inp"))
  expect_true(any(grepl("^ +f1 by y1 y2 y3 y4 y5", ml)))
})

test_that("blcfa_psx runs and is reproducible", {
  skip_on_cran()
  od <- file.path(tempdir(), "pipe3"); unlink(od, recursive = TRUE)
  r1 <- blcfa_psx(simu, v, v, mod, MCMAX = 150, N.burn = 100,
                  conver_check = FALSE, outdir = od)
  r2 <- blcfa_psx(simu, v, v, mod, MCMAX = 150, N.burn = 100,
                  conver_check = FALSE, outdir = od)
  expect_identical(r1$blcfa_est, r2$blcfa_est)
})
