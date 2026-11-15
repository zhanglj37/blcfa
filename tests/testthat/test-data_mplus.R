simu <- system.file("extdata", "simu_data.txt", package = "blcfa")
if (!nzchar(simu)) simu <- testthat::test_path("../../inst/extdata/simu_data.txt")

test_that("read_data validates names and columns", {
  v <- paste0("y", 1:10)
  d <- read_data(simu, v, c("y3", "y1"))
  expect_equal(names(d), c("y3", "y1"))
  expect_error(read_data(simu, v[1:9], v[1:9]), "first line")
  expect_error(read_data(simu, v, c("y1", "z")), "not found in 'varnames'")
  expect_error(read_data("nofile.txt", v, v), "not found")
})

test_that("mark_na handles numeric codes together with existing NA", {
  d <- data.frame(a = c(1, -9, NA), b = c(-9, 2, 3))
  out <- mark_na(d, -9)
  expect_equal(is.na(as.matrix(out)), cbind(a = c(FALSE, TRUE, TRUE), b = c(TRUE, FALSE, FALSE)))
  expect_identical(mark_na(d, NA), d)
})

test_that("Mplus name lists respect the 90-character limit and keep all names", {
  vn <- c("gender", paste0("item", 1:40))
  lines <- mplus_wrap("    NAMES = ", vn)
  expect_true(all(nchar(lines) <= 90))
  toks <- strsplit(sub(";$", "", sub("^ *NAMES = ", "", paste(lines, collapse = " "))), " +")[[1]]
  expect_equal(toks[toks != ""], vn)
  expect_match(lines[length(lines)], ";$")
})

test_that("USEV is complete when length(varnames) is a multiple of 10 (old bug)", {
  vn <- paste0("v", 1:20); us <- paste0("v", 1:17)
  lines <- mplus_wrap("    USEV = ", us)
  toks <- strsplit(sub(";$", "", sub("^ *USEV = ", "", paste(lines, collapse = " "))), " +")[[1]]
  expect_equal(toks[toks != ""], us)
})

test_that("long data paths are copied next to the Mplus input", {
  long_dir <- file.path(tempdir(), paste(rep("very long directory name", 4), collapse = " "))
  dir.create(long_dir, showWarnings = FALSE)
  src <- file.path(long_dir, "simu_data.txt")
  file.copy(simu, src, overwrite = TRUE)
  out <- file.path(tempdir(), "out_long"); dir.create(out, showWarnings = FALSE)
  expect_message(line <- mplus_data_line(src, out), "copied")
  expect_equal(line, "DATA: FILE = simu_data.txt;")
  expect_true(file.exists(file.path(out, "simu_data.txt")))
  # file in the output directory -> plain file name
  expect_equal(mplus_data_line(file.path(out, "simu_data.txt"), out),
               "DATA: FILE = simu_data.txt;")
})

test_that("second-step input files have the paper's structure", {
  od <- file.path(tempdir(), "inp"); dir.create(od, showWarnings = FALSE)
  old <- setwd(od); on.exit(setwd(old))
  v <- paste0("y", 1:10)
  data <- list(line = "DATA: FILE = simu_data.txt;", names = v)
  ind <- list(paste0("y", 1:5), paste0("y", 6:10))
  sp <- c("y6 with y1", "y7 with y2", "y10 with y9")
  write_mplus_second("ml", data, v, c("f1", "f2"), ind, sp)
  write_mplus_second("bayes", data, v, c("f1", "f2"), ind, sp)
  ml <- readLines("blcfa_ml.inp"); by <- readLines("blcfa_bayes.inp")
  expect_true(all(nchar(c(ml, by)) <= 90))
  expect_true(any(grepl("^ +f1 by y1 y2 y3 y4 y5;$", ml)))
  expect_true(any(grepl("^ +y10 with y9;$", ml)))
  expect_true(any(grepl("ESTIMATOR = ML;", ml)))
  expect_true(any(grepl("OUTPUT: TECH1 STDYX;", ml)))
  expect_true(all(c("    ESTIMATOR = BAYES;", "    ALGORITHM = GIBBS(RW);",
                    "    PROC = 2;", "    BITERATIONS = (10000);",
                    "OUTPUT: TECH1 TECH8 STDYX;", "PLOT: TYPE = PLOT2;") %in% by))
})

test_that("find_mplus honours an explicit path and errors on a wrong one", {
  expect_error(find_mplus("/no/such/mplus"), "not found")
  f <- tempfile(); file.create(f)
  expect_equal(find_mplus(f), normalizePath(f, winslash = "/"))
})
