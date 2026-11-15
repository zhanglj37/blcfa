u <- paste0("y", 1:10)

test_that("paper syntax is parsed into the expected IDY matrix", {
  m <- parse_model("
f1 =~ y1 + y2 + y3 + y4 + y5
f2 =~ y6 + y7 + y8 + y9 + y10
", u)
  expect_equal(m$factor_names, c("f1", "f2"))
  expect_equal(m$IDY0[, 1], c(9, 1, 1, 1, 1, 0, 0, 0, 0, 0))
  expect_equal(m$IDY0[, 2], c(0, 0, 0, 0, 0, 9, 1, 1, 1, 1))
  expect_equal(m$indicators[[2]], paste0("y", 6:10))
})

test_that("compact, commented, CRLF and multi-line syntax give the same model", {
  ref <- parse_model("f1 =~ y1 + y2 + y3\nf2 =~ y4 + y5 + y6", u[1:6])
  variants <- c(
    "f1=~y1+y2+y3\nf2=~y4+y5+y6",
    "# 1. CFA\r\nf1 =~ y1 + y2 + y3 # comment\r\n\r\nf2 =~ y4 + y5 + y6\r\n",
    "f1 =~ y1 +\n  y2\n  + y3\n\tf2 =~\n y4 + y5 + y6",
    "f1 =~ y1 + y2 + y3; f2 =~ y4 + y5 + y6 ! Mplus-style comment",
    "f1 =~ y1 + y2\nf2 =~ y4 + y5 + y6\nf1 =~ y3"
  )
  for (v in variants) {
    m <- parse_model(v, u[1:6])
    expect_equal(m$IDY0, ref$IDY0, info = v)
    expect_equal(m$factor_names, ref$factor_names, info = v)
  }
})

test_that("user factor names are kept and cross-loadings allowed", {
  m <- parse_model("PAP =~ y1 + y2 + y3\nPAB =~ y4 + y5 + y3", u[1:5])
  expect_equal(m$factor_names, c("PAP", "PAB"))
  expect_equal(m$IDY0[3, ], c(1, 1))
})

test_that("invalid syntax gives informative errors", {
  expect_error(parse_model("f1 =~ y1 + y2 + z3", u), "not found in 'usevar'")
  expect_error(parse_model("f1 =~ y1", u), "fewer than 2")
  expect_error(parse_model("f1 =~ y1 + y2\ny1 ~~ y2", u), "~~")
  expect_error(parse_model("f1 =~ 1*y1 + y2", u), "modifiers")
  expect_error(parse_model("f1 =~ y1 + + y2", u), "Empty indicator")
  expect_error(parse_model("f1 =~ y1 + y2 + y1", u), "more than once")
  expect_error(parse_model("y1 =~ y2 + y3", u), "also observed")
  expect_error(parse_model("", u), "No factor")
})

test_that("IDY matrix input is validated", {
  IDY <- matrix(c(9, 1, 1, -1, -1, -1, -1, 9, 1, 1), 5, 2)
  m <- parse_model(IDY, u[1:5])
  expect_equal(m$source, "matrix")
  expect_equal(m$factor_names, c("f1", "f2"))
  expect_equal(m$indicators[[2]], c("y3", "y1", "y2", "y4", "y5"))
  expect_error(parse_model(IDY, u[1:4]), "rows")
  bad <- IDY; bad[1, 1] <- 2
  expect_error(parse_model(bad, u[1:5]), "entries must be")
  bad <- IDY; bad[2, 1] <- 9
  expect_error(parse_model(bad, u[1:5]), "exactly one 9")
})

test_that("unused variables trigger a warning", {
  expect_warning(parse_model("f1 =~ y1 + y2", u[1:3]), "y3")
})

test_that("IDY row names are matched to usevar; case-insensitive duplicates rejected", {
  IDY <- matrix(c(9, 1, 1, 0, 0, 0, 0, 0, 9, 1), 5, 2,
                dimnames = list(c("y5", "y4", "y3", "y2", "y1"), NULL))
  m <- parse_model(IDY, paste0("y", 1:5))
  expect_equal(m$IDY0[, 1], c(0, 0, 1, 1, 9))
  rownames(IDY)[1] <- "z"
  expect_error(parse_model(IDY, paste0("y", 1:5)), "row names")
  expect_error(parse_model("F1 =~ y1 + y2\nf1 =~ y3 + y4", paste0("y", 1:4)), "case")
})

test_that("ms accepts numbers and numeric strings", {
  d <- data.frame(a = c(1, -99))
  expect_true(is.na(mark_na(d, "-99")$a[2]))
  expect_error(mark_na(d, "x"), "numeric")
  expect_equal(mplus_missing_line("-99"), "    MISSING = ALL(-99);")
})
