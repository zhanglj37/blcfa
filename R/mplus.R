# Writing and running Mplus input files.
# Mplus limits input lines to 90 characters and variable names to 8 characters.

MPLUS_MAX_LINE <- 90L
MPLUS_INDENT <- "    "

# Wrap `tokens` after `prefix` into lines of at most `width` characters,
# terminated by ";".
mplus_wrap <- function(prefix, tokens, width = MPLUS_MAX_LINE,
                       indent = MPLUS_INDENT) {
  tokens <- c(tokens[-length(tokens)], paste0(tokens[length(tokens)], ";"))
  lines <- character(0)
  cur <- prefix
  first <- TRUE
  for (tk in tokens) {
    cand <- if (first) paste0(cur, tk) else paste(cur, tk)
    if (!first && nchar(cand) > width) {
      lines <- c(lines, cur)
      cur <- paste0(indent, tk)
    } else {
      cur <- cand
    }
    first <- FALSE
  }
  c(lines, cur)
}

warn_long_names <- function(names) {
  long <- names[nchar(names) > 8]
  if (length(long) > 0) {
    warning("Mplus allows at most 8 characters in variable names; these will ",
            "fail in Mplus: ", paste(long, collapse = ", "), call. = FALSE)
  }
}

# DATA: FILE line. The data path is kept if the line fits Mplus' 90-character
# limit; otherwise the (unchanged) data file is copied into the output
# directory and referenced by its file name.
mplus_data_line <- function(datafile, outdir = getwd()) {
  datafile <- normalizePath(datafile, winslash = "/", mustWork = TRUE)
  outdir <- normalizePath(outdir, winslash = "/", mustWork = TRUE)
  make_line <- function(path) {
    if (grepl(" ", path, fixed = TRUE)) path <- paste0("\"", path, "\"")
    paste0("DATA: FILE = ", path, ";")
  }
  if (identical(dirname(datafile), outdir)) {
    return(make_line(basename(datafile)))
  }
  line <- make_line(datafile)
  if (nchar(line) <= MPLUS_MAX_LINE) return(line)

  target <- file.path(outdir, basename(datafile))
  if (file.exists(target) &&
      !identical(unname(tools::md5sum(target)), unname(tools::md5sum(datafile)))) {
    target <- file.path(outdir, paste0("blcfa_data_", basename(datafile)))
  }
  file.copy(datafile, target, overwrite = TRUE)
  message("Note: the data path exceeds Mplus' 90-character line limit; the data ",
          "file was copied to ", target, " for Mplus.")
  line <- make_line(basename(target))
  if (nchar(line) > MPLUS_MAX_LINE) {
    stop("The data file name is too long for Mplus (90-character limit).",
         call. = FALSE)
  }
  line
}

# MODEL statements: "f1 by y1 y2 ...;" and "y6 with y1;"
mplus_model_lines <- function(factor_names, indicators, sigpsxname = character(0)) {
  by <- unlist(lapply(seq_along(factor_names), function(k) {
    mplus_wrap(paste0(MPLUS_INDENT, factor_names[k], " by "), indicators[[k]])
  }))
  with <- if (length(sigpsxname) > 0) paste0(MPLUS_INDENT, sigpsxname, ";") else character(0)
  c(by, with)
}

# Assemble an input file from its parts and write it.
write_inp <- function(file, data_line, names, usev, model_lines,
                      analysis_lines, output_lines, missing_line = NULL,
                      title = "Bayesian Lasso CFA") {
  lines <- c(
    paste0("TITLE: ", title),
    data_line,
    "VARIABLE:",
    mplus_wrap(paste0(MPLUS_INDENT, "NAMES = "), names),
    mplus_wrap(paste0(MPLUS_INDENT, "USEV = "), usev),
    missing_line,
    "ANALYSIS:",
    paste0(MPLUS_INDENT, analysis_lines),
    "MODEL:",
    model_lines,
    "",
    output_lines
  )
  too_long <- lines[nchar(lines) > MPLUS_MAX_LINE]
  if (length(too_long) > 0) {
    warning("Mplus input line(s) longer than 90 characters in ", file, ": ",
            paste(too_long, collapse = " | "), call. = FALSE)
  }
  writeLines(lines, file)
  invisible(file)
}

mplus_missing_line <- function(ms) {
  ms <- check_ms(ms)
  if (!is.na(ms)) {
    paste0(MPLUS_INDENT, "MISSING = ALL(", format(ms, scientific = FALSE), ");")
  } else NULL
}

# ---- locating and running Mplus ----------------------------------------------

#' Locate the Mplus executable
#'
#' @param mplus_path Optional path to the Mplus executable. If \code{NULL},
#'   the option \code{blcfa.mplus_path} and then the system \code{PATH} are
#'   searched.
#' @return The path of the executable, or \code{""} if Mplus is not found.
#' @keywords internal
find_mplus <- function(mplus_path = NULL) {
  if (is.null(mplus_path)) mplus_path <- getOption("blcfa.mplus_path")
  if (!is.null(mplus_path)) {
    if (!file.exists(mplus_path)) {
      stop("Mplus executable not found at '", mplus_path, "'.", call. = FALSE)
    }
    return(normalizePath(mplus_path, winslash = "/"))
  }
  for (cmd in c("mplus", "Mplus")) {
    p <- Sys.which(cmd)
    if (nzchar(p)) return(unname(p))
  }
  ""
}

mplus_not_found_msg <- function() {
  cat("Error: Failed to run the Mplus software, check whether the Mplus has",
      "been installed or the environment variable of your computer includes",
      "the path of Mplus (or set the 'mplus_path' argument).\n")
}

run_mplus <- function(inp, mplus_cmd) {
  if (!nzchar(mplus_cmd)) {
    mplus_not_found_msg()
    return(invisible(FALSE))
  }
  MplusAutomation::runModels(inp, Mplus_command = mplus_cmd)
  invisible(file.exists(sub("\\.inp$", ".out", inp)))
}

# Data file for Mplus runs on the raw (unstandardised) data.
# If the file contains NA (which Mplus cannot read), the analysis variables are
# written to `fallback_file` with "." for missing values (format conversion
# only, values unchanged). Otherwise the original file is used with
# MISSING = ALL(ms) when a numeric missing code is given.
mplus_data_source <- function(filename, varnames, dataset, ms, outdir, fallback_file) {
  # any NA in the file (also in columns not analysed) makes Mplus fail
  raw_has_na <- anyNA(read_data(filename, varnames, varnames))
  if (raw_has_na) {
    utils::write.table(dataset, file.path(outdir, fallback_file), col.names = FALSE,
                       row.names = FALSE, na = ".", sep = "\t")
    list(data_line = paste0("DATA: FILE = ", fallback_file, ";"),
         names = names(dataset), missing_line = paste0(MPLUS_INDENT, "MISSING = .;"))
  } else {
    list(data_line = suppressMessages(mplus_data_line(filename, outdir)),
         names = varnames, missing_line = mplus_missing_line(ms))
  }
}

# ---- non-normality check -----------------------------------------------------

#' Detect non-normality with Mplus (TECH12 univariate skewness and kurtosis)
#'
#' Returns \code{TRUE} if any |skewness| > 2 or |kurtosis| > 7
#' (West, Finch & Curran, 1995). If Mplus is not available, \code{FALSE}
#' (ML) is returned.
#' @keywords internal
normal_detect <- function(filename, varnames, usevar, ms = NA, mplus_cmd = "",
                          outdir = getwd(), dataset = NULL) {
  if (!nzchar(mplus_cmd)) {
    cat("Error: Failed to run Mplus to detect non-normality, choose the",
        "maximum likelihood estimator as default\n")
    return(FALSE)
  }
  if (is.null(dataset)) dataset <- mark_na(read_data(filename, varnames, usevar), ms)
  src <- mplus_data_source(filename, varnames, dataset, ms, outdir, "normal_data.dat")
  write_inp("normal.inp", src$data_line, src$names, usevar,
            model_lines = c(paste0(MPLUS_INDENT, "%OVERALL%"),
                            mplus_wrap(MPLUS_INDENT, usevar)),
            analysis_lines = "TYPE = MIXTURE;",
            output_lines = "OUTPUT: TECH12;",
            missing_line = c(src$missing_line, paste0(MPLUS_INDENT, "CLASSES = C(1);")))
  run_mplus("normal.inp", mplus_cmd)
  tech12 <- tryCatch(MplusAutomation::readModels("normal.out", what = "tech12",
                                                  quiet = TRUE)$tech12,
                     error = function(e) NULL)
  sk <- suppressWarnings(as.numeric(tech12$obsSkewness))
  ku <- suppressWarnings(as.numeric(tech12$obsKurtosis))
  if (length(sk) == 0 || all(is.na(sk))) {
    warning("Could not read skewness/kurtosis from normal.out; using ML.",
            call. = FALSE)
    return(FALSE)
  }
  isTRUE(max(abs(sk), na.rm = TRUE) > 2 || max(abs(ku), na.rm = TRUE) > 7)
}

# ---- second-step input files -------------------------------------------------

# type: "ml" or "bayes". `data` is list(line, names) describing the data file.
write_mplus_second <- function(type, data, usevar, factor_names, indicators,
                               sigpsxname, nonnormal = FALSE, missing_line = NULL) {
  model_lines <- mplus_model_lines(factor_names, indicators, sigpsxname)
  if (type == "ml") {
    file <- "blcfa_ml.inp"
    analysis <- if (nonnormal) "ESTIMATOR = MLM;" else "ESTIMATOR = ML;"
    output <- "OUTPUT: TECH1 STDYX;"
  } else {
    file <- "blcfa_bayes.inp"
    analysis <- c("ESTIMATOR = BAYES;", "ALGORITHM = GIBBS(RW);", "PROC = 2;",
                  "BITERATIONS = (10000);")
    output <- c("OUTPUT: TECH1 TECH8 STDYX;", "PLOT: TYPE = PLOT2;")
  }
  write_inp(file, data$line, data$names, usevar, model_lines, analysis, output,
            missing_line = missing_line)
}

#' Generate and run a traditional CFA in Mplus
#'
#' Checks the normality of the items (Mplus TECH12) and writes \code{ml.inp}
#' with ML (or MLM if non-normal) estimation for the model in \code{myModel},
#' without the Bayesian lasso step.
#'
#' @inheritParams blcfa
#' @return Invisibly, the path of the Mplus input file.
#' @examples
#' \dontrun{
#' filename <- system.file("extdata", "simu_data.txt", package = "blcfa")
#' varnames <- paste0("y", 1:10)
#' myModel <- "
#' f1 =~ y1 + y2 + y3 + y4 + y5
#' f2 =~ y6 + y7 + y8 + y9 + y10
#' "
#' mplus_ml(filename, varnames, varnames, myModel)
#' }
#' @export
mplus_ml <- function(filename, varnames, usevar, myModel, ms = NA,
                     outdir = getwd(), mplus_path = NULL) {
  filename <- normalizePath(filename, winslash = "/", mustWork = TRUE)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  outdir <- normalizePath(outdir, winslash = "/")
  old_wd <- setwd(outdir)
  on.exit(setwd(old_wd), add = TRUE)

  dataset <- mark_na(read_data(filename, varnames, usevar), ms)
  model <- parse_model(myModel, usevar)
  warn_long_names(c(varnames, model$factor_names))
  mplus_cmd <- find_mplus(mplus_path)
  nonnormal <- normal_detect(filename, varnames, usevar, ms, mplus_cmd, outdir,
                             dataset = dataset)

  src <- mplus_data_source(filename, varnames, dataset, ms, outdir, "ml_data.dat")
  data_line <- src$data_line
  if (anyNA(dataset)) data_line <- c(data_line, paste0(MPLUS_INDENT, "LISTWISE = ON;"))
  write_inp("ml.inp", data_line, src$names, usevar,
            model_lines = mplus_model_lines(model$factor_names, model$indicators),
            analysis_lines = if (nonnormal) "ESTIMATOR = MLM;" else "ESTIMATOR = ML;",
            output_lines = "OUTPUT: TECH1 STDY;",
            missing_line = src$missing_line)
  run_mplus("ml.inp", mplus_cmd)
  invisible(file.path(outdir, "ml.inp"))
}
