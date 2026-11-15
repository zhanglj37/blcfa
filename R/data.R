# Data input: read the raw ASCII file, select usevar, mark missing values,
# standardise.

#' Read the analysis data
#'
#' @param filename Path to a whitespace-delimited ASCII data file without a
#'   header row.
#' @param varnames Names of all columns in the file.
#' @param usevar Names of the variables used in the analysis.
#' @return A data.frame with the \code{usevar} columns in \code{usevar} order.
#' @keywords internal
read_data <- function(filename, varnames, usevar) {
  if (!is.character(filename) || length(filename) != 1 || !file.exists(filename)) {
    stop("Data file not found: ", filename, call. = FALSE)
  }
  if (anyDuplicated(varnames)) {
    stop("Duplicated names in 'varnames': ",
         paste(unique(varnames[duplicated(varnames)]), collapse = ", "),
         call. = FALSE)
  }
  if (anyDuplicated(usevar)) {
    stop("Duplicated names in 'usevar'.", call. = FALSE)
  }
  missing_vars <- setdiff(usevar, varnames)
  if (length(missing_vars) > 0) {
    stop("Variable(s) in 'usevar' not found in 'varnames': ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  first <- scan(filename, what = "", nlines = 1, quiet = TRUE)
  if (length(first) != length(varnames)) {
    stop("The first line of the data file has ", length(first),
         " values but 'varnames' has ", length(varnames), " names. ",
         "The data file must not contain a header row.", call. = FALSE)
  }

  dataset_origin <- utils::read.table(file = filename, col.names = varnames,
                                      check.names = FALSE)
  dataset <- dataset_origin[usevar]
  not_num <- usevar[!vapply(dataset, is.numeric, logical(1))]
  if (length(not_num) > 0) {
    stop("Non-numeric values found in variable(s): ",
         paste(not_num, collapse = ", "),
         ". Missing values must be NA or a numeric code given in 'ms'.",
         call. = FALSE)
  }
  dataset
}

# Replace the missing-value code `ms` by NA (NA in the file is kept as NA).
mark_na <- function(dataset, ms) {
  ms <- check_ms(ms)
  if (is.na(ms)) return(dataset)
  dataset[!is.na(dataset) & dataset == ms] <- NA
  dataset
}

# Standardise each variable (mean 0, sd 1, computed from the observed values)
# and return the NY x N matrix expected by the Gibbs samplers.
standardize_data <- function(dataset) {
  t(scale(as.matrix(dataset)))
}

# 'ms' must be NA or a single number (a numeric string such as "-99" is
# accepted and converted).
check_ms <- function(ms) {
  if (length(ms) != 1) stop("'ms' must be a single value.", call. = FALSE)
  if (is.na(ms)) return(NA)
  if (is.character(ms)) {
    v <- suppressWarnings(as.numeric(ms))
    if (is.na(v)) stop("'ms' must be numeric (or NA).", call. = FALSE)
    return(v)
  }
  if (!is.numeric(ms)) stop("'ms' must be numeric (or NA).", call. = FALSE)
  ms
}
