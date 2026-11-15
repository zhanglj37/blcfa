# Model specification: lavaan-style CFA syntax or an IDY loading-pattern matrix.
#
# Both input forms are converted to the same internal representation:
#   list(IDY0 = NY x NZ matrix with codes 9 / 1 / -1 / 0 (rows follow usevar),
#        factor_names = character(NZ),
#        indicators = named list, one character vector of indicators per factor
#                     (in the order given by the user; first = marker),
#        source = "syntax" | "matrix")
#
# IDY0 codes: 9 = loading fixed at 1 (marker), 1 = free loading without shrinkage,
#            -1 = free loading with lasso shrinkage, 0 = loading fixed at 0.

#' Parse a blcfa model specification
#'
#' @param myModel Either a character string in lavaan-style syntax
#'   (\code{f1 =~ y1 + y2 + y3}) or a numeric loading-pattern matrix (IDY) with
#'   one row per variable in \code{usevar} and one column per factor.
#' @param usevar Character vector of the variables used in the analysis.
#' @return A list with elements \code{IDY0}, \code{factor_names},
#'   \code{indicators} and \code{source}.
#' @keywords internal
parse_model <- function(myModel, usevar) {
  if (is.matrix(myModel) || is.data.frame(myModel)) {
    return(parse_model_matrix(as.matrix(myModel), usevar))
  }
  if (!is.character(myModel) || length(myModel) < 1) {
    stop("'myModel' must be a model string (e.g. \"f1 =~ y1 + y2 + y3\") ",
         "or an IDY loading-pattern matrix.", call. = FALSE)
  }
  parse_model_syntax(paste(myModel, collapse = "\n"), usevar)
}

# ---- syntax ------------------------------------------------------------------

parse_model_syntax <- function(txt, usevar) {
  statements <- split_statements(txt)
  if (length(statements) == 0) {
    stop("No factor definition ('=~') found in 'myModel'.", call. = FALSE)
  }

  indicators <- list()
  for (st in statements) {
    if (grepl("~~", st, fixed = TRUE)) {
      stop("Residual covariances ('~~') cannot be specified in 'myModel'; ",
           "they are detected by blcfa. Offending statement: '", st, "'",
           call. = FALSE)
    }
    parts <- strsplit(st, "=~", fixed = TRUE)[[1]]
    if (length(parts) != 2) {
      stop("Cannot parse statement '", st, "'. Each statement must have the ",
           "form 'factor =~ item1 + item2 + ...'.", call. = FALSE)
    }
    lhs <- trimws(parts[1])
    rhs <- trimws(strsplit(parts[2], "+", fixed = TRUE)[[1]])

    if (grepl("~", parts[2], fixed = TRUE) || grepl("~", lhs, fixed = TRUE)) {
      stop("Only measurement statements ('=~') are supported. Offending ",
           "statement: '", st, "'", call. = FALSE)
    }
    if (!is_valid_name(lhs)) {
      stop("Invalid factor name '", lhs, "' in statement '", st, "'.",
           call. = FALSE)
    }
    if (any(rhs == "")) {
      stop("Empty indicator in statement '", st, "' (check the '+' signs).",
           call. = FALSE)
    }
    if (any(grepl("*", rhs, fixed = TRUE))) {
      stop("Parameter modifiers such as '1*y1' or 'NA*y1' are not supported. ",
           "The first indicator of each factor is fixed at 1. Offending ",
           "statement: '", st, "'", call. = FALSE)
    }
    bad <- rhs[!vapply(rhs, is_valid_name, logical(1))]
    if (length(bad) > 0) {
      stop("Invalid indicator name(s) ", paste0("'", bad, "'", collapse = ", "),
           " in statement '", st, "'.", call. = FALSE)
    }
    # lavaan semantics: repeated definitions of the same factor are merged
    indicators[[lhs]] <- c(indicators[[lhs]], rhs)
  }

  factor_names <- names(indicators)
  if (anyDuplicated(tolower(factor_names))) {
    stop("Factor names must differ by more than upper/lower case (Mplus is ",
         "not case-sensitive).", call. = FALSE)
  }
  for (f in factor_names) {
    dup <- unique(indicators[[f]][duplicated(indicators[[f]])])
    if (length(dup) > 0) {
      stop("Indicator(s) ", paste(dup, collapse = ", "),
           " listed more than once for factor '", f, "'.", call. = FALSE)
    }
    if (length(indicators[[f]]) < 2) {
      stop("Factor '", f, "' has fewer than 2 indicators; at least 2 are ",
           "required.", call. = FALSE)
    }
  }
  clash <- factor_names[tolower(factor_names) %in% tolower(usevar)]
  if (length(clash) > 0) {
    stop("Factor name(s) ", paste(clash, collapse = ", "),
         " are also observed variable names.", call. = FALSE)
  }
  unknown <- setdiff(unlist(indicators), usevar)
  if (length(unknown) > 0) {
    stop("Indicator(s) not found in 'usevar': ",
         paste(unique(unknown), collapse = ", "), call. = FALSE)
  }

  IDY0 <- matrix(0, nrow = length(usevar), ncol = length(factor_names))
  for (k in seq_along(factor_names)) {
    loc <- match(indicators[[k]], usevar)
    IDY0[loc[-1], k] <- 1
    IDY0[loc[1], k] <- 9
  }
  warn_unused_items(IDY0, usevar)

  list(IDY0 = IDY0, factor_names = factor_names,
       indicators = unname(indicators), source = "syntax")
}

# Split model text into statements.
# - comments start with '#' or '!' and run to the end of the line
# - statements are separated by newlines or ';'
# - a line that starts with '+' or follows a line ending in '+' / '=~'
#   continues the previous statement
split_statements <- function(txt) {
  txt <- gsub("\r\n?", "\n", txt)
  txt <- gsub("\t", " ", txt, fixed = TRUE)
  lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  lines <- sub("[#!].*$", "", lines)
  lines <- unlist(strsplit(lines, ";", fixed = TRUE))
  lines <- trimws(lines)
  lines <- lines[lines != ""]

  out <- character(0)
  for (ln in lines) {
    n <- length(out)
    continues <- n > 0 && (startsWith(ln, "+") ||
                             grepl("(\\+|=~)$", out[n]))
    if (continues) {
      out[n] <- paste(out[n], ln)
    } else {
      out <- c(out, ln)
    }
  }
  # normalise spacing (purely cosmetic, used in error messages)
  gsub("\\s+", " ", out)
}

is_valid_name <- function(x) {
  grepl("^[A-Za-z][A-Za-z0-9_.]*$", x)
}

# ---- matrix ------------------------------------------------------------------

parse_model_matrix <- function(IDY0, usevar) {
  if (!is.numeric(IDY0)) {
    stop("The IDY matrix must be numeric.", call. = FALSE)
  }
  rn <- rownames(IDY0)
  if (!is.null(rn)) {
    if (!setequal(rn, usevar) || anyDuplicated(rn)) {
      stop("The row names of the IDY matrix must be the variables in 'usevar'.",
           call. = FALSE)
    }
    IDY0 <- IDY0[match(usevar, rn), , drop = FALSE]
  }
  if (nrow(IDY0) != length(usevar)) {
    stop("The IDY matrix has ", nrow(IDY0), " rows but 'usevar' has ",
         length(usevar), " variables; rows must correspond to 'usevar'.",
         call. = FALSE)
  }
  if (anyNA(IDY0) || !all(IDY0 %in% c(9, 1, -1, 0))) {
    stop("IDY entries must be 9 (fixed at 1), 1 (free), -1 (free with lasso ",
         "shrinkage) or 0 (fixed at 0).", call. = FALSE)
  }
  n9 <- colSums(IDY0 == 9)
  if (any(n9 != 1)) {
    stop("Each column (factor) of the IDY matrix must contain exactly one 9 ",
         "(the marker indicator). Problem in column(s): ",
         paste(which(n9 != 1), collapse = ", "), call. = FALSE)
  }
  NZ <- ncol(IDY0)
  cn <- colnames(IDY0)
  factor_names <- if (!is.null(cn) && all(vapply(cn, is_valid_name, logical(1))) &&
                      !anyDuplicated(cn)) cn else paste0("f", seq_len(NZ))
  indicators <- lapply(seq_len(NZ), function(k) {
    c(usevar[IDY0[, k] == 9], usevar[IDY0[, k] %in% c(1, -1)])
  })
  IDY0 <- unname(IDY0)
  storage.mode(IDY0) <- "double"
  warn_unused_items(IDY0, usevar)
  list(IDY0 = IDY0, factor_names = factor_names,
       indicators = indicators, source = "matrix")
}

warn_unused_items <- function(IDY0, usevar) {
  unused <- usevar[rowSums(IDY0 != 0) == 0]
  if (length(unused) > 0) {
    warning("Variable(s) in 'usevar' not loading on any factor: ",
            paste(unused, collapse = ", "), call. = FALSE)
  }
}
