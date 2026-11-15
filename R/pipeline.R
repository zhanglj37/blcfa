# Shared pipeline of blcfa(), blcfa_ly() and blcfa_psx():
#   1. read data + model          4. posterior summaries, significance
#   2. run MCMC chain(s)          5. results tables / files
#   3. convergence check (EPSR)   6. second step: Mplus input files (+ run)
#
# The Gibbs samplers (gibbs_fun, gibbs_ly_fun, gibbs_psx_fun) and the initial
# values (set_ly_int) are used unchanged.

normalize_estimation <- function(estimation) {
  e <- tolower(as.character(estimation))
  if (length(e) != 1 || !e %in% c("ml", "bayes", "bayesian", "both")) {
    stop("'estimation' must be \"ml\", \"bayes\" or \"both\".", call. = FALSE)
  }
  if (e == "bayesian") "bayes" else e
}

# One MCMC chain. Runs inside a worker process when chains run in parallel,
# so it sets the working directory itself (the samplers write log.txt there).
run_one_chain <- function(init_cir, log_cir, sampler, wd, MCMAX, N.burn, NZ, NY,
                          N, Y, IDY0) {
  old_wd <- setwd(wd)
  on.exit(setwd(old_wd), add = TRUE)
  sink("log.txt", append = TRUE)
  on.exit(sink(), add = TRUE)

  LY_int <- set_ly_int(init_cir, IDY0)
  IDY <- IDY0
  IDY[which(IDY0 == 9)] <- 0
  # gibbs_fun and gibbs_ly_fun call set.seed(1) themselves; gibbs_psx_fun
  # does not, so seed here to make every sampler reproducible
  set.seed(1)
  sampler_fun <- get(sampler, mode = "function")
  chain <- sampler_fun(MCMAX, NZ, NY, N, Y, LY_int, IDY0, IDY,
                       1, N.burn, log_cir)

  # imputation: posterior mean of the stored draws; the full NY x N x MCMAX
  # array is dropped to save memory
  keep <- (N.burn + 1):MCMAX
  chain$Y_missing_mean <- if (sum(chain$missing_ind) > 0) {
    apply(chain$Y_missing[, , keep, drop = FALSE], c(1, 2), mean)
  } else NULL
  chain$Y_missing <- NULL
  chain
}

# Can PSOCK workers load the same blcfa that is running here?
# (Not the case under devtools::load_all().)
workers_can_load_blcfa <- function() {
  ns_path <- tryCatch(getNamespaceInfo("blcfa", "path"), error = function(e) "")
  inst <- tryCatch(find.package("blcfa", lib.loc = .libPaths(), quiet = TRUE),
                   error = function(e) character(0))
  length(inst) == 1 && nzchar(ns_path) &&
    identical(normalizePath(ns_path, winslash = "/"),
              normalizePath(inst, winslash = "/")) &&
    file.exists(file.path(inst, "Meta", "package.rds"))
}

run_chains <- function(chains, sampler, wd, MCMAX, N.burn, NZ, NY, N, Y, IDY0) {
  args <- list(sampler = sampler, wd = wd, MCMAX = MCMAX, N.burn = N.burn,
               NZ = NZ, NY = NY, N = N, Y = Y, IDY0 = IDY0)
  call_chain <- function(ch) {
    do.call(run_one_chain, c(list(init_cir = ch$init, log_cir = ch$log), args))
  }
  use_parallel <- length(chains) > 1 &&
    isTRUE(parallel::detectCores() - 1 > 1) && workers_can_load_blcfa()

  if (use_parallel) {
    worker_fun <- function(ch, args) {
      do.call(getFromNamespace("run_one_chain", "blcfa"),
              c(list(init_cir = ch$init, log_cir = ch$log), args))
    }
    environment(worker_fun) <- asNamespace("utils")   # avoid serialising this frame
    set_lib <- function(p) invisible(.libPaths(p))
    environment(set_lib) <- globalenv()
    load_pkg <- function() requireNamespace("blcfa", quietly = TRUE)
    environment(load_pkg) <- globalenv()
    # only the cluster set-up may fall back to serial computation; errors
    # inside the samplers are passed on
    cl <- tryCatch({
      cl <- parallel::makeCluster(length(chains))
      # workers must find blcfa in the same libraries as this session
      parallel::clusterCall(cl, set_lib, .libPaths())
      ok <- unlist(parallel::clusterCall(cl, load_pkg))
      if (!all(ok)) {
        parallel::stopCluster(cl)
        stop("blcfa cannot be loaded on the workers")
      }
      cl
    }, error = function(e) {
      message("Parallel computation not available (", conditionMessage(e),
              "); running the chains serially.")
      NULL
    })
    if (!is.null(cl)) {
      on.exit(parallel::stopCluster(cl), add = TRUE)
      return(parallel::parLapply(cl, chains, worker_fun, args = args))
    }
  }
  lapply(chains, call_chain)
}

run_blcfa <- function(mode, filename, varnames, usevar, myModel, estimation,
                      ms, MCMAX, N.burn, bloutput, interval, conver_check,
                      criterion, outdir, mplus_path) {
  # ---- arguments ---------------------------------------------------------------
  estimation <- normalize_estimation(estimation)
  criterion <- resolve_criterion(criterion, interval)
  if (length(MCMAX) != 1 || length(N.burn) != 1 ||
      !is.numeric(MCMAX) || !is.numeric(N.burn) ||
      !is.finite(MCMAX) || !is.finite(N.burn) || N.burn < 1) {
    stop("'MCMAX' and 'N.burn' must be positive numbers.", call. = FALSE)
  }
  MCMAX <- as.integer(MCMAX)
  N.burn <- as.integer(N.burn)
  if (MCMAX <= N.burn) stop("'MCMAX' must be larger than 'N.burn'.", call. = FALSE)
  if (conver_check && N.burn < 100) {
    stop("'N.burn' must be at least 100 for the convergence check ",
         "(EPSR of the last 100 burn-in iterations).", call. = FALSE)
  }
  filename <- normalizePath(filename, winslash = "/", mustWork = FALSE)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  outdir <- normalizePath(outdir, winslash = "/")
  old_wd <- setwd(outdir)
  on.exit(setwd(old_wd), add = TRUE)
  # results of an earlier run in the same folder must not be mistaken for
  # this run's results (e.g. by sum_second() after non-convergence)
  stale <- c(outer(c("blcfa_ml", "blcfa_bayes", "normal"),
                   c(".inp", ".out", ".gh5"), paste0),
             "data_imputed.txt", "normal_data.dat", "EPSR.png",
             file.path("results", c("ppp.csv", "ly.csv", "mu.csv", "phi.csv",
                                    "phi_cormatrix.csv", "psx.csv",
                                    "psx_sig.csv", "ly_sig.csv", "EPSR.png")))
  stale <- normalizePath(stale[file.exists(stale)], winslash = "/")
  unlink(setdiff(stale, filename))          # never delete the input data

  # ---- data and model ----------------------------------------------------------
  dataset <- mark_na(read_data(filename, varnames, usevar), ms)
  model <- parse_model(myModel, usevar)
  IDY0 <- model$IDY0
  if (mode == "psx" && any(IDY0 == -1)) {
    stop("blcfa_psx() does not support lasso loadings (-1); use blcfa() or ",
         "blcfa_ly().", call. = FALSE)
  }
  if (ncol(IDY0) == 1 && anyNA(dataset)) {
    stop("Missing values with a one-factor model are not supported by the ",
         "sampler; please remove or impute the missing values first.",
         call. = FALSE)
  }
  # NAMES lists all varnames when the original file is passed to Mplus
  warn_long_names(c(varnames, model$factor_names))
  Y <- standardize_data(dataset)
  N <- nrow(dataset)
  NY <- ncol(dataset)
  NZ <- ncol(IDY0)

  # ---- MCMC --------------------------------------------------------------------
  sampler <- switch(mode, full = "gibbs_fun", ly = "gibbs_ly_fun",
                    psx = "gibbs_psx_fun")
  cat("The program is running. See 'log.txt' for details.  \n")
  writeLines("", "log.txt")
  # chain 2 provides the estimates; chain 1 only enters the EPSR check.
  # Without the check a single chain with chain-2 initial values is run.
  chains <- if (conver_check) {
    list(list(init = 1, log = 1), list(init = 2, log = 2))
  } else {
    list(list(init = 2, log = 1))
  }
  out <- run_chains(chains, sampler, outdir, MCMAX, N.burn, NZ, NY, N, Y, IDY0)
  chain2 <- out[[length(out)]]

  # ---- posterior summaries ------------------------------------------------------
  res <- caculate_results(chain2, MCMAX, N.burn)
  hpd <- hpd_fun(chain2, MCMAX, N.burn)
  sigpsx <- if (mode == "ly") sig_psx_diag(res, hpd, usevar) else
    sig_psx_fun(res, hpd, usevar, criterion)
  sigly <- sig_ly_fun(res, hpd, IDY0, usevar, model$factor_names, criterion,
                      keep_all = (mode == "psx"))

  if (conver_check) {
    epsrlist <- caculate_epsr(N.burn, out[[1]], chain2)
    convergence <- epsrlist$convergence
    epsr <- epsrlist$epsr
  } else {
    convergence <- NA
    epsr <- NULL
  }
  cat("Gibbs sampling ended up, specific results are being calculated.  \n")

  tables <- write_results(res, hpd, sigpsx, sigly, usevar, model$factor_names,
                          epsr, N.burn, bloutput)

  # ---- second step -------------------------------------------------------------
  if (isFALSE(convergence)) {
    cat("Error: Failed to satisfy the convergence criterion. Check the epsr",
        "graph and increase the values of N.burn and MCMAX.  \n")
    EPSR_figure(epsr, N.burn, "EPSR.png")
  } else {
    ismissing <- impute_ms(Y, chain2, file = "data_imputed.txt")
    data <- if (ismissing) {
      list(line = "DATA: FILE = data_imputed.txt;", names = usevar)
    } else {
      list(line = mplus_data_line(filename, outdir), names = varnames)
    }
    indicators <- retained_indicators(IDY0, usevar, sigly)
    mplus_cmd <- find_mplus(mplus_path)
    if (estimation %in% c("bayes", "both")) {
      write_mplus_second("bayes", data, usevar, model$factor_names, indicators,
                         sigpsx$sigpsxname)
      run_mplus("blcfa_bayes.inp", mplus_cmd)
    }
    if (estimation %in% c("ml", "both")) {
      nonnormal <- normal_detect(filename, varnames, usevar, ms, mplus_cmd,
                                 outdir, dataset = dataset)
      write_mplus_second("ml", data, usevar, model$factor_names, indicators,
                         sigpsx$sigpsxname, nonnormal = nonnormal)
      run_mplus("blcfa_ml.inp", mplus_cmd)
    }
  }

  list(blcfa_est = tables, estimation = estimation,
       convergence = convergence, criterion = criterion,
       sig_psx = sigpsx$SIGPSX, sig_ly = sigly$SIGLY, outdir = outdir)
}
