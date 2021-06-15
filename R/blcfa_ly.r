## process:
## prepare model and data   
## prior + init + data -> posterior (gibbs sampling, parallel compuation)
## caculate_results
## generate_output

blcfa_ly<-function(filename, varnames, usevar, IDY0, estimation = 'ml', ms = -999999, 
	MCMAX = 10000, N.burn = 5000, bloutput = FALSE,  interval = TRUE)
	## MCMAX: Total number of iterations;  N.burn: Discard the previous N.burn iteration sample
	## estimation = 'ml' / 'bayes'
	## bloutput: Output detailed results (xlsx file);
	## interval: Detect significant residual correlation based on HPD interval or p-value
	## category & point: used for category data (under development)
{

	nthin<-1  ## MCMC algorithm sampling interval
	CNUM<-2   ## number of chain

	dataset <- read_data(filename, varnames, usevar)

	### source("read_model.r")
	N <- nrow(dataset)    # Sample size (N)
	NY <- ncol(dataset)	 # Number of items (p)
	NZ <- ncol(IDY0)


	## record ms values as NA for standarizing data
	dataset_noms <- mark_na(N, NY, dataset, ms)
	Y <- read_data2(dataset_noms)  # standarized
		
	
	###  prior + init + data -> posterior (gibbs sampling) ############################

	### source("ind.R")
	### source("EPSR_set_int.R")
	### source("read_data.R")
	### source("Gibbs.R")
	cat("The program is running. See 'log.txt' for details.  \n")
	set.seed(1)
	
	#**************** Parallel computation ********************
	ncores <- 1
	if(detectCores()-1 > 1)	{
	  ncores <- 2
	}

	#switch between %do% (serial) and %dopar% (parallel)
	if (ncores == 1){  #serial
	  `%is_par%` <- `%do%`
	}else if(Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
			Sys.info()["sysname"] == "Darwin" && getRversion() >= "4.0.0" ) {
			
			if(versionInfo()$version < "1.3.1056"){
				`%is_par%` <- `%do%`
			ncores <- 1
			}
			
	}else{  #parallel
	  `%is_par%` <- `%dopar%`
		cl <- makeCluster(ncores)
		registerDoParallel(cores = ncores)

	}

	writeLines(c(""), "log.txt") #create or clear the log file recording the ouput of foreach loop

	parList <- foreach (CIR = 1:CNUM,
	                      .packages = c("MASS", "statmod", "MCMCpack"),
	                      .export = c("IDY_matrix", "set_ly_int", "gibbs_fun")) %is_par%
	{
		## Calculate the epsr value by running two chains with two kind of initial values

		IDMU<-rep(1,NY)  
		LY_int <- set_ly_int(CIR, IDY0)
		IDY = IDY0
		IDY[which(IDY0==9)] = 0

		sink("log.txt", append=TRUE) # divert the output to the log file

		chainlist <- gibbs_ly_fun(MCMAX, NZ, NY, N, Y, LY_int, IDY0, IDY, 
						  nthin, N.burn, CIR)
		sink() #revert output back to the console

		list(chainlist, IDY, IDMU) #return chainlist, IDY, IDMU to parList
	}
	
	if(ncores > 1) stopCluster(cl)

	#************************Stop Parallel******************************


	#***********access the parameter************
	chain1 <- parList[[1]][[1]]
	chain2 <- parList[[2]][[1]]
	IDY <- parList[[1]][[2]]
	IDMU <- parList[[1]][[3]]
	#*******************************************

	#***********caculate_results************
    ### source("caculate_results.r")
    ### source("HPD.R")
	### source("sigpsx.r")
	### source("EPSR_caculate.r")
	### source("write_mplus.r")

	resultlist <- caculate_results(chain2, CNUM, MCMAX, NY, NZ, N.burn, nthin, IDMU, IDY)
	hpdlist <- hpd_fun(chain2,NZ,NY,N,IDY,N.burn,MCMAX)
	#sigpsx_list <- sig_psx_fun(NZ, NY, dataset, resultlist, hpdlist, interval)
	sigpsx_list<-list(SIGPSX=0,OUTPSX=0)
	sigly_list <- sig_ly_fun(dataset, resultlist, hpdlist, IDY, interval)

	epsrlist <- caculate_epsr(MCMAX, N.burn, CNUM, NY, NZ, chain1, chain2)
	convergence = epsrlist$convergence
	epsr = epsrlist$epsr
	cat("Gibbs sampling ended up, specific results are being calculated.  \n")
	
	#***********generate_output ************
	if (convergence)
	{
		if (bloutput)
		{
			write_results(MCMAX,N.burn,NZ,NY,resultlist,hpdlist,sigpsx_list,sigly_list,epsr,usevar,IDMU,IDY,bloutput)
		}
		ismissing <- impute_ms(Y, NY, N, chain2, N.burn, MCMAX)
		estimation = tolower(estimation)
		if (estimation == 'bayes' || estimation == 'bayesian')
		{
			write_mplus_bayes(varnames,usevar,filename,sigpsx_list,sigly_list,IDY0,ismissing,IDY0)
		}else if (estimation == 'both'){
			write_mplus_bayes(varnames,usevar,filename,sigpsx_list,sigly_list,IDY0,ismissing,IDY0)
			write_mplus_ml(varnames,usevar,filename,sigpsx_list,sigly_list,IDY0,ismissing,IDY0)
		}else{
			write_mplus_ml(varnames,usevar,filename,sigpsx_list,sigly_list,IDY0,ismissing,IDY0)
		}
    }else{
		cat('Error: Failed to satisfy the convergence criterion. Check the epsr graph and increase the values of N.burn and MCMAX.  \n')

		EPSR_figure(epsr, N.burn)
    }
	

}
