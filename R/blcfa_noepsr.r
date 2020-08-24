## process:
## prepare model and data   
## prior + init + data -> posterior (gibbs sampling, parallel compuation)
## caculate_results
## generate_output

blcfa_noepsr<-function(filename, varnames, usevar, myModel, estimation = 'ml', ms = -999999, 
	MCMAX = 10000, N.burn = 5000, bloutput = FALSE,  interval = TRUE, Conver_check = FALSE)
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
	if (is.matrix(myModel))
	{
		IDY0 = myModel
		NZ <- ncol(IDY0)
	}else{
		### source("read_model.r")
		mmvarorigin<-read_model(myModel) ## IDY0 = myModel
		mmvar<-mmvarorigin[2:length(mmvarorigin)] # List: includes factors and variables under each factor
		factorname<-mmvarorigin[[1]]   # names of factors
		numw<-length(mmvar)   # num of factors
		mmvar_loc<-cfa_loc(mmvar,dataset)  # location of indicators
		NZ<-numw  # Number of factors (q)
		IDY0<-IDY_matrix(dataset,mmvar,mmvar_loc) 
	}


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
	

	CIR = 2
	writeLines(c(""), "log.txt") #create or clear the log file recording the ouput of foreach loop

		
	IDMU<-rep(1,NY)  
	LY_int <- set_ly_int(CIR, IDY0)
	IDY = IDY0
	IDY[which(IDY0==9)] = 0 

	sink("log.txt", append=TRUE) # divert the output to the log file

	chain2 <- gibbs_fun(MCMAX, NZ, NY, N, Y, LY_int, IDY0, IDY, 
					  nthin, N.burn, CIR)
	sink() #revert output back to the console


	#*******************************************

	#***********caculate_results************
    ### source("caculate_results.r")
    ### source("HPD.R")
	### source("sigpsx.r")
	### source("EPSR_caculate.r")
	### source("write_mplus.r")

	resultlist <- caculate_results(chain2, CNUM, MCMAX, NY, NZ, N.burn, nthin, IDMU, IDY)
	hpdlist <- hpd_fun(chain2, NZ, NY, N, IDY)
	sigpsx_list <- sig_psx_fun(NZ, NY, dataset, resultlist, hpdlist, interval)
	sigly_list <- sig_ly_fun(dataset, resultlist, hpdlist, IDY, interval)

	cat("Gibbs sampling ended up, specific results are being calculated.  \n")
	
	#***********generate_output ************


	epsr = 'noepsr'
	resultlist2<-write_results(MCMAX,N.burn,NZ,NY,resultlist,hpdlist,sigpsx_list,sigly_list,
				epsr,usevar,IDMU,IDY,bloutput)
	
	ismissing <- impute_ms(Y, NY, N, chain2, N.burn, MCMAX)
	estimation = tolower(estimation)
	if (estimation == 'bayes' || estimation == 'bayesian')
	{
		write_mplus_bayes(varnames,usevar,filename,sigpsx_list,sigly_list,IDY0,ismissing)
	}else if (estimation == 'both'){
		write_mplus_bayes(varnames,usevar,filename,sigpsx_list,sigly_list,IDY0,ismissing)
		write_mplus_ml(varnames,usevar,filename,sigpsx_list,sigly_list,IDY0,ismissing)
	}else{
		write_mplus_ml(varnames,usevar,filename,sigpsx_list,sigly_list,IDY0,ismissing)
	}

	blcfa_results = list(blcfa_est = resultlist2, estimation = estimation)
	return(invisible(blcfa_results))
	

}

