			
	
blcfa<-function(filename, varnames, usevar, model, ms = -999, MCMAX = 15000, N.burn = 5000,   
			 bloutput = FALSE,  interval_psx = TRUE)
			## MCMAX: Total number of iterations;  N.burn: Discard the previous N.burn iteration sample
			## bloutput: Output detailed results (xlsx file);  
			## interval_psx: Detect significant residual correlation based on HPD interval or p-value
			## category & point: used for category data (wait for development)
{

	nthin<-1  ## MCMC algorithm sampling interval
	CNUM<-2
	ms<-ms


	### prepare model and data  #######################################################
	### source("read_data.r")
	dataset<-read_data(filename,varnames)
	
	### source("read_model.r")
	mmvarorigin<-read_model(myModel)
	mmvar<-mmvarorigin[2:length(mmvarorigin)] # List: includes factors and variables under each factor
	factorname<-mmvarorigin[[1]]   # names of factors
	numw<-length(mmvar)   # num of factors
	mmvar_loc<-cfa_loc(mmvar,dataset)  # location of indicators
	N<-nrow(dataset)    # Sample size (N)          
	NY<-ncol(dataset)	 # Number of items (p)      
	NZ<-numw  # Number of factors (q)
	
	###  prior + init + data = posterior (gibbs sampling) ##############################################
	set.seed(1)	
	for(CIR in 1:CNUM){
		## Calculate the epsr value by running two chains with two kind of initial values

		### source("ind.R")   
		IDY<-IDY_matrix_fun(dataset,mmvar,mmvar_loc)
			#vector indicating whether MU should be added to each measurement equation or not. 
			#Added(1), remove(0).
		IDMU<-rep(1,NY)  # now MU is estimated
		IDMUA<-any(as.logical(IDMU))
		
	    
	    ### source("EPSR_1_set_int.R")
		LY_int<-set_int_fun(CIR,dataset,mmvar,mmvar_loc)
		
		### source("read_observed.R")  
		Y<-read_dataset(dataset)
		
	    
		cat(paste("num of chain: ",CIR,"\n"))
		
		chainlist<-gibbs_fun(MCMAX,NZ,NY,N,Y,LY_int,IDMU,IDMUA,IDY,
							  nthin, mmvar, mmvar_loc, N.burn, ms) 
	    assign(paste0("chain",CIR),chainlist)
				
	}#end of GIBs
	chain1<-get(paste0("chain",1))
	chain2<-get(paste0("chain",2))
	resultlist<-caculate_results(chain2,CNUM,MCMAX,NY,NZ,N.burn,nthin,IDMU,IDY)
	hpdlist<-hpd_fun(chain2,NZ,NY,N,IDY)
	sigpsx_list<-sig_psx_fun(NZ,NY,dataset,resultlist,hpdlist,interval_psx)
	

	epsrlist<-caculate_epsr(MCMAX,N.burn,CNUM,NY,NZ,chain1,chain2)
	convergence=epsrlist$convergence
	epsr=epsrlist$epsr
	
	if (convergence)
	{ 
		write_mplus(varnames,usevar,myModel,filename,sigpsx_list)
		if (bloutput)
		{
			write_results(MCMAX,NZ,NY,NLY,resultlist,hpdlist,
							sigpsx_list,epsr,mmvar,factorname,IDMU,IDY)  
		}
    }else{
		cat('Error: The convergence criterion is not satisfied.  \n')
		cat('Please refer to the epsr graph and increase the value of N.burn and MCMAX.')
	
		mycolsi <- rainbow(ncol(epsr), s = 1, v = 1, start = 0, 
				end = max(1, ncol(epsr) - 1)/ncol(epsr), alpha = 1) 
					#One color per parameter

		repxlim<-c(1:(MCMAX-1))
		plot(x = repxlim , y = epsr[,1], type="l", 
			xlab = "iterations", ylab = "EPSR", ylim = c(0,3), col = mycolsi[1])
		for (i in 2:ncol(epsr))
		{
			lines(x = repxlim, y = epsr[,i], col = mycolsi[i])
		}

		savePlot(filename = "EPSR",
				type ="png",
				device = dev.cur(),
				restoreConsole = TRUE)
	}
	
	
  
}