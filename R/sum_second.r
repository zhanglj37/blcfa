sum_second <- function(blcfa_results){

estimation = blcfa_results$estimation
estimation = tolower(estimation)
if (estimation == 'bayes' || estimation == 'bayesian')
{
	mplus_file = 'blcfa_bayes.out'
}else if (estimation == 'both'){
	mplus_file = c('blcfa_bayes.out', 'blcfa_ml.out')
}else{
	mplus_file = 'blcfa_ml.out'
#	fit_index = c('ChiSqM_Value', 'ChiSqM_DF', 'ChiSqM_PValue',
#		'CFI', 'TLI', 'AIC', 'BIC', 'aBIC', 'RMSEA_Estimate', 'RMSEA_90CI_LB', 
#		'RMSEA_90CI_UB', 'SRMR')
}

## fit index




## run
path_detect = Sys.getenv("PATH")
path_detect = tolower(path_detect)
if (str_detect(path_detect,'mplus')){

	if (length(mplus_file) == 1){
	mplus_results= readModels(mplus_file)
	fit_sum = mplus_results$summaries
	par_est = mplus_results$parameters[[1]]
	par_est_std = mplus_results$parameters[[3]]
	
	results_second = list(#fit_summary = fit_sum,
		estimates = par_est,
		std_estimates = par_est_std)
		
	}else{
	bayes_results = readModels(mplus_file[1])
	ml_results = readModels(mplus_file[2])
	
	bayes_fit = bayes_results$summaries
	bayes_par_est = bayes_results$parameters[[1]]
	bayes_par_est_std = bayes_results$parameters[[3]]
	
	ml_fit = ml_results$summaries
	ml_par_est = ml_results$parameters[[1]]
	ml_par_est_std = ml_results$parameters[[3]]
	
	results_second = list(bayes_fit = bayes_fit,
		bayes_par_est = bayes_par_est,
		bayes_par_est_std = bayes_par_est_std,
		ml_fit = ml_fit,
		ml_par_est = ml_par_est,
		ml_par_est_std = ml_par_est_std)	
	
	
	}

	return(results_second)
}else{
	print('Error: Failed to run the Mplus software and get the results of the second-step analysis')
}


}


