\name{blcfa}
\alias{blcfa}
\title{Bayesian Covariance Lasso Prior Confirmatory Factor Analysis}
\description{Uses Bayesian covariance lasso Prior confirmatory factor analysis to 
             detect significant residual correlations and generate the corresponding mplus file.
			 After running this function, you can get an mplus file that includes your model and 
			 the significant residual correlations detected by Bayesian Covariance Lasso Prior Confirmatory Factor Analysis.}
\usage{
blcfa(filename, varnames, usevar, IDY0, MCMAX = 15000, N.burn = 5000,   
			ms = -999, bloutput = FALSE,  interval = TRUE)

}			 
\arguments{
\item{filename}{Name of the data file (eg. "Y.txt"). 
			Make sure the data file is in the format of dat or txt, 
			and the variable name is not included in the file.}
\item{varnames}{Colnames of your dataset (eg. c("gender",paste("y", 1:19, sep = ""))).}
\item{usevar}{select a subset of variables for analysis (eg. c(paste("y", 1:19, sep = ""))))}
\item{IDY0}{define your model by a matrix (eg. 
			IDY0<-matrix(c(
			  9,0,0,
			  1,0,0,
			  1,0,0,
			  1,0,0,
			  1,0,0,
			  -1,9,0,
			  -1,1,0,
			  -1,1,0,
			  0,1,0,
			  0,1,0,
			  0,1,0,
			  0,0,9,
			  0,0,1,
			  0,0,1,
			  0,0,1,
			  0,0,1,
			  0,0,1
			),ncol=NZ,byr=T).
			9:fixed at one for identifing the model, 1:estimate this parameter without shrinkage, -1:estimate this parameter using lasso shrinkage, 0:fixed at zero.}
\item{MCMAX}{Total number of MCMC samples for inference (the default value is 15000).}
\item{N.burn}{Number of burn-in MCMC samples. Discarded (the default value is 5000).
			Besides, if the model does not converge in the N.burn iteration, 
			the result and the mplus file will not be presented, 
			and you need to increase the value of N.burn and MCMAX.}
\item{ms}{define missing value (the default value is NA, which means null value in your dataset).}
\item{bloutput}{Results of bayesian covariance lasso prior cfa, incluse: ppp, 
			estimated value, standard error and hpd interval of ly, mu, phi and psx.
			The default setting is not output these results, if you just want
			the mplus input file then you don't need to change it.}
\item{interval}{Detect the significant residual correlations by hpd interval or p-value,
			the default setting is using hpd interval.}
}
 
\references{
Pan, J., Ip, E. H., & Dubé, L. (2017). An alternative to post hoc model modification in confirmatory factor analysis: the Bayesian lasso. *Psychological Methods, 22*(4), 687–704.
Chen, J.S.*, Guo, Z.H., Zhang, L.J., Pan, J.H.* (accepted). A Partially Confirmatory Approach to Scale Development with the Bayesian Lasso. Psychological Methods.
}
