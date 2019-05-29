\name{blcfa}
\alias{blcfa}
\title{Bayesian Covariance Lasso Prior Confirmatory Factor Analysis}
\description{Uses Bayesian covariance lasso Prior confirmatory factor analysis to 
             detect significant residual correlations and generate the corresponding mplus file.
			 After running this function, you can get an mplus file that includes your model and 
			 the significant residual correlations detected by Bayesian Covariance Lasso Prior Confirmatory Factor Analysis.}
\usage{
blcfa(filename, varnames, usevar, model, MCMAX = 15000, N.burn = 5000,   
			ms = -999, bloutput = FALSE,  interval_psx = FALSE)

}			 
\arguments{
\item{filename}{Name of the data file (eg. "Y.txt"). 
			Make sure the data file is in the format of dat or txt, 
			and the variable name is not included in the file.}
\item{varnames}{Colnames of your dataset (eg. c("gender",paste("y", 1:19, sep = ""))).}
\item{usevar}{select a subset of variables for analysis (eg. c(paste("y", 1:19, sep = ""))))}
\item{model}{define your model by a string (eg. 
			myModel<-'
			# 1 ss model
			f1 =~ y1 + y2 + y3 + y4 + y5 + y6 + y7
			f2 =~ y8 + y9 + y10 + y11 + y12 + y13 + y14
			f3 =~ y15 + y16 + y17 + y18 + y19').
			Wrap when you define a new factor, 
			define he relationship between factors and indicators by the =~ symbol,
			make sure there is a space between each variable or symbol
			(eg. y1 + y2 (right), y1+y2(wrong)).}
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
\item{interval_psx}{Detect the significant residual correlations by hpd interval or p-value,
			the default setting is using hpd interval.}
}			 
\references{
Yves Rosseel (2012). lavaan: An R Package for Structural Equation
Modeling. Journal of Statistical Software, 48(2), 1-36. URL
http://www.jstatsoft.org/v48/i02/.}
\seealso{
\code{\link{lavaan}}
}