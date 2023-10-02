\name{blcfa}
\alias{blcfa}
\title{Bayesian Covariance Lasso Prior Confirmatory Factor Analysis}
\description{Uses Bayesian covariance lasso Prior confirmatory factor analysis to
             detect significant corss-loadings and residual correlations and generate the corresponding mplus file.
			 After running this function, you can get an mplus file that includes your model and
			 the significant cross-loadings and residual correlations detected by Bayesian Lasso Prior Confirmatory Factor Analysis.}
\usage{
blcfa(filename, varnames, usevar, myModel, estimation = "ml", ms = -999999, 
      MCMAX = 10000, N.burn = 5000, bloutput = FALSE,  
	  interval = TRUE, conver_check = TRUE)

}
\arguments{
\item{filename}{Name of the data file (eg. "Y.txt").
			Make sure the data file is in the format of dat or txt,
			and the variable name is not included in the file.}
\item{varnames}{Colnames of your dataset (eg. c("gender",paste("y", 1:19, sep = ""))).}
\item{usevar}{select a subset of variables for analysis (eg. c(paste("y", 1:19, sep = ""))))}
\item{myModel}{define your model by a string }
\item{estimation}{Estimator in the Mplus: 'ML' or 'Bayes'.}
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
\item{conver_check}{TRUE: use two MCMC chains to caculate the EPSR valus and check whether the model converge.
	FALSE: use one MCMC chain to get the estimates without convergence check}
}

\examples{
# SimuExample
## run model
library(blcfa)
filename = system.file("extdata", "simu_data.txt", package = "blcfa")
varnames<-c(paste("y", 1:10, sep = ""))
usevar <- varnames
myModel<-"
f1 =~ y1 + y2 + y3 + y4 + y5
f2 =~ y6 + y7 + y8 + y9 + y10
"
set.seed(1)
results <- blcfa(filename, varnames, usevar, myModel, estimation = "both", MCMAX = 5000, N.burn = 2500, bloutput = TRUE, interval = TRUE)

## Mplus input for SimuExample
TITLE: Bayesian Lasso CFA
DATA: FILE =  D:/Software/R/R-4.0.1/library/blcfa/extdata/simu_data.txt ;
VARIABLE:
 NAMES = y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 ;
 USEV = y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 ;
ANALYSIS:
	ESTIMATOR = ML;
MODEL:
	f1 by 	y1  y2  y3  y4  y5  ;
	f2 by 	y6  y7  y8  y9  y10  ;
	y6  with  y1 ;
	y7  with  y2 ;
	y10  with  y9 ;
	
OUTPUT: TECH1  STDYX;

## results of the first-step analysis
> results
$blcfa_est
$blcfa_est$ppp
[1] 0.5156

$blcfa_est$ly
              est    se p-value HPD_lower HPD_upper
f1  by  y2  0.936 0.043       0     0.850     1.013
...

## results of the second-step analysis
> sum_second(results)
Reading model:  blcfa_bayes.out
Reading model:  blcfa_ml.out
$bayes_fit
...
  Parameters   CFI   TLI      BIC      DIC     pD RMSEA_Estimate RMSEA_90CI_LB
1         34   1   0.999  10991.69 10846.96 33.226          0.001            0
  RMSEA_90CI_UB RMSEA_pLT05 ObsRepChiSqDiff_95CI_LB ObsRepChiSqDiff_95CI_UB
1         0.032           1                 -28.816                    29.6
  PostPred_PValue        Filename
1           0.465 blcfa_bayes.out

$bayes_par_est
          paramHeader param   est posterior_sd pval lower_2.5ci upper_2.5ci   sig
1               F1.BY    Y1 1.000        0.000    0       1.000       1.000 FALSE
2               F1.BY    Y2 0.777        0.026    0       0.728       0.831  TRUE
...

$bayes_par_est_std
          paramHeader param   est posterior_sd pval lower_2.5ci upper_2.5ci   sig
1               F1.BY    Y1 0.863        0.013    0       0.835       0.888  TRUE
2               F1.BY    Y2 0.811        0.017    0       0.776       0.842  TRUE
...

$ml_fit
...
  ChiSqBaseline_PValue        LL UnrestrictedLL   CFI   TLI      AIC      BIC
1                    0  -5389.93       -5373.401 0.999 0.999 10847.86 10991.16
...

$ml_par_est
          paramHeader param   est    se  est_se pval
1               F1.BY    Y1 1.000 0.000 999.000  999
...

$ml_par_est_std
          paramHeader param   est    se  est_se pval
1               F1.BY    Y1 0.864 0.013  64.475    0
...

}


\references{
Pan, J., Ip, E. H., & Dubé, L. (2017). An alternative to post hoc model modification in confirmatory factor analysis: the Bayesian lasso. Psychological Methods, 22(4), 687-704.

Chen, J.S.*, Guo, Z.H., Zhang, L.J., Pan, J.H.* (2020). A Partially Confirmatory Approach to Scale Development with the Bayesian Lasso. Psychological Methods. Advance online publication.

Zhang, L., Pan, J.*, Dubé, L., Ip, E.H. (accepted). blcfa: An R Package for Bayesian Model Modification in Confirmatory Factor Analysis. Structural Equation Modeling: A Multidisciplinary Journal.
}
