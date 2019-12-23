
# blcfa
author: Junhao Pan, Lijin Zhang

published: August 7, 2019

[![Build Status](https://travis-ci.org/zhanglj37/blcfa.svg)](https://travis-ci.org/zhanglj37/blcfa)
[![Downloads from the RStudio CRAN mirror](http://cranlogs.r-pkg.org/badges/blcfa)](https://cran.r-project.org/package=blcfa)
[![](https://cranlogs.r-pkg.org/badges/grand-total/blcfa)](https://cran.r-project.org/package=blcfa)

## Description
The 'blcfa' package uses Bayesian covariance lasso prior confirmatory factor analysis to detect significant residual covariances and generate the corresponding mplus file.

If you would like to know the details about Bayesian covariance lasso prior confirmatory factor analysis, please refer to Pan, Ip and Dubé(2017).

(Pan, J., Ip, E. H., & Dubé, L. (2017). An alternative to post hoc model modification in confirmatory factor analysis: the Bayesian lasso. *Psychological Methods, 22*(4), 687–704.)

## Installation
```r
install.packages("blcfa")  #It hasn't been released yet. 2019-08-07
```

If you want to try out the latest development 'blcfa' code, you can install it  from github using Hadley Wickham's 'devtools' package. 

```r
install.packages("devtools")
library(devtools)

install_github("zhanglj37/blcfa")
```


## Examples

### ex1

```r
library(blcfa)

filename <- "ss.txt"  
varnames <- c("gender",paste("y", 1:17, sep = ""))
usevar <- c(paste("y", 1:17, sep = ""))
myModel<-'   
# 1. CFA
f1 =~ y1 + y2 + y3 + y4 + y5 
f2 =~ y6 + y7 + y8 + y9 + y10 + y11
f3 =~ y12 + y13 + y14 + y15 + y16 + y17  
'
# make sure there is a space between each variable or symbol

blcfa(filename, varnames, usevar, myModel, estimation = 'Bayes', ms = -9)
# estimation ( = 'ML' / 'Bayes', the default value is 'Bayes') denotes the estimation method in Mplus file
# ms represents missing value (you don't need to define it if -999 or NA represents missing value in the dataset).
```

After running this function:
```r
The program is running. See 'log.txt' for details.

Gibbs sampling ended up, specific results are being calculated.
```
('log.txt' records the process of parallel computing of two MCMC chains)

You will get Mplus input file and output file  that include significant residual correlations detected by Bayesian covariance lasso prior confirmatory factor analysis. For example:
```
TITLE: Bayesian Lasso CFA
DATA: FILE =  ss.txt ; 
VARIABLE:
NAMES = gender y1 y2 y3 y4 y5 y6 y7 y8 y9 
	y10 y11 y12 y13 y14 y15 y16 y17 ;
USEV = y1 y2 y3 y4 y5 y6 y7 y8 y9 
	y10 y11 y12 y13 y14 y15 y16 y17 ;
ANALYSIS:
	 ESTIMATOR = BAYES;
	 PROC = 2;
	 BITERATIONS = (10000);
MODEL:
	! 1. CFA;
	f1  BY  y1  y2  y3  y4  y5 ;
	f2  BY  y6  y7  y8  y9  y10  y11;
	f3  BY  y12  y13  y14  y15  y16  y17  ;
	 
	y11  with  y13 ;
	y11  with  y14 ;
	y13  with  y14 ;
	
 OUTPUT: TECH1  TECH8  STDY;
 PLOT: TYPE= PLOT2;

```

### ex2
The convergence criterion is epsr value < 1.2. If the model does not converge within the number of burn-in MCMC samples(N.burn) (the default value = 5000), you will get an epsr graph (for reference) and the warnings:
```r
Error: The convergence criterion is not satisfied.
Please refer to the epsr graph and increase the MCMAX.
```


Then you should increase the value of N.burn and MCMAX (Total number of MCMC samples for inference, the default value = 15000).
```r
blcfa(filename,varnames,usevar,myModel,ms=-9,MCMAX=30000,N.burn=15000)
```

### ex3
If you want to get the detailed results of the Bayesian covariance lasso prior confirmatory factor analysis:
```r
blcfa(filename,varnames,usevar,myModel,ms=-9,MCMAX = 10000, N.burn = 5000,bloutput = TRUE)
```

Then you will get the results folder includes: ppp, epsr graph,
			estimated value, standard error and hpd interval of parameters (ly, mu, phi and psx).

### ex4
Detect significant residual correlations by p-value rather than Highest Posterior Density (HPD)  interval.
```r
blcfa(filename,varnames,usevar,myModel,ms=-9,MCMAX = 10000, N.burn = 5000,bloutput = TRUE,interval_psx = FALSE)
```

## BugsReports

https://github.com/zhanglj37/blcfa/issues

or contact with me: zhanglj37@mail2.sysu.edu.cn.

## Functions under development

Bayesian lasso partial confirmatory factor analysis models: detect cross-loadings and residual correlations simultaneously (expected in 2020/05).

Bayesian lasso confirmatory factor analysis models with ordered categorical data (expected in 2021).

If you have any suggestions or are willing to join in the improvement of this package, please contact with me.  I really hope that we can jointly promote the improvement of this package.

## Acknowledgement

Thanks to YaTing Deng for adding parallel computing to this package.
