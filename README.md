
# blcfa
author: Junhao Pan, Lijin Zhang

date: May 21, 2019


[![Build Status](https://travis-ci.org/zhanglj37/blcfa.svg)](https://travis-ci.org/zhanglj37/blcfa)
[![Downloads from the RStudio CRAN mirror](http://cranlogs.r-pkg.org/badges/blcfa)](https://cran.r-project.org/package=blcfa)

## Description
The 'blcfa' package Uses Bayesian covariance lasso prior confirmatory factor analysis to detect significant residual correlations and generate the corresponding mplus file.

If you want to know more about Bayesian covariance lasso prior confirmatory factor analysis, please refer to Pan, Ip and Dubé(2017).

(Pan, J., Ip, E. H., & Dubé, L. (2017). An alternative to post hoc model modification in confirmatory factor analysis: the Bayesian lasso. *Psychological Methods, 22*(4), 687–704.)

## Installation
```r
install.packages("blcfa")
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

filename <- "ss.txt"  # use null value to represent missing value in the dataset
varnames <- c("gender",paste("y", 1:17, sep = ""))
usevar <- c(paste("y", 1:17, sep = ""))

myModel<-'   
# 1. CFA
f1 =~ y1 + y2 + y3 + y4 + y5 
f2 =~ y6 + y7 + y8 + y9 + y10 + y11
f3 =~ y12 + y13 + y14 + y15 + y16 + y17  
'
# make sure there is a space between each variable or symbol

blcfa(filename,varnames,usevar,myModel)
```

After running this function(two chain):
```r
num of chain:1
num of iteration:100
num of iteration:200
....
num of chain:2
...
```

you will get Mplus input file and output file  that include significant residual correlations detected by Bayesian covariance lasso Prior confirmatory factor analysis. for example:
```
TITLE: Bayesian Lasso CFA
 DATA: FILE =  ss.txt ; 
 VARIABLE:
 NAMES = gender y1 y2 y3 y4 y5 y6 y7 y8 y9 
	y10 y11 y12 y13 y14 y15 y16 y17 ;
 USEV = gender y1 y2 y3 y4 y5 y6 y7 y8 y9 
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
The convergence criterion is epsr value < 1.2. If the model does not converge within the number of burn-in MCMC samples(N.burn) (the default num of burn-in MCMC samples=5000), you will get an epsr graph and the warnings:
```r
Error: The convergence criterion is not satisfied.
Please refer to the epsr graph and increase the MCMAX.
```


Then you should increase the value of N.burn and MCMAX(Total number of MCMC samples for inference).
```r
blcfa(filename,varnames,usevar,myModel,MCMAX=30000,N.burn=15000)
```

### ex3
If you want to get the detailed results of the Bayesian covariance lasso prior confirmatory factor analysis:
```r
blcfa(filename,varnames,usevar,myModel,MCMAX = 10000, N.burn = 5000,bloutput = TRUE)
```

Then you will get the results folder includes: ppp, epsr graph,
			estimated value, standard error and hpd interval of ly, mu, phi and psx.

### ex4
Detect significant residual correlations by p-value rather than Highest Posterior Density (HPD)  interval.
```r
blcfa(filename,varnames,usevar,myModel,MCMAX = 10000, N.burn = 5000,bloutput = TRUE,interval_psx = FALSE)
```

## BugsReports

https://github.com/zhanglj37/blcfa/issues

or contact with me: zhanglj37@mail2.sysu.edu.cn. I really hope that we can jointly promote the improvement of this package.
