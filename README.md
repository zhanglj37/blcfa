
# blcfa
author: Junhao Pan, Lijin Zhang

published: August 7, 2019

[![Build Status](https://travis-ci.org/zhanglj37/blcfa.svg)](https://travis-ci.org/zhanglj37/blcfa)
[![Downloads from the RStudio CRAN mirror](http://cranlogs.r-pkg.org/badges/blcfa)](https://cran.r-project.org/package=blcfa)
[![](https://cranlogs.r-pkg.org/badges/grand-total/blcfa)](https://cran.r-project.org/package=blcfa)

## Description
The 'blcfa' package was built to conduct Bayesian model modification in confirmatory factor analysis. It can shrink weak residual correlations toward zero and detect significant residual correlations by Bayesian lasso method.

This package aims to: (1) detect significant cross-loadings and/or residual covariances different from zero by Bayesian covariance Lasso CFA; (2.1) free the identified significant parameters; (2.2) automatically feed the output from (2.1) into M*plus* to obtain an appropriately modified CFA model using Maximum likelihood (ML) or Bayesian estimation. 

If you would like to know the details about Bayesian covariance lasso prior confirmatory factor analysis, please refer to 1. Pan, Ip and Dubé(2017), 2. Chen, Guo, Zhang and Pan (accepted).

1. Pan, J., Ip, E. H.\*, & Dubé, L. (2017). An alternative to post hoc model modification in confirmatory factor analysis: the Bayesian lasso. *Psychological Methods, 22*(4), 687–704.
2. Chen, J.S.\*, Guo, Z.H., Zhang, L.J., Pan, J.H.\* (accepted). A Partially Confirmatory Approach to Scale Development with the Bayesian Lasso. *Psychological Methods.* 

We are also preparing a paper to introduce this package:

Pan, J.H., Zhang, L.J. (co-first author), Ip, E.H.\* (manuscript drafted). BLCFA: An R Package for Bayesian Model Modification in Confirmatory Factor Analysis. 


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

### ex1: detect cross-loadings and residual covariances simultaneously

```r
library(blcfa)

filename <- "ss.txt"  
varnames <- c("gender",paste("y", 1:17, sep = ""))
usevar <- c(paste("y", 1:17, sep = ""))
NZ=3
IDY<-matrix(c(
  9,-1,-1,
  1,-1,-1,
  1,-1,-1,
  1,-1,-1,
  1,-1,-1,
  -1,9,-1,
  -1,1,-1,
  -1,1,-1,
  -1,1,-1,
  -1,1,-1,
  -1,1,-1,
  -1,-1,9,
  -1,-1,1,
  -1,-1,1,
  -1,-1,1,
  -1,-1,1,
  -1,-1,1
),ncol=NZ,byr=T)
# NZ: number of factors
# 9: fixed at one for identifing the factor
# 1: estimate this parameter without shrinkage
# -1: estimate this parameter using lasso shrinkage
# 0: fixed at zero.

# To illustrate this model structure, the corresponding relationships between factors and loadings were listed as follows:
# f1: y1@1 y2-y5 y6-y17(eatimate with lasso shrinkage)
# f2: y1-y5(eatimate with lasso shrinkage) y6@1 y7-y11 y12-y17(eatimate with lasso shrinkage)
# f3: y1-y11(eatimate with lasso shrinkage) y12@1 y13-y17

blcfa(filename, varnames, usevar, IDY, estimation = 'Bayes', ms = -9)
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
	f1 by 	y1  y2  y3  y4  y5  y17  ;
	f2 by 	y6  y7  y8  y9  y10  y11  y13  y14  ;
	f3 by 	y12  y5  y13  y14  y15  y16  y17  ;
	 
	y11  with  y13 ;
	y11  with  y14 ;
	y13  with  y14 ;
	
 OUTPUT: TECH1  TECH8  STDY;
 PLOT: TYPE= PLOT2;

```

### ex2: detect cross-loadings 

```r
# same as ex1

blcfa_ly(filename, varnames, usevar, IDY, estimation = 'Bayes', ms = -9)

```


### ex3: detect residual covariances

```r
# same as ex1
NZ=3
IDY0<-matrix(c(
  9,0,0,
  1,0,0,
  1,0,0,
  1,0,0,
  1,0,0,
  0,9,0,
  0,1,0,
  0,1,0,
  0,1,0,
  0,1,0,
  0,1,0,
  0,0,9,
  0,0,1,
  0,0,1,
  0,0,1,
  0,0,1,
  0,0,1
),ncol=NZ,byr=T)
# 0: fixed at zero.

# To illustrate this model structure, the corresponding relationships between factors and loadings were listed as follows:
# f1: y1@1 y2-y5
# f2: y6@1 y7-y11 
# f3: y12@1 y13-y17

blcfa_psx(filename, varnames, usevar, IDY, estimation = 'Bayes', ms = -9)

```

### Tips 1

The convergence criterion is epsr value < 1.2. If the model does not converge within the number of burn-in MCMC samples(N.burn) (the default value = 5000), you will get an epsr graph (for reference) and the warnings:
```r
Error: The convergence criterion is not satisfied.
Please refer to the epsr graph and increase the MCMAX.
```


Then you should increase the value of N.burn and MCMAX (Total number of MCMC samples for inference, the default value = 15000).
```r
blcfa(filename,varnames,usevar,myModel,ms=-9,MCMAX=30000,N.burn=15000)
```

### Tips 2
If you want to get the detailed results of the Bayesian covariance lasso prior confirmatory factor analysis:
```r
blcfa(filename,varnames,usevar,myModel,ms=-9,MCMAX = 10000, N.burn = 5000,bloutput = TRUE)
```

Then you will get the results folder includes: ppp, epsr graph,
			estimated value, standard error and hpd interval of parameters (ly, mu, phi and psx).

### Tips 3
Detect significant residual correlations by p-value rather than Highest Posterior Density (HPD)  interval.
```r
blcfa(filename,varnames,usevar,myModel,ms=-9,MCMAX = 10000, N.burn = 5000,bloutput = TRUE,interval_psx = FALSE)
```

## BugsReports

https://github.com/zhanglj37/blcfa/issues

or contact with me: zhanglj37@mail2.sysu.edu.cn.

## Functions under development

Bayesian lasso confirmatory factor analysis models with ordered categorical data (expected in 2021).

If you have any suggestions or are willing to join in the improvement of this package, please contact with me.  I really hope that we can jointly promote the improvement of this package.

## Acknowledgement

Thanks to YaTing Deng for adding parallel computing to this package.
