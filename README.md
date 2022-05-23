
# blcfa

update: April 27, 2020

[![](https://img.shields.io/github/stars/zhanglj37/blcfa?style=social)](https://github.com/zhanglj37/blcfa/stargazers)


* [Description](#Description)
* [Installation](#Installation)
* [Examples](#Examples)
  * [Example code in the paper.](#example-code-in-the-paper)
  * [EX1: Detect Cross-loadings and Residual Covariances Simultaneously.](#ex1-detect-cross-loadings-and-residual-covariances-simultaneously)
  * [EX2: Detect Residual Covariances.](#ex2-detect-residual-covariances)
  * [EX3: Detect Cross-Loadings.](#ex3-detect-cross-loadings)
  * [EX4: Detect Significant Loadings to Explore the Model Structure.](#ex4-detect-significant-loadings-to-explore-the-model-structure)
  * [Tips1: Model Convergence](#Tips1)
  * [Tips2: Detailed Results](#Tips2)
  * [Tips3: Threshold and HPD interval](#Tips3)
  * [Tips4: Normality Test](#Tips4)
  * [Tips5: Missing Data](#Tips5)
* [BugsReports](#BugsReports)
* [Functions under development](#functions-under-development)
* [Acknowledgement](#Acknowledgement)

## Description
This 'blcfa' package aims to: (1) detect significant cross-loadings and/or residual covariances different from zero by Bayesian covariance Lasso CFA; (2.1) free the identified significant parameters; (2.2) automatically feed the output from (2.1) into M*plus* to obtain an appropriately modified CFA model using Maximum likelihood (ML) or Bayesian estimation. 

For the details about Bayesian lasso confirmatory factor analysis, please refer to 1. Pan, Ip and Dubé(2017), 2. Chen, Guo, Zhang and Pan (2020).

1. Pan, J., Ip, E. H.\*, & Dubé, L. (2017). An alternative to post hoc model modification in confirmatory factor analysis: the Bayesian lasso. *Psychological Methods, 22*(4), 687–704. 
2. Chen, J.\*, Guo, Z., Zhang, L., Pan, J.\* (2020). A Partially Confirmatory Approach to Scale Development with the Bayesian Lasso. *Psychological Methods.* Advance Online Publication.

Citation of this package:

Zhang, L., Pan, J., Dubé, L., Ip, E.H. (2021). blcfa: An R Package for Bayesian Model Modification in Confirmatory Factor Analysis. *Structural Equation Modeling: A Multidisciplinary Journal*. Advance Online Publication 

## Installation

If you want to try out the latest development 'blcfa' code, you can install it  from github using Hadley Wickham's 'devtools' package. 

```r
install.packages("devtools")
library(devtools)

install_github("zhanglj37/blcfa")
```

## Examples

### Example code in the paper

```R
library(blcfa)
setwd("C:/Users/Desktop/SimuExample/")
filename = system.file("extdata", "simu_data.txt", package = "blcfa")
varnames<-c(paste("y", 1:10, sep = ""))
usevar <- varnames
myModel<-"
f1 =~ y1 + y2 + y3 + y4 + y5
f2 =~ y6 + y7 + y8 + y9 + y10
"
set.seed(1)
results <- blcfa(filename, varnames, usevar, myModel, estimation = "both", MCMAX = 5000, N.burn = 2500, bloutput = TRUE, interval = TRUE)

```

```

## Mplus input 
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
```

```R

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


```



### EX1: Detect Cross-loadings and Residual Covariances Simultaneously.

```r
library(blcfa)

filename = "ss.txt"  
varnames = c("gender",paste("y", 1:17, sep = ""))  # variables in dataset
usevar = c(paste("y", 1:17, sep = ""))  # variables used in the analysis
NZ = 3  # number of factors
IDY = matrix(c(
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
# ms represents missing value (you don't need to define it if NA represents missing value in the dataset).
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

### EX2: Detect Residual Covariances.

```r
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

blcfa(filename, varnames, usevar, IDY, estimation = 'Bayes', ms = -9)

```

### EX3: Detect Cross-Loadings. 

```r
# The model structure is same as EX1

blcfa_ly(filename, varnames, usevar, IDY, estimation = 'Bayes', ms = -9)
# residual covariances were set at zero in "blcfa_ly" function

```


### EX4: Detect Significant Loadings to Explore the Model Structure. 

```r
NZ=3
IDY<-matrix(c(
  9,-1,-1,
  1,-1,-1,
  -1,-1,-1,
  -1,-1,-1,
  -1,-1,-1,
  -1,9,-1,
  -1,1,-1,
  -1,-1,-1,
  -1,-1,-1,
  -1,-1,-1,
  -1,-1,-1,
  -1,-1,9,
  -1,-1,1,
  -1,-1,-1,
  -1,-1,-1,
  -1,-1,-1,
  -1,-1,-1
),ncol=NZ,byr=T)

# This is an extra application of Bayesian Lasso CFA, Check Chen et al (accepted) for the details of this method

# To illustrate this model structure, the corresponding relationships between factors and loadings were listed as follows:
# f1: y1@1 y2 y3-y17(eatimate with lasso shrinkage)
# f2: y1-y5(eatimate with lasso shrinkage) y6@1 y7 y8-y17(eatimate with lasso shrinkage)
# f3: y1-y11(eatimate with lasso shrinkage) y12@1 y13 y14-y17(eatimate with lasso shrinkage)

# make sure there are at least two identified loadings per factor (Chen et al., accepted)

blcfa_ly(filename, varnames, usevar, IDY, estimation = 'Bayes', ms = -9)

```


### Tips1

The convergence criterion is epsr value < 1.2. If the model does not converge within the number of burn-in MCMC samples(N.burn) (the default value = 5000), you will get an epsr graph (for reference) and the warnings:
```r
Error: The convergence criterion is not satisfied.
Please refer to the epsr graph and increase the MCMAX.
```


Then you should increase the value of N.burn and MCMAX (Total number of MCMC samples for inference, the default value = 15000).
```r
blcfa(filename,varnames,usevar,myModel,ms=-9,MCMAX=30000,N.burn=15000)
```

### Tips2
If you want to get the detailed results of the Bayesian covariance lasso prior confirmatory factor analysis:
```r
blcfa(filename,varnames,usevar,myModel,ms=-9,MCMAX = 10000, N.burn = 5000, bloutput = TRUE)
```

Then you will get the results folder includes: ppp, epsr graph,
			estimated value, standard error and hpd interval of parameters (ly, mu, phi and psx).

### Tips3
Detect significant cross-loadings and residual correlations by threshold (0.1 for cross-loadings, 0.15 for residual correlations) rather than 95% Highest Posterior Density (HPD)  interval.
```r
blcfa(filename,varnames,usevar,myModel,ms=-9,MCMAX = 10000, N.burn = 5000, bloutput = TRUE, interval = FALSE)
```

### Tips4
When the estimation is set at ml, the 'blcfa' package will select a maximum likelihood estimator based on the normality of items.
MLM method will be specified if the dataset did not satisfy multivariate normal distribution (Skewness > 2 or Kurtosis > 7; West, Finch, & Curran, 1995). 

West, S. G., Finch, J. F., & Curran, P. J. (1995). Structural equation models with nonnormal variables. *Structural equation modeling: Concepts, issues, and applications*, 56-75. 

By the way, the mplus_ml() function in this package can detect the non-normality of data and generate the Mplus file (traditional CFA model) without doing Bayesian Lasso CFA analysis.
```r
myModel<-'   
# 1. CFA
f1 =~ y1 + y2 + y3 + y4 + y5 
f2 =~ y6 + y7 + y8 + y9 + y10 + y11
f3 =~ y12 + y13 + y14 + y15 + y16 + y17  
'

mplus_ml(filename, varnames, usevar, myModel, ms = -9)
```

### Tips5
For missing data that is assumed missing at random (MAR), this package can automatically impute these missing values with Gibbs sampling method and generated a new data set. 

## BugsReports

https://github.com/zhanglj37/blcfa/issues

or contact with me: zhanglj37@mail2.sysu.edu.cn.

## Functions under development

Bayesian lasso confirmatory factor analysis models with ordered categorical data (expected in 2022).

If you have any suggestions or are willing to join in the improvement of this package, please contact with me. 

## Acknowledgement

Thanks to YaTing Deng for adding parallel computing to this package.
