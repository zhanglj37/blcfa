

library('blcfa')
library(doParallel)
						
setwd("E:19ZLJ/3_lasso_package/Lasso_users' guide/blcfa_renew")
source("read_model.r")
source("EPSR_set_int.R")
source("read_data.R")
source("Gibbs.R")
source("Gibbs_psx.R")
source("Gibbs_ly.R")
source("caculate_results.r")
source("HPD.R")
source("sigpsx.r")
source("sigly.r")
source("EPSR_caculate.r")
source("write_mplus_bayes.r")
source('idy.r')
source('impute_ms.r')
source("write_results.r")
source("EPSR_figure.r")

filename <- "ss.txt" 
varnames <- c("gender",paste("y", 1:17, sep = ""))
usevar <- c(paste("y", 1:17, sep = ""))

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

IDY0<-matrix(c(
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


ms=999
estimation = 'Bayes'
MCMAX = 15000
N.burn = 10000
bloutput = FALSE
interval_psx = TRUE
CIR=1
nthin=1
CNUM=2
NY=14
NZ=3
N=218
IDY0<-matrix(c(
  9,-1,-1,
  -1,9,-1,
  1,-1,-1,
  -1,1,-1,
  1,-1,-1,
  1,-1,-1,
  -1,1,-1,
  -1,1,-1,
  -1,1,-1,
  -1,-1,9,
  -1,-1,1,
  -1,-1,1,
  -1,-1,1,
  -1,1,1
),ncol=NZ,byr=T)
fname<-"Math_learning_J14.csv"
Y<-t(read.csv(fname,header = FALSE))
