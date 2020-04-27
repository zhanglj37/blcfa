#setwd("E:/R code/")

if(!require(MCMCpack)) install.packages(pkgs="MCMCpack",repos="http://cran.r-project.org")
if(!require(msm)) install.packages(pkgs="msm",repos="http://cran.r-project.org")
if(!require(statmod)) install.packages(pkgs="statmod",repos="http://cran.r-project.org")
if(!require(psychometric)) install.packages(pkgs="psychometric",repos="http://cran.r-project.org")

library(stats)
library(MASS)
library(MCMCpack)
library(msm)
library(statmod)
#library(psychometric)

rm(list=ls())
#fname<-"CON_N250_J10.dat"
#source("def_con.R")
#Sample size (N)
N<-250             
#Number of items (p)
NY<-18	       
#Number of factors (q)
NZ<-3		       
#Number of replication; for real data analysis, CNUM should be 1
CNUM<-2  
#Number of burn-in MCMC samples. Discarded          
N.burn<-2000 
#Total number of MCMC samples for inference.      
MCMAX<-4000     
#Save one sample for inference every nthin samples.
nthin<-1  

# ecr=0 #0: no error covariance
#source("ind.R")
#########Revise the parts below according to dimension of measurement equation.

IDY<-matrix(c(
  1,-1,-1,
  1,-1,-1,
  -1,-1,-1,
  -1,-1,-1,
  -1,-1,-1,
  -1,-1,-1,
  -1,-1,-1,
  -1,1,-1,
  -1,1,-1,
  -1,-1,-1,
  -1,-1,-1,
  -1,-1,-1,
  -1,-1,-1,
  -1,-1,-1,
  -1,-1,1,
  -1,-1,1,
  -1,-1,-1,
  -1,-1,-1
),ncol=NZ,byr=T)


#vector indicating whether MU should be added to each measurement equation or not. 
#Added(1), remove(0).
IDMU<-rep(0,NY)  # now MU is estimated
IDMUA<-any(as.logical(IDMU))

#source("def_rec.R")
NM<-0	             #dimension of eta (q_1)
NK<-NM+NZ	       #dimension of latent variables (eta+xi);  number of factors

############################## Automatic done
#Y<-array(0,dim=c(NY,N))			#observed data Y

xi<-array(dim=c(NZ,N))			#independent latent variable xi
#eta<-array(dim=c(NM,N))			#dependent latent variable eta
TUE<-Omega<-array(dim=c(NK,N))	#latent variable omega

NMU<-sum(IDMU)			      #number of Mu in measurement equation.
NLY<-sum(IDY!=0)				#number of free lambda need to be estimated in Lambda.
Nrec<-(MCMAX-N.burn)/nthin		#number of samples after burn-in.

EMU<-array(0,dim=c(Nrec,NMU))		#Store retained trace of MU
ELY<-array(0,dim=c(Nrec,NLY))		#Store retained trace of Lambda
#ELY<-array(0,dim=c(Nrec,NY*NZ))

EPSX<-array(0,dim=c(Nrec,NY,NY))	#Store retained trace of PSX
EinvPSX<-array(0,dim=c(Nrec,NY,NY)) #Store retained trace of inv(PSX)
EPHI<-array(0,dim=c(Nrec,NZ,NZ))	#Store retained trace of PHI
EXI<-array(0,dim=c(Nrec,NZ,N))	#Store retained trace of xi
Elambda<-array(0,dim=c(Nrec,1))     #Store retained trace of shrinkage paraemter lambda
Elamsq<-array(0,dim=c(Nrec,1))          #Store retained trace of shrinkage paraemter lamsq 

indmx<-matrix(1:NY^2, nrow=NY, ncol=NY)
temp<-indmx[upper.tri(indmx)]
upperind<-temp[temp>0]

indmx_t<-t(indmx)
temp<-indmx_t[upper.tri(indmx_t)]
lowerind<-temp[temp>0]

ind_noi_all<-array(0,dim=c(NY-1,NY))

for(i in 1:NY){
  
  if(i==1) {ind_noi<-2:NY}
  else if(i==NY) {ind_noi<-1:(NY-1)}
  else ind_noi<-c(1:(i-1),(i+1):NY)
  
  ind_noi_all[,i]<-ind_noi
  
} # end of i


# Empostp<-numeric(CNUM)
# chainpsx<-array(0,dim=c(Nrec,NY*(NY+1)/2))
# chainphi<-array(0,dim=c(Nrec,NZ*(NZ-1)/2))

ptm<-proc.time()
for(rep in 1:CNUM){
#rep=1
set.seed(12345+N+NY+rep)
  
  tausq<-array(0,dim=c(NY,NZ))
  tausq[IDY==-1]=1

  #source("Prior.R")
  ######Hyperparameter of prior distribution.
  
  #Measurement equation.
  
  #Prior mean of Lambda, NY*NZ
  PLY<-matrix(0.0,ncol=NZ,nrow=NY)
  
  #Prior mean of MU, NY*1
  PMU<-rep(0.0,NY)		    
  
  #Inverse prior variance of unknown parameters in factor loading matrix
  sigly<-0.25
  
  #Inverse prior variance of unknown parameters in intercept
  sigmu<-0.25
  
  rou.scale<-1.0
  
  #rho_0, hyperparameters of Wishart distribution			
  rou.zero<-rou.scale+NZ+1	
  
  #Matrix R_0, hyperparameters of Wishart distribution
  #R.zero<-rou.scale*diag(1,NZ)
  R.zero<-matrix(.1,nrow=NZ,ncol=NZ)
  diag(R.zero)<-rou.scale
  
  #Hyperparameters of Gamma distribution for the shrinkage parameter
  a_lamsq<-1
  b_lamsq<-0.01
  
  NLY1<-sum(IDY==-1)				#number of lasso lambda.
  a_psx<-1
  b_psx<-.01    
  #stau<-0
  LY_eps<- 0 #threhods for LY sign change
  
  
  #source("read_observed.R") 
  #Read the observed data from a text file named data.txt with dimension N*NY
  #Y[,]<-matrix(scan("data.txt",skip=(CIR-1)*N,nlines=N,sep="\t"),nrow=NY,ncol=N) 
  fname<-paste("CON_N",N,"_J",NY,"_i",rep,".dat",sep="")
  # Y<-t(read.csv(fname,header = FALSE))
   Y<-t(read.table(fname,header = FALSE))
  # Y[Y==-1]=NA
  #Standardized the data
   Y.temp<-t(Y)
   Y<-t(scale(Y.temp))
  
  
  #Creat the matrix of missing indicators where 1 represents missing
  # missing_ind<-array(0, dim=c(NY, N))
  # for(i in 1:NY)
  #    for(j in 1:N)
  #       if(is.na(Y[i,j])) missing_ind[i,j]<-1
  #       #if(Y[i,j]==-1) missing_ind[i,j]<-1
  # 
  # missingness=(sum(missing_ind)>0)
  
  
  
  #source("init1.R")
  #####Parameter in measurement equation
  
  #initial value of Lambda, 
  LY<-matrix(0,ncol=NK,nrow=NY)
  LY[IDY==1]=.7
  LY[IDY==-1]=.1
  #initial value of MU
  MU<-rep(1.0,NY)
  
  #initial value of PHI
  PHI<-matrix(0.3,nrow=NZ,ncol=NZ)
  diag(PHI[,])<-1.0								
  CPH<-PHI
  
  #initial value of PSX
  xi<-t(mvrnorm(N,mu=rep(0,NZ),Sigma=CPH)) # NZ*N
  #xi<-matrix(runif(NZ*N),nrow=NZ)+2# NZ*N
  
  PSX<-matrix(0.0,nrow=NY,ncol=NY)
  diag(PSX)<-0.5
  #invD_tau<-0
  
  
  #initial value of PHI^(-1)
  #inv.PHI<-chol2inv(chol(PHI))
  inv.CPH<-chol2inv(chol(CPH))
  
  #initial value of PHI^(-1/2)
  #c.inv.PHI<-chol(inv.PHI)			
  
  #initial value of PSX^(-1)
  inv.PSX<-chol2inv(chol(PSX))
  
  #initial value of PSX^(-1/2)
  inv.sqrt.PSX<-chol(inv.PSX)
  
  
  if(IDMUA==F) MU<-rep(0,NY)
  
  # initial values for missing data in Y
  for(i in 1:NY)
    for(j in 1:N)
      if(is.na(Y[i,j])) Y[i,j]<-rnorm(1)
  
  
  #source("Gibbs.R")
  Epostp<-array(0, dim=c(Nrec,1))
  
  
  for(g in 1:MCMAX){
    #g=1
    gm<-g-N.burn
    
    #Generate the latent factors from its conditinal distribution
    #source("Gibbs_Omega.R") 
    ##################  update Omega ##############################################################
    
    ISG<-crossprod(inv.sqrt.PSX%*%LY)+inv.CPH 
    SIG<-chol2inv(chol(ISG)) 
    Ycen<-Y   
    if(IDMUA==T) Ycen<-Y-MU   
    Mean<-SIG%*%t(LY)%*%inv.PSX%*%Ycen  
    for(i in 1:N) Omega[,i]<-xi[,i]<-mvrnorm(1, Mean[,i], Sigma=SIG) 
    
    ##################  end of update Omega #######################################################
    
    
    #Generate the unknown parameter of intercept in CFA from its conditinal distribution
    #source("Gibbs_MU.R") 
    ###################    update MU  ################################################################# 
    if(IDMUA==T){
      
      calsm<-chol2inv(chol(N*inv.PSX+diag(rep(sigmu,NY)))) # inv[sigma0^(-1)+N*inv.PSX]
      Ycen<-Y-LY%*%Omega		       
      temp<-rowSums(Ycen)
      mumu<-calsm%*%(inv.PSX%*%temp+rep(sigmu,NY)*PMU)
      MU<-mvrnorm(1,mumu,Sig=calsm)
      
    } 
    ###################    end of update MU  ##########################################################	
    
    # if (ecr==0){ #no error covariance
    # for (j in 1:NY){
    #  Ptemp<-1-t(LY[j,])%*%CPH%*%LY[j,]
    #  PSX[j,j]<-pmax(Ptemp,10^(-6))
    # }
    # inv.PSX<-chol2inv(chol(PSX))
    # inv.sqrt.PSX<-chol(inv.PSX)
    # }else{
    #Generate the unknown parameter of covariance matrix of measurement errors in CFA from its conditinal distribution
    #source("Gibbs_PSX.R")
    ###################    update LY  #################################################################   
    
    temp<-Y-MU-LY%*%Omega  # NY*N
    S<-temp%*%t(temp)      # NY*NY 
    
    apost1<-a_lamsq+NLY1;
    stau<-sum(tausq)
    
    #sample lambda
    bpost1<-b_lamsq + stau/2  # C is the presicion matrix
    lamsq<- rgamma(1, shape=apost1, rate=bpost1)
    
    
    lambda<-sqrt(lamsq)
    
    #sample PSX and inv(PSX)
    for(i in 1:NY){
      #i=1
      #subs<-(IDY[i,]==-1)
      ind<-which(IDY[i,]==-1)
      #len<-length(LY[i,subs])
      len<-length(ind)
      if(len>0){
        #gam<-rgamma(1, shape=N/2+1, rate=(S[i,i]+lambda)/2)
        #invD_tau<-diag(1/tausq[i,ind])
        if(len==1){
          invD_tau<-1/tausq[i,ind]
        }else{
          invD_tau<-diag(1/tausq[i,ind])
        }
        
        tmp<-t(LY[i,ind])%*%invD_tau%*%LY[i,ind]
        gam<-rgamma(1, shape=a_psx+(N+len)/2-1, rate=b_psx+(S[i,i]+tmp)/2)
        inv.PSX[i,i]<-gam
        PSX[i,i]<-1/gam
      }
    } # end of i, sample Sig and C=inv(Sig)
    
    inv.sqrt.PSX<-chol(inv.PSX)
    
    #count.n<-1 
    for(j in 1:NY){
      #j=1        
      temp1<-chol2inv(chol(PSX[-j,-j]))
      convar<-PSX[j,j]-PSX[j,-j]%*%temp1%*%PSX[-j,j]
      #invconvar<-chol2inv(chol(convar))
      invconvar<-1/convar[1,1] 
      subs<-(IDY[j,]==1)
      len<-length(LY[j,subs])    
      if(len>0){
        
        Ycen<-Y[j,]-MU[j]  # 1*N	  
        #Ycen<-Ycen-matrix(LY[j,(!subs),drop=F],nrow=1)%*%matrix(Omega[(!subs),,drop=F],ncol=N) # 1*N
        
        Ycen<-Ycen-matrix(LY[j,(!subs)],nrow=1)%*%matrix(Omega[(!subs),],ncol=N)-PSX[j,-j]%*%temp1%*%(Y[-j,]-MU[-j]-LY[-j,]%*%Omega) # 1*N 
        Ycen<-as.vector(Ycen) # vector      			     
        if(len==1){omesub<-matrix(Omega[subs,],nrow=1)
        }else{omesub<-Omega[subs,]}
        PSiginv<-diag(len)
        diag(PSiginv)<-rep(sigly,len)
        Pmean<-PLY[j,subs]
        #calsmnpsx<-chol2inv(chol(invconvar%*%tcrossprod(omesub)+PSiginv))
        #temp<-(omesub%*%Ycen%*%invconvar+PSiginv*Pmean)
        calsmnpsx<-chol2inv(chol(invconvar*tcrossprod(omesub)+PSiginv))
        temp<-(omesub%*%Ycen*invconvar+PSiginv%*%Pmean)
        LYnpsx<-calsmnpsx%*%temp
        LY[j,subs]<-mvrnorm(1,LYnpsx,Sig=(calsmnpsx))
        #if((gm>0)&&(gm%%nthin==0)){ELY[gm/nthin,count.n:(count.n+len-1)]<-LY[j,subs]}
        #    count.n<-count.n+len
      } # end len>0
      
      subs<-(IDY[j,]==-1)
      ind<-which(IDY[j,]==-1)
      len<-length(ind)
      
      if(len>0){
        
        Ycen<-Y[j,]-MU[j]  # 1*N	  
        #Ycen<-Ycen-matrix(LY[j,(!subs),drop=F],nrow=1)%*%matrix(Omega[(!subs),,drop=F],ncol=N) # 1*N
        #temp1<-chol2inv(chol(PSX[-j,-j]))  
        Ycen<-Ycen-matrix(LY[j,(!subs)],nrow=1)%*%matrix(Omega[(!subs),],ncol=N)-PSX[j,-j]%*%temp1%*%(Y[-j,]-MU[-j]-LY[-j,]%*%Omega) # 1*N 
        #Ycen<-Ycen-matrix(LY[j,(!subs)],nrow=1)%*%matrix(Omega[(!subs),],ncol=N)
        
        Ycen<-as.vector(Ycen) # vector
        #convar<-PSX[j,j]
        #invconvar<-1/convar
        Cadj1<-pmax((LY[j,subs])^2,10^(-6))
        #mu_p1<-pmin(sqrt(lamsq/Cadj1), 10^12)
        mu_p1<-pmin(sqrt(lamsq*PSX[j,j]/Cadj1), 10^12)
        #lambda_prime<-lambda^2
        #tausq<-rep(0,length(mu_p1))
        #for(i in 1:length(mu_p1)){
        for(i in 1:len){
          tausq[j,ind[i]]<-1/rinvgauss(1, mean=mu_p1[i], dispersion=1/lamsq)
        }        
        #invD<-diag(1/tausq[j,ind])
        
        if(len==1){
          omesub<-matrix(Omega[subs,],nrow=1)
          invD_tau<-1/tausq[j,ind]
        }else{
          omesub<-Omega[subs,]
          invD_tau<-diag(1/tausq[j,ind])
        }
        
        #stau<-stau+sum(tausq)
        #calsmnpsx<-chol2inv(chol(invconvar*tcrossprod(omesub)+invD_tau))
        #temp<-(omesub%*%Ycen*invconvar)
        calsmnpsx<-chol2inv(chol(tcrossprod(omesub)+invD_tau))
        temp<-(omesub%*%Ycen)
        LYnpsx<-calsmnpsx%*%temp
        LY[j,subs]<-mvrnorm(1,LYnpsx,Sig=(PSX[j,j]*calsmnpsx))
        
      } # end len>0
      
      
      for(k in 1:NZ){
        #k=2
        ind=!(IDY[,k]==0)
        if(mean(LY[ind,k])<LY_eps){
          LY[ind,k]<--LY[ind,k]
          Omega[k,]<--Omega[k,]
        }
      }    
      
      
    } # end of NY
    
    ##################  end of update LY ########################################################
    
    
    #Generate the unknown parameter of covariance matrix of latent factors in CFA from its conditinal distribution
    #source("Gibbs_PHI.R")
    ########  update PHI ########################################################
    
    inv.PHI<-rwish(rou.zero+N, solve(tcrossprod(Omega)+R.zero)) 
    PHI<-chol2inv(chol(inv.PHI))
    #c.inv.PHI<-chol(inv.PHI)
    
    tmp<-chol2inv(chol(sqrt(diag(diag(PHI)))))
    CPH0<-tmp%*%PHI%*%tmp
    acc<-exp((NZ+1)/2*(log(det(CPH0))-log(det(CPH))))
    acid<-acc>runif(1)
    CPH<-CPH0*acid+CPH*(1-acid)
    inv.CPH<-chol2inv(chol(CPH))
    ########  end of update PHI #################################################
    
    
    #Generate the missing reponse in CFA from its conditinal distribution
    #source("Gibbs_MISY.R")
    # if(missingness){
    # for(j in 1:NY)
    #    for(i in 1:N)
    #       if(missing_ind[j,i]==1){
    #         mean<-MU[j]+LY[j,]%*%Omega[,i]+PSX[j,-j]%*%chol2inv(chol(PSX[-j,-j]))%*%(Y[-j,i]-MU[-j]-LY[-j,]%*%Omega[,i])
    #         var<-PSX[j,j]-PSX[j,-j]%*%chol2inv(chol(PSX[-j,-j]))%*%PSX[-j,j]
    #         Y[j,i]<-rnorm(1, mean, var)
    # 
    #       }
    # }#end missingness
    
    #Save results
    if((gm>0)&&(gm%%nthin==0)){
      gm<-gm/nthin
      EPHI[gm,,]<-CPH[,]
      EPSX[gm,,]<-PSX
      #EinvPSX[gm,,]<-inv.PSX
      #EMU[gm,]<-MU
      Elambda[gm,]<-lambda
      Elamsq[gm,]<-lamsq
      ELY[gm,]<-LY[IDY!=0]
      
      
    }#end save results
    
    if(g%%1000==0)print(g)
    
    
  }#end of g MCMAX

rm=apply(ELY,FUN=mean,MAR=c(2))
rsd=apply(ELY,FUN=sd,MAR=c(2))
hpd.95=HPDinterval(mcmc(ELY),prob=.95)
hpd.90=HPDinterval(mcmc(ELY),prob=.90)   
res=cbind(rm,rsd,hpd.95,hpd.90)
fn=paste("ld_necr2_N",N,"_J",NY,".out",sep="")
write.table(res,fn,row.names=F,col.names=F,append = T)

chainpsx<-array(0,dim=c(Nrec,NY))
for(i in 1:NY){
  chainpsx[,i]<-EPSX[,i,i]  
}
rm=apply(chainpsx,FUN=mean,MAR=c(2))
rsd=apply(chainpsx,FUN=sd,MAR=c(2))
hpd.95=HPDinterval(mcmc(chainpsx),prob=.95)          
hpd.90=HPDinterval(mcmc(chainpsx),prob=.90)   
res=cbind(rm,rsd,hpd.95,hpd.90)         
fn=paste("psx_necr2_N",N,"_J",NY,".out",sep="")
write.table(res,fn,row.names=F,col.names=F,append = T)

chainphi<-array(0,dim=c(Nrec,NZ*(NZ-1)/2))
k<-1
for(i in 1:NZ){
  for(j in 1:NZ){
    if(i>j) {chainphi[,k]<-EPHI[,i,j];k<-k+1}
  }
}
rm=apply(chainphi,FUN=mean,MAR=c(2))
rsd=apply(chainphi,FUN=sd,MAR=c(2))
hpd.95=HPDinterval(mcmc(chainphi),prob=.95)          
hpd.90=HPDinterval(mcmc(chainphi),prob=.90)   
res=cbind(rm,rsd,hpd.95,hpd.90)
fn=paste("phi_necr2_N",N,"_J",NY,".out",sep="")
write.table(res,fn,row.names=F,col.names=F,append = T)

tlam=Elamsq
rm=apply(tlam,FUN=mean,MAR=c(2))
rsd=apply(tlam,FUN=sd,MAR=c(2))
hpd.95=HPDinterval(mcmc(tlam),prob=.95)          
hpd.90=HPDinterval(mcmc(tlam),prob=.90)
res=cbind(rm,rsd,hpd.95,hpd.90)   
fn=paste("lam_necr2_N",N,"_J",NY,".out",sep="")
write.table(res,fn,row.names=F,col.names=F,append = T)

#traceplot(mcmc(ELY[,1]))
#densplot(mcmc(EPHI[,2]))
print(rep)
} #end replication
proc.time()-ptm
