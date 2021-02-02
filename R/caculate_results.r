caculate_results<-function(chain2,CNUM,MCMAX,NY,NZ,N.burn,nthin,IDMU,IDY)
{
	### def_rec  ##########################################################
	EMU=chain2$EMU
	ELY=chain2$ELY
	EPHI=chain2$EPHI
	EPSX=chain2$EPSX
	Epostp=chain2$Epostp
#	EinvPSX=chain2$EinvPSX
	
	NM<-0	             #dimension of eta (q_1)
	NK<-NM+NZ	       #dimension of latent variables (eta+xi);  number of factors

	NMU<-sum(IDMU)			      #number of Mu in measurement equation.
	NLY<-sum(IDY!=0)				#number of free lambda need to be estimated in Lambda.
	#Nrec<-(MCMAX-N.burn)/nthin		#number of samples after burn-in.
	Nrec<-MCMAX/nthin #save all then extract

	EmMU<-array(0,dim=c(1,NMU))		#Store estimates of MU 
	EmLY<-array(0,dim=c(1,NLY))		#Store estimates of Lambda 
	EmPSX<-array(0,dim=c(1,NY,NY))		#Store estimates of PSX
#	EminvPSX<-array(0,dim=c(1,NY,NY))	#Store estimates of inv(PSX)
	EmPHI<-array(0,dim=c(1,(NZ*NZ)))		#Store estimates of PHI
	Emlambda<-array(0,dim=c(1,1))          #Store eetimates of shrinkage paraemter lambda
	#if (category)
	#{
	#	EALPHA<-array(0,dim=c(Nrec,NS,NH+1))#Store retained trace of thresholds
	#	EmALPHA<-array(0,dim=c(CNUM,NS,NH+1))     #Store estimates of thresholds
	#	SEALPHA<-array(0,dim=c(CNUM,NS,NH+1))     #Store standard error of thresholds
	#}

	SEMU<-array(0,dim=c(1,NMU))		#Store standard error of estimates of MU 
	SELY<-array(0,dim=c(1,NLY))		#Store standard error of estimates of Lambda
	SEPSX<-array(0,dim=c(1,NY,NY))		#Store standard error of estimates of PSX
#	SEinvPSX<-array(0,dim=c(1,NY,NY))	#Store standard error of estimates of inv(PSX)
	SEPHI<-array(0,dim=c(1,(NZ*NZ)))       #Store standard error of estimates of PHI
	SElambda<-array(0,dim=c(1,1))          #Store standard error of estimates of shrinkage paraemter lambda
	PMU<-array(0,dim=c(1,NMU))		#Store p-value of MU 
	PLY<-array(0,dim=c(1,NLY))		#Store p-value of Lambda 
	PPSX<-array(0,dim=c(1,NY,NY))		#Store p-value of PSX
#	EminvPSX<-array(0,dim=c(1,NY,NY))	#Store p-value of inv(PSX)
	PPHI<-array(0,dim=c(1,(NZ*NZ)))		#Store p-value of PHI

	Empostp<-0 #numeric(CNUM)
	
	
	
	### Record the result of the last chain for output ###########################################
	EmLY[1,]<-apply(ELY[(N.burn+1):MCMAX,],FUN=mean,MARGIN=c(2))
	EmMU[1,]<-apply(EMU[(N.burn+1):MCMAX,],FUN=mean,MARGIN=c(2))
    EmPSX[1,,]<-apply(EPSX[(N.burn+1):MCMAX,,],FUN=mean,MARGIN=c(2,3))
#    EminvPSX[1,,]<-apply(EinvPSX[(N.burn+1):MCMAX,,],FUN=mean,MARGIN=c(2,3))
	if(dim(EPHI)[2] > 1){
    	EmPHI[1,]<-apply(EPHI[(N.burn+1):MCMAX,],FUN=mean,MARGIN=c(2))
	}else{
		EmPHI[1,]<-mean(EPHI[(N.burn+1):MCMAX,])
	}
#    Emlambda[1]<-mean(Elambda[(N.burn+1):MCMAX,])  
	      
    SELY[1,]<-apply(ELY[(N.burn+1):MCMAX,],FUN=sd,MARGIN=c(2))
    SEMU[1,]<-apply(EMU[(N.burn+1):MCMAX,],FUN=sd,MARGIN=c(2))
    SEPSX[1,,]<-apply(EPSX[(N.burn+1):MCMAX,,],FUN=sd,MARGIN=c(2,3))
#    SEinvPSX[1,,]<-apply(EinvPSX[(N.burn+1):MCMAX,,],FUN=sd,MARGIN=c(2,3))
	if(dim(EPHI)[2] > 1){
	    SEPHI[1,]<-apply(EPHI[(N.burn+1):MCMAX,],FUN=sd,MARGIN=c(2))
	}else{
		SEPHI[1,]<-sd(EPHI[(N.burn+1):MCMAX,]) 
	}
#    SElambda[1]<-sd(Elambda[(N.burn+1):MCMAX,])

	for (nlyi in 1:NLY){
		if(EmLY[1,nlyi]>0){
			PLY[1,nlyi]=length(which(ELY[(N.burn+1):MCMAX,nlyi]<0))/length(ELY[(N.burn+1):MCMAX,nlyi])
		}
		else{
			PLY[1,nlyi]=length(which(ELY[(N.burn+1):MCMAX,nlyi]>0))/length(ELY[(N.burn+1):MCMAX,nlyi])
		}
	}
	for (nmui in 1:NMU){
		if(EmMU[1,nmui]>0){
			PMU[1,nmui]=length(which(EMU[(N.burn+1):MCMAX,nmui]<0))/length(EMU[(N.burn+1):MCMAX,nmui])
		}
		else{
			PMU[1,nmui]=length(which(EMU[(N.burn+1):MCMAX,nmui]>0))/length(EMU[(N.burn+1):MCMAX,nmui])
		}
	}
	for (nphii in 1:dim(PPHI)[2]){
		if(EmPHI[1,nphii]>0){
			PPHI[1,nphii]=length(which(EPHI[(N.burn+1):MCMAX,nphii]<0))/length(EPHI[(N.burn+1):MCMAX,nphii])
		}
		else{
			PPHI[1,nphii]=length(which(EPHI[(N.burn+1):MCMAX,nphii]>0))/length(EPHI[(N.burn+1):MCMAX,nphii])
		}
	}
	for (npsxi1 in 1:NY){
		for (npsxi2 in 1:NY){
		if(EmPSX[1,npsxi1,npsxi2]>0){
			PPSX[1,npsxi1,npsxi2]=length(which(EPSX[(N.burn+1):MCMAX,npsxi1,npsxi2]<0))/length(EPSX[(N.burn+1):MCMAX,npsxi1,npsxi2])
		}
		else{
			PPSX[1,npsxi1,npsxi2]=length(which(EPSX[(N.burn+1):MCMAX,npsxi1,npsxi2]>0))/length(EPSX[(N.burn+1):MCMAX,npsxi1,npsxi2])
		}
		}
	}
			
	
#	if (category)
#	{
#	    EmALPHA[1,,]=apply(EALPHA[(N.burn+1):MCMAX,],FUN=mean,MARGIN=c(2,3))  
#	    SEALPHA[1,,]=apply(EALPHA[(N.burn+1):MCMAX,],FUN=sd,MARGIN=c(2,3)) 
#	}
	
    Empostp[1]<-mean(Epostp[,])
	

		#MARGINk significant residual correlation
		
	resultlist<-list(EmLY=EmLY,EmMU=EmMU,EmPHI=EmPHI,EmPSX=EmPSX,#EminvPSX=EminvPSX,
					SELY=SELY,SEMU=SEMU,SEPHI=SEPHI,SEPSX=SEPSX,
					PMU=PMU, PLY=PLY, PPHI=PPHI, PPSX=PPSX,
					#SEinvPSX=SEinvPSX,
					Empostp=Empostp)
	return(resultlist)
}	