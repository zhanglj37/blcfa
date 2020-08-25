write_results<-function(MCMAX,N.burn,NZ,NY,resultlist,hpdlist,sigpsx_list,sigly_list,
						epsr,usevar,IDMU,IDY,bloutput)
{	

EmLY=resultlist$EmLY
EmMU=resultlist$EmMU
EmPHI=resultlist$EmPHI
EmPSX=resultlist$EmPSX
SELY=resultlist$SELY
SEMU=resultlist$SEMU
SEPHI=resultlist$SEPHI
SEPSX=resultlist$SEPSX
Empostp=resultlist$Empostp
HPD_LY1=hpdlist$HPD_LY1
HPD_MU1=hpdlist$HPD_MU1
HPD_PHI1=hpdlist$HPD_PHI1
HPD_PSX1=hpdlist$HPD_PSX1
OUTPSX=sigpsx_list$OUTPSX
SIGPSX=sigpsx_list$SIGPSX

NMU<-sum(IDMU)			      #number of Mu in measurement equation.
NLY<-sum(IDY!=0)				#number of free lambda need to be estimated in Lambda.

### combine lambda results---------------------------------------
### sigly.r

### combine mu results ----------------------------------------------------------
EmMU1<-matrix(EmMU[1,],NMU,1)
SEMU1<-matrix(SEMU[1,],NMU,1)
ZMU=PMU=EmMU1/SEMU1
for (i in 1:NMU)
{
	if (ZMU[i]>0)
	{
		PMU[i]<-2*(1-pnorm(ZMU[i]))
	}
	else 
	{
		PMU[i]<-2*pnorm(ZMU[i])
	}
}
MU_matrix<-cbind(EmMU1,SEMU1,PMU,HPD_MU1)
murownames<-usevar
rownames(MU_matrix)<-murownames
colnames(MU_matrix)<-c("EstMU","SeMU","p-value","HPD_lower","HPD_upper")


### combine phi results-----------------------------------------------------
EmPHI1<-matrix(EmPHI[1,],NZ,NZ)
SEPHI1<-matrix(SEPHI[1,],NZ,NZ)
ZPHI=PPHI=EmPHI[1,]/SEPHI[1,]
numphi = NZ*NZ
for (i in 1:numphi)
{
	if (ZPHI[i]>0)
	{
		PPHI[i]<-2*(1-pnorm(ZPHI[i]))
	}
	else{
		PPHI[i]<-2*pnorm(ZPHI[i])
	}
}
PPHI1=CORPHI=matrix(PPHI,NZ,NZ)
for (i in 1:NZ)
	for (j in 1:NZ)
		CORPHI[i,j] = EmPHI1[i,j]/( sqrt(EmPHI1[i,i])* sqrt(EmPHI1[j,j]) )
fname = c(paste("f", 1:NZ, sep = ""))
colnames(CORPHI) = fname


philoc=phiest=phise=phicor=phip=hpd_low=hpd_up=rep(0,(NZ*(NZ+1))/2)
k<-1
for(i in 1:NZ)
{
	for(j in i:NZ)
	{
			philoc[k]<-paste(paste0('f',i)," with ", paste0('f',j))
			phiest[k]<-EmPHI1[i,j]
			phise[k]<-SEPHI1[i,j]
			phicor[k]<-CORPHI[i,j]
			phip[k]<-PPHI1[i,j]
			hpd_low[k]<-HPD_PHI1[(j-1)*NZ + i,1]  ###取决于上下三角
			hpd_up[k]<-HPD_PHI1[(j-1)*NZ + i,2]
			k<-k+1
	}
}
OUTPHI<-cbind(phiest,phise,phicor,phip,hpd_low,hpd_up)   #### error
colnames(OUTPHI)<-c("est","se","cor","p-value","HPD_lower","HPD_upper")
rownames(OUTPHI)<-philoc


if (bloutput){
### write results----------------------------------------------------
wd_origin<-getwd()
if (file.exists('results')){
    setwd('results')
} else {
    dir.create('results')
    setwd('results')
}

write.csv(Empostp[1],file = paste('ppp.csv', sep = ''), row.names = FALSE)
 

#if (category)
#{
#    write.xlsx(t(NUM2/MCMAX),xlsxname,"acc_threshold", row.names = FALSE, col.names = FALSE)		
#	EmALPHA1<-matrix(EmALPHA[1,,],NY,piont-3)
#	SEALPHA1<-matrix(SEALPHA[1,,],NY,piont-3)
#	ALPHA_matrix2<-cbind(EmALPHA1,SEALPHA1)
#	rownames(ALPHA_matrix2)<-murownames
#	write.xlsx(ALPHA_matrix2,xlsxname,"alpha", row.names = FALSE, col.names = FALSE)		
#}
write.csv(sigly_list$OUTLY,file = paste('ly.csv', sep = ''), row.names = TRUE)
write.csv(MU_matrix,file = paste('mu.csv', sep = ''), row.names = TRUE)



write.csv(OUTPHI,file = paste('phi.csv', sep = ''), row.names = TRUE)

write.csv(CORPHI,file = paste('phi_cormatrix.csv', sep = ''), row.names = FALSE)

write.csv(OUTPSX,file = paste('psx.csv', sep = ''), row.names = TRUE)
if (SIGPSX[1] != 0 )
{
	write.csv(SIGPSX,file = paste('psx_sig.csv', sep = ''), row.names = TRUE)
}else{
	write.csv( "no sig residual correlation",file = paste('psx_sig.csv', sep = ''), row.names = TRUE)
}




# epsr graph
if (epsr[1]!='noepsr')
{
	EPSR_figure(epsr, N.burn)
}
setwd('..')
}

resultlist2 = list(ppp = Empostp, ly = round(sigly_list$OUTLY,3),
	mu = round(MU_matrix,3), phi = round(OUTPHI,3),
	psx = round(OUTPSX,3))
return(resultlist2)
}

