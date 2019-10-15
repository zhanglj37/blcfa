write_results<-function(MCMAX,NZ,NY,NLY,resultlist,hpdlist,sigpsx_list,
						epsr,mmvar,factorname,IDMU,IDY)
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
NLY<-sum(IDY)				#number of free lambda need to be estimated in Lambda.

### combine lambda results---------------------------------------
EmLY1<-matrix(EmLY[1,],NLY,1)
SELY1<-matrix(SELY[1,],NLY,1)
ZLY=PPLY=EmLY1/SELY1
for (i in 1:NLY)
{
	if (ZLY[i]>0)
	{
		PPLY[i]<-2*(1-pnorm(ZLY[i]))
	}
	else 
	{
		PPLY[i]<-2*pnorm(ZLY[i])
	}
}
LY_matrix<-cbind(EmLY1,SELY1,PPLY,HPD_LY1)
lyrownames<-c(mmvar[[1]][2])
for (i in 1:NZ)
{
	for (j in 2:length(mmvar[[i]]))
	{
		lyrownames<-c(lyrownames,mmvar[[i]][j])
	}
}
lyrownames<-lyrownames[2:(NLY+1)]
rownames(LY_matrix)<-lyrownames
colnames(LY_matrix)<-c("EstLY","SeLY","p-value","HPD_lower","HPD_upper")

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
murownames<-c(mmvar[[1]][1])
for (i in 1:NZ)
{
	for (j in 1:length(mmvar[[i]]))
	{
		murownames<-c(murownames,mmvar[[i]][j])
	}
}
murownames<-murownames[2:(NMU+1)]
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



philoc=phiest=phise=phicor=phip=hpd_low=hpd_up=rep(0,(NZ*(NZ+1))/2)
k<-1
for(i in 1:NZ)
{
	for(j in i:NZ)
	{
			philoc[k]<-paste(factorname[i]," with ", factorname[j])
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



### write results----------------------------------------------------
wd_origin<-getwd()
if (file.exists('results')){
    setwd('results')
} else {
    dir.create('results')
    setwd('results')
}

write.table(Empostp[1],file = paste('ppp.dat', sep = ''), sep = '\t', row.names = FALSE, col.names = FALSE)
 

#if (category)
#{
#    write.xlsx(t(NUM2/MCMAX),xlsxname,"acc_threshold", row.names = FALSE, col.names = FALSE)		
#	EmALPHA1<-matrix(EmALPHA[1,,],NY,piont-3)
#	SEALPHA1<-matrix(SEALPHA[1,,],NY,piont-3)
#	ALPHA_matrix2<-cbind(EmALPHA1,SEALPHA1)
#	rownames(ALPHA_matrix2)<-murownames
#	write.xlsx(ALPHA_matrix2,xlsxname,"alpha", row.names = FALSE, col.names = FALSE)		
#}
write.table(LY_matrix,file = paste('ly.dat', sep = ''), sep = '\t', row.names = TRUE, col.names = TRUE)
write.table(MU_matrix,file = paste('mu.dat', sep = ''), sep = '\t', row.names = TRUE, col.names = TRUE)



write.table(OUTPHI,file = paste('phi.dat', sep = ''), sep = '\t', row.names = TRUE, col.names = TRUE)
rownames(CORPHI)<-factorname
colnames(CORPHI)<-factorname
write.table(CORPHI,file = paste('phi_cormatrix.dat', sep = ''), sep = '\t', row.names = TRUE, col.names = TRUE)

write.table(OUTPSX,file = paste('psx.dat', sep = ''), sep = '\t', row.names = TRUE, col.names = TRUE)
if (SIGPSX != 0 )
{
	write.table(SIGPSX,file = paste('psx_sig.dat', sep = ''), sep = '\t', row.names = TRUE, col.names = TRUE)
}else{
	write.table( "no sig residual correlation",file = paste('psx_sig.dat', sep = ''), sep = '\t', row.names = TRUE, col.names = TRUE)
}


## x,file,sheet
## names需要改

# 4. 绘图并输出
mycolsi <- rainbow(ncol(epsr), s = 1, v = 1, start = 0, 
			end = max(1, ncol(epsr) - 1)/ncol(epsr), alpha = 1) 
	#每个参数一种颜色

repxlim<-c(1:(MCMAX-1))
plot(x = repxlim , y = epsr[,1], type="l", 
	xlab = "iterations", ylab = "EPSR", ylim = c(0,3), col = mycolsi[1])
for (i in 2:ncol(epsr))
{
	lines(x = repxlim, y = epsr[,i], col = mycolsi[i])
}

savePlot(filename = "EPSR",
         type ="png",
         device = dev.cur(),
         restoreConsole = TRUE)
setwd('..')

}

