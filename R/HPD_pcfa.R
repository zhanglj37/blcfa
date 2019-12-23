

hpd_pcfa_fun<-function(chain2,NZ,NY,N,IDY)
{
	NLY<-sum(IDY!=0)				#number of free lambda need to be estimated in Lambda.
	EMU=chain2$EMU
	ELY=chain2$ELY
	EPHI=chain2$EPHI
	EPSX=chain2$EPSX
	Epostp=chain2$Epostp
	chainpsx=chain2$chainpsx
	
c<-numeric(NY*(NY+1)/2)
temp.sig<-array(0,dim=c(NY,NY))

inter<-HPDinterval(mcmc(chainpsx))   # 95% HPD interval for PSX
#for(i in 1:(NY*(NY+1)/2)){
#    if(inter[i,1]<0 && inter[i,2]>0) c[i]<-0
#    else c[i]<-1;
#}  # 1 indicates significance  

#k<-1
#for(i in 1:NY){
#   for(j in 1:NY){
#      if(i>=j){temp.sig[i,j]<-c[k];k<-k+1}
#   }
#}     
#position.sig<-which(temp.sig==1,arr.ind=T)
HPD_PSX1<-inter


temp.interval<-array(0,dim=c(NY,NY,2))
k<-1
for(i in 1:NY){
   for(j in 1:NY){
      if(i>=j){temp.interval[i,j,]<-inter[k,];k<-k+1}
   }
}     

#temp.length<-length(position.sig[,1])
#for(i in 1:temp.length){
#   a<-position.sig[i,]
#   write(t(temp.interval[a[1],a[2],]),"Result/Est/HPD_PSX_sig.txt", ncolumns=2,append=TRUE,sep="\t")             
#}


c<-numeric(NLY)
temp.sig<-array(0,dim=c(NLY))
HPD_LY1 = array(0,dim=c(2,NY,NZ))
for (i in 1:NY)
{
	for (j in 1:NZ)
	{
		inter<-HPDinterval(mcmc(ELY[,i,j])) 
		HPD_LY1[,i,j] = inter
	}
}

#inter<-HPDinterval(mcmc(ELY))   # 95% HPD interval for LY
#HPD_LY1<-inter
#write(t(inter),"Result/Est/HPD_LY.txt", ncolumns=2,append=TRUE,sep="\t")


c<-numeric(NY)
temp.sig<-array(0,dim=c(NY))
inter<-HPDinterval(mcmc(EMU))   # 95% HPD interval for MU
HPD_MU1<-inter
#write(t(inter),"Result/Est/HPD_MU.txt", ncolumns=2,append=TRUE,sep="\t")

c<-numeric(NZ*NZ)
temp.sig<-array(0,dim=c(NZ*NZ))
inter<-HPDinterval(mcmc(EPHI))   # 95% HPD interval for PHI
HPD_PHI1<-inter
#write(t(inter),"Result/Est/HPD_PHI.txt", ncolumns=2,append=TRUE,sep="\t")

hpdlist<-list(HPD_LY1 = HPD_LY1, HPD_MU1 = HPD_MU1, HPD_PHI1 = HPD_PHI1, HPD_PSX1 = HPD_PSX1)

}


