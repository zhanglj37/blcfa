## Calculate the epsr chain of all parameters (lambda,mu,phi,sig_psx)
## Pseudo code is as follows:
## 1. (epsr_2_read_chain.r)
## 2. Generate a variable for each parameter, including CNUM chains of the parameter
##		such as the ly1 variable including n chains of load parameter 1.
## 3. Calculate the epsr chain for each parameter
##	3.1 Starting from the second iteration, for the kth iteration, calculate the epsr_k value of the n chains of each parameter in the k iterations
##	3.2 Save the calculated epsr_k value to the epsr chain of the parameter
##	3.3 Integrate the epsr chains of all parameters together and output
## 4. drawgraph


caculate_epsr<-function(MCMAX,N.burn,CNUM,NY,NZ,chain1,chain2)
{
	names<-c("EALPHA","ELY","EMU","EPHI","EPSX")
	names2<-c("alpha","ly","mu","phi","psx")
	names3<-c("epsralpha","epsrly","epsrmu","epsrphi","epsrpsx")



	for (i in 1:CNUM)
	{
		for (j in 2:length(names))
		{
			
			if(j == length(names))
			{
				temp<-matrix(get(paste0('chain',i))[[j-1]],MCMAX,NY*NY)
				num=(NY*(NY+1))/2
				assign(paste0(names[j],i),temp[,1:num])
			}else{
				assign(paste0(names[j],i),get(paste0('chain',i))[[j-1]])
			}
		}
	}
		
	

	
for (namesi in 2:length(names))
{

	
## 2. Generate a variable for each parameter, including CNUM chains of the parameter
	nparc<-ncol(get(paste0(names[namesi],1)))
	npar<-1:nparc
	
	for (pari in 1:nparc)
	{
		parchain<-cbind(get(paste0(names[namesi],1))[ ,pari],
						get(paste0(names[namesi],2))[ ,pari])
		
		assign(paste0(names2[namesi],pari),parchain)
	}

## 3. Calculate the epsr chain for each parameter
	## To save time, the epsr values are caculated based on the c(1:N.burn+100) iterations rather than MCMAX iterations
	assign(paste0(names3[namesi]),array(0, dim=c(N.burn,length(npar) ) ) )
	epsrpar<-get(paste0(names3[namesi]))
	for (burni in 2:(N.burn+1))
	{
		for(npari in 1:length(npar))
		{
			pari<-npar[npari]
			epsrpar[burni-1,npari]<-potscalered.mcmc(
							get(paste0(names2[namesi],pari))[1:burni,])
		}
	}
	assign(paste0(names3[namesi]),epsrpar)
}



epsr<-cbind(	get(paste0("epsrly")),
		get(paste0("epsrmu")),
		get(paste0("epsrphi")),
		get(paste0("epsrpsx")) )



## Determine whether the model converges in the N.burn iterations
epsr_reserve<-epsr[(N.burn-99):(N.burn),]
tempindex<-arrayInd(sort.list(epsr_reserve,decreasing=T)[1],
			dim(epsr_reserve))

if (epsr_reserve[tempindex]<1.2)
{
	convergence<-TRUE
}else
{
	convergence<-FALSE
}

epsrlist<-list(epsr=epsr,convergence=convergence)

return(epsrlist)
}
