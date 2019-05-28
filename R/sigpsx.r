
sig_psx_fun<-function(NZ,NY,dataset,resultlist,hpdlist,interval_psx)
{ 

EmPSX=resultlist$EmPSX
SEPSX=resultlist$SEPSX
HPD_PSX1=hpdlist$HPD_PSX1

#### PSX-----------------------------------------------------------
EmPSX1<-matrix(as.vector(EmPSX[1,,]),NY,NY)
SEPSX1<-matrix(as.vector(SEPSX[1,,]),NY,NY)
ZPSX=PPSX= (as.vector(EmPSX[1,,])/as.vector(SEPSX[1,,]))

## caculate p-value
for (i in 1:length(ZPSX))
{
	if (ZPSX[i]>0)	{ PPSX[i]<-2*(1-pnorm(ZPSX[i])) }
	else{ PPSX[i]<-2*pnorm(ZPSX[i])	}
}
PPSX1=CORPSX=matrix(PPSX,NY,NY)

## caculate correlation
for (i in 1:NY)
{
	for (j in 1:NY)
	{
		CORPSX[i,j] = EmPSX1[i,j]/( sqrt(EmPSX1[i,i])* sqrt(EmPSX1[j,j]) )
	}
}

## combine psx output
# psxloc<-array("na", dim=c((NY*(NY+1))/2,2))
psxloc=psxest=psxse=psxcor=psxp=rep(0,(NY*(NY+1))/2)
k<-1
for(i in 1:NY)
{
	for(j in i:NY)
	{
			psxloc[k]<-paste(colnames(dataset)[i]," with ", colnames(dataset)[j])
			psxest[k]<-EmPSX1[i,j]
			psxse[k]<-SEPSX1[i,j]
			psxcor[k]<-CORPSX[i,j]
			psxp[k]<-PPSX1[i,j]
			k<-k+1
	}
}
OUTPSX<-cbind(psxest,psxse,psxcor,psxp,HPD_PSX1[1:((NY*(NY+1))/2),])
colnames(OUTPSX)<-c("est","se","cor","p-value","HPD_lower","HPD_upper")
rownames(OUTPSX)<-c(psxloc)



### sig psx-----------------------------------------------------------
HPD_PSX2<-array(HPD_PSX1,dim=c(NY,NY,2))

## count num of sig psx(two ways:interval or p-value)
count_sig_interval<-function(NY,matrix)
{
	count=0
	for(i in 1:(NY-1))
	{
		ii<-i+1
		for (j in ii:NY)
		{
			if(matrix[i,j,1]*matrix[i,j,2]>0)
			{
				count=count+1
			}
		}
	}
return(count)
}

count_sig_p<-function(NY,matrix)
{
	count=0
	for(i in 1:(NY-1))
	{
		ii<-i+1
		for (j in ii:NY)
		{
			if(matrix[i,j] < 0.05)
			count=count+1
		}
	}
return(count)
}

## mark sig psx 
out_sig<-function(NY,NZ,PPSX1,CORPSX,HPD_PSX2,dataset)
{
	if (interval_psx)
	{
		count<-count_sig_interval(NY,HPD_PSX2)
	}else{
		count<-count_sig_p(NY,PPSX1)
	}
	
	sigloc=sighpdpsx=matrix(0,count,2)
	signame=sigpsxcor=sigpsxp=rep(0,count)
	k=m=1
	for(i in 1:(NY-1))
	{
		ii<-i+1
		for(j in ii:NY)
		{
			if (interval_psx)
			{
				if (matrix[i,j,1]*matrix[i,j,2]>0)
				{
					signame[k]<-paste(colnames(dataset)[i]," with ",colnames(dataset)[j])
					sigpsxcor[k]<-CORPSX[i,j]
					sigpsxp[k]<-PPSX1[i,j]
					sigloc[k,1]<-i
					sigloc[k,2]<-j
					sighpdpsx[k,]<-HPD_PSX2[i,j,]
					
					k<-k+1
				}
			}else{
				if ( PPSX1[i,j] <0.05)
				{
					signame[k]<-paste(colnames(dataset)[i]," with ",colnames(dataset)[j])
					sigpsxcor[k]<-CORPSX[i,j]
					sigpsxp[k]<-PPSX1[i,j]
					sigloc[k,1]<-i
					sigloc[k,2]<-j
					sighpdpsx[k,]<-HPD_PSX2[i,j,]
				
					k<-k+1
				}			
			}
		}
	}

	
	SIGPSX<-cbind(sigloc,signame,sigpsxcor,sigpsxp,sighpdpsx)
	return(SIGPSX)

}
SIGPSX<-apply(out_sig(NY,NZ,PPSX1,CORPSX,HPD_PSX2,dataset)[,4:7],2,as.numeric)
sigloc<-apply(out_sig(NY,NZ,PPSX1,CORPSX,HPD_PSX2,dataset)[,1:2],2,as.numeric)
sigpsxname<-out_sig(NY,NZ,PPSX1,CORPSX,HPD_PSX2,dataset)[,3]
colnames(SIGPSX)<-c("cor","p-value","HPD_lower","HPD_upper")
rownames(SIGPSX)<-sigpsxname

sigpsx_list<-list(SIGPSX=SIGPSX,sigloc=sigloc)
return(sigpsx_list)
}