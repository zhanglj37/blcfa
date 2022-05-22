
sig_ly_fun<-function(dataset,resultlist,hpdlist,IDY,interval)
{ 

EmLY=resultlist$EmLY
SELY=resultlist$SELY
PLY=resultlist$PLY[1,]
HPD_LY1=hpdlist$HPD_LY1
NLY=sum(IDY!=0)
NZ=ncol(IDY)
NY=nrow(IDY)

#### LY-----------------------------------------------------------

#ZLY = PLY = EmLY[1,]/SELY[1,]

## caculate p-value
#for (i in 1:length(ZLY))
#{
#	if (ZLY[i]>0)	{ PLY[i]<-2*(1-pnorm(ZLY[i])) }
#	else{ PLY[i]<-2*pnorm(ZLY[i])	}
#}


## combine LY output

lyloc=lyest=lyse=lycor=lyp=rep(0,NLY)
k<-1
for(j in 1:NZ)
{
	for(i in 1:NY)
	{
		if(IDY[i,j]!=0)
		{
			lyloc[k]<-paste(paste0('f',j)," by ", colnames(dataset)[i])
			lyest[k]<-EmLY[1,k]
			lyse[k]<-SELY[1,k]
			lyp[k]<-PLY[k]
			k<-k+1
		}
	}
}
OUTLY<-cbind(lyest,lyse,lyp,HPD_LY1)
colnames(OUTLY)<-c("est","se","p-value","HPD_lower","HPD_upper")
rownames(OUTLY)<-c(lyloc)



### sig LY-----------------------------------------------------------
  

## count num of sig LY(two ways:interval or p-value)
count_sig_interval<-function(NLY,matrix)
{
	count=0
	for(i in 1:NLY)
	{
		if((matrix[i,1]*matrix[i,2])>0)
		{
			count=count+1
		}
	}
return(count)
}

count_sig_p<-function(NLY,matrix)
{
	count=0
	for(i in 1:NLY)
	{
		if(matrix[i]<0.05)
		{
			count=count+1
		}
	}
return(count)
}

count_sig_threshold<-function(NLY,matrix)
{
	count=0
	for(i in 1:NLY)
	{
		if(abs(matrix[i])>0.1) #should be stand
		{
			count=count+1
		}
	}
return(count)
}

if (interval)
{
	count<-count_sig_interval(NLY,HPD_LY1)
}else{
	count<-count_sig_threshold(NLY,lyest)
}


## mark sig LY 
out_sig<-function(NZ,NY,NLY,PLY,EmLY,HPD_LY1,IDY,dataset,count)
{

	sigloc=sighpdly=matrix(0,count,2)
	signame=siglyp=sigest=rep(0,count)
	k=m=1
	for(j in 1:NZ)
	{
		for(i in 1:NY)
		{
		if(IDY[i,j]!=0)
		{
			if (interval)
			{
				if ((HPD_LY1[m,1]*HPD_LY1[m,2])>0)
				{
					signame[k]<-paste(paste0('f',j)," by ", colnames(dataset)[i])
					siglyp[k]<-PLY[m]
					sigloc[k,1]<-j
					sigloc[k,2]<-i
					sighpdly[k,]<-HPD_LY1[m,]
					sigest[k]<-EmLY[1,m]
					k<-k+1
				}
			}else{
				if ( abs(EmLY[1,m]) >0.1)
				{
					signame[k]<-paste(paste0('f',i)," by ", colnames(dataset)[j])
					siglyp[k]<-PLY[m]
					sigloc[k,1]<-j
					sigloc[k,2]<-i
					sighpdly[k,]<-HPD_LY1[m,]
					sigest[k]<-EmLY[1,m]
					k<-k+1
				}			
			}
			m=m+1
		}
	}
	}

	
	SIGLY<-cbind(sigloc,signame,sigest,siglyp,sighpdly)
	return(SIGLY)

}

if(count > 0)
{

	if (count == 1)
	{
		SIGLY<-as.numeric(out_sig(NZ,NY,NLY,PLY,EmLY,HPD_LY1,IDY,dataset,count)[,4:7])
		sigloc<-as.numeric(out_sig(NZ,NY,NLY,PLY,EmLY,HPD_LY1,IDY,dataset,count)[,1:2])
		siglyname<-out_sig(NZ,NY,NLY,PLY,EmLY,HPD_LY1,IDY,dataset,count)[,3]
		names(SIGLY)<-siglyname

	}else{
		SIGLY<-apply(out_sig(NZ,NY,NLY,PLY,EmLY,HPD_LY1,IDY,dataset,count)[,4:6],2,as.numeric)
		sigloc<-apply(out_sig(NZ,NY,NLY,PLY,EmLY,HPD_LY1,IDY,dataset,count)[,1:2],2,as.numeric)
		siglyname<-out_sig(NZ,NY,NLY,PLY,EmLY,HPD_LY1,IDY,dataset,count)[,3]
		colnames(SIGLY)<-c("p-value","HPD_lower","HPD_upper")
		rownames(SIGLY)<-siglyname
	}
	
	sigly_list<-list(SIGLY=SIGLY,sigloc=sigloc,OUTLY=OUTLY,siglyname=siglyname)
}else{
	sigly_list<-list(SIGLY=count,OUTLY=OUTLY)
}




return(sigly_list)

}


