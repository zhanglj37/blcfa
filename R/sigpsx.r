
sig_psx_fun<-function(NZ,NY,dataset,resultlist,hpdlist,interval,step='c')
{ 

EmPSX=resultlist$EmPSX
SEPSX=resultlist$SEPSX
PPSX1=resultlist$PPSX[1,,]
HPD_PSX1=hpdlist$HPD_PSX1


if (step == 'e')
{
	EmPSX_diagonal = SEPSX_diagonal = PPSX_diagonal = rep(0,NY)
	HPD_PSX<-array(0,dim=c(NY,2))
	k1<-1
	k2<-1
	for(i in 1:NY)
	{
		for(j in 1:NY)
		{
			if(i>=j)
			{
				if (i == j)
				{
					HPD_PSX[k2,] = HPD_PSX1[k1,]
					EmPSX_diagonal[k2] = EmPSX[1,i,j]
					SEPSX_diagonal[k2] = SEPSX[1,i,j]
					PPSX_diagonal[k2] = 2*(1-pnorm(EmPSX_diagonal[k2]/SEPSX_diagonal[k2]))
					k2 = k2 +1
				}
				k1 = k1 + 1
			}
		}
	}

	psxloc=psxest=psxse=psxcor=psxp=rep(0,NY)
	k<-1
	for(i in 1:NY)
	{
		for(j in 1:NY)
		{
			if(i==j)
			{
				psxloc[k]<-paste(colnames(dataset)[i])
				psxest[k]<-EmPSX_diagonal[k]
				psxse[k]<-SEPSX_diagonal[k]
				psxp[k]<-PPSX_diagonal[k]
				k<-k+1
			}
		}
	}
	OUTPSX<-cbind(psxest,psxse,rep(1,NY),psxp,HPD_PSX)
	colnames(OUTPSX)<-c("est","se","cor","p-value","HPD_lower","HPD_upper")
	rownames(OUTPSX)<-c(psxloc)
	sigpsx_list<-list(SIGPSX=0,OUTPSX=OUTPSX)

}else if (step == 'c')
{
	HPD_PSX3<-array(0,dim=c(NY,NY,2))
	k=1
	for(i in 1:NY)
	{
		for(j in 1:NY)
		{
			if(i>=j)
			{
				HPD_PSX3[i,j,]<-HPD_PSX1[k,]
				k=k+1
			}
		}
	}	

	#### PSX-----------------------------------------------------------
	EmPSX1<-matrix(as.vector(EmPSX[1,,]),NY,NY)
	SEPSX1<-matrix(as.vector(SEPSX[1,,]),NY,NY)
#	ZPSX=PPSX= (as.vector(EmPSX[1,,])/as.vector(SEPSX[1,,]))

	## caculate p-value
#	for (i in 1:length(ZPSX))
#	{
#		if (ZPSX[i]>0)	{ PPSX[i]<-2*(1-pnorm(ZPSX[i])) }
#		else{ PPSX[i]<-2*pnorm(ZPSX[i])	}
#	}
#	PPSX1=CORPSX=matrix(PPSX,NY,NY)
	CORPSX=PPSX1
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
		for(j in 1:NY)
		{
			if(i>=j)
			{
				psxloc[k]<-paste(colnames(dataset)[i]," with ", colnames(dataset)[j])
				psxest[k]<-EmPSX1[i,j]
				psxse[k]<-SEPSX1[i,j]
				psxcor[k]<-CORPSX[i,j]
				psxp[k]<-PPSX1[i,j]
				k<-k+1
			}
		}
	}
	OUTPSX<-cbind(psxest,psxse,psxcor,psxp,HPD_PSX1)
	colnames(OUTPSX)<-c("est","se","cor","p-value","HPD_lower","HPD_upper")
	rownames(OUTPSX)<-c(psxloc)



	### sig psx-----------------------------------------------------------
	  

	## count num of sig psx(two ways:interval or p-value)
	count_sig_interval<-function(NY,matrix)
	{
		count=0
		for(i in 1:NY)
		{
		for(j in 1:NY)
		{
			if(i>j)
			{
				if((matrix[i,j,1]*matrix[i,j,2])>0)
				{
					count=count+1
				}
			}
		}
		}
	return(count)
	}

	count_sig_p<-function(NY,matrix)
	{
		count=0
		for(i in 1:NY)
		{
		for(j in 1:NY)
		{
			if(i>j)
			{
				if(matrix[i,j] < 0.05)
				count=count+1
			}
		}
		}
	return(count)
	}
	count_sig_threshold<-function(NY,matrix)
	{
		count=0
		for(i in 1:NY)
		{
		for(j in 1:NY)
		{
			if(i>j)
			{
				if(abs(matrix[i,j]) > 0.1)
				count=count+1
			}
		}
		}
	return(count)
	}

	if (interval)
	{
		count<-count_sig_interval(NY,HPD_PSX3)
	}else{
		count<-count_sig_threshold(NY,CORPSX)
	}


	## mark sig psx 
	out_sig<-function(NY,NZ,PPSX1,CORPSX,HPD_PSX3,dataset,count)
	{

		
		sigloc=sighpdpsx=matrix(0,count,2)
		signame=sigpsxcor=sigpsxp=rep(0,count)
		k=m=1
		for(i in 1:NY)
		{
		for(j in 1:NY)
		{
			if(i>j)
			{
				if (interval)
				{
					if ((HPD_PSX3[i,j,1]*HPD_PSX3[i,j,2])>0)
					{
						signame[k]<-paste(colnames(dataset)[i]," with ",colnames(dataset)[j])
						sigpsxcor[k]<-CORPSX[i,j]
						sigpsxp[k]<-PPSX1[i,j]
						sigloc[k,1]<-i
						sigloc[k,2]<-j
						sighpdpsx[k,]<-HPD_PSX3[i,j,]
						
						k<-k+1
					}
				}else{
					if ( abs(CORPSX[i,j]) >0.15)
					{
						signame[k]<-paste(colnames(dataset)[i]," with ",colnames(dataset)[j])
						sigpsxcor[k]<-CORPSX[i,j]
						sigpsxp[k]<-PPSX1[i,j]
						sigloc[k,1]<-i
						sigloc[k,2]<-j
						sighpdpsx[k,]<-HPD_PSX3[i,j,]
					
						k<-k+1
					}			
				}
			}
		}
		}

		
		SIGPSX<-cbind(sigloc,signame,sigpsxcor,sigpsxp,sighpdpsx)
		return(SIGPSX)

	}

	if(count > 0)
	{

		if (count == 1)
		{
			SIGPSX<-as.numeric(out_sig(NY,NZ,PPSX1,CORPSX,HPD_PSX3,dataset,count)[,4:7])
			sigloc<-as.numeric(out_sig(NY,NZ,PPSX1,CORPSX,HPD_PSX3,dataset,count)[,1:2])
			sigpsxname<-out_sig(NY,NZ,PPSX1,CORPSX,HPD_PSX3,dataset,count)[,3]
			names(SIGPSX)<-sigpsxname

		}else{
			SIGPSX<-apply(out_sig(NY,NZ,PPSX1,CORPSX,HPD_PSX3,dataset,count)[,4:7],2,as.numeric)
			sigloc<-apply(out_sig(NY,NZ,PPSX1,CORPSX,HPD_PSX3,dataset,count)[,1:2],2,as.numeric)
			sigpsxname<-out_sig(NY,NZ,PPSX1,CORPSX,HPD_PSX3,dataset,count)[,3]
			colnames(SIGPSX)<-c("cor","p-value","HPD_lower","HPD_upper")
			rownames(SIGPSX)<-sigpsxname
		}
		
		sigpsx_list<-list(SIGPSX=SIGPSX,sigloc=sigloc,OUTPSX=OUTPSX,sigpsxname=sigpsxname)
	}else{
		sigpsx_list<-list(SIGPSX=count,OUTPSX=OUTPSX)
	}


}

return(sigpsx_list)

}


