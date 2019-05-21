### function---read_model----
### SYSU Lijin Zhang,190312
	
## 	Output factors and indicators
read_model<-function(myModel)
{
	mmsplit<-strsplit(myModel,'\n')
	mmsplit<-mmsplit[[1]]
	mmlength<-length(mmsplit)
	
	## Delete blank line
	for (i in 1:mmlength)
	{
		tempstr<-mmsplit[i]
		tempsplit<-strsplit(tempstr,' ')
		tempsplit<-tempsplit[[1]]
		templength<-length(tempsplit)
	
		## Delete blank line
		tempempty<-NULL
		for (j in 1:length(tempsplit))
		{
			tempempty=paste(tempempty,"")
		}
		

		if(tempstr == tempempty)
		{
			mm2split<-mmsplit[-i]
		}
	}	
	
	
	mm2length<-length(mm2split)
	mmvar<-NULL
	numw<-1  ## num of latent variables(w)
	factorname<-NULL
	
	## Process each line of characters separately
	for (i in 1:mm2length)
	{
		tempstr<-mm2split[i]
		tempsplit<-strsplit(tempstr,' ')
		tempsplit<-tempsplit[[1]]
		templength<-length(tempsplit)
	
	## fault tolerance
	## Remove multi-entered spaces
		temp2split<-tempsplit
		loc<-0
		for (j in 1:templength)
		{
			if (tempsplit[j] == "")
			{
				temp2split<-temp2split[-j+loc]
				loc<-loc+1
			}
		}	
		
	## Locate the mathematical symbols in the model to locate each variable	
		templength<-length(temp2split)
		locj<-1
		cfa_loc<-rep(0,templength)	
	
		for (j in 1:templength)
		{
	############## Need to be modified after the subsequent SEM is included here.
			if (temp2split[j] == "=~" || temp2split[j] == "+")
			{
				cfa_loc[locj]<-j
				locj=locj+1
			}
		}
		cfa_loc<-cfa_loc[1:(locj-1)]	
		
		
	## If there is a =~ symbol, the number of latent variables +1.	
		if (cfa_loc[1] != 0)
		{
			numw=numw+1
			cfa_loc1<-cfa_loc+1
			mmvar[[numw]]<-temp2split[cfa_loc1]				
		}
		factorname[(numw-1)]<-temp2split[1]
	}
	mmvar[[1]]<-factorname
	return(mmvar)
}


## Calculate the location of each indicator in the corresponding dataset
cfa_loc<-function(mmvar,data)
{
	
	mmvar_loc<-NULL
	for (i in 1:length(mmvar))
	{
		tempsplit<-strsplit(mmvar[[i]]," ")
		col_loc<-rep(0,length(mmvar[[i]]))
		for (j in 1:length(mmvar[[i]]))
		{
			tempstr<-as.character(tempsplit[j])
			for (q in 1:ncol(data))
			{		
			## Fault tolerance: capitalization
				if (tempstr ==  colnames(data)[q])
				{
					col_loc[j]<-q
				}
			}			
		}
		mmvar_loc[[i]]<-col_loc
	}
	return(mmvar_loc)
}




		
	

	
	
	
	
	
	
	

	
	