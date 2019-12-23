## Set different initial values for each MCMC chain

set_int_fun<-function(CIR,dataset,mmvar,mmvar_loc)
{
if (CIR == 1)
{
	LY_int<-matrix(0.0,nrow=ncol(dataset),ncol=length(mmvar))
	for (i in 1:length(mmvar))
	{
		LY_int[mmvar_loc[[i]][1],i]<-1.0
		LY_int[mmvar_loc[[i]][2:length(mmvar_loc[[i]])],i]<-1.3
	}
	LY_int[mmvar_loc[[1]][2:length(mmvar_loc[[1]])],1]<-1.1
}


if (CIR == 2)
{
	LY_int<-matrix(0.0,nrow=ncol(dataset),ncol=length(mmvar))
	for (i in 1:length(mmvar))
	{
		LY_int[mmvar_loc[[i]][1],i]<-1.0
		LY_int[mmvar_loc[[i]][2:length(mmvar_loc[[i]])],i]<-1.3
	}
}

if (CIR == 3)
{
	LY_int<-matrix(0.0,nrow=ncol(dataset),ncol=length(mmvar))
	for (i in 1:length(mmvar))
	{
		LY_int[mmvar_loc[[i]][1],i]<-1.0
		LY_int[mmvar_loc[[i]][2:length(mmvar_loc[[i]])],i]<-1.3
	}
	LY_int[mmvar_loc[[1]][2:length(mmvar_loc[[1]])],1]<-1.5
}
return(LY_int)
}


set_int_pcfa_fun<-function(CIR,dataset,mmvar,mmvar_loc,IDY)
{

	if (CIR == 1)
	{
		LY_int<-matrix(0.1,nrow=ncol(dataset),ncol=length(mmvar))
		for (i in (dim(IDY)[1]))
		{
			for (j in (dim(IDY)[2]))
			{
				if (IDY[i,j]==1)
				{
					LY_int[i,j]=0.7
				}
			}
		}
		LY_int[1,1]=1.0
	}
	if (CIR == 2)
	{
		LY_int<-matrix(0.1,nrow=ncol(dataset),ncol=length(mmvar))
		for (i in (dim(IDY)[1]))
		{
			for (j in (dim(IDY)[2]))
			{
				if (IDY[i,j]==1)
				{
					LY_int[i,j]=0.7
				}
			}
		}
	}
	if (CIR == 3)
	{
		LY_int<-matrix(0.1,nrow=ncol(dataset),ncol=length(mmvar))
		for (i in (dim(IDY)[1]))
		{
			for (j in (dim(IDY)[2]))
			{
				if (IDY[i,j]==1)
				{
					LY_int[i,j]=0.7
				}
			}
		}
		LY_int[1,1]=1.3
	}

return(LY_int)
}