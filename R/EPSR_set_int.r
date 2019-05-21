## Set different initial values for each MCMC chain

set_int_fun<-function(CIR,dataset,mmvar,mmvar_loc)
{
if (CIR == 1)
{
	LY_int<-matrix(0.0,nrow=ncol(dataset),ncol=length(mmvar))
	for (i in 1:length(mmvar))
	{
		LY_int[mmvar_loc[[i]][1],i]<-1.0
		LY_int[mmvar_loc[[i]][2:length(mmvar_loc[[i]])],i]<-1.5
	}
	LY_int[mmvar_loc[[1]][2:length(mmvar_loc[[1]])],1]<-1.1
}


if (CIR == 2)
{
	LY_int<-matrix(0.0,nrow=ncol(dataset),ncol=length(mmvar))
	for (i in 1:length(mmvar))
	{
		LY_int[mmvar_loc[[i]][1],i]<-1.0
		LY_int[mmvar_loc[[i]][2:length(mmvar_loc[[i]])],i]<-1.5
	}
}

if (CIR == 3)
{
	LY_int<-matrix(0.0,nrow=ncol(dataset),ncol=length(mmvar))
	for (i in 1:length(mmvar))
	{
		LY_int[mmvar_loc[[i]][1],i]<-1.0
		LY_int[mmvar_loc[[i]][2:length(mmvar_loc[[i]])],i]<-1.5
	}
	LY_int[mmvar_loc[[1]][2:length(mmvar_loc[[1]])],1]<-1.3
}
return(LY_int)
}