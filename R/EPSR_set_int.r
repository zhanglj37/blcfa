## Set different initial values for each MCMC chain


set_ly_int<-function(CIR,IDY0)
{

	if (CIR == 1)
	{
		LY_int<-matrix(0.1,nrow=nrow(IDY0),ncol=ncol(IDY0))
		for (i in 1:(dim(IDY0)[1]))
		{
			for (j in 1:(dim(IDY0)[2]))
			{
				if (IDY0[i,j]==9)
				{
					LY_int[i,j]=1.0
				}			
				if (IDY0[i,j]==1)
				{
					LY_int[i,j]=1.0
				}
				if (IDY0[i,j]==-1)
				{
					LY_int[i,j]=0.1
				}
				if (IDY0[i,j]==0)
				{
					LY_int[i,j]=0.0
				}
			}
		}
	}
	if (CIR == 2)
	{
		LY_int<-matrix(0.1,nrow=nrow(IDY0),ncol=ncol(IDY0))
		for (i in 1:(dim(IDY0)[1]))
		{
			for (j in 1:(dim(IDY0)[2]))
			{
				if (IDY0[i,j]==9)
				{
					LY_int[i,j]=1.0
				}	
				if (IDY0[i,j]==1)
				{
					LY_int[i,j]=0.9
				}
				if (IDY0[i,j]==-1)
				{
					LY_int[i,j]=0.1
				}
				if (IDY0[i,j]==0)
				{
					LY_int[i,j]=0.0
				}
			}
		}
	}
	if (CIR == 3)
	{
		LY_int<-matrix(0.1,nrow=nrow(IDY0),ncol=ncol(IDY0))
		for (i in 1:(dim(IDY0)[1]))
		{
			for (j in 1:(dim(IDY0)[2]))
			{
				if (IDY0[i,j]==9)
				{
					LY_int[i,j]=1.0
				}	
				if (IDY0[i,j]==1)
				{
					LY_int[i,j]=1.1
				}
				if (IDY0[i,j]==-1)
				{
					LY_int[i,j]=0.1
				}
				if (IDY0[i,j]==0)
				{
					LY_int[i,j]=0.0
				}
			}
		}
	}

return(LY_int)
}


set_ly_int_psx <-function (CIR,dataset,mmvar,mmvar_loc)
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

