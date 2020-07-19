######### Revise the parts below according to dimension of measurement equation.

IDY_matrix <- function(dataset,mmvar,mmvar_loc)
{
	#Matrix indicating which element in Lambda is free to estmate (1) or fixed(0) 
	IDY_matrix<-matrix(0.0,nrow=ncol(dataset),ncol=length(mmvar))
	for (i in 1:length(mmvar))
	{
		IDY_matrix[mmvar_loc[[i]][1],i]<-9
		IDY_matrix[mmvar_loc[[i]][2:length(mmvar_loc[[i]])],i]<-1.0
	}

	IDY<-IDY_matrix

	return(IDY)
}




#if(category)
#{
#	IDD<-1:NY                         # define the position of ordinal data
	#ROAL<-c(NH-1, NH-1, NH-1, NH-1, NH-1, NH-1, NH-1, NH-1, NH-1, NH-1)      # define the last threshold to estimate
                                                 # the first and last thresholds are fixed to identification
#	ROAL<-rep(NH-1,NY)
#}
