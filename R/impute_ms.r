
impute_ms <- function(Y, NY, N, chain2, N.burn, MCMAX){

	Y_missing = chain2$Y_missing
	missing_ind = chain2$missing_ind
	
	dataset_new = t(Y)
	if(sum(missing_ind) > 0)
	{
		ismissing = 1
		missing_mean = c((N.burn+1):MCMAX)
		for(i in 1:NY)
		{
			for (j in 1:N)
			{
				if(missing_ind[i,j]==1)
				{					
					dataset_new[j,i] = mean(Y_missing[i,j,missing_mean])
				}
			}
		}
		write.table(dataset_new,'data_imputed.txt',col.names=FALSE,row.names=FALSE,sep='\t')
	}else{
		ismissing = 0
	}
	return(ismissing)
}