
read_data<-function(filename,varnames,usevar)
{	
	dataset_origin <- read.table(file = filename, col.names = varnames)
		## original data
	loc<-c(1:length(usevar))
	for (i in 1:length(usevar))
	{
		loc[i]<-which(varnames==usevar[i])
	} 
	## Select the data you need to use
	dataset<-dataset_origin[loc]

	return(dataset)
}

## Read the observed data from a text file named data.txt with dimension N*NY

read_data2<-function(dataset)
{
	Y<-t(dataset)
	#Standardized the data
	Y.temp<-t(Y)   
	Y<-t(scale(Y.temp)) 

	return(Y)
}


 
	  
#if (category)
#{
#	Z<-Y[IDD,]
	#read Z. Dimension is NS*N, each row contains n samples of one response variable.

	#Calculate the total number of observation for ordered categorical variable in each categories
#	NAZ<-array(0,dim=c(NS, NH))
#	for(j in 1:NS)
#		for(i in 1:N){
#			k<-Z[j,i]
#			NAZ[j,k]<-NAZ[j,k]+1
#		}
#}

mark_na <- function(N, NY, dataset, ms){
	dataset_noms = dataset
	for(i in 1:N)
	{
		for(j in 1:NY)
		{
			if(dataset[i,j] == ms) 
			{
				dataset_noms[i,j] = NA
			}
		}
	}
	return(dataset_noms)
}