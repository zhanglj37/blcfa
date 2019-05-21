
read_data<-function(filename,varnames)
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
