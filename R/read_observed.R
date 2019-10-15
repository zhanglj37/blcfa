## Read the observed data from a text file named data.txt with dimension N*NY

read_dataset<-function(dataset)
{
Y<-t(dataset)
#Standardized the data
#Y.temp<-t(Y)   
#Y<-t(scale(Y.temp)) 

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