
write_mplus_ml<-function(varnames,usevar,filename,sigpsx_list,sigly_list,IDY0,ismissing,myModel)
{
if(file.exists("blcfa_ml.inp"))
{
 file.remove("blcfa_ml.inp")
}

SIGPSX=sigpsx_list$SIGPSX
SIGLY=sigly_list$SIGLY

## The variable name part is split into multiple lines in order to be smaller than the 90 characters constrained by Mplus.
if (length(varnames)%%10 == 0)
{
	var_row_num<- length(varnames)%/%10
}else{
	var_row_num<- length(varnames)%/%10+1
}

for (i in 1:var_row_num)
{
	if (i == var_row_num)
	{
		tempnum = length(varnames)-10*(i-1)
		if (tempnum == 1)
		{
			assign(paste0("combine_names",i),varnames[(i-1)*10+1])		
		}else if (tempnum ==2)
		{
			assign(paste0("combine_names",i),paste(varnames[(i-1)*10+1],varnames[(i-1)*10+2]))		
		}else{
			assign(paste0("combine_names",i),paste(varnames[(i-1)*10+1],varnames[(i-1)*10+2]))
			for (j in 3:tempnum)
			{
				assign(paste0("combine_names",i),paste(get(paste0("combine_names",i)),varnames[(i-1)*10+j]))
			}	
		}			
	}else{
		assign(paste0("combine_names",i),paste(varnames[(i-1)*10+1],varnames[(i-1)*10+2]))
		for (j in 3:10)
		{
			assign(paste0("combine_names",i),paste(get(paste0("combine_names",i)),varnames[(i-1)*10+j]))
		}
	}
	
}


if (length(varnames)%%10 == 0)
{
	usevar_row_num<- length(usevar)%/%10
}else{
	usevar_row_num<- length(usevar)%/%10+1
}

for (i in 1:usevar_row_num)
{
	if (i == usevar_row_num)
	{
		tempnum = length(usevar)-10*(i-1)
		if (tempnum == 1)
		{
			assign(paste0("combine_usenames",i),usevar[(i-1)*10+1])		
		}else if (tempnum ==2)
		{
			assign(paste0("combine_usenames",i),paste(usevar[(i-1)*10+1],usevar[(i-1)*10+2]))		
		}else{
			assign(paste0("combine_usenames",i),paste(usevar[(i-1)*10+1],usevar[(i-1)*10+2]))
			for (j in 3:tempnum)
			{
				assign(paste0("combine_usenames",i),paste(get(paste0("combine_usenames",i)),usevar[(i-1)*10+j]))
			}	
		}			
	}else{
		assign(paste0("combine_usenames",i),paste(usevar[(i-1)*10+1],usevar[(i-1)*10+2]))
		for (j in 3:10)
		{
			assign(paste0("combine_usenames",i),paste(get(paste0("combine_usenames",i)),usevar[(i-1)*10+j]))
		}
	}
	
}



## generate input file
cat(
	"TITLE: Bayesian Lasso CFA\n", #_fix211_
	file = paste("blcfa_ml.inp", sep = ''), append = T)

if(ismissing == 1)
{
	var_row_num<-usevar_row_num
	cat(
		"DATA: FILE = data_imputed.txt ;",
		file = paste("blcfa_ml.inp", sep = ''), append = T)
}else{
	cat(
		"DATA: FILE = ", filename, ";",
		file = paste("blcfa_ml.inp", sep = ''), append = T)
}
cat(
	"\n",
	"VARIABLE:\n",
	"NAMES = ",
	file = paste("blcfa_ml.inp", sep = ''), append = T)


if (var_row_num == 1)
{  
	cat(get(paste0("combine_names",var_row_num)),";\n",
		file = paste("blcfa_ml.inp", sep = ''), append = T)
	
}else{
	for (i in 1:(var_row_num-1))
	{
		cat(get(paste0("combine_names",i)),"\n\t",
			file = paste("blcfa_ml.inp", sep = ''), append = T)
	}
		cat(get(paste0("combine_names",var_row_num)),";\n",
			file = paste("blcfa_ml.inp", sep = ''), append = T)
}	
	
cat(
	"USEV = ",
	file = paste("blcfa_ml.inp", sep = ''), append = T)



if (usevar_row_num == 1)
{
	cat(get(paste0("combine_usenames",usevar_row_num)),";\n",
		file = paste("blcfa_ml.inp", sep = ''), append = T)
}else{
	for (i in 1:(usevar_row_num-1))
	{
		cat(get(paste0("combine_usenames",i)),"\n\t",
			file = paste("blcfa_ml.inp", sep = ''), append = T)
	}
		cat(get(paste0("combine_usenames",usevar_row_num)),";\n",
			file = paste("blcfa_ml.inp", sep = ''), append = T)
} 

	
#cat(
#	"Define:\n",
#	file = paste("blcfa_ml.inp", sep = ''), append = T)

### 
### if (usevar_row_num == 1)
### {
### 	cat(get(paste0("combine_usenames",usevar_row_num)),";\n",
### 		file = paste("blcfa_ml.inp", sep = ''), append = T)
### }else{
### 	for (i in 1:(usevar_row_num-1))
### 	{
### 		cat(get(paste0("combine_usenames",i)),"\n\t",
### 			file = paste("blcfa_ml.inp", sep = ''), append = T)
### 	}
### 		cat(get(paste0("combine_usenames",usevar_row_num)),";\n",
### 			file = paste("blcfa_ml.inp", sep = ''), append = T)
### } 
### 


#if (is.numeric(ms))
#{
#	cat(
#		"missing: ALL(",
#		ms,
#		")\n\n\t",
#		file = paste("blcfa_ml.inp", sep = ''), append = T)
#}

cat(
	"ANALYSIS:\n\t",
	file = paste("blcfa_ml.inp", sep = ''), append = T)

nonnormal = normal_detect(filename,varnames,usevar)
if (nonnormal == 1)
{
	cat(
		"ESTIMATOR = MLM;\n\t",
		file = paste("blcfa_ml.inp", sep = ''), append = T)
}else{
	cat(
		"ESTIMATOR = ML;\n",
		file = paste("blcfa_ml.inp", sep = ''), append = T)
}

cat(
	"MODEL:\n\t",
	file = paste("blcfa_ml.inp", sep = ''), append = T)


if (is.matrix(myModel)){  #length(which(IDY0[,]==-1))>0
siglyname<-sigly_list$siglyname
siglyloc<-sigly_list$sigloc
NZ=dim(IDY0)[2]
NY=dim(siglyloc)[1]	
NLY=NZ+NY
if (SIGLY[1] != 0)
{
	
	for (i in 1:NZ)
	{
		temp = c(which(IDY0[,i]==9), siglyloc[which(siglyloc[,1]==i),2])
		cat(paste(paste0('f',i), "by"),"\t",
			file = paste("blcfa_ml.inp", sep = ''), append = T)
		for (j in 1:length(temp))
		{
			cat(usevar[temp[j]]," ",
				file = paste("blcfa_ml.inp", sep = ''), append = T)
		}

		cat(";\n\t",
			file = paste("blcfa_ml.inp", sep = ''), append = T)
	}
}else{

	for (i in 1:NZ)
	{
		temp = which(IDY0[,i]==9)
		cat(paste(paste0('f',i), "by"),"\t",
			file = paste("blcfa_ml.inp", sep = ''), append = T)

		cat(usevar[temp]," ",
				file = paste("blcfa_ml.inp", sep = ''), append = T)

		cat(";\n\t",
			file = paste("blcfa_ml.inp", sep = ''), append = T)
	}

}
}else{

model1<-gsub("#","!",myModel)
model2<-gsub("\\+","",model1)
model3<-gsub("\\=~"," BY ",model2)
model4<-gsub("\n",";\n\t",model3)

cat(
	model4,
	"\n\n\t",
	file = paste("blcfa_ml.inp", sep = ''), append = T)
 

}


if (SIGPSX[1] != 0)
{
	sigpsxname<-sigpsx_list$sigpsxname
	for (i in 1:length(sigpsxname))
	{
		cat(sigpsxname[i],";\n\t",
			file = paste("blcfa_ml.inp", sep = ''), append = T)
	}
}

cat(
	"\n",
	"OUTPUT: TECH1  STDYX;\n",
	file = paste("blcfa_ml.inp", sep = ''), append = T)
 
## run
path_detect = Sys.getenv("PATH")
path_detect = tolower(path_detect)
if (str_detect(path_detect,'mplus')){
	runModels("blcfa_ml.inp")
}else{
	print('Error: Failed to run the Mplus software, check whether the Mplus has been installed or the environment variable of your computer includes the path of Mplus')
}

#if(file.exists("Mplus Run Models.log"))
#{
# file.remove("Mplus Run Models.log")
#}

}
