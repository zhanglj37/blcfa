 
normal_detect<-function(filename,varnames,usevar)
{
if(file.exists("normal.inp"))
{
 file.remove("normal.inp")
}



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
	file = paste("normal.inp", sep = ''), append = T)
cat(
	"DATA: FILE = ", filename, ";",
	file = paste("normal.inp", sep = ''), append = T)
cat(
	"\n",
	"VARIABLE:\n",
	"NAMES = ",
	file = paste("normal.inp", sep = ''), append = T)

if (var_row_num == 1)
{  
	cat(get(paste0("combine_names",var_row_num)),";\n",
		file = paste("normal.inp", sep = ''), append = T)
	
}else{
	for (i in 1:(var_row_num-1))
	{
		cat(get(paste0("combine_names",i)),"\n\t",
			file = paste("normal.inp", sep = ''), append = T)
	}
		cat(get(paste0("combine_names",var_row_num)),";\n",
			file = paste("normal.inp", sep = ''), append = T)
}	
	
cat(
	"USEV = ",
	file = paste("normal.inp", sep = ''), append = T)

if (usevar_row_num == 1)
{
	cat(get(paste0("combine_usenames",usevar_row_num)),";\n",
		file = paste("normal.inp", sep = ''), append = T)
}else{
	for (i in 1:(usevar_row_num-1))
	{
		cat(get(paste0("combine_usenames",i)),"\n\t",
			file = paste("normal.inp", sep = ''), append = T)
	}
		cat(get(paste0("combine_usenames",usevar_row_num)),";\n",
			file = paste("normal.inp", sep = ''), append = T)
} 
#if (is.numeric(ms))
#{
#	cat(
#		"missing: ALL(",
#		ms,
#		")\n\n\t",
#		file = paste("normal.inp", sep = ''), append = T)
#}
 
cat(
	"CLASSES = C(1);\n",
	"ANALYSIS:\n\t",
	"TYPE = MIXTURE;\n",
	"MODEL:\n\t",
	file = paste("normal.inp", sep = ''), append = T)

cat(
	"%OVERALL%\n\t",
	"\n\n\t",
	file = paste("normal.inp", sep = ''), append = T)
 

if (usevar_row_num == 1)
{
	cat(get(paste0("combine_usenames",usevar_row_num)),";\n",
		file = paste("normal.inp", sep = ''), append = T)
}else{
	for (i in 1:(usevar_row_num-1))
	{
		cat(get(paste0("combine_usenames",i)),"\n\t",
			file = paste("normal.inp", sep = ''), append = T)
	}
		cat(get(paste0("combine_usenames",usevar_row_num)),";\n",
			file = paste("normal.inp", sep = ''), append = T)
} 




cat(
	"\n",
	"OUTPUT: TECH12;\n",
	file = paste("normal.inp", sep = ''), append = T)
 
## run
path_detect = Sys.getenv("PATH")
path_detect = tolower(path_detect)
if (str_detect(path_detect,'mplus')){
runModels("normal.inp")

normality = try(readModels('normal.out')$tech12)
obsSkewness = normality$obsSkewness
obsKurtosis = normality$obsKurtosis

	if (max(abs(obsSkewness)) > 2 || max(abs(obsKurtosis)) > 7)
	{
		nonnormal = 1
	}else{
		nonnormal = 0
	}


}else{
	nonnormal = 0
	print('Error: Failed to run Mplus to detect non-normality, choose the maximum likelihood estimator as default')
}


return(nonnormal)

}
